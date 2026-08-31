#!/usr/bin/env ruby

require "digest"
require "fileutils"
require "json"
require "optparse"
require "time"

BASELINE_LABEL = "direct_injection"
VARIANT_LABEL = "internal_current_variables"

MODEL_CLASS_BY_METRIC = {
  "branches" => "Branch",
  "regca" => "Regca",
  "loadz" => "LoadZ",
  "loadzip" => "LoadZIP",
  "genrou" => "Genrou",
  "gensal" => "Gensal",
  "genclassical" => "GenClassical"
}.freeze

PROFILE_FAMILIES = %w[regca branch loadz loadzip genrou gensal genclassical].freeze

EXPECTED_STATE_INCREASES = {
  "ACTIVSg10k" => 63_180,
  "ACTIVSg2000" => 16_702,
  "WECC240" => 2_276
}.freeze

options = {
  runs: 5,
  cpu: 4,
  warmup: true,
  output: File.expand_path("../../build/benchmark/network-current/#{Time.now.strftime("%Y-%m-%d")}", __dir__)
}

OptionParser.new do |parser|
  parser.banner = "Usage: run_bench.rb --baseline-binary PATH --variant-binary PATH [options]"
  parser.on("--baseline-binary PATH", "Direct-injection DynamicSimulation binary") { |value| options[:baseline_binary] = File.expand_path(value) }
  parser.on("--variant-binary PATH", "Internal-current-variable DynamicSimulation binary") { |value| options[:variant_binary] = File.expand_path(value) }
  parser.on("--baseline-ref REF", "Baseline git ref recorded in metadata") { |value| options[:baseline_ref] = value }
  parser.on("--variant-ref REF", "Variant git ref recorded in metadata") { |value| options[:variant_ref] = value }
  parser.on("--output PATH", "Generated input and result directory") { |value| options[:output] = File.expand_path(value) }
  parser.on("--runs N", Integer, "Measured trials per cell (default: 5)") { |value| options[:runs] = value }
  parser.on("--cpu N", Integer, "Logical CPU used by taskset (default: 4)") { |value| options[:cpu] = value }
  parser.on("--no-warmup", "Do not run discarded warm-ups") { options[:warmup] = false }
end.parse!

%i[baseline_binary variant_binary].each do |key|
  path = options[key]
  raise "missing --#{key.to_s.tr("_", "-")}" unless path
  raise "binary is not executable: #{path}" unless File.executable?(path)
end
raise "--runs must be positive" unless options[:runs].positive?
raise "--cpu must be nonnegative" if options[:cpu].negative?
raise "/usr/bin/taskset is required" unless File.executable?("/usr/bin/taskset")
raise "/usr/bin/time is required" unless File.executable?("/usr/bin/time")

baseline_sha256 = Digest::SHA256.file(options[:baseline_binary]).hexdigest
variant_sha256 = Digest::SHA256.file(options[:variant_binary]).hexdigest
raise "baseline and variant binaries are identical" if baseline_sha256 == variant_sha256

repo_root = File.expand_path("../..", __dir__)
output_dir = options[:output]
raise "refusing to overwrite existing trials: #{output_dir}" if File.exist?(File.join(output_dir, "trials.json"))

inputs_dir = File.join(output_dir, "inputs")
raw_dir = File.join(output_dir, "raw")
FileUtils.mkdir_p(inputs_dir)
FileUtils.mkdir_p(raw_dir)

def remove_monitoring!(value)
  case value
  when Hash
    value.delete("mon")
    value.delete("monitors")
    value.each_value { |child| remove_monitoring!(child) }
  when Array
    value.each { |child| remove_monitoring!(child) }
  end
end

def device_counts(data)
  classes = data.fetch("devices").each_with_object(Hash.new(0)) do |device, counts|
    counts[device.fetch("class")] += 1
  end
  MODEL_CLASS_BY_METRIC.transform_values { |class_name| classes[class_name] }
end

case_sources = {
  "ACTIVSg10k" => File.join(repo_root, "cases/PhasorDynamics/ACTIVSg10k/ACTIVSg10k.case.json"),
  "ACTIVSg2000" => File.join(repo_root, "examples/PhasorDynamics/Large/Texas/texas.case.json"),
  "WECC240" => File.join(repo_root, "cases/PhasorDynamics/WECC240/WECC240.case.json")
}

solver_paths = {}
case_sha256s = {}
solver_sha256s = {}
case_device_counts = {}
case_sources.each do |case_name, source|
  raise "case source does not exist: #{source}" unless File.file?(source)

  data = JSON.parse(File.read(source))
  counts = device_counts(data)
  computed_increase = 4 * counts.fetch("branches") + 2 * %w[loadz loadzip genrou gensal].sum { |key| counts.fetch(key) }
  expected_increase = EXPECTED_STATE_INCREASES.fetch(case_name)
  unless computed_increase == expected_increase
    raise "#{case_name} device counts imply a #{computed_increase}-state increase, expected #{expected_increase}"
  end

  case_device_counts[case_name] = counts
  remove_monitoring!(data)
  generated_case = File.join(inputs_dir, "#{case_name}.case.json")
  File.write(generated_case, JSON.pretty_generate(data))

  solver = {
    "system_model_file" => generated_case,
    "dt_monitor" => 0.0,
    "tmax" => 10.0,
    "rel_tol" => 1.0e-5,
    "abs_tol" => 1.0e-7,
    "dt_fixed" => 0.0,
    "max_steps" => 1_000_000,
    "consistent_ic_type" => "ya_ydp",
    "events" => [
      {"time" => 1.0, "type" => "fault_on", "element_id" => 0},
      {"time" => 1.1, "type" => "fault_off", "element_id" => 0}
    ]
  }
  solver_path = File.join(inputs_dir, "#{case_name}.solver.json")
  File.write(solver_path, JSON.pretty_generate(solver))
  solver_paths[case_name] = solver_path
  case_sha256s[case_name] = Digest::SHA256.file(generated_case).hexdigest
  solver_sha256s[case_name] = Digest::SHA256.file(solver_path).hexdigest
end

def numeric(value)
  return value.to_i if value.match?(/\A[-+]?\d+\z/)
  Float(value)
end

def parse_output(path)
  metrics = {}
  block_markers = Hash.new(0)
  record_metric = lambda do |key, value|
    raise "duplicate metric #{key} in #{path}" if metrics.key?(key)
    metrics[key] = value
  end
  labels = {
    "Steps" => "steps",
    "Residual evals" => "residual_evals",
    "Jacobian evals" => "jacobian_evals",
    "Linear solver setups" => "linear_solver_setups",
    "Error test failures" => "error_test_failures",
    "Nonlinear iterations" => "nonlinear_iterations",
    "Nonlinear convergence failures" => "nonlinear_convergence_failures"
  }

  File.foreach(path) do |line|
    if (match = line.match(/\A(GRIDKIT_SYSTEM|APPLICATION_PROFILE|GRIDKIT_PROFILE|SYSTEM_PROFILE)_(BEGIN|END)\s*\z/))
      block_markers[match[0].strip] += 1
    elsif (match = line.match(/Complete in ([0-9.eE+-]+) seconds/))
      record_metric.call("complete_cpu_seconds", match[1].to_f)
    elsif (match = line.match(/\A([a-z][a-z0-9_]*)=([-+0-9.eE]+)\s*\z/))
      record_metric.call(match[1], numeric(match[2]))
    elsif (match = line.match(/\A\s*(.+?)\s*:\s*(\d+)\s*\z/))
      key = labels[match[1]]
      record_metric.call(key, match[2].to_i) if key
    end
  end

  %w[GRIDKIT_SYSTEM APPLICATION_PROFILE GRIDKIT_PROFILE SYSTEM_PROFILE].each do |block|
    %w[BEGIN END].each do |boundary|
      marker = "#{block}_#{boundary}"
      raise "expected one #{marker} in #{path}" unless block_markers[marker] == 1
    end
  end

  required = %w[
    buses branches regca loadz loadzip genrou gensal genclassical
    states differential_variables algebraic_variables jacobian_nnz
    complete_cpu_seconds simulation_wall_seconds application_wall_seconds
    steps residual_evals jacobian_evals linear_solver_setups error_test_failures
    nonlinear_iterations nonlinear_convergence_failures
    residual_calls residual_seconds residual_model_seconds
    jacobian_calls jacobian_seconds jacobian_model_seconds
    linear_setup_calls linear_setup_seconds linear_solve_calls linear_solve_seconds
    system_residual_calls system_residual_seconds
    system_jacobian_calls system_jacobian_seconds
  ]
  PROFILE_FAMILIES.each do |family|
    required << "#{family}_residual_seconds"
    required << "#{family}_jacobian_seconds"
  end
  missing = required.reject { |key| metrics.key?(key) }
  raise "missing required metrics in #{path}: #{missing.join(", ")}" unless missing.empty?

  metrics
end

def parse_time(path)
  metrics = File.readlines(path, chomp: true).each_with_object({}) do |line, values|
    key, value = line.split("=", 2)
    values[key] = numeric(value) if key && value
  end
  required = %w[
    time_wall_seconds user_seconds sys_seconds maxrss_kb major_page_faults
    minor_page_faults voluntary_context_switches involuntary_context_switches
  ]
  missing = required.reject { |key| metrics.key?(key) }
  raise "missing GNU time metrics in #{path}: #{missing.join(", ")}" unless missing.empty?
  metrics
end

def safe_ratio(numerator, denominator, scale = 1.0)
  return nil unless numerator && denominator && denominator.to_f.nonzero?
  numerator.to_f / denominator.to_f * scale
end

def add_derived_metrics!(metrics)
  metrics["nnz_per_state"] = safe_ratio(metrics["jacobian_nnz"], metrics["states"])
  metrics["residual_microseconds_per_call"] = safe_ratio(metrics["residual_seconds"], metrics["residual_calls"], 1.0e6)
  metrics["jacobian_microseconds_per_call"] = safe_ratio(metrics["jacobian_seconds"], metrics["jacobian_calls"], 1.0e6)
  metrics["linear_setup_microseconds_per_call"] = safe_ratio(metrics["linear_setup_seconds"], metrics["linear_setup_calls"], 1.0e6)
  metrics["linear_solve_microseconds_per_call"] = safe_ratio(metrics["linear_solve_seconds"], metrics["linear_solve_calls"], 1.0e6)

  if metrics.values_at("consistent_ic_seconds", "ida_solve_seconds").all?
    metrics["solver_envelope_seconds"] = metrics["consistent_ic_seconds"] + metrics["ida_solve_seconds"]
    nested = %w[residual_seconds jacobian_seconds linear_setup_seconds linear_solve_seconds].sum { |key| metrics.fetch(key, 0.0) }
    metrics["sundials_internal_seconds"] = metrics["solver_envelope_seconds"] - nested
  end

  residual_parts = %w[residual_input_seconds residual_model_seconds residual_output_seconds].sum { |key| metrics.fetch(key, 0.0) }
  metrics["residual_callback_other_seconds"] = metrics["residual_seconds"] - residual_parts if metrics["residual_seconds"]

  jacobian_parts = %w[jacobian_input_seconds jacobian_model_seconds jacobian_zero_seconds jacobian_structure_seconds jacobian_values_seconds].sum { |key| metrics.fetch(key, 0.0) }
  metrics["jacobian_callback_other_seconds"] = metrics["jacobian_seconds"] - jacobian_parts if metrics["jacobian_seconds"]

  metrics.delete_if { |_key, value| value.nil? }
end

def validate_device_counts!(case_name, metrics, expected_counts)
  mismatches = expected_counts.filter_map do |key, expected|
    actual = metrics[key]
    "#{key}=#{actual.inspect} (expected #{expected})" unless actual == expected
  end
  return if mismatches.empty?
  raise "#{case_name} model-count mismatch: #{mismatches.join(", ")}"
end

def validate_state_deltas!(records)
  records.group_by { |record| [record.fetch("trial"), record.fetch("case")] }.each do |(trial, case_name), paired|
    by_label = paired.to_h { |record| [record.fetch("label"), record.fetch("metrics")] }
    missing = [BASELINE_LABEL, VARIANT_LABEL].reject { |label| by_label.key?(label) }
    raise "trial #{trial} #{case_name} is missing formulations: #{missing.join(", ")}" unless missing.empty?

    baseline_states = by_label.fetch(BASELINE_LABEL).fetch("states")
    variant_states = by_label.fetch(VARIANT_LABEL).fetch("states")
    actual = variant_states - baseline_states
    expected = EXPECTED_STATE_INCREASES.fetch(case_name)
    next if actual == expected
    raise "trial #{trial} #{case_name} has a #{actual}-state increase (#{baseline_states} -> #{variant_states}), expected #{expected}"
  end
end

time_format = [
  "time_wall_seconds=%e",
  "user_seconds=%U",
  "sys_seconds=%S",
  "maxrss_kb=%M",
  "major_page_faults=%F",
  "minor_page_faults=%R",
  "voluntary_context_switches=%w",
  "involuntary_context_switches=%c"
].join("\\n")

taskset = ["/usr/bin/taskset", "-c", options[:cpu].to_s]
cells = case_sources.keys.flat_map do |case_name|
  [
    {"label" => BASELINE_LABEL, "case" => case_name, "binary" => options[:baseline_binary]},
    {"label" => VARIANT_LABEL, "case" => case_name, "binary" => options[:variant_binary]}
  ]
end

def run_cell(cell, solver_path, raw_dir, name, time_format, taskset, expected_counts)
  stdout_path = File.join(raw_dir, "#{name}.stdout.txt")
  stderr_path = File.join(raw_dir, "#{name}.stderr.txt")
  time_path = File.join(raw_dir, "#{name}.time.txt")
  command = ["/usr/bin/time", "-f", time_format, "-o", time_path, "--", *taskset, cell.fetch("binary"), solver_path]

  start = Process.clock_gettime(Process::CLOCK_MONOTONIC)
  success = system({"LC_ALL" => "C"}, *command, out: stdout_path, err: stderr_path, chdir: File.dirname(solver_path))
  runner_wall = Process.clock_gettime(Process::CLOCK_MONOTONIC) - start
  raise "benchmark failed (#{name}); see #{stderr_path} and #{stdout_path}" unless success

  metrics = parse_output(stdout_path).merge(parse_time(time_path))
  validate_device_counts!(cell.fetch("case"), metrics, expected_counts)
  metrics["runner_wall_seconds"] = runner_wall
  add_derived_metrics!(metrics)
  metrics
end

if options[:warmup]
  cells.each_with_index do |cell, index|
    name = format("warmup-%02d-%s-%s", index + 1, cell.fetch("label"), cell.fetch("case"))
    warn "running #{name}"
    run_cell(cell, solver_paths.fetch(cell.fetch("case")), raw_dir, name, time_format, taskset, case_device_counts.fetch(cell.fetch("case")))
  end
end

records = []
options[:runs].times do |round|
  cells.rotate(round).each do |cell|
    name = format("trial-%02d-%s-%s", round + 1, cell.fetch("label"), cell.fetch("case"))
    warn "running #{name}"
    metrics = run_cell(cell, solver_paths.fetch(cell.fetch("case")), raw_dir, name, time_format, taskset, case_device_counts.fetch(cell.fetch("case")))
    records << {
      "trial" => round + 1,
      "label" => cell.fetch("label"),
      "case" => cell.fetch("case"),
      "metrics" => metrics
    }
  end
end
validate_state_deltas!(records)

metadata = {
  "created_at" => Time.now.iso8601,
  "runs" => options[:runs],
  "warmup" => options[:warmup],
  "cpu" => options[:cpu],
  "baseline_label" => BASELINE_LABEL,
  "variant_label" => VARIANT_LABEL,
  "baseline_ref" => options[:baseline_ref],
  "variant_ref" => options[:variant_ref],
  "baseline_binary" => options[:baseline_binary],
  "variant_binary" => options[:variant_binary],
  "baseline_sha256" => baseline_sha256,
  "variant_sha256" => variant_sha256,
  "case_sources" => case_sources,
  "case_sha256s" => case_sha256s,
  "solver_paths" => solver_paths,
  "solver_sha256s" => solver_sha256s,
  "device_counts" => case_device_counts,
  "expected_state_increases" => EXPECTED_STATE_INCREASES
}

File.write(File.join(output_dir, "trials.json"), JSON.pretty_generate({"metadata" => metadata, "records" => records}))

metric_keys = records.flat_map { |record| record.fetch("metrics").keys }.uniq.sort
File.open(File.join(output_dir, "trials.tsv"), "w") do |file|
  file.puts (["label", "case", "trial"] + metric_keys).join("\t")
  records.each do |record|
    metrics = record.fetch("metrics")
    file.puts ([record.fetch("label"), record.fetch("case"), record.fetch("trial")] + metric_keys.map { |key| metrics[key] }).join("\t")
  end
end

def median(values)
  sorted = values.sort
  middle = sorted.length / 2
  sorted.length.odd? ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2.0
end

File.open(File.join(output_dir, "summary.tsv"), "w") do |file|
  file.puts %w[label case metric median min max mad samples].join("\t")
  records.group_by { |record| [record.fetch("label"), record.fetch("case")] }.sort.each do |(label, case_name), group|
    metric_keys.each do |key|
      values = group.filter_map { |record| record.fetch("metrics")[key] }.map(&:to_f)
      next if values.empty?
      center = median(values)
      mad = median(values.map { |value| (value - center).abs })
      file.puts [label, case_name, key, center, values.min, values.max, mad, values.length].join("\t")
    end
  end
end

puts "wrote #{records.length} measured trials to #{output_dir}"
