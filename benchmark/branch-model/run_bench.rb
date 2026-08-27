#!/usr/bin/env ruby

require "digest"
require "fileutils"
require "json"
require "optparse"
require "time"

options = {
  runs: 5,
  cpu: 4,
  warmup: true,
  output: File.expand_path("../../build/benchmark/branch-model/#{Time.now.strftime("%Y-%m-%d")}", __dir__)
}

OptionParser.new do |parser|
  parser.banner = "Usage: run_bench.rb --baseline-binary PATH --variant-binary PATH [options]"
  parser.on("--baseline-binary PATH", "Direct-coupling DynamicSimulation binary") { |value| options[:baseline_binary] = File.expand_path(value) }
  parser.on("--variant-binary PATH", "Four-current-variable DynamicSimulation binary") { |value| options[:variant_binary] = File.expand_path(value) }
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

case_sources = {
  "ACTIVSg10k" => File.join(repo_root, "cases/PhasorDynamics/ACTIVSg10k/ACTIVSg10k.case.json"),
  "ACTIVSg2000" => File.join(repo_root, "examples/PhasorDynamics/Large/Texas/texas.case.json")
}

solver_paths = {}
case_sources.each do |case_name, source|
  data = JSON.parse(File.read(source))
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
    branches states jacobian_nnz simulation_wall_seconds
    steps residual_evals jacobian_evals linear_solver_setups
    residual_calls residual_seconds residual_model_seconds
    jacobian_calls jacobian_seconds jacobian_model_seconds
    linear_setup_calls linear_setup_seconds linear_solve_calls linear_solve_seconds
    system_residual_calls system_residual_seconds branch_residual_seconds
    system_jacobian_calls system_jacobian_seconds branch_jacobian_seconds
  ]
  missing = required.reject { |key| metrics.key?(key) }
  raise "missing required metrics in #{path}: #{missing.join(", ")}" unless missing.empty?

  metrics
end

def parse_time(path)
  File.readlines(path, chomp: true).each_with_object({}) do |line, metrics|
    key, value = line.split("=", 2)
    metrics[key] = numeric(value) if key && value
  end
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
cells = [
  {"label" => "direct", "case" => "ACTIVSg10k", "binary" => options[:baseline_binary]},
  {"label" => "four_currents", "case" => "ACTIVSg10k", "binary" => options[:variant_binary]},
  {"label" => "direct", "case" => "ACTIVSg2000", "binary" => options[:baseline_binary]},
  {"label" => "four_currents", "case" => "ACTIVSg2000", "binary" => options[:variant_binary]}
]

def run_cell(cell, solver_path, raw_dir, name, time_format, taskset)
  stdout_path = File.join(raw_dir, "#{name}.stdout.txt")
  stderr_path = File.join(raw_dir, "#{name}.stderr.txt")
  time_path = File.join(raw_dir, "#{name}.time.txt")
  command = ["/usr/bin/time", "-f", time_format, "-o", time_path, "--", *taskset, cell.fetch("binary"), solver_path]

  start = Process.clock_gettime(Process::CLOCK_MONOTONIC)
  success = system({"LC_ALL" => "C"}, *command, out: stdout_path, err: stderr_path)
  runner_wall = Process.clock_gettime(Process::CLOCK_MONOTONIC) - start
  raise "benchmark failed (#{name}); see #{stderr_path} and #{stdout_path}" unless success

  metrics = parse_output(stdout_path).merge(parse_time(time_path))
  metrics["runner_wall_seconds"] = runner_wall
  add_derived_metrics!(metrics)
  metrics
end

if options[:warmup]
  cells.each_with_index do |cell, index|
    name = format("warmup-%02d-%s-%s", index + 1, cell.fetch("label"), cell.fetch("case"))
    warn "running #{name}"
    run_cell(cell, solver_paths.fetch(cell.fetch("case")), raw_dir, name, time_format, taskset)
  end
end

records = []
options[:runs].times do |round|
  cells.rotate(round).each do |cell|
    name = format("trial-%02d-%s-%s", round + 1, cell.fetch("label"), cell.fetch("case"))
    warn "running #{name}"
    metrics = run_cell(cell, solver_paths.fetch(cell.fetch("case")), raw_dir, name, time_format, taskset)
    records << {
      "trial" => round + 1,
      "label" => cell.fetch("label"),
      "case" => cell.fetch("case"),
      "metrics" => metrics
    }
  end
end

metadata = {
  "created_at" => Time.now.iso8601,
  "runs" => options[:runs],
  "warmup" => options[:warmup],
  "cpu" => options[:cpu],
  "baseline_ref" => options[:baseline_ref],
  "variant_ref" => options[:variant_ref],
  "baseline_binary" => options[:baseline_binary],
  "variant_binary" => options[:variant_binary],
  "baseline_sha256" => baseline_sha256,
  "variant_sha256" => variant_sha256,
  "case_sources" => case_sources,
  "solver_paths" => solver_paths
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
