#!/usr/bin/env ruby

require "digest"
require "fileutils"
require "json"
require "open3"
require "optparse"
require "rbconfig"
require "time"

FLAMEGRAPH_COMMIT = "41fee1f99f9276008b7cd112fca19dc3ea84ac32"
BASELINE_LABEL = "direct_injection"
VARIANT_LABEL = "internal_current_variables"

def find_executable(name)
  ENV.fetch("PATH", "").split(File::PATH_SEPARATOR).each do |directory|
    candidate = File.join(directory, name)
    return candidate if File.executable?(candidate) && !File.directory?(candidate)
  end
  nil
end

options = {
  captures: 3,
  cases: [],
  max_unknown_percent: 5.0
}

OptionParser.new do |parser|
  parser.banner = "Usage: run_flamegraphs.rb --trials PATH --flamegraph-dir PATH [options]"
  parser.on("--trials PATH", "trials.json produced by run_bench.rb") { |value| options[:trials] = File.expand_path(value) }
  parser.on("--flamegraph-dir PATH", "External FlameGraph checkout at the required commit") { |value| options[:flamegraph_dir] = File.expand_path(value) }
  parser.on("--perf PATH", "perf executable (default: search PATH)") { |value| options[:perf] = File.expand_path(value) }
  parser.on("--output PATH", "Profile output directory (default: beside trials.json)") { |value| options[:output] = File.expand_path(value) }
  parser.on("--cpu N", Integer, "Logical CPU (default: value recorded in trials.json)") { |value| options[:cpu] = value }
  parser.on("--captures N", Integer, "Captured runs per formulation/case (default: 3)") { |value| options[:captures] = value }
  parser.on("--case NAME", "Profile one recorded case; may be repeated") { |value| options[:cases] << value }
  parser.on("--max-unknown-percent N", Float, "Maximum sample share containing an unknown frame (default: 5.0)") { |value| options[:max_unknown_percent] = value }
end.parse!

raise "missing --trials" unless options[:trials]
raise "trials file does not exist: #{options[:trials]}" unless File.file?(options[:trials])
raise "missing --flamegraph-dir" unless options[:flamegraph_dir]
raise "FlameGraph directory does not exist: #{options[:flamegraph_dir]}" unless File.directory?(options[:flamegraph_dir])
raise "--captures must be positive" unless options[:captures].positive?
unless options[:max_unknown_percent].between?(0.0, 100.0)
  raise "--max-unknown-percent must be between 0 and 100"
end
raise "/usr/bin/taskset is required" unless File.executable?("/usr/bin/taskset")

perf = options[:perf] || find_executable("perf")
unless perf && File.executable?(perf)
  raise "perf executable not found; install a perf build matching this Linux/WSL kernel or pass --perf PATH. No profile captures were started."
end

repo_root = File.expand_path("../..", __dir__)
allowed_output_root = File.join(repo_root, "build", "benchmark", "network-current")
options[:output] ||= File.join(File.dirname(options[:trials]), "flamegraphs")
output_dir = File.expand_path(options[:output])
unless output_dir == allowed_output_root || output_dir.start_with?(allowed_output_root + File::SEPARATOR)
  raise "profile output must be under #{allowed_output_root}: #{output_dir}"
end
raise "refusing to overwrite existing profile metadata: #{output_dir}" if File.exist?(File.join(output_dir, "profiles.json"))

trial_document = JSON.parse(File.read(options[:trials]))
metadata = trial_document.fetch("metadata")
unless metadata.fetch("baseline_label") == BASELINE_LABEL && metadata.fetch("variant_label") == VARIANT_LABEL
  raise "trials.json does not use the expected #{BASELINE_LABEL}/#{VARIANT_LABEL} formulation labels"
end

binaries = {
  BASELINE_LABEL => {
    "path" => File.expand_path(metadata.fetch("baseline_binary")),
    "sha256" => metadata.fetch("baseline_sha256")
  },
  VARIANT_LABEL => {
    "path" => File.expand_path(metadata.fetch("variant_binary")),
    "sha256" => metadata.fetch("variant_sha256")
  }
}
binaries.each do |label, binary|
  path = binary.fetch("path")
  raise "recorded #{label} binary is not executable: #{path}" unless File.executable?(path)
  actual = Digest::SHA256.file(path).hexdigest
  expected = binary.fetch("sha256")
  raise "recorded #{label} binary hash changed: #{actual}, expected #{expected}" unless actual == expected
end

solver_paths = metadata.fetch("solver_paths").transform_values { |path| File.expand_path(path) }
solver_sha256s = metadata.fetch("solver_sha256s")
case_sha256s = metadata.fetch("case_sha256s")
recorded_cases = solver_paths.keys
selected_cases = options[:cases].empty? ? recorded_cases : options[:cases]
unknown_cases = selected_cases - recorded_cases
raise "cases not present in trials.json: #{unknown_cases.join(", ")}" unless unknown_cases.empty?
raise "no cases selected" if selected_cases.empty?
selected_cases.each do |case_name|
  raise "unsafe case name in trials.json: #{case_name.inspect}" unless case_name.match?(/\A[A-Za-z0-9_.-]+\z/)

  solver_path = solver_paths.fetch(case_name)
  raise "recorded solver input does not exist: #{solver_path}" unless File.file?(solver_path)
  solver_sha256 = Digest::SHA256.file(solver_path).hexdigest
  expected_solver_sha256 = solver_sha256s.fetch(case_name)
  unless solver_sha256 == expected_solver_sha256
    raise "recorded #{case_name} solver input hash changed: #{solver_sha256}, expected #{expected_solver_sha256}"
  end

  solver = JSON.parse(File.read(solver_path))
  case_path = File.expand_path(solver.fetch("system_model_file"), File.dirname(solver_path))
  raise "recorded #{case_name} generated case does not exist: #{case_path}" unless File.file?(case_path)
  case_sha256 = Digest::SHA256.file(case_path).hexdigest
  expected_case_sha256 = case_sha256s.fetch(case_name)
  unless case_sha256 == expected_case_sha256
    raise "recorded #{case_name} generated case hash changed: #{case_sha256}, expected #{expected_case_sha256}"
  end
end

recorded_cells = trial_document.fetch("records").map { |record| [record.fetch("label"), record.fetch("case")] }.uniq
selected_cases.each do |case_name|
  [BASELINE_LABEL, VARIANT_LABEL].each do |label|
    raise "trials.json has no measured #{label}/#{case_name} cell" unless recorded_cells.include?([label, case_name])
  end
end

def capture3!(*command)
  stdout, stderr, status = Open3.capture3({"LC_ALL" => "C"}, *command)
  return stdout if status.success?
  raise "command failed (#{command.join(" ")}): #{stderr.strip}"
end

git_head = capture3!("git", "-C", options[:flamegraph_dir], "rev-parse", "HEAD").strip
unless git_head == FLAMEGRAPH_COMMIT
  raise "FlameGraph checkout is at #{git_head}, required #{FLAMEGRAPH_COMMIT}"
end

flamegraph_tools = {
  "stackcollapse" => File.join(options[:flamegraph_dir], "stackcollapse-perf.pl"),
  "render" => File.join(options[:flamegraph_dir], "flamegraph.pl"),
  "diff" => File.join(options[:flamegraph_dir], "difffolded.pl")
}
flamegraph_tools.each_value do |path|
  raise "FlameGraph tool is not executable: #{path}" unless File.executable?(path)
end

perf_version = capture3!(perf, "version").strip

options[:cpu] ||= Integer(metadata.fetch("cpu"))
raise "--cpu must be nonnegative" if options[:cpu].negative?

FileUtils.mkdir_p(output_dir)
raw_dir = File.join(output_dir, "raw")
folded_dir = File.join(output_dir, "folded")
svg_dir = File.join(output_dir, "svg")
FileUtils.mkdir_p(raw_dir)
FileUtils.mkdir_p(folded_dir)
FileUtils.mkdir_p(svg_dir)

def run_to_files!(command, stdout_path, stderr_path, failure_message, chdir: nil)
  spawn_options = {out: stdout_path, err: stderr_path}
  spawn_options[:chdir] = chdir if chdir
  success = system({"LC_ALL" => "C"}, *command, **spawn_options)
  return if success
  raise "#{failure_message}; see #{stderr_path} and #{stdout_path}"
end

smoke_data = File.join(output_dir, "perf-smoke.data")
smoke_stdout = File.join(output_dir, "perf-smoke.stdout.txt")
smoke_stderr = File.join(output_dir, "perf-smoke.stderr.txt")
smoke_script = File.join(output_dir, "perf-smoke.script.txt")
smoke_command = [
  "/usr/bin/taskset", "-c", options[:cpu].to_s,
  perf, "record", "-F", "99", "-e", "cpu-clock:u", "--call-graph", "dwarf,16384", "-o", smoke_data, "--",
  RbConfig.ruby, "-e", "deadline = Process.clock_gettime(Process::CLOCK_MONOTONIC) + 0.1; value = 0; value += 1 while Process.clock_gettime(Process::CLOCK_MONOTONIC) < deadline"
]
begin
  run_to_files!(smoke_command, smoke_stdout, smoke_stderr, "perf prerequisite smoke test failed before any profile capture")
  run_to_files!([perf, "script", "-i", smoke_data], smoke_script, smoke_stderr, "perf could not read its prerequisite smoke capture")
  if !File.file?(smoke_script) || File.zero?(smoke_script)
    raise "perf prerequisite smoke capture contained no samples; no profile captures were started"
  end
ensure
  FileUtils.rm_f([smoke_data, smoke_stdout, smoke_stderr, smoke_script]) unless $!
end

cells = selected_cases.flat_map do |case_name|
  [
    {"label" => BASELINE_LABEL, "case" => case_name, "binary" => binaries.fetch(BASELINE_LABEL).fetch("path")},
    {"label" => VARIANT_LABEL, "case" => case_name, "binary" => binaries.fetch(VARIANT_LABEL).fetch("path")}
  ]
end
taskset = ["/usr/bin/taskset", "-c", options[:cpu].to_s]

cells.each_with_index do |cell, index|
  name = format("warmup-%02d-%s-%s", index + 1, cell.fetch("label"), cell.fetch("case"))
  warn "running #{name}"
  solver_path = solver_paths.fetch(cell.fetch("case"))
  run_to_files!(
    [*taskset, cell.fetch("binary"), solver_path],
    File.join(raw_dir, "#{name}.stdout.txt"),
    File.join(raw_dir, "#{name}.stderr.txt"),
    "profile warm-up failed (#{name})",
    chdir: File.dirname(solver_path)
  )
end

def folded_stats(path)
  stacks = 0
  samples = 0.0
  unknown_stack_samples = 0.0
  frames = +""
  File.foreach(path) do |line|
    match = line.match(/\A(.*)\s+(\d+(?:\.\d+)?)\s*\z/)
    raise "malformed folded stack in #{path}: #{line.inspect}" unless match
    stack = match[1]
    count = match[2].to_f
    stacks += 1
    samples += count
    unresolved = stack.split(";").any? do |frame|
      frame.match?(/\A(?:\[unknown\]|unknown|0x[0-9a-f]+)(?:_\[[^\]]+\])?\z/i)
    end
    unknown_stack_samples += count if unresolved
    frames << stack << "\n"
  end
  {
    "stacks" => stacks,
    "samples" => samples,
    "unknown_stack_samples" => unknown_stack_samples,
    "unknown_stack_percent" => samples.zero? ? 100.0 : unknown_stack_samples / samples * 100.0,
    "frames" => frames
  }
end

capture_records = []
options[:captures].times do |round|
  cells.rotate(round).each do |cell|
    name = format("capture-%02d-%s-%s", round + 1, cell.fetch("label"), cell.fetch("case"))
    warn "running #{name}"
    solver_path = solver_paths.fetch(cell.fetch("case"))
    data_path = File.join(raw_dir, "#{name}.perf.data")
    stdout_path = File.join(raw_dir, "#{name}.stdout.txt")
    stderr_path = File.join(raw_dir, "#{name}.stderr.txt")
    unfolded_path = File.join(raw_dir, "#{name}.unfolded.txt")
    script_stderr_path = File.join(raw_dir, "#{name}.perf-script.stderr.txt")
    folded_path = File.join(folded_dir, "#{name}.folded.txt")
    collapse_stderr_path = File.join(raw_dir, "#{name}.stackcollapse.stderr.txt")
    command = [
      *taskset,
      perf, "record", "-F", "99", "-e", "cpu-clock:u", "--call-graph", "dwarf,16384", "-o", data_path, "--",
      cell.fetch("binary"), solver_path
    ]

    started_at = Process.clock_gettime(Process::CLOCK_MONOTONIC)
    run_to_files!(command, stdout_path, stderr_path, "perf capture failed (#{name})", chdir: File.dirname(solver_path))
    elapsed_seconds = Process.clock_gettime(Process::CLOCK_MONOTONIC) - started_at
    raise "perf produced an empty data file for #{name}" unless File.file?(data_path) && !File.zero?(data_path)

    run_to_files!([perf, "script", "-i", data_path], unfolded_path, script_stderr_path, "perf script failed (#{name})")
    raise "perf script produced no stack samples for #{name}" unless File.file?(unfolded_path) && !File.zero?(unfolded_path)

    run_to_files!([flamegraph_tools.fetch("stackcollapse"), "--all", unfolded_path], folded_path, collapse_stderr_path, "stack collapse failed (#{name})")
    stats = folded_stats(folded_path)
    raise "collapsed stack file contains no samples for #{name}" if stats.fetch("stacks").zero? || stats.fetch("samples").zero?

    capture_records << {
      "capture" => round + 1,
      "label" => cell.fetch("label"),
      "case" => cell.fetch("case"),
      "elapsed_seconds" => elapsed_seconds,
      "perf_data" => data_path,
      "unfolded" => unfolded_path,
      "folded" => folded_path,
      "samples" => stats.fetch("samples"),
      "unknown_stack_percent" => stats.fetch("unknown_stack_percent")
    }
  end
end

def aggregate_folded(inputs, output)
  counts = Hash.new(0.0)
  inputs.each do |path|
    File.foreach(path) do |line|
      match = line.match(/\A(.*)\s+(\d+(?:\.\d+)?)\s*\z/)
      raise "malformed folded stack in #{path}: #{line.inspect}" unless match
      counts[match[1]] += match[2].to_f
    end
  end
  File.open(output, "w") do |file|
    counts.sort.each do |stack, count|
      formatted = count == count.to_i ? count.to_i : count
      file.puts "#{stack} #{formatted}"
    end
  end
end

def validate_svg!(path)
  raise "FlameGraph SVG is empty: #{path}" unless File.file?(path) && !File.zero?(path)
  svg = File.read(path)
  required = ["<svg", "<script", "id=\"frames\""]
  missing = required.reject { |marker| svg.include?(marker) }
  raise "FlameGraph SVG is not interactive (missing #{missing.join(", ")}): #{path}" unless missing.empty?
end

aggregate_records = []
selected_cases.each do |case_name|
  [BASELINE_LABEL, VARIANT_LABEL].each do |label|
    captures = capture_records.select { |record| record.fetch("case") == case_name && record.fetch("label") == label }
    if captures.length != options[:captures]
      raise "expected #{options[:captures]} captures for #{label}/#{case_name}, found #{captures.length}"
    end

    aggregate_path = File.join(folded_dir, "#{label}-#{case_name}.folded.txt")
    aggregate_folded(captures.map { |record| record.fetch("folded") }, aggregate_path)
    stats = folded_stats(aggregate_path)
    raise "aggregate folded stack contains no samples for #{label}/#{case_name}" if stats.fetch("stacks").zero? || stats.fetch("samples").zero?

    frame_patterns = {
      "GridKit" => /(?:GridKit|DynamicSimulation|SystemModel)/,
      "IDA" => /(?:\bIDA[A-Za-z0-9_]*\b|\bida[A-Za-z0-9_]*\b)/,
      "KLU" => /(?:\bKLU[A-Za-z0-9_]*\b|\bklu[A-Za-z0-9_]*\b)/
    }
    missing_frames = frame_patterns.reject { |_name, pattern| stats.fetch("frames").match?(pattern) }.keys
    unless missing_frames.empty?
      raise "#{label}/#{case_name} aggregate lacks resolvable #{missing_frames.join(", ")} frames"
    end
    if stats.fetch("unknown_stack_percent") > options[:max_unknown_percent]
      raise format(
        "%s/%s has %.3f%% of samples containing an unknown frame, above the %.3f%% limit",
        label,
        case_name,
        stats.fetch("unknown_stack_percent"),
        options[:max_unknown_percent]
      )
    end

    svg_path = File.join(svg_dir, "#{label}-#{case_name}.svg")
    render_stderr_path = File.join(raw_dir, "#{label}-#{case_name}.flamegraph.stderr.txt")
    title = "#{case_name}: #{label.tr("_", " ")}"
    run_to_files!([flamegraph_tools.fetch("render"), "--hash", "--title", title, aggregate_path], svg_path, render_stderr_path, "FlameGraph render failed (#{label}/#{case_name})")
    validate_svg!(svg_path)

    aggregate_records << {
      "label" => label,
      "case" => case_name,
      "folded" => aggregate_path,
      "svg" => svg_path,
      "stacks" => stats.fetch("stacks"),
      "samples" => stats.fetch("samples"),
      "unknown_stack_percent" => stats.fetch("unknown_stack_percent")
    }
  end
end

differential_records = selected_cases.map do |case_name|
  baseline = aggregate_records.find { |record| record.fetch("case") == case_name && record.fetch("label") == BASELINE_LABEL }
  variant = aggregate_records.find { |record| record.fetch("case") == case_name && record.fetch("label") == VARIANT_LABEL }
  diff_path = File.join(folded_dir, "#{VARIANT_LABEL}-vs-#{BASELINE_LABEL}-#{case_name}.folded.txt")
  diff_stderr_path = File.join(raw_dir, "#{VARIANT_LABEL}-vs-#{BASELINE_LABEL}-#{case_name}.difffolded.stderr.txt")
  run_to_files!(
    [flamegraph_tools.fetch("diff"), "-n", baseline.fetch("folded"), variant.fetch("folded")],
    diff_path,
    diff_stderr_path,
    "normalized folded diff failed (#{case_name})"
  )
  raise "normalized folded diff is empty for #{case_name}" unless File.file?(diff_path) && !File.zero?(diff_path)

  svg_path = File.join(svg_dir, "#{VARIANT_LABEL}-vs-#{BASELINE_LABEL}-#{case_name}.svg")
  render_stderr_path = File.join(raw_dir, "#{VARIANT_LABEL}-vs-#{BASELINE_LABEL}-#{case_name}.flamegraph.stderr.txt")
  title = "#{case_name}: internal current variables vs direct injection"
  run_to_files!([flamegraph_tools.fetch("render"), "--title", title, diff_path], svg_path, render_stderr_path, "differential FlameGraph render failed (#{case_name})")
  validate_svg!(svg_path)
  {"case" => case_name, "folded" => diff_path, "svg" => svg_path, "widths" => VARIANT_LABEL, "normalized" => true}
end

profile_metadata = {
  "created_at" => Time.now.iso8601,
  "trials" => options[:trials],
  "trials_sha256" => Digest::SHA256.file(options[:trials]).hexdigest,
  "cpu" => options[:cpu],
  "captures" => options[:captures],
  "cases" => selected_cases,
  "perf" => perf,
  "perf_version" => perf_version,
  "perf_event" => "cpu-clock:u",
  "perf_frequency_hz" => 99,
  "perf_call_graph" => "dwarf,16384",
  "flamegraph_dir" => options[:flamegraph_dir],
  "flamegraph_commit" => git_head,
  "max_unknown_percent" => options[:max_unknown_percent],
  "binaries" => binaries,
  "solver_paths" => selected_cases.to_h { |case_name| [case_name, solver_paths.fetch(case_name)] },
  "solver_sha256s" => selected_cases.to_h { |case_name| [case_name, solver_sha256s.fetch(case_name)] },
  "capture_records" => capture_records,
  "aggregate_records" => aggregate_records,
  "differential_records" => differential_records
}
File.write(File.join(output_dir, "profiles.json"), JSON.pretty_generate(profile_metadata))

puts "wrote #{capture_records.length} captures and #{differential_records.length} differential FlameGraphs to #{output_dir}"
