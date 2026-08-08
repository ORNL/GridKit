require "open3"

directory = File.expand_path(__dir__)
repo_root = File.expand_path("../../..", directory)
application = File.join(repo_root, "build/application/PhasorDynamics/DynamicSimulation")
cases = %w[ACTIVSg200 ACTIVSg500 WECC240 ACTIVSg10k]
runs = Hash.new { |hash, key| hash[key] = [] }

monitor_sinks = %w[
  ACTIVSg200.omega.csv
  ACTIVSg500.omega.csv
  WECC240.omega.csv
].map { |name| File.join(directory, name) }
created_monitor_sinks = []
monitor_sinks.each do |path|
  if File.exist?(path) || File.symlink?(path)
    raise "monitor sink already exists: #{path}" unless File.symlink?(path) && File.readlink(path) == "/dev/null"
  else
    File.symlink("/dev/null", path)
    created_monitor_sinks << path
  end
end
at_exit do
  created_monitor_sinks.each do |path|
    File.unlink(path) if File.symlink?(path) && File.readlink(path) == "/dev/null"
  end
end

results_file = File.open(File.join(directory, "results.tsv"), "w")
results_file.puts "type\ttrial\tcase\tcomplete_seconds\tsteps\tresidual_seconds\tjacobian_seconds\tlinear_solve_seconds\tmonitor_calls"

def capture_metric(output, pattern, name)
  match = output.match(pattern)
  raise "missing #{name}" unless match
  match[1]
end

3.times do |round|
  cases.rotate(round).each do |case_name|
    solver = File.join(directory, "#{case_name}.solver.json")
    output, status = Open3.capture2e(application, solver, chdir: directory)
    unless status.success?
      warn "#{case_name} failed:\n#{output}"
      exit(status.exitstatus || 1)
    end

    result = {
      complete: capture_metric(output, /Complete in ([0-9.eE+-]+) seconds/, "Complete in").to_f,
      buses: capture_metric(output, /^buses=(\d+)$/, "buses").to_i,
      states: capture_metric(output, /^states=(\d+)$/, "states").to_i,
      nnz: capture_metric(output, /^jacobian_nnz=(\d+)$/, "jacobian_nnz").to_i,
      steps: capture_metric(output, /^\s*Steps\s*:\s*(\d+)/, "Steps").to_i,
      residual_seconds: capture_metric(output, /^residual_seconds=([0-9.eE+-]+)$/, "residual_seconds").to_f,
      jacobian_seconds: capture_metric(output, /^jacobian_seconds=([0-9.eE+-]+)$/, "jacobian_seconds").to_f,
      linear_solve_seconds: capture_metric(output, /^linear_solve_seconds=([0-9.eE+-]+)$/, "linear_solve_seconds").to_f,
      monitor_calls: capture_metric(output, /^monitor_print_calls=(\d+)$/, "monitor_print_calls").to_i
    }
    runs[case_name] << result

    line = ["trial", round + 1, case_name, result[:complete], result[:steps],
            result[:residual_seconds], result[:jacobian_seconds],
            result[:linear_solve_seconds], result[:monitor_calls]].join("\t")
    puts line
    results_file.puts line
    results_file.flush
    STDOUT.flush
  end
end

summary_header = "summary\tcase\tbuses\tstates\tnnz\tmedian_complete\tmedian_steps\tmedian_residual_seconds\tmedian_jacobian_seconds\tmedian_linear_solve_seconds\tmonitor_calls"
puts summary_header
results_file.puts summary_header
cases.each do |case_name|
  values = runs.fetch(case_name)
  median = lambda do |key|
    values.map { |result| result.fetch(key) }.sort[1]
  end
  first = values.first
  line = ["summary", case_name, first[:buses], first[:states], first[:nnz],
          median.call(:complete), median.call(:steps), median.call(:residual_seconds),
          median.call(:jacobian_seconds), median.call(:linear_solve_seconds),
          first[:monitor_calls]].join("\t")
  puts line
  results_file.puts line
end
results_file.close
