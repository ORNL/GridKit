# ContingencyAnalysis

Source: `application/PhasorDynamics/ContingencyAnalysis.cpp`

Input format: [Phasor Dynamics](README.md)

Set `contingency_stats_file` in the solver input to write a JSON document with
one record per supplied bus fault. Each record identifies the bus and fault
index, success/failure status, diagnostic text, total IDA counters, and counters
for each event-delimited segment. The output path is relative to the working
directory, like monitor output paths; its parent directory must exist.

Counters are captured before every event restart and after the final segment.
They include IDA's consistent-initialization work at the fault events. The
initial solve keeps the application's existing `findConsistent=false` policy.
Jacobian evaluations come from `IDAGetNumJacEvals`; linear-solver setups are
reported separately. Failed records have null totals and are not valid samples
for a Jacobian-evaluations-per-accepted-step plot. Results are written after
workers finish, including when some faults fail; the application still returns
a failure exit status in that case.
