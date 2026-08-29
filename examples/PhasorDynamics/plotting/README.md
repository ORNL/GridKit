# PhasorDynamics plotting

These scripts plot standard `DynamicSimulation` outputs. They accept any
solver JSON that specifies `output_file` or `step_trace_file`; no case-specific
preparation or simulation wrapper is required.

From the repository root, for example:

```bash
python3 examples/PhasorDynamics/plotting/plot_bus_frequency.py \
  build/examples/PhasorDynamics/Huge/ACTIVSg25k/ACTIVSg25k.plot.solver.json \
  --output examples/PhasorDynamics/Huge/ACTIVSg25k/ACTIVSg25k_bus_frequency.png

python3 examples/PhasorDynamics/plotting/plot_adaptive_step.py \
  build/examples/PhasorDynamics/Large/ACTIVSg10k/ACTIVSg10k.plot.solver.json \
  build/examples/PhasorDynamics/Huge/ACTIVSg25k/ACTIVSg25k.plot.solver.json \
  build/examples/PhasorDynamics/Huge/ACTIVSg70k/ACTIVSg70k.plot.solver.json \
  --output examples/PhasorDynamics/ACTIVSg_adaptive_step.png
```

The frequency plot calculates

\[
f(t) = f_0 + \frac{1}{2\pi}\frac{d\theta(t)}{dt}
\]

from monitored bus voltage angle. Samples adjacent to switching events remain
visible but do not set the automatic frequency-axis limits or color scale.

The adaptive-step script overlays every supplied study in one comparison
figure. It does not generate separate case plots.
