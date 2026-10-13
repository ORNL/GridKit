# Delayed signal

This example delays a cosine source by 1 ms and 150 ms. Both outputs share one
adaptive source history. The expected output is 1 for $t<\tau$ and
$\cos(t-\tau)$ afterward, using constant equilibrium prehistory.

From the configured example directory, run:

```sh
/path/to/build/application/PhasorDynamics/DynamicSimulation DelayedSignal.solver.json
```

`mon.csv` contains both delayed outputs. Set `signal_history.compress` to `false`
in the solver file to retain every accepted history segment within the delay
window. The 1 ms delay does not impose a 1 ms maximum step.

See the [model contract](../../../../../GridKit/Model/PhasorDynamics/SignalDelay/README.md).

See the [signal delay gallery](gallery/index.html) for smooth inputs, different
delay durations, adaptive steps, and history compression.
