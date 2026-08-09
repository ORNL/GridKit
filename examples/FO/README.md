# WECC240 forced-oscillation response gallery

This gallery contains 100 independent ten-second WECC240 studies: 63
synchronous-machine/control injections and one injection at each of the 37
renewable units. Together they exercise all 17 connected signal-input types in
the validated case.

Every forcing starts at 1 s, remains active through the end of the simulation,
and is monitored at 0.05 s. The set covers sine, square, triangle, and sawtooth
carriers together with frequency chirps and decaying envelopes.

Each numbered directory contains exactly four system-response figures:

1. bus-voltage magnitude;
2. bus-voltage angle;
3. synchronous-machine and renewable-converter active power; and
4. synchronous-machine and renewable-converter reactive power.

The narrow upper panel in every figure shows the applied forcing. Responses are
changes from the final sample before activation; all system traces are shown and
the associated bus or resource is highlighted. Power is on system base, voltage
angle is shown in degrees, and voltage magnitude is per unit.

[`manifest.csv`](manifest.csv) records the target, waveform, parameters, and
response extrema for each study. Numerical monitor output is temporary and is
not retained.

Regenerate the gallery from the repository root after building
`DynamicSimulation`:

```sh
python3 examples/FO/generate_gallery.py
```
