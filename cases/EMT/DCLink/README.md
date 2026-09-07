# DC-link charge and discharge

This case connects a 20 mF DC-link capacitor to a Converter with fixed bridge
state `(1, 0, 0)`. Each AC phase has a 1 ohm source resistance and a 9 ohm load.
The source current steps from 80 A to 20 A at 0.5 s. The capacitor starts at 600 V.
This is a fixed-state bridge test, with no PWM or voltage controller.

The bridge projection is `(2/3, -1/3, -1/3)`, so its input current is
`idc = vdc / 15 ohm`. The exact capacitor time constant is 0.3 s:

```math
v_{\mathrm{dc}}(t)=
\begin{cases}
1200-600e^{-t/0.3}, & 0\le t\le0.5,\\
300+[v_{\mathrm{dc}}(0.5)-300]e^{-(t-0.5)/0.3}, & t>0.5.
\end{cases}
```

The normal EMT application uses IDA with the sparse KLU solver. From the
repository root:

```sh
cmake --build build --target EMTDynamicSimulation -j 10
python3 cases/EMT/DCLink/validate.py \
  --exe build/application/EMT/EMTDynamicSimulation \
  --results build/validation/DCLink --plot
```

Plotting requires matplotlib. The validator itself uses only the Python standard
library and is registered as `EMTDcLinkTransient` when Enzyme, sparse SUNDIALS,
and Python are enabled. Without `--results`, it uses a temporary directory.

Validation checks all 3,002 samples, the exact voltage, inferred differential
classification, current signs, converter power balance, voltage continuity and
derivative jump at the event, and capacitor energy against integrated net power.
Both sides of the event are retained for quadrature. The energy-integral tolerance
accounts for trapezoidal quadrature at the 0.5 ms monitor interval.

Outputs are `mon.csv`, `state.csv` and its layout, `simulation.log`, `metrics.json`,
and, with `--plot`, `dc_link.png` and `dc_link.svg`.
