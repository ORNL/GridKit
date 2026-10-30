# SixBusGas Supply Pressure Drop

This study uses the reusable
[SixBusGas case](../../../../../cases/PhasorDynamics/Toy/README.md): a six-bus
system whose two GENROU machines take their mechanical power from a gas network
through `GasInputPm`, with no governors. Gas bus 101 is a `GasBusInfinite`
supply, and pipes lead through junction 102 to the machine gas buses 103 (G2)
and 104 (G3).

The case starts in steady state. [`steady_state.py`](steady_state.py) solves
the power flow for each machine's mechanical power, converts it to a fuel
withdrawal with the specific fuel consumption, solves the gas network for the
bus pressures, and writes the case file.

At 5 s the supply pressure drops from 5.0 to 4.5 MPa. The pressure wave reaches
the junction about 13 s later and then the machine buses, and each machine's
mechanical power falls with its gas-bus pressure.

![SixBusGas response](SixBusGas.png)

To regenerate the case, run `python3 steady_state.py`. To regenerate the figure,
run the study in its build directory and pass the output to
[`plot.py`](plot.py):

```bash
DynamicSimulation SixBusGas.solver.json
python3 plot.py SixBusGas.csv
```
