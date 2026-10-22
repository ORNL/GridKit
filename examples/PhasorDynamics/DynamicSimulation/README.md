# DynamicSimulation

These examples demonstrate running phasor-dynamics simulations using JSON
inputs.

| Example | Description |
| --- | --- |
| [ThreeBusBasic](Toy/ThreeBusBasic/README.md) | A three-bus line opened and reclosed. |
| [ThreeBusConstantSource](Toy/ThreeBusConstantSource/README.md) | A three-bus constant signal source example. |
| [ThreeBusRelay](Toy/ThreeBusRelay/README.md) | A mid-line fault cleared by definite-time and inverse-time relays. |
| [ACTIVSg10k](ACTIVSg10k/README.md) | A short simulation without disturbances using the reusable ACTIVSg10k case. |
| [WECC240Relay](WECC240Relay/README.md) | A Southern California 500 kV line fault isolated by overcurrent relays. |

The relay examples are plotted with [plot.py](plot.py) after a run:

```bash
python plot.py path/to/build/.../WECC240Relay.solver.json --dir WECC240Relay/figures
```
