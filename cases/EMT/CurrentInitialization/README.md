# Initial currents in an RL circuit

The normal `EMTDynamicSimulation` application runs a three-phase RL branch
and an RL load supplied by a 10 V peak, 20 rad/s ideal voltage source. The
branch has R = 2 ohm and L = 0.05 H; the load has R = 4 ohm and L = 0.1 H.
A 5 ohm resistive load checks algebraic-current initialization. All devices
are inside `network`, exercising qualified paths in the state file.

`Steady.state.json` supplies sinusoidal steady-state currents.
`Perturbed.state.json` adds known initial offsets and supplies deliberately
inconsistent algebraic current guesses for the resistor. `Default.state.json`
omits currents, retaining the de-energized inductive start. The solver files
use `ya_ydp`, which preserves differential currents while finding consistent
algebraic values and derivatives.

Run a case from an output directory:

```sh
EMTDynamicSimulation /path/to/CurrentInitialization/Steady.solver.json
```

Run all three cases and check the independent analytic solution:

```sh
python3 validate.py --exe /path/to/EMTDynamicSimulation
```

For each uncoupled phase, the validation uses

```math
I = \frac{V}{R+j\omega L},\qquad
i(t) = \operatorname{Re}(I e^{j\omega t})
       + [i(0)-\operatorname{Re}(I)]e^{-Rt/L}.
```

For load injection, use `-V` in the expression for `I`. The check covers
all phase currents and their differential-state derivatives, initial-current
preservation, algebraic-current correction, finite output, and sample times.
Bounds are 2e-7 A and 2e-5 A/s for the 1e-10 solver tolerances. Runs use
temporary output directories; only case inputs and the oracle are retained.
