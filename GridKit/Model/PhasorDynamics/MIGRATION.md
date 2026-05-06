# PhasorDynamics Migration Notes

Temporary notes for converting legacy files from [INPUT_FORMAT.md](INPUT_FORMAT.md)
to the canonical input format.


| Legacy | Canonical |
| --- | --- |
| Bus `v_base`, `freq_base`, `va_base` | Bus `params.v_base`, `params.freq_base`, `params.va_base` |
| Bus `init.Vr`, `init.Vi` | Bus `init.v.re`, `init.v.im` |
| Machine `params.p0`, `params.q0` | Machine `init.p`, `init.q` |
| Branch objects in `devices` | Branch objects in root `branches` |
| Branch class `Branch` | Branch class `Line` |
| Bus fault `params.state0` | Bus fault `init.enabled` |


Each component README should identify which input fields are required, optional,
or defaulted for that component. This includes `params`, `ports`, `init`, and
monitor variables when applicable. Keep detailed component contracts in the
component READMEs, not in this migration note or the type directory.
