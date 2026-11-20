# ParaEMT WECC240 Validation

This study uses the [case](WECC240.case.json) and [state](WECC240.state.json)
in this folder, built from the public ParaEMT WECC240 data. The Palo Verde
nuclear unit trips at 1 s, and each phase of its breaker opens at its next
current zero.

## Development

The inverter plants use `Repca`, `Reecb` and `Regca`, each with its own `Pmu`.
