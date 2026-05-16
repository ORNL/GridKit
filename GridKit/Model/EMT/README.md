# Electromagnetic Transients (EMT)

This directory contains design documentation for electromagnetic transient
(EMT) component models in instantaneous abc coordinates.

## Conventions

These conventions reflect the current EMT model draft and may change as the
EMT design and implementation develop.

- Phase order is $a$, $b$, $c$.
- Equations use SI units unless a model states otherwise.
- Current injection terms are written as positive into buses.


## Model Categories

The current EMT documentation is organized into two categories:
- `Bus`
- `Component`

Branch models such as `BranchLumpedConstant` are documented under
`Component/Branch` because they are EMT components connected to buses.

## Open Design Notes

Distributed parameter lines are placeholders until internal signal delay support
is designed.
- Initial electrical wiring will use Delta configuration only
