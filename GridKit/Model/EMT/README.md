# Electromagnetic Transients (EMT)

This directory contains design documentation for electromagnetic transient
(EMT) models and reusable operators in instantaneous phase coordinates.

## Conventions

These conventions reflect the current EMT model draft and may change as the
EMT design and implementation develop.

- EMT model equations are written for $N$ phases unless a model states a
  narrower implementation scope.
- Equations use SI units unless a model states otherwise.
- Time derivatives are written explicitly as $\mathrm{d}(\cdot)/\mathrm{d}t$.
- Current injection terms are written as positive into buses.

## Directories

The current EMT documentation is organized into:

- `Bus` (See [Bus](Bus/README.md))
- `Component` (See [Component](Component/README.md))
- `Operators` (See [Operators](Operators/README.md))
