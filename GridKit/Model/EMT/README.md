# Electromagnetic Transients (EMT)

This directory documents electromagnetic transient (EMT) model specifications
and reusable operators in instantaneous phase coordinates.

## Conventions

- $N$ denotes a phase or output dimension, $K$ a conductor or input dimension,
  $M$ a mode count, and $Q$ a pole count.
- Equations use SI units unless a model states otherwise.
- Time derivatives are written explicitly as $\mathrm{d}(\cdot)/\mathrm{d}t$.
- Current injection terms are written as positive into buses.
- Initialization assignments use $\leftarrow$ and separate input, internal,
  and output initialization. Text on the right-hand side is reserved for
  externally supplied inputs.

## Directories

The EMT documentation is organized into:

- `Bus` (See [Bus](Bus/README.md))
- `Component` (See [Component](Component/README.md))
- `Operators` (See [Operators](Operators/README.md))
- `Parameters` (See [Parameters](Parameters/README.md))

Component models are documented under `Component` because they are EMT models
connected to buses. Reusable transfer and coordinate tools are documented under
`Operators`.
Frequency-domain line-parameter models are documented under `Parameters`
because they generate sampled per-unit-length data rather than EMT bus
residuals.
