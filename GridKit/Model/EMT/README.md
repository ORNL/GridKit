# Electromagnetic Transients (EMT)

This directory documents electromagnetic transient (EMT) model specifications
and reusable operators in instantaneous phase coordinates.

> [!NOTE]
> The formulation supports $N$ phases; initial development targets three
> phases.

## Conventions

- $N$ denotes a phase or output dimension, $K$ a conductor or input dimension,
  $M$ a mode count, and $Q$ a pole count.
- Corresponding calligraphic symbols denote index sets, such as
  $\mathcal{M}$ for modal indices and $\mathcal{Q}$ for pole indices.
- Equations use SI units unless a model states otherwise.
- Time derivatives are written explicitly as $\mathrm{d}(\cdot)/\mathrm{d}t$.
- Current injection terms are written as positive into buses.
- $\mathrm{j}$ denotes the imaginary unit and is not used as an index.
- The Model Variables section lists DAE variables owned directly by the model.
  Each scalar internal variable corresponds to one local residual row;
  external variables are owned elsewhere. Submodel-owned variables are not
  repeated by the parent.
- A terminal connects a model to a bus: the model reads the bus voltage and
  adds its current to the current balance of the bus.
- An input reads a variable. An output is a value that a model computes; it
  is added to the residual of the variable it is connected to.
- An input that names an output reads a signal that follows it, a variable
  with residual $0 = -s + y$ to which the model adds its output $y$.
- An output that no input reads and no case connects is evaluated only for a
  monitor.
- An input that the case leaves out is held. It reads a held signal, a
  variable with residual $0 = -\mathrm{d}s/\mathrm{d}t$, whose value the
  state gives and an event changes.
- A monitor is an internal variable or an output.
- Equation headings follow the assembled DAE row classification; connected
  models may add derivative terms to a residual that contains none locally.

## Assembly

Bus voltage is algebraic unless a connected model contributes a voltage
derivative to the bus current-balance.

A held input and an input that names an output each read a
[Signal](Operators/Reference/Signal/README.md). The system adds one for each
kind.

## Initialization

EMT models are initialized using real-valued, instantaneous phase-coordinate
quantities in $\mathbb{R}^N$. RMS or phasor calculations, when used by an
upstream initialization workflow, are outside the EMT model specification.

Initial differential states, algebraic variables, derivatives, and discrete
inputs must be consistent with the assembled EMT equations at $t_0$. Models
that require additional history document that requirement locally.

## Contents

- [Input Format](INPUT_FORMAT.md)
- [Bus](Bus/README.md)
- [Components](Component/README.md)
- [Operators](Operators/README.md)

## References

- B. Gustavsen and A. Semlyen, [“Rational approximation of frequency domain
  responses by vector fitting,”](https://doi.org/10.1109/61.772353) *IEEE
  Transactions on Power Delivery*, 1999.
- A. Morched, B. Gustavsen, and M. Tartibi, [“A universal model for accurate
  calculation of electromagnetic transients on overhead lines and underground
  cables,”](https://doi.org/10.1109/61.772350) *IEEE Transactions on Power
  Delivery*, 1999.
- B. Gustavsen and A. Semlyen, [“Enforcing passivity for admittance matrices
  approximated by rational functions,”](https://doi.org/10.1109/59.910786)
  *IEEE Transactions on Power Systems*, 2001.
