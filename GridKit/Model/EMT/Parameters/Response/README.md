# Response Models

`Response` contains the propagation and characteristic models derived from
assembled full-conductor overhead-line parameters.

These models consume the per-unit-length matrices from
[`SeriesImpedance`](../Effects/SeriesImpedance/README.md) and
[`ShuntAdmittance`](../Effects/ShuntAdmittance/README.md). They are evaluated
at the current system parameter $\omega$, and every quantity they own is
analytic in $\omega$, so all of them belong in the sweep DAE.

Model | Description
----- | -----------
[`Gamma`](Gamma/README.md) | Computes the propagation constant matrix
[`Yc`](Yc/README.md) | Computes characteristic admittance
[`Zc`](Zc/README.md) | Computes characteristic impedance

Modal quantities (transformations, modal constants, delays, and propagation
factors) are not response models: they are observations of `Gamma` produced
by [`ModalDecomposition`](../Modal/README.md) at monitor emission, outside
the DAE.

## Inputs

Quantity | Source | Used By
-------- | ------ | -------
$\mathbf{R}'$, $\mathbf{L}'$ | [`SeriesImpedance`](../Effects/SeriesImpedance/README.md) | `Gamma`, `Yc`, `Zc`
$\mathbf{G}'$, $\mathbf{C}'$ | [`ShuntAdmittance`](../Effects/ShuntAdmittance/README.md) | `Gamma`, `Yc`, `Zc`

## Outputs

Quantity | Produced By
-------- | -----------
$\boldsymbol{\Gamma}$ | [`Gamma`](Gamma/README.md)
$\mathbf{Y}_c$ | [`Yc`](Yc/README.md)
$\mathbf{Z}_c$ | [`Zc`](Zc/README.md)
