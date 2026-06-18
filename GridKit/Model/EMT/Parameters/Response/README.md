# Response Models

`Response` contains the propagation and characteristic models derived from
assembled full-conductor overhead-line parameters.

These models consume the per-unit-length matrices from
[`SeriesImpedance`](../Effects/SeriesImpedance/README.md) and
[`ShuntAdmittance`](../Effects/ShuntAdmittance/README.md). Finite-line response
models also consume the path length from [`Path`](../Geometry/Path/README.md).
They are evaluated at the current system parameter $\omega$.

Model | Description
----- | -----------
[`Gamma`](Gamma/README.md) | Computes propagation constant matrices and modal propagation constants
[`Tau`](Tau/README.md) | Computes finite-length modal phase delay
[`H`](H/README.md) | Computes modal finite-length propagation function
[`Yc`](Yc/README.md) | Computes characteristic admittance
[`Zc`](Zc/README.md) | Computes characteristic impedance

## Inputs

Quantity | Source | Used By
-------- | ------ | -------
$\mathbf{R}'$, $\mathbf{L}'$ | [`SeriesImpedance`](../Effects/SeriesImpedance/README.md) | `Gamma`, `Yc`, `Zc`
$\mathbf{G}'$, $\mathbf{C}'$ | [`ShuntAdmittance`](../Effects/ShuntAdmittance/README.md) | `Gamma`, `Yc`, `Zc`
$\ell$ | [`Path`](../Geometry/Path/README.md) | `Tau`, `H`

## Outputs

Quantity | Produced By
-------- | -----------
$\boldsymbol{\Gamma}$, $\mathbf{T}_v$, $\mathbf{T}_i$, $\mathbf{a}$, $\mathbf{b}$ | [`Gamma`](Gamma/README.md)
$\boldsymbol{\tau}$ | [`Tau`](Tau/README.md)
$\mathbf{H}$ | [`H`](H/README.md)
$\mathbf{Y}_c$ | [`Yc`](Yc/README.md)
$\mathbf{Z}_c$ | [`Zc`](Zc/README.md)
