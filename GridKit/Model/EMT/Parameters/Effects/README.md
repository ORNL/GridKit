# Effect Models

`Effects` contains the physical effect models that assemble full-conductor
per-unit-length overhead-line data:

```math
\mathbf{R}'(\omega),\quad \mathbf{L}'(\omega),\quad
\mathbf{G}'(\omega),\quad \mathbf{C}'(\omega)
```

These models are frequency-domain parameter models, not EMT network components.
They consume geometry, conductor, and earth data and produce algebraic signals
for the response models. They are evaluated at the current system parameter
$\omega$; sweep endpoints, sampling, and linear or logarithmic grids belong to
the application or study driver.

Model | Description
----- | -----------
[`SeriesImpedance`](SeriesImpedance/README.md) | Builds full-conductor per-unit-length series resistance and inductance
[`GeometricInductance`](GeometricInductance/README.md) | Computes external geometric inductance entries
[`SkinEffect`](SkinEffect/README.md) | Computes conductor-internal skin-effect resistance and inductance entries
[`ShuntAdmittance`](ShuntAdmittance/README.md) | Builds full-conductor per-unit-length shunt conductance and capacitance
[`ShuntPotential`](ShuntPotential/README.md) | Computes potential-derived zero conductance and capacitance
[`InsulatorLeakage`](InsulatorLeakage/README.md) | Computes direct shunt leakage conductance
[`Carson`](Carson/README.md) | Computes Carson earth-return resistance and inductance entries

## Assembly

Output | Built From
------ | ----------
$\mathbf{R}'$ | Internal skin-effect resistance and Carson earth-return resistance
$\mathbf{L}'$ | Internal skin-effect inductance, geometric inductance, and Carson earth-return inductance
$\mathbf{G}'$ | Shunt potential terms and optional direct leakage terms
$\mathbf{C}'$ | Shunt potential terms

The aggregate [`Overhead`](../README.md#overhead-aggregate) model wires these
effect outputs into the response models.
