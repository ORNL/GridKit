# DCLink Model

`DCLink` represents an ideal DC-link capacitor. It owns one differential voltage
and exchanges current and voltage signals with the source and converter.

## Block Diagram

None.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$C$ | [F] | `C` | DC-link capacitance | Required, finite and positive

### Parameter Validation

```math
C > 0
```

The capacitance must be finite.

### Derived Parameters

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$i_{\mathrm{src}}$ | `isrc` | Input | [A] | Current supplied to the DC link | Required
$i_{\mathrm{dc}}$ | `idc` | Input | [A] | Current drawn from the DC link | Required
$v_{\mathrm{dc}}$ | `vdc` | Output | [V] | Capacitor voltage | Owned differential variable

Connect `idc` to the [Converter](../../../Operators/Converter/README.md) DC-current
output and connect `vdc` to its DC-voltage input. Positive source current charges
the capacitor; positive converter current discharges it. Negative current permits
regeneration. Voltage is signed; the model imposes no clamp or protection logic.

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$v_{\mathrm{dc}}$ | [V] | DC-link voltage | One owned state

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$i_{\mathrm{src}}$ | [A] | Source current | Positive into the capacitor
$i_{\mathrm{dc}}$ | [A] | Converter current | Positive out of the capacitor

## Model Equations

### Internal Equations

#### Differential

```math
0 = -C\dfrac{\mathrm{d}v_{\mathrm{dc}}}{\mathrm{d}t}
    + i_{\mathrm{src}} - i_{\mathrm{dc}}
```

#### Algebraic

None.

### External Equations

None.

## Initialization

The state-file key `vdc` sets the initial voltage, defaulting to zero. Only finite
values are accepted. The consistent-initial-condition solve preserves this voltage
and obtains $\mathrm{d}v_{\mathrm{dc}}/\mathrm{d}t=(i_{\mathrm{src}}-i_{\mathrm{dc}})/C$ from the
connected network. Initialization does not depend on the current inputs, so the
converter feedback connection adds no initialization ordering cycle.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`vdc` | [V] | DC-link voltage |
`isrc` | [A] | Source current | Positive into the capacitor
`idc` | [A] | Converter current | Positive out of the capacitor
`energy` | [J] | Stored energy | $E=\tfrac{1}{2}Cv_{\mathrm{dc}}^2$

The energy balance is

```math
\dfrac{\mathrm{d}E}{\mathrm{d}t}=v_{\mathrm{dc}}i_{\mathrm{src}}-v_{\mathrm{dc}}i_{\mathrm{dc}}.
```

See the [charge/discharge case](../../../../../../cases/EMT/DCLink/README.md) for
an exact-solution validation and plotting command.
