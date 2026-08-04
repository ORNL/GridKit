# OPF case format

Version 0.1 uses UTF-8 JSON and the `.opf.json` extension. It stores network
topology and parameters; operating values and statuses belong in the companion
`.state.json` file.

The root object contains exactly four required fields:

| Field | Type | Contents |
|---|---|---|
| `header` | object | Format metadata |
| `params` | object | System bases |
| `buses` | array | Bus records |
| `devices` | array | Device records |

Unknown fields are invalid. Required fields cannot be `null`. All numerical
values must be finite. Electrical quantities are per unit on `va_base` unless
stated otherwise.

## Header and system parameters

| Object | Field | Required | Type |
|---|---|---:|---|
| `header` | `format_version` | yes | non-negative integer; must be `0` |
| `header` | `format_revision` | yes | non-negative integer; must be `1` |
| `header` | `case_name` | yes | string |
| `header` | `case_date_time` | no | string or `null` |
| `header` | `case_description` | no | string or `null` |
| `header` | `case_comments` | no | string or `null` |
| `params` | `freq_base` | yes | hertz |
| `params` | `va_base` | yes | volt-amperes |

## Buses

Each bus record contains exactly `class`, `id`, and `params`. `class` must be
`"Bus"`; `id` is a string.

| Parameter | Required | Description |
|---|---:|---|
| `number` | yes | Non-negative integer used by device bus maps |
| `kv` | yes | Voltage base in kilovolts |
| `vmin`, `vmax` | no | Voltage-magnitude bounds; missing or `null` is unbounded |

The OPF model uses the lowest-numbered bus as its voltage-angle gauge.

## Devices

Each device record contains exactly `class`, `buses`, `id`, and `params`.
`id` is a string. Bus maps and parameter objects accept only the fields listed
below.

| Class | Bus map | Required parameters | Optional parameters |
|---|---|---|---|
| `Branch` | `from`, `to` | `R`, `X`, `G`, `B` | `smax` |
| `Generator` | `bus` | `mva` | `pmin`, `pmax`, `qmin`, `qmax`, `c0`, `c1`, `c2` |
| `Load` | `bus` | none | `pmin`, `pmax`, `qmin`, `qmax` |
| `Shunt` | `bus` | `G`, `B` | none |

`R` and `X` are series impedance. Branch `G` and `B` are total shunt
admittance. `smax` is the apparent-power limit at each terminal. `mva` is the
generator machine base in megavolt-amperes. Generator cost is

\[
c_0+c_1p+c_2p^2,
\]

with `p` in per unit. Missing cost coefficients default to zero; unlike limits,
explicitly `null` cost coefficients are invalid. All optional bounds and
`smax` may be omitted or set to `null`.

Identifiers, topology, limit ordering, and operating-state consistency are
validated when the OPF system is allocated rather than by this syntax parser.
