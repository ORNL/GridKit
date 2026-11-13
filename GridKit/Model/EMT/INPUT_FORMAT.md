# EMT Case Format

EMT cases use the root layout of the
[grid dynamics case format](../PhasorDynamics/INPUT_FORMAT.md). Values are SI,
complex values are `[re, im]`, and matrices are arrays of rows.

## System Parameters

  Name      | Value
  ----------|-------------------------------------------------------
  `va_base` | System power base [VA]

## Buses

  Name     | Value
  ---------|-------------------------------------------------------
  `number` | Bus number
  `class`  | `"Bus"`
  `params` | [Bus](Bus/README.md) parameters
  `mon`    | Monitored variables (optional)

## Devices

  Name     | Value
  ---------|-------------------------------------------------------
  `class`  | Model class
  `id`     | Device id
  `ports`  | Bus number at each port
  `params` | Model parameters by README `JSON` name
  `mon`    | Monitored variables (optional)

## State File

Initial values by bus number and device id, with submodels nested by symbol.
Omitted values are zero.

```json
{
  "buses": {"1": {"v": [1.0]}},
  "devices": {"LD1": {"i": [-0.5], "z": {"w": [0.25]}}}
}
```
