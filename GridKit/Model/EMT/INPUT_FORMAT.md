# EMT Case Format

EMT cases use the root layout of the
[grid dynamics case format](../PhasorDynamics/INPUT_FORMAT.md). Values are SI,
complex values are `[re, im]`, and matrices are arrays of rows.

## System Parameters

  Name        | Value
  ------------|-------------------------------------------------------
  `freq_base` | System frequency base [Hz]
  `va_base`   | System power base [VA]

## Buses

  Name     | Value
  ---------|-------------------------------------------------------
  `number` | Bus number
  `class`  | `"Bus"`
  `params` | [Bus](Bus/README.md) parameters
  `mon`    | Monitors by README name (optional)

## Devices

  Name        | Value
  ------------|-------------------------------------------------------
  `class`     | Model class
  `id`        | Device id
  `terminals` | Bus number at each terminal
  `inputs`    | What each input reads: a bus number for its voltage, or `"owner.name"` for a variable or an output of device or bus `owner`. An input left out is held (optional)
  `outputs`   | Bus number whose current balance each output is added to (optional)
  `params`    | Model parameters by README `JSON` name
  `mon`       | Monitors by README name, with those of a submodel as `symbol.name` (optional)

A reference `"owner.name"` names device `owner`, or bus number `owner` when
no device has that id. A variable is read in place, and an output through a
signal that the system adds.

## State File

Initial values by bus number and device id: internal variables, held inputs,
and submodels nested by symbol. Omitted values are zero.

```json
{
  "buses": {"1": {"v": [1.0]}},
  "devices": {"LD1": {"i": [-0.5], "z": {"w": [0.25]}}, "Y1": {"u": [1.0]}}
}
```
