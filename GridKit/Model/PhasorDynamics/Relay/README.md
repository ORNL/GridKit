# Relay Models

Relays measure network quantities and issue breaker commands.

## Types

- [DefiniteTimeRelay](DefiniteTimeRelay/README.md): Definite-time relay with optional lockout
- [InverseTimeRelay](InverseTimeRelay/README.md): Inverse-time overcurrent relay with optional lockout

## Elements

ID[^device-numbers] | Short Name | Long Name                | Model
--------------------|------------|--------------------------|------------------
27                  | UV         | Undervoltage             | DefiniteTimeRelay
50                  | OC         | Overcurrent              | DefiniteTimeRelay
51                  | TOC        | Inverse-time overcurrent | InverseTimeRelay
59                  | OV         | Overvoltage              | DefiniteTimeRelay
81O                 | OF         | Overfrequency            | DefiniteTimeRelay
81U                 | UF         | Underfrequency           | DefiniteTimeRelay

[^device-numbers]: [IEEE Std C37.2-2022, *IEEE Standard for Electrical Power System Device Function Numbers, Acronyms, and Contact Designations*](https://standards.ieee.org/ieee/C37.2/6043/).
