# Illinois

## One-Line Diagram

![](illinois.png)

Figure 1: Oneline of the ACTIVSg200 Case, from Texas A&M University [Grid Repository](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/activsg200/)

## Case Description

Geographically located in the state of Illinois, the ACTIVSg200 case is a 200 bus power system test case that is entirely synthetic, built from public information and a statistical analysis of real power systems. It bears no relation to the actual grid in this location, except that generation and load profiles are similar. The dynamics of this case are fully modeled in GridKit.

## Model Inventory

Model | Count
---|---
Bus | 200
Branch | 246
GENROU | 40
TGOV1 | 40
SEXS-PTI | 40
LoadZIP | 164
SignalNode | 120

## Data Notes

The reference `.pwb` file had no machine model at `Bus 197 GIBSON CITY 1 2`, so a `GENROU` model was added. The case does not initialize in steady state without it, and a negative-impedance load is not an appropriate substitute.

## Outstanding

### Statics

- Transformers
- Switched Shunts
