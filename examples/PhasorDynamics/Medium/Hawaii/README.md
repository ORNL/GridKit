# **Synthetic Hawaii**

## One-Line Diagram

![](hawaii.png)

Figure 1: Oneline of the synthetic Hawaii case, courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/) (WIP Updated Oneline)

## Case Description

This case is a modified version of the synthetic Hawaii case, available at the Texas A&M University [Electric Grid Test Case Repository](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/). It has been modified to include GridKit-compatible generator models.

## Data Notes

None.

## Outstanding

### Dynamics

The following models are not implemented in GridKit and are represented using surrogate models.

- IEEEST (In GridKit, not yet added to this case)
- REGCA
- IEEEG1
- GGOV1
- ESST4B
- REECA1
- EXST1_PTI
- ESST1A

The following examples needs to be constructed with this case. Currently the only events that have been implemented are faults. This is a practical case to use for modeling outages, or even forced oscillations.
- Line Outage
- Generator Outage
- Forced Oscillations

### Statics

GridKit models Switched Shunts as [LoadZ](../../../../GridKit/Model/PhasorDynamics/Load/LoadZ/README.md) with constant impedance,but this erases information and our ability to statically re-initialze the model if initalizing at a different operating point.

- Transformers
- Switched Shunts
