# **Synthetic Hawaii Case**

## One-Line Diagram

<div align="center">
   <img align="center" src="newengland-oneline.png">
   
  Figure 1: Oneline of the New England IEEE 39-bus case, courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Case Description

This case is a modified version of the NE 39-bus case, available at the Texas A&M University [Electric Grid Test Case Repository](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/). It has been modified to include GridKit-compatible generator models.

Model       | Count  
------------|--------
[Bus](../../../../GridKit/Model/PhasorDynamics/Bus/README.md)         | 39
[Branch](../../../../GridKit/Model/PhasorDynamics/Branch/README.md)     | 46
[GENROU](../../../../GridKit/Model/PhasorDynamics/SynchronousMachine/GENROUwS/README.md)       | 10
[TGOV1](../../../../GridKit/Model/PhasorDynamics/Governor/Tgov1/README.md)        | 10
[IEEET1](../../../../GridKit/Model/PhasorDynamics/Exciter/IEEET1/README.md)  | 10

## Case Events

The following examples needs to be constructed with this case. Currently the only events that have been implemented are faults. This is a practical case to use for modeling outages, or even forced oscillations.
- Bus Fault
- Line Outage
- Generator Outage
- Forced Oscillations
