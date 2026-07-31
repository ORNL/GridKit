# **Two-Area Kundur Case**

## One-Line Diagram

![](twoarea_oneline.png)

Figure 1: Oneline of the two-area Case, courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Case Description

This is a model derived from the well-known Kundur case in Literature. It has been modified to include exciters for the sake of implementation within GridKit.

Model       | Count  
------------|--------
[Bus](../../../../GridKit/Model/PhasorDynamics/Bus/README.md)         | 10
[Branch](../../../../GridKit/Model/PhasorDynamics/Branch/README.md)     | 15 
[GENROU](../../../../GridKit/Model/PhasorDynamics/Machine/GENROU/README.md)       | 3
[TGOV1](../../../../GridKit/Model/PhasorDynamics/Governor/Tgov1/README.md)        | 3
[IEEET1](../../../../GridKit/Model/PhasorDynamics/Exciter/IEEET1/README.md)  | 3

## Case Events

The following examples needs to be constructed with this case. Currently the only events that have been implemented are faults. This is a practical case to use for modeling a line or generator outage.
- Bus Fault
- Line Outage
- Generator Outage
