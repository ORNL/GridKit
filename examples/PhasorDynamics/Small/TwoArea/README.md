# **Two-Area Kundur Case**

## One-Line Diagram

<div align="center">
   <img align="center" src="twoarea_oneline.png">
   
  Figure 1: Oneline of the two-area Case, courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Case Description

This is a model derived from the well-known Kundur case in Literature. It has been modified to include exciters for the sake of implementation within GridKit.

Model       | Count  
------------|--------
[Bus](../../../../src/Model/PhasorDynamics/Bus/README.md)         | 10
[Branch](../../../../src/Model/PhasorDynamics/Branch/README.md)     | 15 
[GENROU](../../../../src/Model/PhasorDynamics/SynchronousMachine/GENROUwS/README.md)       | 3
[TGOV1](../../../../src/Model/PhasorDynamics/Governor/Tgov1/README.md)        | 3
[IEEET1](../../../../src/Model/PhasorDynamics/Exciter/IEEET1/README.md)  | 3

## Case Events

The following examples needs to be constructed with this case. Currently the only events that have been implemented are faults. This is a practical case to use for modeling a line or generator outage.
- Bus Fault
- Line Outage
- Generator Outage
