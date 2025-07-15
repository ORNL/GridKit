
The Load used for the Microgrid found in references 1 and 2.

Parameters:
 + $R$ - Resistance
 + $L$ - Inductance

Variables (External):
 + $\omega_{ref}$ - Reference Rotor Angle
 + $v_{D}$   - Incoming Bus Voltage (D)
 + $v_{Q}$   - Incoming Bus Voltage (Q)
 + $i_{D}$    - Line Current (D)
 + $i_{Q}$    - Line Current (Q)

Equations (External, Residuals):
 + $0$ &emsp;&emsp; (Reference Rotor Residual)
 + $-i_{D}$ &emsp;&emsp; (Bus Residuals)
 + $-i_{Q}$

Equations (Internal):
 + $\frac{di_{D}}{dt} = -(\frac{R}{L}) i_{D} + \omega_{ref} i_{Q} + \frac{v_{D1} - v_{D2}}{L}$
 + $\frac{di_{Q}}{dt} = -(\frac{R}{L}) i_{Q} + \omega_{ref} i_{D} + \frac{v_{Q1} - v_{Q2}}{L}$



1. Pogaku, Nagaraju, Milan Prodanovic, and Timothy C. Green. "Modeling, analysis and testing of autonomous operation of an inverter-based microgrid." IEEE Transactions on power electronics 22.2 (2007): 613-625.
2. Bidram, Ali, Frank L. Lewis, and Ali Davoudi. "Distributed control systems for small-scale power networks: Using multiagent cooperative control theory." IEEE Control systems magazine 34.6 (2014): 56-77.
