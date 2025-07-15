
The Bus used for the Microgrid found in references 1 and 2.

Parameters:
 + $RN$ - Virtual Resistance 

Variables (External):
 + $v_{D}$          - Incoming Bus Voltage (D)
 + $v_{Q}$          - Incoming Bus Voltage (Q)

Equations (External, Residuals):
 + $\frac{-v_D}{RN}$
 + $\frac{-v_Q}{RN}$

There are no internal variables to this system. Only residuals to be added from existing externals. As $RN \rightarrow \infty$ then the bus represent Kirchhoff's current law.


1. Pogaku, Nagaraju, Milan Prodanovic, and Timothy C. Green. "Modeling, analysis and testing of autonomous operation of an inverter-based microgrid." IEEE Transactions on power electronics 22.2 (2007): 613-625.
2. Bidram, Ali, Frank L. Lewis, and Ali Davoudi. "Distributed control systems for small-scale power networks: Using multiagent cooperative control theory." IEEE Control systems magazine 34.6 (2014): 56-77.
