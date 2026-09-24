
The bus model used for the microgrid found in Pogaku et al.[^1] and Bidram et al.[^2]

Parameters:
 + $R_N$ - Virtual Resistance 

Variables (External):
 + $v_{D}$          - Bus voltage along direct axis
 + $v_{Q}$          - Bus voltage along quadrature axis

Component external contributions:
 + $\frac{-v_D}{R_N}$
 + $\frac{-v_Q}{R_N}$

There are no internal variables to this system. Only residuals to be added from existing externals. As $RN \rightarrow \infty$ then the bus represent Kirchhoff's current law.


[^1]: Pogaku, Nagaraju, Milan Prodanovic, and Timothy C. Green. "Modeling, analysis and testing of autonomous operation of an inverter-based microgrid." IEEE Transactions on power electronics 22.2 (2007): 613-625.

[^2]: Bidram, Ali, Frank L. Lewis, and Ali Davoudi. "Distributed control systems for small-scale power networks: Using multiagent cooperative control theory." IEEE Control systems magazine 34.6 (2014): 56-77.
