

The model is an expanded form of the model in reference 1 and 2. Utilizes the PowerElectronicsModel composer framework. Given a free parameter $N \ge 1$ then the amount of components the system have are as follows.
 + Distributed Generators : $`2N`$
 + Microgrid Lines : $`2N-1`$
 + Microgrid Loads : $`N`$
 + Microgrid Buses : $`N`$

When $N=2$ the original case in reference 1 is constructed.
All parameters are given by duplicating the parameters in reference 1 over the extended system. The system extends by connecting new systems over the even indexed buses.
Reference solutions are provided for the cases of when $N=\{2, 4, 8\}$. These are generated from a MATLAB based code with tight tolerances.


1. Pogaku, Nagaraju, Milan Prodanovic, and Timothy C. Green. "Modeling, analysis and testing of autonomous operation of an inverter-based microgrid." IEEE Transactions on power electronics 22.2 (2007): 613-625.
2. Bidram, Ali, Frank L. Lewis, and Ali Davoudi. "Distributed control systems for small-scale power networks: Using multiagent cooperative control theory." IEEE Control systems magazine 34.6 (2014): 56-77.