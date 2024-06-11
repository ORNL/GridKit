## Component Models

GridKit™ provides component models for power flow and electromechanical transient simulations, as well as experimental component models for dynamic constrained optimal power flow analysis. GridKit™ assembles components into a grid model using power flow equations. 

## Network Equations

The relation between all the bus current injections and bus voltage is given by the following node equation: 
```math
I_{i}=\sum_{j=1}^{n} Y_{ij}V_{j} ~~~ i=1,2,...,n
```
where  $`n`$ is the number of buses in the network. $`I_{i}`$ and $`V_{j}`$ are injected current at bus $`i`$ and voltage at bus $`j`$. $`Y_{ij}`$ are the elements of the admittance matrix **Y**. Diagonal elements $`Y_{ii}`$ are equal to the sum of all admittances of all devices incident to the bus $`i`$. Off-diagonal elements $`Y_{ij}`$ are equal to the **negative** of the sum of the admittances that are joining buses $`i`$ and $`j`$. In case that there is shift transformer at the bus, $`Y_{ij}`$ should be calculated as explained in the branch section.

In the power system, complex voltage and current values are unknown, but rather real power injections at the generator buses and voltage magnitude setpoint as well as complex power (S) consumed by the load.
The relation between injected current and power at the node is given as:
```math
S_{i}=P_{i}+jQ_{i}=V_{i}I^*_{i}
```
Mode information on the sign convention used can be found [here](./Bus/README.md).

Now the complex power flow equations are given as:

```math
P_{i}+jQ_{i}=V_{i}\sum_{j=1}^{n} Y_{ij}^*V_{j}^* \;\;\;\;\;\; i=1,2,...,n
```

Considering:
```math
V_{i}=\vert V_{i} \vert e^{j\theta_{i}}
```

```math
Y_{ij}=\vert Y_{ij} \vert e^{j\psi_{ij}}=G_{ij}+jB_{ij}
```

```math
\theta_{ij}=\theta_{i}-\theta_{j}
```

The complex set can be rewritten as two sets of real power balance equations:

```math
P_{i}= \vert V_{i} \vert \sum_{j=1}^{n}\vert Y_{ij} \vert \vert V_{j} \vert \cos (\theta_{ij}-\psi_{ij})
```

```math
Q_{i}= \vert V_{i} \vert \sum_{j=1}^{n} \vert Y_{ij} \vert \vert V_{j} \vert \sin (\theta_{ij}-\psi_{ij})
```

or

```math
P_{i}= \vert V_{i} \vert \sum_{j=1}^{n}\vert V_{j} \vert (G_{ij}\cos\theta_{ij}+B_{ij}\sin\theta_{ij})
```

```math
Q_{i}= \vert V_{i} \vert \sum_{j=1}^{n} \vert V_{j} \vert (G_{ij}\sin\theta_{ij}-B_{ij}\cos\theta_{ij})
```
