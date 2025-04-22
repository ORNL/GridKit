# **General Synchronous Machine Model**

> [!NOTE]  
> Only the GENROU model has been implemented.

## Convention


<div align="center">
   <img align="center" src="../../../../docs/Figures/SM1.JPG">
   
   
  Figure 1: Synchronous Machine. Figure courtesy of [PowerWorld](https://www.powerworld.com/files/Synchronous-Machines.pdf)
</div>

The following conventians are used for the d-q reference frame.
- The q-axis leads the d-axis
- The Rotor angle is w.r.t. to q-axis

## Types

There are two main variations 

- Round Rotor (See [GENROU](GENROUwS/README.md))
- Salient Rotor/Pole (See [GENSAL](GENSALwS/README.md))
- GENPWS
- GENTPF
- GENTPJ
- GENQEC

### Per-Unit Basis

In relevant models, the terminal impedences are on the generator impedance base. To convert to network base, the following must be performed.
``` math
\begin{aligned}
Z_{term} &\mapsto Z_{term}\dfrac{S_{base,sys}}{S_{base,machine}}
\end{aligned}
```

#### Saturation

Saturation means increasingly large amounts of current are needed to increase the flux density. There are various methods to include the saturation (it is not standardized yet). We are going to use the approach implemented in PTI PSS/E and PowerWorld Simulator (scaled quadratic). 

``` math
\begin{aligned}
  k_{sat} = 
  \begin{cases}
    S_B(\psi''-S_A)^2 &\text{if } \psi''>S_A \\
    0 &\text{if } \psi''\leq S_A
  \end{cases}
\end{aligned}
```

