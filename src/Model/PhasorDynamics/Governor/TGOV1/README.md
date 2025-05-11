# **Steam Turbine-Governor Model (TGOV1)**

## Control Diagram

Standard model of the stream turbine

<div align="center">
   <img align="center" src="../../../../../docs/Figures/TGOV1.JPG">
   
   
  Figure 1: Governor TGOV1 model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Nomenclature

### Algebraic Variables
- $\Delta\omega$ - Per-unit Speed Deviation (p.u.)
- $P_{ref}$ - Reference Power (p.u.)
- $P_{mech}$ - Mechnical Power to Generator (p.u.)

### Differential Variables
- $P_{tx}$ - Turbine Power (p.u.) (State 1 in Fig. 1)
- $P_{v}$ - Valve Position (State 2 in Fig. 1)

### Parameters
- $R$ - Droop Constant (p.u.)
- $T_1$ - Valve Time Delay (sec)
- $T_2$ - Turbine Numerator Time Constat (sec)
- $T_3$ - Turbine Delay (sec)
- $P_{vmax}$ - Max Valve Position
- $P_{vmin}$ -  Min Valve Position
- $D_t$ - Turbine Damping Coefficient (p.u.)

## Equations


### Algebraic Equations
The algebraic equation dictating the mechnical power output.
```math
\begin{aligned}
   P_{mech} &= \dfrac{1}{T_3}(P_{tx}+T_2P_v) - D_t \Delta\omega \\
\end{aligned}
```

### Differential Equations
The TGOV1 differential equations, as derived from the model diagram.
```math
\begin{aligned}
   \dot{P}_{tx}   &= P_v - \dfrac{1}{T_3}(P_{tx}+T_2P_v) \\
   \dot{P}_{v}    &= 
      \dfrac{1}{T_1}
   \begin{cases}
      -P_{v} + \dfrac{1}{R}(P_{ref}-\Delta\omega)
         &  \text{if } P_v \in (P_{vmin}, P_{vmax})\\
      0  
         &  \text{else } \\
   \end{cases}
\end{aligned}
```
The domain of the state variable $P_{v}\in(P_{vmin}, P_{vmax})$ is enforced
through the piece-wise definition above. This may need to be expressed as a
smooth approximation (smooth indicator $\phi$) expressed generically as follows.
```math
\begin{aligned}
   \dot{P}_{v}   
            &= 
            \phi(P_v) \cdot \dfrac{1}{T_1}
            \left[
               -P_{v} + \dfrac{1}{R}(P_{ref}-\Delta\omega)
            \right] \\
\end{aligned}
```
