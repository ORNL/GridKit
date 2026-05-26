# **Steam Turbine-Governor Model (TGOV1)**

## Block Diagram

Standard model of the stream turbine

<div align="center">
   <img align="center" src="../../../../../docs/Figures/TGOV1.JPG">
   
  Figure 1: Governor TGOV1 model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Model Parameters

Symbol      | Units  | Description                       | Typical Value | Note
------------|--------|-----------------------------------|---------------| ------
$R$         | [p.u.] | Droop Constant                    | 0.05 |
$T_1$       | [sec]  | Valve Time Delay                  | 0.5  |
$T_2$       | [sec]  | Turbine Numerator Time Constant   | 2.5  |
$T_3$       | [sec]  | Turbine Delay                     | 7.5  |
$P_v^{\max}$ | [p.u.] | Max Valve Position                | 1    |
$P_v^{\min}$ | [p.u.] | Min Valve Position                | 0    |
$D_t$       | [p.u.] | Turbine Damping Coefficient       | 0    | 

## Model Variables 

### Internal Variables

#### Differential

Symbol    | Units  | Description                       | Note
----------|--------|-----------------------------------|-------
$P_{tx}$  | [p.u.] | Turbine Power (State 1 in Fig. 1) |
$P_v$     | [p.u.] | Valve Position (State 2 in Fig. 1)|

#### Algebraic
Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$P_m$           | [p.u.] | Mechnical Power to Generator      | Read by a Machine Model

### External Variables

#### Differential
Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$\omega$  | [p.u.] | Machine Speed Deviation                   | Read from a Machine Model

#### Algebraic
Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$P_{ref}$       | [p.u.] | Reference Power                   | Either a constant parameter or external variable

## Model Equations

### Differential Equations
The TGOV1 differential equations, as derived from the model diagram. Define the pre-limit derivative of $P_v$

```math
f = \dfrac{1}{T_1}\left[-P_v + \dfrac{1}{R}(P_{ref} - \omega)\right]
```

so that $\dot P_v$ can be written in piecewise form compactly.
```math
\begin{aligned}
   \dot P_{tx}   &= P_v - \dfrac{1}{T_3}(P_{tx}+T_2P_v) \\
   \dot P_v      &=
   \begin{cases}
      f
         &  \text{if } (P_v^{\min} < P_v < P_v^{\max}) & \lor \\
         &  \quad (P_v \leq P_v^{\min} \land f>0)       & \lor \\
         &  \quad(P_v \geq P_v^{\max} \land f<0)            \\
      0  
         &  \text{else}
   \end{cases}
\end{aligned}
```

### Algebraic Equations
The algebraic equation dictating the mechnical power output.
```math
\begin{aligned}
   P_m &= \dfrac{1}{T_3}(P_{tx}+T_2P_v) - D_t \omega \\
\end{aligned}
```

In simulation the piecewise form above is replaced with a smooth approximation where $\phi$ is GridKit's smooth anti-windup indicator. See [CommonMath: Anti-Windup Indicator](../../../../CommonMath.md#anti-windup-indicator) for its definition, behavior, and design rationale.

## Initialization
At steady state we assume that $P_v$ is at or within its limits. This implies the initial conditions are a function of $P_m$ which is equal to the electric torque.
```math
\begin{aligned}
   P_{tx}  &= (T_3-T_2) P_m\\
   P_v     &= P_m\\
   \dot P_{tx} &=0\\
   \dot P_v    &=0\\
\end{aligned}
```

And if the reference power is a constant parameter, we can determine the value by solving the steady state equations.
```math
\begin{aligned}
   P_{ref}  &= R P_m\\
\end{aligned}
```
