# **Steam Turbine-Governor Model (TGOV1)**

## Control Diagram

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
$P_{vmax}$  | [p.u.] | Stator leakage reactance          | 1    | 
$P_{vmin}$  | [p.u.] | Max Valve Position                | 0    | 
$D_t$       | [p.u.] | Turbine Damping Coefficient       | 0    | 

## Model Variables 

### Internal Variables

#### Differential

Symbol    | Units  | Description                       | Note
----------|--------|-----------------------------------|-------
$P_{tx}$  | [p.u.] | Turbine Power (State 1 in Fig. 1) |
$P_{v}$   | [p.u.] | Valve Position (State 2 in Fig. 1)| 

#### Algebraic
Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$P_{mech}$      | [p.u.] | Mechnical Power to Generator      | Read by a Machine Model

### External Variables

#### Differential
Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$\Delta\omega$  | [p.u.] | Speed Deviation                   |

#### Algebraic
Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$P_{ref}$       | [p.u.] | Reference Power                   | 

## Model Equations

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

### Algebraic Equations
The algebraic equation dictating the mechnical power output.
```math
\begin{aligned}
   P_{mech} &= \dfrac{1}{T_3}(P_{tx}+T_2P_v) - D_t \Delta\omega \\
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

## Initialization
Assuming no limits are reached, the initial conditions are a function of $P_{mech}$ which is equal to the electric torque.
```math
\begin{aligned}
   P_{tx}  &= (T_3-T_2) P_{mech}\\
   P_{v}   &= P_{mech}\\
   \dot{P}_{tx} &=0\\
   \dot{P}_{v}  &=0\\
\end{aligned}
```

And if the reference power is a parameter, we can deduce
```math
\begin{aligned}
   P_{ref}  &= R P_{mech}\\
\end{aligned}
```
