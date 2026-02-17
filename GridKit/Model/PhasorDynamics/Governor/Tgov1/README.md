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
$P_{vmax}$  | [p.u.] | Max Valve Position                | 1    | 
$P_{vmin}$  | [p.u.] | Min Valve Position                | 0    | 
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
$P_{m}$         | [p.u.] | Mechnical Power to Generator      | Read by a Machine Model

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
The TGOV1 differential equations, as derived from the model diagram. By defining an auxillary valve function $f=-P_{v} + (P_{ref}-\omega)/R$ one can write the piecewise definition of $\dot P_v$ compactly.
```math
\begin{aligned}
   \dot{P}_{tx}   &= P_v - \dfrac{1}{T_3}(P_{tx}+T_2P_v) \\
   \dot{P}_{v}    &= \dfrac{1}{T_1}
   \begin{cases}
      -P_{v} + \dfrac{1}{R}(P_{ref}-\omega)
         &  \text{if } (P_{vmin} < P_v < P_{vmax}) & \lor \\
         &  \quad (P_v \leq P_{vmin} \land f>0)    & \lor \\
         &  \quad(P_v \geq P_{vmax} \land f<0)            \\
      0  
         &  \text{else}
   \end{cases}
\end{aligned}
```

### Algebraic Equations
The algebraic equation dictating the mechnical power output.
```math
\begin{aligned}
   P_{m} &= \dfrac{1}{T_3}(P_{tx}+T_2P_v) - D_t \omega \\
\end{aligned}
```

#### Smooth Piecewise Approximation
The domain of the state variable is enforced
through the piece-wise definition above. This may need to be expressed as a
smooth approximation (smooth indicator $\phi$) expressed generically as follows.
```math
\begin{aligned}
   f(P_v)      &= \dfrac{1}{T_1}
      \left[
         -P_{v} + \dfrac{1}{R}(P_{ref}-\omega)
      \right] \\
   \dot{P}_{v} &= 
            \phi(P_v, f) \cdot f(P_v)
\end{aligned}
```

The indicator function $\phi$ can be defined in terms of a scaled activation function ($\sigma$, sigmoid) and the $P_v$ limits $(P_{vmin}, P_{vmax})$.
```math
\begin{aligned}
   \phi_L(P_v)&= \sigma(P_v-V_{rmin}) \\
   \phi_U(P_v)&= \sigma(V_{rmax}-P_v) \\
   \phi_0(P_v)&= \phi_L + \phi_U - 1  \\
   \phi(P_v,f)&= \phi_0 \delta(f) + \left[ \phi_L \sigma(-f) + \phi_U \sigma(f)\right] \left[ 1-\delta(f) \right] \\
\end{aligned}
```

The scale of the sigmoid function ($\alpha$ on the order of $10^2$) should be chosen so that for all practical parameters of the TGOV1 model, the sigmoid acts as a step function.
```math
\begin{aligned}
   \sigma(x) =\dfrac{1}{1+\exp(-\alpha x)}
\end{aligned}
```

The derivative of the sigmoid then approximates the delta function.
```math
\begin{aligned}
   \delta(x) = 4 \sigma(x) \left[ 1- \sigma(x) \right]
\end{aligned}
```

## Initialization
At steady state we assume that $P_v$ is at or within its limits. This implies the initial conditions are a function of $P_{m}$ which is equal to the electric torque.
```math
\begin{aligned}
   P_{tx}  &= (T_3-T_2) P_{m}\\
   P_{v}   &= P_{m}\\
   \dot{P}_{tx} &=0\\
   \dot{P}_{v}  &=0\\
\end{aligned}
```

And if the reference power is a constant parameter, we can determine the value by solving the steady state equations.
```math
\begin{aligned}
   P_{ref}  &= R P_{m}\\
\end{aligned}
```
