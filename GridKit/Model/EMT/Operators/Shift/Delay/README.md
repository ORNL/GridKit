# Delay Model

`Delay` represents a smooth approximation of a transport delay on a vector input
signal. The approximation uses a chain of $n$ identical first-order lag stages.

The Laplace domain representation is:

```math
e^{-s\tau}\mathbf{U}(s) =
\lim_{n \to \infty}\left(\dfrac{1}{1 + s\tau/n}\right)^n \mathbf{U}(s)
```

The time domain convolutional form is:

```math
\mathbf{u}(t-\tau) = \delta(t-\tau) * \mathbf{u}(t)
```

## Block Diagram

<div align="center">
   <img align="center" src="../../../../../../docs/Figures/EMT/Delay/diagram.png">

  Figure 1: Lag-chain approximation of a pure delay
</div>

## Model Parameters

Symbol            | Units | JSON     | Description          | Typical Value | Note
------------------|-------|----------|----------------------|---------------|-----
$\tau$            | [s]   | `delay`  | Delay to approximate | --            | Required, positive
$\Delta t_{\min}$ | [s]   | `dt_min` | Block resolution     | --            | Required, positive

### Parameter Validation

```math
\begin{aligned}
\tau &> 0 \\
\Delta t_{\min} &> 0
\end{aligned}
```

### Model Derived Parameters

```math
n = \text{floor}\left(\dfrac{\tau}{\Delta t_{\min}}\right)
```

## Model Variables

### Internal Variables

#### Differential

Symbol                            | Units | Description                                               | Note
----------------------------------|-------|-----------------------------------------------------------|-----
$\mathbf{x}_1,\dots,\mathbf{x}_n$ | [-]   | Lag-block states; $\mathbf{x}_n$ is the delayed signal    | $nK$ states

#### Algebraic

None.

### External Variables

#### Differential

Symbol       | Units | Description | Note
-------------|-------|-------------|-----
$\mathbf{u}$ | [-]   | Input vector | $\mathbf{u} \in \mathbb{R}^K$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{u}$ | `input` | Input | [-] | Input vector port | $\mathbf{u} \in \mathbb{R}^K$
$\mathbf{y}$ | `out` | Output | [-] | Output contribution port | $\mathbf{y} \in \mathbb{R}^K$

## Model Equations

### Differential Equations

The lag-chain residuals are:

```math
\begin{aligned}
0 &= -\tau\,\dot{\mathbf{x}}_1 + n\,(\mathbf{u} - \mathbf{x}_1) \\
0 &= -\tau\,\dot{\mathbf{x}}_2 + n\,(\mathbf{x}_1 - \mathbf{x}_2) \\
&\vdots \\
0 &= -\tau\,\dot{\mathbf{x}}_n + n\,(\mathbf{x}_{n-1} - \mathbf{x}_n)
\end{aligned}
```

### Algebraic Equations

None.

### Port Equations

```math
\mathbf{y} = \mathbf{x}_n
```

## Initialization

For a constant input $\mathbf{u}_0$ at $t_0$, the chain is at rest:

```math
\begin{aligned}
\mathbf{x}_1(t_0) = \mathbf{x}_2(t_0) = \cdots = \mathbf{x}_n(t_0) &= \mathbf{u}_0 \\
\dot{\mathbf{x}}_1(t_0) = \dot{\mathbf{x}}_2(t_0) = \cdots = \dot{\mathbf{x}}_n(t_0) &= \mathbf{0}
\end{aligned}
```

A steady input therefore passes through unchanged at $t_0$ and downstream
consumers initialize consistently:

```math
\mathbf{y}_0 = \mathbf{x}_n(t_0)
```

## Monitors

None.
