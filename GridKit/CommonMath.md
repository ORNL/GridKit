# CommonMath

Smooth, autodiff-friendly replacements for piecewise functions used across GridKit component models. See `CommonMath.hpp` for implementation details.

## Primitives

The scale $\mu=4\cdot f_{\text{sync}}=240$ is chosen so $\sigma$ behaves like a step on inputs while keeping derivatives finite. As $\mu \to \infty$, these functions approach their exact targets.

| Name | Description | Usage |
|------|-------------|-------|
| `sigmoid` | Step function | `GENSAL`, `GENROU`, `REECA` |
| `ramp` | Smooth one-sided ramp | `REGCA`, `REECA`, `REPCA` |
| `qramp` | Exact one-sided quadratic ramp | `IEEET1` |

### $\sigma$ - `sigmoid`

The sigmoid is also known as the logistic function. The equivalent `tanh` form is used for numerical stability because the exponential form can divide by a very large value.

```math
\begin{aligned}
  \sigma(x)
    &=
      \begin{cases}
        0 & x\le 0 \\[0pt]
        1 & x\gt 0
      \end{cases} \\[0pt]
    &\approx \dfrac{1}{2}\left(1+\tanh\left(\dfrac{\mu x}{2}\right)\right)
\end{aligned}
```

![](../docs/Figures/CommonMath/sigmoid.svg)

### $\rho$ - `ramp`

`ramp` is the softplus approximation to the one-sided ramp. We do not use $x\sigma(x)$ directly because it introduces a negative tail for $x \lt 0$, while softplus stays nonnegative and approaches $\max(x, 0)$ as the smoothing becomes sharp.

```math
\begin{aligned}
  \rho(x)
    &= x\,\sigma(x) \\[0pt]
    &\approx \dfrac{x+\lvert x\rvert}{2}+\dfrac{\ln\!\left(1+e^{-\mu\lvert x\rvert}\right)}{\mu}
\end{aligned}
```

![](../docs/Figures/CommonMath/ramp.svg)

### $q$ - `qramp`

*Note*: the implementation of the quadratic ramp `q(x)` could be optimized with Enzyme features down the road so that we don't need the smooth approximation.

```math
q(x)=x^2\,\sigma(x)
```

![](../docs/Figures/CommonMath/qramp.svg)

## Derived Functions

| Name | Description | Usage |
|------|-------------|-------|
| `max` | Smooth binary maximum | `REGCA`, `REECA`, `REECB` |
| `min` | Smooth binary minimum | `REGCA`, `REECA` |
| `clamp` | Bounded saturation | `IEEEST`, `REGCA`, `REECA`, `REECB`, `REPCA` |
| `deadband1` | Type 1 no-offset signed two-sided deadband | - |
| `deadband2` | Type 2 offset signed two-sided deadband | `REECA`, `REECB`, `REPCA` |
| `slew` | Symmetric slew-rate limiter | - |
| `linseg` | Saturated linear segment contribution | `REGCA`, `REECA` |
| `above` | Above-lower-limit indicator | `REPCA` |
| `below` | Below-upper-limit indicator | - |
| `inside` | Interior pulse indicator | `REECB` |
| `outside` | Outside-band indicator | `REECA` |
| `antiwindup` | Anti-windup limited derivative | `IEEET1`, `SEXS-PTI`, `TGOV1`, `REECA`, `REECB`, `REPCA` |

### `max`

```math
\begin{aligned}
  \text{max}(x,y)
    &=
      \begin{cases}
        x & x\gt y \\[0pt]
        y & x\le y
      \end{cases} \\[0pt]
    &\approx y+\rho(x-y)
\end{aligned}
```

![](../docs/Figures/CommonMath/max.svg)

### `min`

```math
\begin{aligned}
  \text{min}(x,y)
    &=
      \begin{cases}
        x & x\lt y \\[0pt]
        y & x\ge y
      \end{cases} \\[0pt]
    &\approx x-\rho(x-y)
\end{aligned}
```

![](../docs/Figures/CommonMath/min.svg)

### `clamp`

```math
\begin{aligned}
  \text{clamp}(x;\ell,u)
    &=
      \begin{cases}
        \ell & x\lt \ell \\[0pt]
        x    & \ell\le x\le u \\[0pt]
        u    & x\gt u
      \end{cases} \\[0pt]
    &\approx \ell+\rho(x-\ell)-\rho(x-u)
\end{aligned}
```

![](../docs/Figures/CommonMath/clamp.svg)

### `deadband1`

```math
\begin{aligned}
  \text{deadband1}(x;\ell,u)
    &=
      \begin{cases}
        x & x\lt \ell \\[0pt]
        0 & \ell\le x\le u \\[0pt]
        x & x\gt u
      \end{cases} \\[0pt]
    &\approx x\left[\sigma(\ell-x)+\sigma(x-u)\right]
\end{aligned}
```

![](../docs/Figures/CommonMath/deadband1.svg)

### `deadband2`

```math
\begin{aligned}
  \text{deadband2}(x;\ell,u)
    &=
      \begin{cases}
        x-\ell & x\lt \ell \\[0pt]
        0      & \ell\le x\le u \\[0pt]
        x-u    & x\gt u
      \end{cases} \\[0pt]
    &\approx \rho(x-u)-\rho(\ell-x)
\end{aligned}
```

![](../docs/Figures/CommonMath/deadband2.svg)

### `slew`

```math
\begin{aligned}
  \text{slew}(f;r)
    &=
      \begin{cases}
        -r & f\lt -r \\[0pt]
        f  & -r\le f\le r \\[0pt]
        r  & f\gt r
      \end{cases} \\[0pt]
    &\approx -r+\rho(f+r)-\rho(f-r)
\end{aligned}
```

![](../docs/Figures/CommonMath/slew.svg)

### `linseg`

```math
\begin{aligned}
  \text{linseg}(x;a,b,h)
    &=
      \begin{cases}
        0                       & x\lt a \\[0pt]
        \dfrac{h}{b-a}(x-a)     & a\le x\le b \\[0pt]
        h                       & x\gt b
      \end{cases} \\[0pt]
    &\approx \dfrac{h}{b-a}\left[\rho(x-a)-\rho(x-b)\right]
\end{aligned}
```

![](../docs/Figures/CommonMath/linseg.svg)

### `above`

```math
\begin{aligned}
  \text{above}(x;\ell)
    &=
      \begin{cases}
        0 & x\le \ell \\[0pt]
        1 & x\gt \ell
      \end{cases} \\[0pt]
    &\approx \sigma(x-\ell)
\end{aligned}
```

![](../docs/Figures/CommonMath/above.svg)

### `below`

```math
\begin{aligned}
  \text{below}(x;u)
    &=
      \begin{cases}
        1 & x\lt u \\[0pt]
        0 & x\ge u
      \end{cases} \\[0pt]
    &\approx \sigma(u-x)
\end{aligned}
```

![](../docs/Figures/CommonMath/below.svg)

### `inside`

```math
\begin{aligned}
  \text{inside}(x;\ell,u)
    &=
      \begin{cases}
        1 & \ell\lt x\lt u \\[0pt]
        0 & \text{else}
      \end{cases} \\[0pt]
    &\approx \sigma(x-\ell)+\sigma(u-x)-1
\end{aligned}
```

![](../docs/Figures/CommonMath/inside.svg)

### `outside`

```math
\begin{aligned}
  \text{outside}(x;\ell,u)
    &=
      \begin{cases}
        1 & x\lt \ell\ \lor\ x\gt u \\[0pt]
        0 & \text{else}
      \end{cases} \\[0pt]
    &\approx \sigma(\ell-x)+\sigma(x-u)
\end{aligned}
```

![](../docs/Figures/CommonMath/outside.svg)

### `antiwindup`

<a id="anti-windup-indicator"></a>

```math
\begin{aligned}
  \text{antiwindup}(x,f;\ell,u)
    &=
      \begin{cases}
        f & \ell\lt x\lt u \\[0pt]
        f & x\le\ell\ \land\ f\gt 0 \\[0pt]
        f & x\ge u\ \land\ f\lt 0 \\[0pt]
        0 & \text{otherwise}
      \end{cases} \\[0pt]
  \phi_L
    &= \text{above}(x;\ell) \\[0pt]
  \phi_U
    &= \text{below}(x;u) \\[0pt]
  \phi(x,f)
    &= \phi_L\phi_U+(1-\phi_U)\sigma(-f)+(1-\phi_L)\sigma(f) \\[0pt]
  \text{antiwindup}(x,f;\ell,u)
    &\approx \phi(x,f)\,f
\end{aligned}
```
