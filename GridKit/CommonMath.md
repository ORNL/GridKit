# CommonMath

Smooth, autodiff-friendly replacements for piecewise functions used across GridKit component models. See [CommonMath.hpp](CommonMath.hpp) for implementation details.

## Primitives

The scale $\mu=4\cdot f_{\text{sync}}=240$ is chosen so $\sigma$ behaves like a step on inputs while keeping derivatives finite. As $\mu \to \infty$, these functions approach their exact targets.

| Name | Description | Usage |
|------|-------------|-------|
| `sigmoid` | Step function | `GENSAL`, `GENROU`, `REGCA`, `REECA` |
| `ramp` | Smooth one-sided ramp | `REGCA`, `REECA`, `REPCA` |
| `qramp` | Exact one-sided quadratic ramp | `IEEET1` |

### $\sigma$ - `sigmoid`

The sigmoid is also known as the logistic function. The equivalent `tanh` form is used for numerical stability because the exponential form can divide by a very large value.

```math
\begin{aligned}
\sigma(x) &= \begin{cases}0 & x\le 0\\ 1 & x>0\end{cases} \\
          &\approx \dfrac{1}{2}\left(1+\tanh\left(\dfrac{\mu x}{2}\right)\right)
\end{aligned}
```

<div align="center"><img src="../docs/Figures/CommonMath/sigmoid.svg" width="250"></div>

### $\rho$ - `ramp`

`ramp` is the softplus approximation to the one-sided ramp. We do not use $x\sigma(x)$ directly because it introduces a negative tail for $x < 0$, while softplus stays nonnegative and approaches $\max(x, 0)$ as the smoothing becomes sharp.

```math
\begin{aligned}
\rho(x) &= x\,\sigma(x) \\
        &\approx \dfrac{x+\lvert x\rvert}{2}+\dfrac{\ln\!\left(1+e^{-\mu\lvert x\rvert}\right)}{\mu}
\end{aligned}
```

<div align="center"><img src="../docs/Figures/CommonMath/ramp.svg" width="250"></div>

### $q$ - `qramp`

*Note*: the implementation of the quadratic ramp `q(x)` could be optimized with Enzyme features down the road so that we don't need the smooth approximation.

```math
q(x)=x^2\,\sigma(x)
```

<div align="center"><img src="../docs/Figures/CommonMath/qramp.svg" width="250"></div>

## Derived Functions

| Name | Description | Usage |
|------|-------------|-------|
| `max` | Smooth binary maximum | `REECA`, `REECB` |
| `min` | Smooth binary minimum | `REECA` |
| `clamp` | Bounded saturation | `IEEEST`, `REECA`, `REECB`, `REPCA` |
| `deadband1` | Type 1 no-offset signed two-sided deadband | - |
| `deadband2` | Type 2 offset signed two-sided deadband | `REECA`, `REECB`, `REPCA` |
| `slew` | Symmetric slew-rate limiter | - |
| `linseg` | Saturated linear segment contribution | `REGCA`, `REECA` |
| `above` | Above-lower-limit indicator | `REPCA` |
| `below` | Below-upper-limit indicator | - |
| `inside` | Interior pulse indicator | - |
| `outside` | Outside-band indicator | `REECA`, `REECB` |
| `antiwindup` | Anti-windup limited derivative | `IEEET1`, `SEXS-PTI`, `TGOV1`, `REECA`, `REECB`, `REPCA` |

### `max`

```math
\begin{aligned}
\text{max}(x,y) &= \begin{cases}x & x>y\\ y & x\le y\end{cases} \\
                &\approx y+\rho(x-y)
\end{aligned}
```

<div align="center"><img src="../docs/Figures/CommonMath/max.svg" width="250"></div>

### `min`

```math
\begin{aligned}
\text{min}(x,y) &= \begin{cases}x & x<y\\ y & x\ge y\end{cases} \\
                &\approx x-\rho(x-y)
\end{aligned}
```

<div align="center"><img src="../docs/Figures/CommonMath/min.svg" width="250"></div>

### `clamp`

```math
\begin{aligned}
\text{clamp}(x;\ell,u) &= \begin{cases}\ell & x<\ell\\ x & \ell\le x\le u\\ u & x>u\end{cases} \\
                       &\approx \ell+\rho(x-\ell)-\rho(x-u)
\end{aligned}
```

<div align="center"><img src="../docs/Figures/CommonMath/clamp.svg" width="250"></div>

### `deadband1`

```math
\begin{aligned}
\text{deadband1}(x;\ell,u) &= \begin{cases}x & x<\ell\\ 0 & \ell\le x\le u\\ x & x>u\end{cases} \\
                           &\approx x\left[\sigma(\ell-x)+\sigma(x-u)\right]
\end{aligned}
```

<div align="center"><img src="../docs/Figures/CommonMath/deadband1.svg" width="250"></div>

### `deadband2`

```math
\begin{aligned}
\text{deadband2}(x;\ell,u) &= \begin{cases}x-\ell & x<\ell\\ 0 & \ell\le x\le u\\ x-u & x>u\end{cases} \\
                           &\approx \rho(x-u)-\rho(\ell-x)
\end{aligned}
```

<div align="center"><img src="../docs/Figures/CommonMath/deadband2.svg" width="250"></div>

### `slew`

```math
\begin{aligned}
\text{slew}(f;r) &= \begin{cases}-r & f<-r\\ f & -r\le f\le r\\ r & f>r\end{cases} \\
                 &\approx -r+\rho(f+r)-\rho(f-r)
\end{aligned}
```

<div align="center"><img src="../docs/Figures/CommonMath/slew.svg" width="250"></div>

### `linseg`

```math
\begin{aligned}
\text{linseg}(x;a,b,h) &= \begin{cases}0 & x<a\\ \dfrac{h}{b-a}(x-a) & a\le x\le b\\ h & x>b\end{cases} \\
                       &\approx \dfrac{h}{b-a}\left[\rho(x-a)-\rho(x-b)\right]
\end{aligned}
```

<div align="center"><img src="../docs/Figures/CommonMath/linseg.svg" width="250"></div>

### `above`

```math
\begin{aligned}
\text{above}(x;\ell) &= \begin{cases}0 & x\le \ell\\ 1 & x>\ell\end{cases} \\
                     &\approx \sigma(x-\ell)
\end{aligned}
```

<div align="center"><img src="../docs/Figures/CommonMath/above.svg" width="250"></div>

### `below`

```math
\begin{aligned}
\text{below}(x;u) &= \begin{cases}1 & x<u\\ 0 & x\ge u\end{cases} \\
                  &\approx \sigma(u-x)
\end{aligned}
```

<div align="center"><img src="../docs/Figures/CommonMath/below.svg" width="250"></div>

### `inside`

```math
\begin{aligned}
\text{inside}(x;\ell,u) &= \begin{cases}1 & \ell<x<u\\ 0 & \text{else}\end{cases} \\
                        &\approx \sigma(x-\ell)+\sigma(u-x)-1
\end{aligned}
```

<div align="center"><img src="../docs/Figures/CommonMath/inside.svg" width="250"></div>

### `outside`

```math
\begin{aligned}
\text{outside}(x;\ell,u) &= \begin{cases}1 & x<\ell\ \lor\ x>u\\ 0 & \text{else}\end{cases} \\
                         &\approx \sigma(\ell-x)+\sigma(x-u)
\end{aligned}
```

<div align="center"><img src="../docs/Figures/CommonMath/outside.svg" width="250"></div>

### `antiwindup`

<a id="anti-windup-indicator"></a>

```math
\begin{aligned}
\text{antiwindup}(x,f;\ell,u) &= \begin{cases}f & \ell<x<u\\ f & x\le\ell\ \land\ f>0\\ f & x\ge u\ \land\ f<0\\ 0 & \text{otherwise}\end{cases} \\
\phi_L &= \text{above}(x;\ell) \\
\phi_U &= \text{below}(x;u) \\
\phi(x,f) &= \phi_L\phi_U+(1-\phi_U)\sigma(-f)+(1-\phi_L)\sigma(f) \\
\text{antiwindup}(x,f;\ell,u) &\approx \phi(x,f)\,f
\end{aligned}
```
