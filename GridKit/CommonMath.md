# CommonMath

Smooth, autodiff-friendly replacements for piecewise functions used across GridKit component models. See [CommonMath.hpp](CommonMath.hpp) for implementation details.

## Primitives

### Descriptions

| Name | Description | Usage |
|------|-------------|-------|
| `sigmoid` | Step function | `IEEET1` |
| `ramp` | Smooth one-sided ramp | `REGCA` |
| `qramp` | Exact one-sided quadratic ramp | `IEEET1` |

### Exact Equations

```math
\begin{aligned}
\sigma(x)
&=
\begin{cases}
0 & x \le 0 \\
1 & x > 0
\end{cases}
\\
\rho(x) &= x\,\sigma(x) \\
q(x) &= x^2\,\sigma(x)
\end{aligned}
```

### Smooth Equations

```math
\begin{aligned}
\sigma(x) &= \dfrac{1}{1+e^{-\mu x}} \\
\rho(x) &= \dfrac{(\mu x+\lvert\mu x\rvert)/2+\log(1+e^{-\lvert\mu x\rvert})}{\mu} \\
q(x) &= x^2\,\sigma(x)
\end{aligned}
```

The scale $\mu=4\cdot f_{\text{sync}}=240$ is chosen so $\sigma$ behaves like a step on inputs of order 1 while keeping derivatives finite. As $\mu \to \infty$, these functions approach their exact targets. (*Note*: the implementation of the quadratic ramp `q(x)` could be optimized with Enzyme features down the road).

## Derived Functions

### Descriptions

| Name | Description | Usage |
|------|-------------|-------|
| `max` | Smooth binary maximum | `REECA` |
| `min` | Smooth binary minimum | `REECA` |
| `clamp` | Bounded saturation | `IEEEST` |
| `deadband` | Signed two-sided deadband | `REECA` |
| `slew` | Symmetric slew-rate limiter | - |
| `linseg` | Saturated linear segment contribution | `REGCA`, `REECA` |
| `above` | Above-lower-limit indicator | - |
| `below` | Below-upper-limit indicator | - |
| `inside` | Interior pulse indicator | `IEEEST` |
| `outside` | Outside-band indicator | `REECA` |
| `antiwindup` | Anti-windup limited derivative | `IEEET1`, `TGOV1`, `SEXS-PTI`, `REECA` |

### Exact Equations

```math
\begin{aligned}
\text{max}(x,y)
&=
\begin{cases}
x & x > y \\
y & x \le y
\end{cases}
\\
\text{min}(x,y)
&=
\begin{cases}
x & x < y \\
y & x \ge y
\end{cases}
\\
\text{clamp}(x;\ell,u)
&=
\begin{cases}
\ell & x < \ell \\
x & \ell \le x \le u \\
u & x > u
\end{cases}
\\
\text{deadband}(x;\ell,u)
&=
\begin{cases}
x-\ell & x < \ell \\
0 & \ell \le x \le u \\
x-u & x > u
\end{cases}
\\
\text{slew}(f;r)
&=
\begin{cases}
-r & f < -r \\
f & -r \le f \le r \\
r & f > r
\end{cases}
\\
\text{linseg}(x;a,b,h)
&=
\begin{cases}
0 & x < a \\
\dfrac{h}{b-a}(x-a) & a \le x \le b \\
h & x > b
\end{cases}
\\
\text{above}(x;\ell)
&=
\begin{cases}
0 & x \le \ell \\
1 & x > \ell
\end{cases}
\\
\text{below}(x;u)
&=
\begin{cases}
1 & x < u \\
0 & x \ge u
\end{cases}
\\
\text{inside}(x;\ell,u)
&=
\begin{cases}
1 & \ell < x < u \\
0 & \text{else}
\end{cases}
\\
\text{outside}(x;\ell,u)
&=
\begin{cases}
1 & x < \ell \lor x > u \\
0 & \text{else}
\end{cases}
\\
\text{antiwindup}(x,f;\ell,u)
&=
\begin{cases}
f & \ell < x < u \\
f & x \le \ell \land f > 0 \\
f & x \ge u \land f < 0 \\
0 & \text{otherwise}
\end{cases}
\end{aligned}
```

### Smooth Equations

<a id="anti-windup-indicator"></a>

```math
\begin{aligned}
\text{max}(x,y) &= y+\rho(x-y) \\
\text{min}(x,y) &= x-\rho(x-y) \\
\text{clamp}(x;\ell,u) &= \ell+\rho(x-\ell)-\rho(x-u) \\
\text{deadband}(x;\ell,u) &= \rho(x-u)-\rho(\ell-x) \\
\text{slew}(f;r) &= -r+\rho(f+r)-\rho(f-r) \\
\text{linseg}(x;a,b,h) &= \dfrac{h}{b-a}\left[\rho(x-a)-\rho(x-b)\right] \\
\text{above}(x;\ell) &= \sigma(x-\ell) \\
\text{below}(x;u) &= \sigma(u-x) \\
\text{inside}(x;\ell,u) &= \sigma(x-\ell)+\sigma(u-x)-1 \\
\text{outside}(x;\ell,u) &= \sigma(\ell-x)+\sigma(x-u) \\
\phi_L &= \text{above}(x;\ell) \\
\phi_U &= \text{below}(x;u) \\
\phi(x,f) &= \phi_L\phi_U+(1-\phi_U)\sigma(-f)+(1-\phi_L)\sigma(f) \\
\text{antiwindup}(x,f;\ell,u) &= \phi(x,f)f
\end{aligned}
```
