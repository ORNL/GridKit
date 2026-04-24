# CommonMath

Numerical utilities in [CommonMath.hpp](CommonMath.hpp): smooth, autodiff-friendly replacements for piecewise functions used across GridKit component models.

## Sigmoid

Smooth approximation to the step function.

```math
\sigma(x) = \dfrac{1}{1+\exp(-\alpha x)}
```

The scale $\alpha$ (currently $240$) is chosen large enough that $\sigma$ behaves as a step on inputs of order 1.

## Limit Indicators

For a state variable $x$ with limits $(x_{\min}, x_{\max})$:

```math
\begin{aligned}
   \phi_L(x) &= \sigma(x - x_{\min}) \\
   \phi_U(x) &= \sigma(x_{\max} - x) \\
   \phi_0(x) &= \phi_L + \phi_U - 1
\end{aligned}
```

$\phi_L$ and $\phi_U$ are soft indicators for "above lower limit" and "below upper limit"; their product (or $\phi_0$) is an interior indicator.

## Anti-Windup Indicator

Component models with a limited state $x \in (x_{\min}, x_{\max})$ and pre-limit derivative $f$ express the gated dynamics piecewise as

```math
\dot x =
   \begin{cases}
      f
         &  \text{if } (x_{\min} < x < x_{\max}) & \lor \\
         &  \quad (x \leq x_{\min} \land f > 0)  & \lor \\
         &  \quad (x \geq x_{\max} \land f < 0)         \\
      0  &  \text{else}
   \end{cases}
```

In simulation this is replaced with the smooth approximation $\dot x = \phi(x, f) \cdot f$, where

```math
\phi(x, f) = \phi_L \phi_U + (1 - \phi_U)\,\sigma(-f) + (1 - \phi_L)\,\sigma(f).
```

The first term passes the dynamics through in the interior. The second re-admits them when $x$ is above $x_{\max}$ and $f$ is pulling back down; the third re-admits them when $x$ is below $x_{\min}$ and $f$ is driving back up. When $x$ is outside a limit *and* $f$ is still driving further beyond it, $\phi$ vanishes smoothly and blocks the windup.

## Callers

Anti-windup indicator (`Math::indicator`):

- [IEEET1](Model/PhasorDynamics/Exciter/IEEET1/README.md): gates $\dot V_R$ on $V_R \in (V_{rmin}, V_{rmax})$
- [TGOV1](Model/PhasorDynamics/Governor/Tgov1/README.md): gates $\dot P_v$ on $P_v \in (P_{vmin}, P_{vmax})$
- [SEXS-PTI](Model/PhasorDynamics/Exciter/SEXS-PTI/README.md): gates $\dot E_{fd}$ on $E_{fd} \in (E_{fd,\min}, E_{fd,\max})$

Interior indicator (`Math::indicator_zero`):

- [IEEEST](Model/PhasorDynamics/Stabilizer/IEEEST/README.md): clips the stabilizer output $v_7$ to $[L_{s\min}, L_{s\max}]$, passing $v_7$ through in the interior and saturating to the limits outside
