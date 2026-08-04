# Bus

A bus owns the polar voltage variables $v_i,\theta_i$ and the active and
reactive balance rows. Its complex voltage is

$$V_i=v_i e^{\mathrm j\theta_i}.$$

`kv` defines the voltage base; optional `vmin` and `vmax` bound $v_i$. The
lowest-numbered bus fixes $\theta_i$ to its input-state angle. See the
[system formulation](../README.md#formulation).
