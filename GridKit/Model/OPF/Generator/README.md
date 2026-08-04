# Generator

An online generator injects $S_g=p_g+\mathrm jq_g$ and contributes

$$f_g=c_0+c_1p_g+c_2p_g^2.$$

Optional `pmin`, `pmax`, `qmin`, and `qmax` bound its variables. `mva` is
positive metadata; cost coefficients apply directly to per-unit $p_g$. An
offline generator is fixed to zero and contributes no cost. See the
[system formulation](../README.md#formulation).
