# Load

A load is a fixed signed injection read from `.state.json`:

$$S_d^{\mathrm{inj}}=a_d(p_d+\mathrm jq_d).$$

Ordinary demand has negative `p` and `q`; an offline load contributes zero.
Optional power limits validate the fixed state and do not create decision
variables. See the [system formulation](../README.md#formulation).
