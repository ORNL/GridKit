# EXST1-PTI

A lead–lag stage drives the regulator with washout feedback. Both field-voltage limits track terminal voltage and commutation loading, from $V_t V_{rmin}-K_c I_{fd}$ to $V_t V_{rmax}-K_c I_{fd}$.

The `Exst1Pti` device reads voltage magnitude from its `bus` and exciter-base field current from the machine’s `ifd` output. The machine initializes the `efd` output. Optional inputs are `vref`, `vs`, `vuel`, and `voel`.
