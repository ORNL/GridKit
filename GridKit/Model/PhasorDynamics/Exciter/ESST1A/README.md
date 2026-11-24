# ESST1A

Two lead–lag stages drive a limited regulator with rate feedback. `UEL` selects the underexcitation limiter location and `VOS` selects the stabilizer location. The field-current limiter subtracts $K_{lr}\max(I_{fd}-I_{lr},0)$. The final lower and upper field-voltage limits are $V_t V_{rmin}$ and $V_t V_{rmax}-K_c I_{fd}$.

The `Esst1a` device reads voltage magnitude from its `bus` and exciter-base field current from the machine’s `ifd` output. The machine initializes the `efd` output. Optional inputs are `vref`, `vs`, `vuel`, and `voel`.
