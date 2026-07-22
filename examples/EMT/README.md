# Electromagnetic Transients (EMT) Examples

This directory preserves the intended model composition for future EMT
examples. The trees describe ownership and nesting only. They do not define a
case-file schema.

> [!NOTE]
> The formulations support $N$ phases; initial development targets three
> phases.

## LineLumped

```text
LineLumped
  parameters: N, K, conductors, dx
  initial state: i12
  submodels
    Zp: VectorFit
    Yp[1]: VectorFit
    Yp[2]: VectorFit
```

Zp represents per-unit-length series impedance. The two Yp instances represent
terminal shunt admittances. Pole-free coefficients give the RLGC form, while
poles and residues give the frequency-dependent form.

## LineDistributed

```text
LineDistributed
  parameters: N, K, conductors
  submodels
    Yc[1]: VectorFit
    Yc[2]: VectorFit
    H21: Propagation
      input: VectorFit
      delays[M]: Delay
      output: VectorFit
    H12: Propagation
      input: VectorFit
      delays[M]: Delay
      output: VectorFit
```

Yc represents characteristic admittance. H21 and H12 represent directional
propagation. The rational factors have no term linear in $s$. Pole-free
factors give the Bergeron form, while fitted poles and residues give the
frequency-dependent form.
