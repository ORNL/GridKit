# Electromagnetic Transients (EMT) Examples

> [!NOTE]
> The formulations support $N$ phases; initial development targets three
> phases.

| Example | Description |
| --- | --- |
| [ThreeBusBasic](Toy/ThreeBusBasic/README.md) | Line energization. |

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
