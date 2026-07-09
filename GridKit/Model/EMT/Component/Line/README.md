# Line Model

EMT line models represent $N$-phase network connections between buses in
instantaneous phase coordinates.

## Types

Lumped transmission line models approximate the line with finite network
elements, sometimes referred to as the $\pi$-model.

Distributed transmission line models preserve propagation and delay.
The formulation allows this model to be used for either constant or
frequency-dependent line parameters.

- `LineLumped` (See [documentation](LineLumped/README.md))
- `LineDistributed` (See [documentation](LineDistributed/README.md))

## Examples

The snippets below show model composition only. Final `.case.json` syntax may
change when the EMT line models are implemented.

### Nominal Pi Model

The nominal pi model is the constant-parameter lumped case.

```js
{
  "class": "LineLumped",
  "params": {
    "Rp": ...,
    "Lp": ...,
    "Gp": ...,
    "Cp": ...,
    "dx": ...
  }
}
```

### Frequency-Dependent Pi Model

The frequency-dependent pi model keeps the lumped topology and replaces the
constant line matrices with rational parameter models.
Note: using `VectorFit` models for `Zp` and `Yp` with no poles or residues
is an equivalent way to define the nominal pi model: for `Zp`,
$\mathbf{D}=\mathbf{R}'$ and $\mathbf{E}=\mathbf{L}'$. For `Yp`,
$\mathbf{D}=\mathbf{G}'$ and $\mathbf{E}=\mathbf{C}'$.

```js
{
  "class": "LineLumped",
  "params": {
    "Zp": {
      "class": "VectorFit",
      "params": {
        "D": ...,
        "E": ...,
        "poles": [...],
        "residues": [...]
      }
    },
    "Yp": {
      "class": "VectorFit",
      "params": {
        "D": ...,
        "E": ...,
        "poles": [...],
        "residues": [...]
      }
    },
    "dx": ...
  }
}
```

### Bergeron Model

The constant parameter case uses constant characteristic admittance
$f^{\mathbf{y}_c}=\mathbf{Y}_0$ and a lossless transport delay
$f^{\mathbf{h}}(s)=e^{-s\tau}$.

```js
{
  "class": "LineDistributed",
  "params": {
    "conductors": [...],
    "Yc": ...,
    "H": {
      "class": "Delay",
      "params": {
        "tau": ...,
        "fmax": ...
      }
    }
  }
}
```

### Universal Line Model

The ULM is the general case with characteristic admittance `Yc` and
current-form propagation function `H`. In equations these operators are written
as $f^{\mathbf{y}_c}(\cdot)$ and $f^{\mathbf{h}}(\cdot)$.

```js
{
  "class": "LineDistributed",
  "params": {
    "conductors": [...],
    "Yc": {
      "class": "VectorFit",
      "params": {
        "D": ...,
        "E": ...,
        "poles": [...],
        "residues": [...]
      }
    },
    "H": {
      "class": "Propagation",
      "params": {
        "input": {
          "class": "VectorFit",
          "params": {
            "D": ...,
            "E": ...,
            "poles": [...],
            "residues": [...]
          }
        },
        "tau": [...],
        "fmax": ...,
        "output": {
          "class": "VectorFit",
          "params": {
            "D": ...,
            "E": ...,
            "poles": [...],
            "residues": [...]
          }
        }
      }
    }
  }
}
```
