# Line Model

EMT line models represent network connections between buses in instantaneous
phase coordinates.

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
    "R": ...,
    "L": ...,
    "G": ...,
    "C": ...,
    "dx": ...
  }
}
```

### Frequency-Dependent Pi Model

The frequency-dependent pi model keeps the lumped topology and replaces the
constant line matrices with rational parameter models.
Note: using a `VectorFit` model with only $\mathbf{D}$ nonzero is an equivalent
way to define the nominal pi model.

```js
{
  "class": "LineLumped",
  "params": {
    "Z": {
      "class": "VectorFit",
      "params": {
        "D": ...,
        "E": ...,
        "poles": [...],
        "residues": [...]
      }
    },
    "Y": {
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
$\mathbf{y}_c=\mathbf{Y}_0$ and a lossless transport delay
$\mathbf{H}(s)=e^{-s\tau}$.

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

The ULM is the general case where a `VectorFit` model defines the
characteristic admittance $\mathbf{y}_c$ and a `Propagation` model defines the
current-form propagation function $\mathbf{H}_i$.

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
