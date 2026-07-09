# Line Models

EMT line models represent $N$-phase network connections between buses in
instantaneous phase coordinates.

## Types

Lumped transmission line models approximate the line with finite network
elements, sometimes referred to as the $\pi$-model.

Distributed transmission line models preserve propagation and delay.
The formulation allows this model to be used for either constant or
frequency-dependent line parameters.

- `LineLumped` (See [LineLumped](LineLumped/README.md))
- `LineDistributed` (See [LineDistributed](LineDistributed/README.md))

## Proposed Model Specifications

The snippets below define the proposed model-composition shape. They are not
yet parser-supported `.case.json` entries. Fields under `params` belong to the
named model, while `submodels` contains child-model specifications. A named
`Yp`, `Yc`, or `H` specification is reused when the detailed model requires
separate terminal or directional instances.

### LineLumped

`LineLumped` owns `N`, `K`, `conductors`, and `dx`. Its `Zp` and `Yp`
specifications configure the series-impedance and shunt-admittance `VectorFit`
submodels. The nominal pi specialization uses empty `poles` and `residues`.
The frequency-dependent specialization supplies fitted poles and residues.
The `D` and `E` coefficients remain owned by the submodels in both cases.

```js
{
  "class": "LineLumped",
  "params": {
    "N": ...,
    "K": ...,
    "conductors": [...],
    "dx": ...
  },
  "submodels": {
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
    }
  }
}
```

### LineDistributed

`LineDistributed` owns `N`, `K`, and `conductors`. Its `Yc` specification
configures the characteristic-admittance `VectorFit`, while `H` configures the
`Propagation` model. `Propagation` owns `K`, `tau`, and `fmax`, and its `input`
and `output` specifications configure its two `VectorFit` factors.

The Bergeron specialization uses pole-free `Yc`, `input`, and `output` models
with zero `E` matrices. The `Yc` `D` matrix provides the constant
characteristic admittance, while the `input` and `output` `D` matrices provide
the constant delay-free modal factors. The `tau` values exclusively supply the
modal delays. The universal line model supplies fitted poles and residues for
`Yc`, `input`, and `output`.

```js
{
  "class": "LineDistributed",
  "params": {
    "N": ...,
    "K": ...,
    "conductors": [...]
  },
  "submodels": {
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
        "K": ...,
        "tau": [...],
        "fmax": ...
      },
      "submodels": {
        "input": {
          "class": "VectorFit",
          "params": {
            "D": ...,
            "E": ...,
            "poles": [...],
            "residues": [...]
          }
        },
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
