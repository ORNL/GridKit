# Electromagnetic Transients (EMT) Examples

## Provisional Line Model Composition

The pseudo-JSON below records the intended ownership and composition of the
`LineLumped` and `LineDistributed` models. It is not parser-supported and does
not define a `.case.json` schema. Object names, nesting, and value encodings
remain provisional; `...` denotes an unspecified value.

In these sketches, `params` groups parameters owned by the named model, and
`submodels` groups child-model specifications. A named `Yp`, `Yc`, or `H`
specification is reused where the model requires separate terminal or
directional instances.

### LineLumped

`LineLumped` owns `N`, `K`, `conductors`, and `dx`. Its `Zp` and `Yp`
specifications configure the series-impedance and shunt-admittance `VectorFit`
submodels. The pole-free RLGC specialization uses empty `poles` and `residues`;
frequency-dependent forms supply fitted values. The `D` and `E` coefficients
remain submodel parameters in both cases.

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
`Propagation` model. `Propagation` owns `K`, `tau`, and `fmax`; its `input` and
`output` specifications configure the two `VectorFit` factors. The `input`
factor's `E` coefficient is zero in every realization.

The Bergeron specialization uses pole-free `Yc`, `input`, and `output` models
with zero `E` matrices. The `Yc` `D` matrix provides the constant
characteristic admittance, while the `input` and `output` `D` matrices provide
the constant delay-free modal factors. The `tau` values supply the modal
delays. The universal line model supplies fitted poles and residues for `Yc`,
`input`, and `output`.

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
