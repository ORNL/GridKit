(catalog-conductors)=
# Conductors

Conductor types carry the dimensions, material, and weight of one physical
conductor. Phase, circuit, tension, and the attachment point are line
data, so the same type serves phase conductors, neutrals, and shield wires
alike. Permeability is relative to vacuum and defaults to 1; the inner
radius is omitted for solid conductors.

```{catalog-conductors} examples/EMT/Lines/north-american.catalog.yaml
```

To add a type, append an entry to a catalog file or define it in the
document's own `catalog` section; the
[input format](#models-emt-parameters-schema) describes the fields and the
validation applied at load.
