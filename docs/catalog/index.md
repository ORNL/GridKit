(catalog)=
# Catalog

The catalog names the reusable types that overhead-line input documents
reference: conductors (dimensions, material, and weight) and towers
(attachment geometry of representative structures). Types live in YAML
catalog files that line documents pull in through `include`, so one
library serves many lines and users extend it by adding entries or files.
The [overhead line input format](#models-emt-parameters-schema) documents
the file formats and how references resolve.

The pages below render the shipped catalog,
[`north-american.catalog.yaml`](../../examples/EMT/Lines/north-american.catalog.yaml).
Values are illustrative, not manufacturer data.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} Conductors
:link: catalog-conductors
:link-type: ref

Conductor types: radii, conductivity, permeability, and weight.
:::

:::{grid-item-card} Towers
:link: catalog-towers
:link-type: ref

Tower types: named attachment points with cross-section figures.
:::

::::

```{toctree}
:maxdepth: 1
:titlesonly:
:hidden:

Conductors <conductors>
Towers <towers>
```
