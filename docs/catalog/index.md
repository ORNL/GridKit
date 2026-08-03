(catalog)=
# Catalog

Catalogs name the reusable conductor and tower types that
[overhead-line input documents](#models-emt-parameters-schema) reference,
in YAML files that stay readable, commentable, and easy to extend. The
pages below render the shipped catalog,
[`north-american.catalog.yaml`](../../examples/EMT/Lines/north-american.catalog.yaml);
values are illustrative, not manufacturer data.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} Conductors
:link: catalog-conductors
:link-type: ref

Dimensions, material, and weight of each conductor type.
:::

:::{grid-item-card} Towers
:link: catalog-towers
:link-type: ref

Named attachment points of each tower type, with cross-sections.
:::

::::

```{toctree}
:maxdepth: 1
:titlesonly:
:hidden:

Conductors <conductors>
Towers <towers>
```
