(gridkit)=
# GridKit

GridKit is a modeling framework for power-system simulation and analysis. Use
the sections below to install GridKit, run studies, inspect supported models,
or contribute to the project.

::::{grid} 1 2 3 3
:gutter: 2

:::{grid-item-card} Installation
:link: INSTALL
:link-type: doc

Build GridKit and configure its optional dependencies.
:::

:::{grid-item-card} Applications
:link: application/README
:link-type: doc

Runnable applications and end-to-end workflows.
:::

:::{grid-item-card} Models
:link: GridKit/Model/README
:link-type: doc

Equations, parameters, and implementation documentation.
:::

:::{grid-item-card} Cases
:link: cases/index
:link-type: doc

Machine-readable datasets for GridKit applications.
:::

:::{grid-item-card} Examples
:link: examples/README
:link-type: doc

Worked examples with solver configurations and validation results.
:::

:::{grid-item-card} API Reference
:link: api
:link-type: doc

Generated reference documentation for the C++ API.
:::

:::{grid-item-card} Development
:link: development/README
:link-type: doc

Contributor guidance, shared utilities, and build documentation.
:::

::::

```{toctree}
:caption: Start Here
:maxdepth: 1
:titlesonly:
:hidden:

Overview <self>
Installation <INSTALL>
```

```{toctree}
:caption: Using GridKit
:maxdepth: 1
:titlesonly:
:hidden:

Examples <examples/README>
Apps <application/README>
```

```{toctree}
:caption: Reference
:maxdepth: 1
:titlesonly:
:hidden:

Cases <cases/index>
Models <GridKit/Model/README>
API <api>
```

```{toctree}
:maxdepth: 1
:titlesonly:
:hidden:

Development <development/README>
```
