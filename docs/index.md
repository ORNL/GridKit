(gridkit)=
# GridKit

GridKit is a modeling framework for power-system simulation and analysis. Use
the sections below to install GridKit, run studies, inspect supported models,
or contribute to the project.

::::{grid} 1 2 3 3
:gutter: 2

:::{grid-item-card} Installation
:link: installation
:link-type: ref

Build GridKit and configure its optional dependencies.
:::

:::{grid-item-card} Applications
:link: apps
:link-type: ref

Command line programs that read a case and run a study.
:::

:::{grid-item-card} Models
:link: models
:link-type: ref

Equations, parameters, and ports for every model.
:::

:::{grid-item-card} Cases
:link: cases
:link-type: ref

Machine-readable datasets for GridKit applications.
:::

:::{grid-item-card} Examples
:link: examples
:link-type: ref

Worked examples with solver configurations and validation results.
:::

:::{grid-item-card} API Reference
:link: api
:link-type: ref

Generated reference documentation for the C++ API.
:::

:::{grid-item-card} Development
:link: development
:link-type: ref

Contributor guidance, build system, and documentation build.
:::

::::

```{toctree}
:caption: Start Here
:maxdepth: 1
:titlesonly:
:hidden:

Overview <self>
Installation <installation>
```

```{toctree}
:caption: Using GridKit
:maxdepth: 2
:titlesonly:
:hidden:

Examples <examples/index>
Applications <applications/index>
```

```{toctree}
:caption: Reference
:maxdepth: 4
:titlesonly:
:hidden:

Models <models/index>
Cases <cases/index>
API <api/index>
```

```{toctree}
:caption: Concepts
:maxdepth: 1
:titlesonly:
:hidden:

Architecture <concepts/architecture>
Initialization <concepts/initialization>
Numerical Methods <concepts/numerical-methods>
```

```{toctree}
:maxdepth: 1
:titlesonly:
:hidden:

Development <development/index>
```
