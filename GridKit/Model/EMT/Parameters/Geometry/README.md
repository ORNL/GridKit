# Geometry Models

`Geometry` contains static parameter models that describe the line layout and
conductor data used by overhead-line calculations.

These models do not own EMT network ports. They provide static inputs or
derived geometry quantities to the effect and response models.

Model | Description
----- | -----------
[`Path`](Path/README.md) | Computes path length for finite-line response models
[`Tower`](Tower/README.md) | Stores conductor attachment coordinates and span-scale geometry
[`Conductor`](Conductor/README.md) | Stores conductor dimensions, material properties, weight, and phase labels

## Shared Inputs

Quantity | Description
-------- | -----------
`tower.x` | Horizontal conductor attachment coordinates
`tower.height` | Conductor attachment heights above earth
`tower.span` | Support-to-support span length
`tower.tension` | Optional span tension parameter
`conductors[*].radius` | Outer conductor radius
`conductors[*].inner_radius` | Inner conductor radius
`conductors[*].conductivity` | Conductor conductivity
`conductors[*].permeability` | Conductor permeability
`conductors[*].weight` | Conductor weight per unit length
`conductors[*].phase` | Conductor phase label
`length`, `path` | Explicit length or GIS path data

The number of physical conductors $K$ is inferred from the conductor list.
Child models own their detailed validation for geometry, material, and path
inputs.
