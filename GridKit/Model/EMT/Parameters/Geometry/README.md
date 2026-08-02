# Geometry Models

`Geometry` contains static parameter models that describe the line layout and
conductor data used by overhead-line calculations.

These models do not own EMT network ports. They provide static inputs or
derived geometry quantities to the effect and response models.

Model | Description
----- | -----------
[`Path`](Path/README.md) | Computes path length for finite-line response models
[`Tower`](Tower/README.md) | Stores conductor attachment coordinates and derives span-scale geometry
[`Conductor`](Conductor/README.md) | Stores conductor dimensions, material properties, weight, and phase labels

## Shared Inputs

Quantity | Description
-------- | -----------
`conductors[*].x` | Horizontal conductor attachment coordinates
`conductors[*].height` | Conductor attachment heights above earth
`path.span` | Support-to-support span length
`conductors[*].tension` | Optional per-conductor tension
`conductors[*].radius` | Outer conductor radius
`conductors[*].inner_radius` | Inner conductor radius
`conductors[*].conductivity` | Conductor conductivity
`conductors[*].permeability` | Conductor permeability
`conductors[*].weight` | Conductor weight per unit length
`conductors[*].phase` | Conductor phase label
`path.length`, `path.points` | Explicit length or GIS route data

Quantities are those of the resolved line document: the loader replaces the
tower type and conductor type references with their catalog data, giving one
flat record per physical conductor.

The number of physical conductors $K$ is inferred from the conductor list.
Child models own their detailed validation for geometry, material, and path
inputs.
