# Partitions in Phasor Dynamics

## Overview and Purpose

The concept of "Partition" is used in phasor dynamics in GridKit to refer
to dividing up a system model into smaller parts
**for the purpose of parallel simulation** and other
solver acceleration techniques. 
Dividing up a system for _modeling convenience_ (like having physically
meaningful areas corresponding to particular utilities, 
containerized subcomponents, template aggregate models
for multi-instantiation, etc.) 
is a fundamentally different goal. It is not
appropriate to use divisions driven by modeling for solver acceleration.

Partitions for a system are defined in a JSON file, usually called
`[model_name].partition.json`, which is separate from the system model 
definition
but refers to buses by number from the system model definition.
Defining a partition file for a system is optional, but needed to enable
parallel model evaluation for some solvers.

This file defines a nested, hierarchical structure of containers
of type `Partition` and `PartitionGroup`.
The format is general and could accommodate any combination of nested
partitions, including the common structure of a single level of partitions
dividing the system for a single level of parallel evaluation.
Any more complex structure can also be easily flattened 
into a single level by solvers that require only a single level of partitioning.
The JSON file schema is defined in `partition.schema.json`, and is also
described below.

In the future it is envisioned that users would not need to think about
partitioning a system. Rather, an automated partitioning routine 
would pick
the best partitioning for maximum solver efficiency on a given computational
environment,
and perhaps write the `partition.json` file as an intermediate step. 
For now, partitions need to be defined by the user-written JSON file.

## Partition Definitions and JSON Structure

As mentioned above, there are two types of objects in a `partition.json` file:
partitions and partition groups.

1. A `Partition` is a container for buses that is a leaf node of the partition
hierarchy. Each partition has the following properties:
    - `name` of type `string` identifying the partition
    - `external` of type `array of integers` containing the 
bus numbers external to the partition. These bus numbers must correspond to 
buses in the model JSON file. Any bus connecting to this partition and also
to other partitions is considered external.
    - `internal` of type `array of integers` containing the bus numbers internal 
to the partition. These bus numbers must correspond to 
buses in the model JSON file. Any bus completely unique to this partition
should go in this list.

2. A `PartitionGroup` is a container for partitions and other 
partition groups. It functions as an intermediate node 
in the partition hierarchy.
Partition groups do not contain any buses directly.
The root-level object of the JSON file is a partition group.
Each partition group has the following properties:
    - `name` of type `string` identifying the partition group.
    - `partitions` of type `array of (Partition or PartitionGroup)` 
    containing any combination of partitions and partition groups.

## Validation Rules and Component Assignment

Any violation of the rules below results in a `BadPartition` exception.
Component assignment is done automatically. This avoids the error-prone 
requirement for the user to manually assign both buses and components and
ensure they are compatible.

1. The file must follow the `partition.schema.json` structure.
2. The root object must be a `PartitionGroup`.
3. Every bus must appear somewhere,
either as `external` or `internal`, but not both.
4. Internal buses must only appear uniquely in one `Partition`.
5. Each `Component` that has bus ports shall be assigned to the first 
`Partition` containing all its bus ports. If none exists,
raise a `BadPartition` exception. 
6. Some components are connected to each other by signals.
If a `Component` assigned to a partition is 
connected to an unassigned `Component`, 
assign to the same partition. This process needs to be repeated until no further
assignments occur.
If a `Component` is left unassigned or if two connected
components are assigned to different partitions,
raise a `BadPartition` exception.

In the future it is envisioned that signals may also be part of the 
partition definition, allowing connected components to be in separate 
partitions. For now, only buses can divide partitions.

## Examples

### Ten Generator Case - Divide into Two Partitions

File name: `TenGen.partition.json`

```JSON
{
  "name": "TenGen",
  "partitions": [
    {
      "name": "A",
      "internal": [
        1,
        2,
        3,
        4,
        5
      ],
      "external": [
        6
      ]
    },
    {
      "name": "B",
      "internal": [
        7,
        8,
        9,
        10
      ],
      "external": [
        6
      ]
    }
  ]
}
```

![TenGen partition](TenGen.png)

*Figure 1: TenGen split into two partitions, A and B. Buses 1–5 are internal to A and buses 7–10 are internal to B (red). Bus 6 is external for both (blue); the dashed line indicates the partitions couple at Bus 6.*

### WECC 240 Case - Divide into Four Partitions

File name: `WECC240.partition.json`

(This file can be found under `examples/PhasorDynamics/Partitions/WECC240`.)

![WECC240 partition](WECC240.png)

*Figure 2: WECC240 split into four regional partitions: Northwest, Southwest, NorCal and SoCal. Internal buses (red) belong to exactly one partition. Boundary buses (blue) are external and appear in every partition that touches them; dashed lines show where
buses couple between partitions.*


### Sample Multi-Level Partition

This example shows how multi-level hierarchical partitions can be implemented
with partition groups.

File name: `deep.partition.json`

```JSON
{
  "name": "root",
  "partitions": [
    {
      "name": "A",
      "internal": [
        1
      ],
      "external": [
        5
      ]
    },
    {
      "name": "B",
      "partitions": [
        {
          "name": "B1",
          "internal": [
            2
          ],
          "external": [
            5,
            6
          ]
        },
        {
          "name": "C",
          "partitions": [
            {
              "name": "C1",
              "internal": [
                3
              ],
              "external": [
                6,
                7
              ]
            },
            {
              "name": "C2",
              "internal": [
                4
              ],
              "external": [
                7
              ]
            }
          ]
        }
      ]
    }
  ]
}
```

![Sample multi-Level partition](deep.png)

*Figure 3: Sample multi-level partition scheme. The partitions are
A, B1, C1, and C2, whereas root, B and C are partition groups.*
