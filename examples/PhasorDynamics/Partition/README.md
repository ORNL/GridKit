# Partition Schema & Rules

## TYPES

A `Partition` has the following properties:
- A `name` of type `str` identifying the partition
- A `external` of type `List<int>` containing the bus numbers external to the partition.
- A `internal` of type `List<int>` containing the bus numbers internal to the partition.

A `PartitionGroup` has the following properties:
- A `name` of type `str` identifying the partition group.
- A `partitions` list of type `List<Partition|PartitionGroup>` containing any combination of partitions and partition groups.

## RULES


- The root object is a `PartitionGroup`.
- Internal buses are unique to one `Partition`.
- Every bus must be exlusively `external` or `internal`.
- A `Component`s `bus` ports must all reside in the same partition

## PROPERTIES

- A `PartitionGroup` can always be written as a flat set of `Partitions`


## Regarding Components and SignalNodes

Until we support `internal`/`external` signals, we should enforce:
- If a `Component` is connected to an `external` bus, it will be assigned to the first partition where that bus appears.
- `SignalNodes` must only connect `Component`s of the same `Partition`
