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
- Assign each `Component` with bus ports to the first `Partition` containing all its bus ports. Reject partition if none exists. 
- If a `Component` assigned to a partition is connected to an unassigned `Component`, assign to the same partition. Repeat until no further assignments occur, then reject the partition if a `Component` is left unassigned or conflicts occur.
