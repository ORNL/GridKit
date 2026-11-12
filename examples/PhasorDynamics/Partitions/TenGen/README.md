# TenGen Partitioned

This example runs the ten bus, nine generator system of the `TenGenGenrou`
integration test through `PartitionedDynamicSimulation`. `TenGen.partition.json`
splits it into two partitions that meet at bus 6.

The study runs for one second without events and only checks that the run
completes.
