# ACTIVSg2000 PostgreSQL

This example runs the reusable
[ACTIVSg2000 case](../../../../cases/PhasorDynamics/ACTIVSg2000/README.md)
through a bus fault and loads its long monitor output into PostgreSQL with
[psql](https://www.postgresql.org/docs/current/app-psql.html).

```bash
DynamicSimulation ACTIVSg2000.solver.json
psql -f ACTIVSg2000.sql
```

The script loads 802,602 rows and lists the generators with the largest speed
deviation:

```text
 class  |      id       | max_speed_deviation
--------+---------------+---------------------
 Genrou | 7357_1_genrou |            0.003144
 Genrou | 7349_1_genrou |            0.002489
 Genrou | 7359_1_genrou |            0.002405
 Genrou | 7360_1_genrou |            0.002402
 Genrou | 7358_1_genrou |            0.002395
(5 rows)
```
