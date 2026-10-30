# ACTIVSg2000 PostgreSQL

This example runs the reusable
[ACTIVSg2000 case](../../../../cases/PhasorDynamics/ACTIVSg2000/README.md)
through a bus fault with `"output_layout": "long"` and loads the output into
PostgreSQL with [psql](https://www.postgresql.org/docs/current/app-psql.html).
The [Docker Compose](https://docs.docker.com/compose/) file runs both in
containers. Run from the example's build or install directory (see
[Running a Study](../../README.md#running-a-study)).

```bash
DynamicSimulation ACTIVSg2000.solver.json
docker compose run --rm psql
docker compose down
```

With Podman, use `podman compose` (or `podman-compose`) in place of
`docker compose`.

`ACTIVSg2000.sql` loads 802,602 rows and lists the generators with the largest
speed deviation:

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
