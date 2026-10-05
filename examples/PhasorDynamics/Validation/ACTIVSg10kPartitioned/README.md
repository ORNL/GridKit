# Partitioned ACTIVSg10k

This study uses eight regional IDA solvers with multicolor Gauss-Seidel
coupling. The eight regions form three colors of three, two, and three regions,
so up to three OpenMP workers advance regions at once. It includes the same
fault at 1 s and clearing at 1.15 s as the monolithic validation study. Build
with `GridKit_ENABLE_OPENMP=ON`.

From the repository root, build and run the monolithic reference followed by
the partitioned study:

```bash
cmake --build build --target DynamicSimulation -j 10
(cd build/examples/PhasorDynamics/Validation/ACTIVSg10k && ../../../../application/PhasorDynamics/DynamicSimulation ACTIVSg10k.solver.json)
(cd build/examples/PhasorDynamics/Validation/ACTIVSg10kPartitioned && ../../../../application/PhasorDynamics/DynamicSimulation ACTIVSg10kPartitioned.solver.json)
```

Set `partition.threads` to 1 to run the same method without concurrency; the
results are identical. The four- and sixteen-region partition files are also
copied here; select one with `partition.file`.

`partition.tol` controls boundary-input mismatch, not trajectory error. The
application reports aggregate speed error against the monolithic output;
assess individual channels and refine the coupling tolerance for accuracy
studies. This full-size study is run manually because of its cost.
