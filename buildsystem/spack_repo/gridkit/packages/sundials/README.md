# IDAS step validation

`idas-step-check.patch` adds a pre-commit accuracy check to SUNDIALS 7.8.0 IDAS.
The recipe applies it automatically. This is a local extension, not an upstream
SUNDIALS API. It contains no GridKit or delay-specific logic.

```c
typedef int (*IDAStepCheckFn)(sunrealtype t0, sunrealtype t1,
                            N_Vector y, N_Vector yp,
                            sunrealtype* hnext, void* user_data);
int IDASetStepCheckFn(void* mem, IDAStepCheckFn check, void* user_data);
```

The callback runs after all native error tests, before `IDACompleteStep`:

- Return zero to accept. A positive `*hnext` caps the next step magnitude;
  zero leaves it to IDAS. The cap is clipped to the configured minimum step.
- Return positive to retry with `0 < *hnext < abs(t1-t0)`. IDAS restores the
  trial state, preserves multistep history and order, and applies its minimum
  step. Repeated rejection is limited by `IDASetMaxErrTestFails`.
- Return negative to stop. Callback errors, invalid caps, an exhausted retry
  budget, or accuracy failure at the floor return `IDA_STEP_CHECK_FAIL`.
- `IDAGetDky` inside the callback evaluates the candidate polynomial on
  `[t0,t1]`. Outside it, the existing accepted-output behavior is unchanged.
  Do not modify the candidate vectors or solver options inside the callback.

Before commitment, IDAS holds scaled divided differences and a correction `ee`.
The extension evaluates the coefficients that commitment would produce,
without copying or modifying the history vectors:

```text
phi_new[j] = sum(phi_star[j..q]) + ee
p^(k)(t)   = sum_j c_j^(k)(t) phi_new[j]
```

The existing dense-output basis recurrence is shared by both paths. No history
snapshots, reinitializations, or additional solution vectors are required.
For fixed-leading-coefficient BDF, the dense polynomial's derivative need not
be identical to the nonlinear solve's `yp`; trial and committed dense output
use exactly the same polynomial.

The patch also makes a changed `IDASetMaxStep` constrain the already selected
next step, and stops redundant native error-test retries once a rejected step is at
`h_min` and order one. Higher orders can still fall back before failure.
Stop-time clipping retains standard IDAS behavior.

GridKit's direct IDAS regression compares trial and accepted polynomials and
all available derivatives, forces a higher-order retry, and lowers an existing
step cap. Adapter tests cover minimum-step failure and C++ exception propagation.
Signal-delay integration tests exercise feedback, continuation and long delays.

For a manual SUNDIALS build, apply the patch in its source directory, rebuild and
install IDAS, then configure GridKit with `SUNDIALS_DIR` pointing to that install's
`lib/cmake/sundials`. Ordinary simulations still build with unpatched SUNDIALS;
signal delays report the missing extension explicitly and their integration test
is marked skipped.
