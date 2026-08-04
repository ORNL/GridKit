# TODO

Earlier findings now resolved and removed from this list: the per-mode
eigen-tracking DAE (replaced by pointwise modal decomposition with
alignment), the scalar Yc/Zc identities (sandwich residuals landed),
missing bundle/shield reduction (K is now the phase count), loose
per-variable absolute tolerances, and silent nonpassive artifacts
(exit code 3 plus the verdict inside yc.model.json).

- `SkinEffect` ignores `inner_radius`: the schema and parser carry it,
  but the conductor model uses only the outer radius, so annular
  catalog conductors (Drake, Cardinal) sit roughly 10-15% low on DC
  resistance. The skin-effect tests force the inner radius to zero.

- `ShuntPotential` guards its capacitance-inversion pivots with assert
  only, so release builds can divide by zero; `Tower` accepts
  overlapping conductors and conductors at or below ground, which feed
  log terms invalid arguments.

- Passivity screening samples a finite grid (`Passivity.cpp`), so a
  violation between screen points can slip through.

- Nearly degenerate aerial modes make individual mode dyads noisy even
  though the noise cancels in the modal sum (69kv-wood-pole worst mode
  2.5e-1 against a 5.6e-2 sum; 345kv-staggered mode 1 fits |D| = 11
  with a sub-ns transient). Fitting delay groups over degenerate
  clusters instead of single modes would remove this.
