# Coupled overhead-line parameters

`model.json` supplies full, reciprocal 3-by-3 phase-domain `Zp` and `Yp`
coefficient sets to the case's `LineLumped` devices. Coefficients use meters:
`Zp` is impedance per meter, `Yp` is admittance per meter, and each device's
`dx` is its length in meters. Both terminal shunts receive the same `Yp`
coefficients; `LineLumped` applies the half-length factor.

## Physical inputs

The example uses a synthetic 13.8 kV, untransposed, three-wire overhead feeder
with catalog conductors and explicitly defined geometry. These are calculated
engineering parameters, not measurements from a surveyed feeder.

| Input | Value |
| --- | --- |
| Conductor | ACSR Linnet, 336.4 kcmil, 26 aluminum / 7 steel strands |
| Outside / core diameter | 18.288 / 6.73608 mm |
| Catalog GMR | 7.415784 mm |
| Catalog DC resistance at 20 C | 0.1660105 ohm/km |
| Phase A coordinates `(x, height)` | `(-1.2, 10.5)` m |
| Phase B coordinates | `(0.0, 11.0)` m |
| Phase C coordinates | `(1.2, 10.5)` m |
| Conductor temperature | 20 C |
| Homogeneous earth resistivity | 100 ohm m |
| Neutral / shield wires | None |
| Transposition | None |

The conductor catalog and parameter calculator are from
[OpenLine revision `56c405add9f082fa0816338978f25a8242c2cc15`](https://github.com/lukelowry/openline/tree/56c405add9f082fa0816338978f25a8242c2cc15).
The bundled Linnet entry attributes its dimensions and resistance to the CME
Wire BareNRG ACSR catalog, 336.4 kcmil 26/7 row. The calculation uses OpenLine's
Schelkunoff conductor internal impedance, perfect-earth external inductance,
Carson earth-return impedance, and Maxwell potential-coefficient capacitance.
Catalog resistance already includes stranding, so no extra stranding multiplier
is added. The 20 C condition matches the catalog's DC resistance measurement;
this OpenLine revision clamps that measurement when queried outside its
temperature range.

All nine matrix entries are retained. For example, near 60 Hz the A-B mutual
impedance is approximately `0.058 + j0.49 ohm/km`, compared with approximately
`0.23 + j0.90 ohm/km` for phase A self impedance. The shunt capacitance matrix is

```text
C [nF/km] = [[ 8.525007, -2.542233, -1.474296],
             [-2.542233,  8.978749, -2.542233],
             [-1.474296, -2.542233,  8.525007]]
```

## Rational realization and validation

`openline_sweep.rs` calculates 801 logarithmically spaced frequencies from
1 Hz to 30 kHz. `generate.py` fits eight shared, negative real poles using real
least squares with pole optimization and relative weighting of every matrix
entry. Even-index samples train the fit; the 400 odd-index samples are held out.
This is a stable real rational fit to OpenLine data; the local `vecfit` checkout
is not required.

The fitted impedance first uses the form

```math
Z'(s)=R_0+sL_\infty+\sum_k A_k\frac{s}{s+a_k},\qquad a_k>0.
```

Small negative eigenvalues in the fitted coefficient matrices are projected to
zero. Consequently `R0`, `Linf`, and every `Ak` are symmetric positive
semidefinite. Each term is positive real for `Re(s) > 0`, giving a passive
matrix impedance. The exported `VectorFit` coefficients are
`D = R0 + sum(Ak)`, `E = Linf`, poles `-ak`, and residues `-ak*Ak`.
The full Maxwell shunt is represented exactly by `Yp(s) = s*C`.

Validation in `fit_validation.json` includes:

- Maximum held-out entrywise relative impedance error: **0.22920%**;
  RMS entrywise error: **0.08875%**.
- Maximum held-out relative Frobenius matrix error: **0.20638%**.
- Shunt relative error below `8e-13`.
- Independent capacitance calculation by inversion of the electrostatic
  image-method potential matrix: maximum relative difference below `5e-10`.
- Independent evaluation of the exported pole-residue model against its
  passive realization.
- Numerical passivity check at 20,001 frequencies from `1e-4` to `1e9` Hz:
  minimum eigenvalue of `Re(Z')` is `1.66033e-4 ohm/m`.

The wide passivity check establishes numerical behavior of the rational model;
physical fit accuracy is assessed only over 1 Hz to 30 kHz. A lumped pi section
also has spatial discretization error, separately from parameter-fit error.
The 300-600 m sections used by the case limit this error near the 900 Hz carrier
and its lower harmonics. They are not a substitute for a distributed traveling
wave model at arbitrary frequency. Air-dielectric loss and corona are omitted.
An independent comparison with the exact uniform-line telegrapher transfer
matrix gives relative 6-by-6 terminal-admittance errors for the longest 600 m
section of **0.00234% at 900 Hz**, **0.0207% at 2.7 kHz**, and **2.53% at 30 kHz**.
Both sides use identical per-meter parameters; these errors isolate the pi
approximation from the rational parameter fit.

## Reproduction

The stored raw sweep allows refitting with Python 3, NumPy, and SciPy:

```bash
python3 cases/EMT/CoupledGrid/lines/generate.py
```

To recalculate physical parameters using a clean OpenLine checkout at the
recorded revision, with Rust/Cargo available:

```bash
python3 cases/EMT/CoupledGrid/lines/generate.py --recompute --openline /home/lukel/openline
```

The driver is built in `/tmp` with the supplied `Cargo.lock`; source repositories
are read only. `CARGO_HOME` and `CARGO_TARGET_DIR` may override the default
temporary dependency/build directories. The validation JSON records source
revision and SHA-256 hashes of the raw sweep and driver.

`openline.csv` is the calculator's original long-form output in ohm/km and S/km.
`response.csv` has the raw and fitted values in ohm/m and S/m for plots:
`freq_hz`, followed by `raw_Z00_re`, `raw_Z00_im`, `fit_Z00_re`, and corresponding
fields for all matrix entries and for `Y`. The header is authoritative for
column order.
