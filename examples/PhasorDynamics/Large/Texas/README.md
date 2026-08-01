# Synthetic Texas Case (ACTIVSg2000)

## Case Description

This case is adapted from the Texas A&M University ACTIVSg2000 synthetic Texas
grid[^activsg2000] and uses GridKit-compatible dynamic models.

## Data Notes

- Generator `6215_1` in source data had `Xl == Xdpp == Xqpp = 0.2011`, which can trigger division-by-zero during GENROU initialization.
- In this case, `Xl` was reduced to `0.1911` to satisfy `Xl < Xdpp`.

## Outstanding

The following source models are not represented directly in this case:

- [ESAC1A](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Exciter%20ESAC1A.htm)
- [ESDC1A](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Exciter%20ESDC1A.htm)
- [ESST4B](../../../../docs/GridKit/Model/PhasorDynamics/Exciter/ESST4B/README.md)
- [EXAC2](../../../../docs/GridKit/Model/PhasorDynamics/Exciter/EXAC2/README.md)
- [EXPIC1](../../../../docs/GridKit/Model/PhasorDynamics/Exciter/EXPIC1/README.md)
- [SCRX](../../../../docs/GridKit/Model/PhasorDynamics/Exciter/SCRX/README.md)
- [GGOV1](../../../../docs/GridKit/Model/PhasorDynamics/Governor/GGOV1/README.md)
- [HYGOV](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Governor%20HYGOV%20and%20HYGOVD.htm)
- [IEEEG1](../../../../docs/GridKit/Model/PhasorDynamics/Governor/IEEEG1/README.md)
- [GENSAL](../../../../docs/GridKit/Model/PhasorDynamics/SynchronousMachine/GENSALwS/README.md)
- [IEEEST](../../../../docs/GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/README.md)

Other outstanding modeling:

- Switched shunts
- LTC transformers and non-unit tap ratios

[^activsg2000]: Texas A&M University,
    “[ACTIVSg2000: 2000-bus synthetic grid on footprint of Texas](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/activsg2000/),”
    *Electric Grid Test Case Repository*.
