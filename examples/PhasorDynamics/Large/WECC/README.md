# Western Electricity Coordinating Council (WECC)

## One-Line Diagram

![](oneline.jpg)

Figure 1: One-line diagram of the WECC case. Source: [National Laboratory of
the Rockies](https://www.nlr.gov/grid/test-case-repository).

## Case Description

This case is adapted from the [National Laboratory of the
Rockies](https://www.nlr.gov/grid/test-case-repository) WECC model and uses
GridKit-compatible dynamic models.

## Outstanding

The following source models are not represented directly in this case:

- [REECB1](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Exciter%20REEC_B.htm)
- [GAST_PTI](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Governor%20GAST_PTI%20and%20GASTD.htm)
- [REGCA](../../../../GridKit/Model/PhasorDynamics/Converter/REGCA/README.md)
- [REPCA1](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Plant%20Controller%20REPC_A.htm)

Other outstanding modeling:

- Switched shunts
- LTC transformers and non-unit tap ratios
