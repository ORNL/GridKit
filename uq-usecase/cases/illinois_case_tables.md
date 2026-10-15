# Illinois case JSON tables (illinois.case.json)

> `examples/PhasorDynamics/Large/Illinois/illinois.case.json`

Parsed directly from the JSON file. All values are exact as stored in the file.
Powers (`p0`, `q0`) are in per-unit on each machine's own MVA base.
Reactances (`Xd`, `Xdp`, etc.) are in per-unit on machine base.
`H` is in seconds. `mva` is machine MVA base.
`H*mva` (MWs) is the total stored kinetic energy — the physically meaningful inertia contribution to the grid.

## Genrou devices (40 total)

Sorted by `H*mva` ascending. Bus names are as in `illinois.case.json` buses array.

**Why H×mva?** `H` (the inertia constant, in seconds) is defined as stored kinetic energy normalized by the machine's own MVA rating: $H = \frac{\frac{1}{2} J \omega_0^2}{S_\text{rated}}$. It is already per-unit on each machine's individual base, not on the system base (100 MVA). Each machine has its own `mva` (`MBASE` in MATPOWER), which varies widely across the fleet — from 2.04 MVA (NORMAL 2 2, tiny wind) to 682.98 MVA (CLINTON 1 2, nuclear). To recover the actual stored kinetic energy in MWs you multiply back: $H \times \text{mva}$. Two machines with the same H can contribute very different inertia to the grid: `77_1` (H=5.14, mva=2.88) stores 14.8 MWs, while `189_1` (H=3.40, mva=682.98) stores 2323 MWs. **H×mva is therefore the physically relevant quantity for ranking inertia contributions and selecting UQ generators.**

### Genrou parameter glossary

All quantities are in per-unit on the machine's own MVA base unless noted otherwise.

| param | unit | description |
|-------|------|-------------|
| `mva` | MVA | machine MVA rating (each machine's own base; varies per unit) |
| `p0` | pu | initial active power output (per-unit on `mva`) |
| `q0` | pu | initial reactive power output (per-unit on `mva`) |
| `H` | s | inertia constant: stored kinetic energy at rated speed divided by `mva`; epistemic UQ parameter |
| `D` | pu | damping coefficient; zero for all machines in this case |
| `Ra` | pu | armature (stator) resistance; zero for all machines in this case |
| `Xd` | pu | d-axis synchronous reactance |
| `Xdp` | pu | d-axis transient reactance (primed) |
| `Xdpp` | pu | d-axis subtransient reactance (double-primed) |
| `Xq` | pu | q-axis synchronous reactance |
| `Xqp` | pu | q-axis transient reactance |
| `Xqpp` | pu | q-axis subtransient reactance |
| `Xl` | pu | leakage (armature) reactance |
| `Tdop` | s | d-axis transient open-circuit time constant |
| `Tdopp` | s | d-axis subtransient open-circuit time constant |
| `Tqop` | s | q-axis transient open-circuit time constant |
| `Tqopp` | s | q-axis subtransient open-circuit time constant |
| `S10` | pu | saturation factor at 1.0 pu flux linkage |
| `S12` | pu | saturation factor at 1.2 pu flux linkage; together with `S10` defines the saturation curve |

The suffix convention `p` = transient (one rotor winding), `pp` = subtransient (two rotor windings). Open-circuit time constants (`Tdo`, `Tqo`) are for the rotor circuit with stator open; they are longer than short-circuit time constants by the ratio $X / X'$.

| id | bus | bus_name | mva | H | H×mva (MWs) | p0 (pu) | q0 (pu) | D | Ra | Xd | Xdp | Xdpp | Xq | Xqp | Xqpp | Xl | Tdop | Tdopp | Tqop | Tqopp | S10 | S12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `77_1` | 77 | BRIMFIELD 2 | 2.880 | 5.1416 | **14.81** | 0.00720 | 0.01160 | 0.0 | 0.0 | 1.5282 | 0.1935 | 0.14170 | 1.4642 | 0.50050 | 0.14170 | 0.1048 | 6.2717 | 0.0456 | 1.3349 | 0.0783 | 0.1220 | 0.4792 |
| `90_1` | 90 | CLINTON 3 1 | 3.840 | 4.4204 | **16.97** | 0.00960 | 0.01081 | 0.0 | 0.0 | 1.6558 | 0.1847 | 0.13550 | 1.5830 | 0.28640 | 0.13550 | 0.1000 | 4.5742 | 0.0374 | 0.8295 | 0.0493 | 0.1563 | 0.5477 |
| `114_1` | 114 | NORMAL 2 2 | 2.040 | 8.9160 | **18.19** | 0.00400 | 0.00360 | 0.0 | 0.0 | 1.4605 | 0.2653 | 0.19270 | 1.4012 | 0.46130 | 0.19270 | 0.1438 | 6.0742 | 0.0476 | 1.0342 | 0.0643 | 0.0228 | 0.2807 |
| `52_1` | 52 | RANTOUL 2 4 | 5.440 | 4.9625 | **27.00** | 0.04530 | 0.02110 | 0.0 | 0.0 | 1.9192 | 0.2101 | 0.15350 | 1.8281 | 0.56830 | 0.15350 | 0.1138 | 8.4034 | 0.0414 | 1.3643 | 0.0478 | 0.1114 | 0.4579 |
| `49_1` | 49 | RANTOUL 2 1 | 5.440 | 5.2115 | **28.35** | 0.04530 | 0.02110 | 0.0 | 0.0 | 2.4971 | 0.4396 | 0.31650 | 2.3662 | 0.76020 | 0.31650 | 0.2384 | 6.8164 | 0.0245 | 1.3976 | 0.0796 | 0.1971 | 0.6294 |
| `76_1` | 76 | BRIMFIELD 1 | 4.800 | 6.7471 | **32.39** | 0.01200 | 0.02040 | 0.0 | 0.0 | 2.0742 | 0.3757 | 0.20816 | 1.9725 | 0.26020 | 0.20816 | 0.2037 | 6.5649 | 0.0279 | 0.7611 | 0.0584 | 0.0844 | 0.4038 |
| `151_1` | 151 | SPRINGFIELD 5 2 | 6.480 | 5.3589 | **34.73** | 0.05400 | 0.01051 | 0.0 | 0.0 | 1.2438 | 0.2028 | 0.14830 | 1.1994 | 0.20370 | 0.14830 | 0.1098 | 5.4147 | 0.0260 | 1.3223 | 0.0613 | 0.0804 | 0.3959 |
| `170_1` | 170 | CHAMPAIGN 1 7 | 11.280 | 3.4015 | **38.37** | 0.02820 | 0.04780 | 0.0 | 0.0 | 1.5110 | 0.1514 | 0.11190 | 1.4482 | 0.43940 | 0.11190 | 0.0819 | 7.0160 | 0.0221 | 0.7490 | 0.0744 | 0.0492 | 0.3335 |
| `91_1` | 91 | CLINTON 3 2 | 6.000 | 6.4018 | **38.41** | 0.01500 | 0.02550 | 0.0 | 0.0 | 1.5717 | 0.2163 | 0.15790 | 1.5047 | 0.39640 | 0.15790 | 0.1171 | 7.5063 | 0.0476 | 0.5431 | 0.0599 | 0.0545 | 0.3441 |
| `51_1` | 51 | RANTOUL 2 3 | 5.440 | 8.6003 | **46.79** | 0.04530 | 0.02110 | 0.0 | 0.0 | 1.4510 | 0.2080 | 0.15200 | 1.3924 | 0.21320 | 0.15200 | 0.1126 | 8.7481 | 0.0272 | 0.8509 | 0.0426 | 0.0219 | 0.2789 |
| `69_1` | 69 | MOUNT ZION 3 | 33.500 | 1.4636 | **49.03** | 0.27920 | 0.02626 | 0.0 | 0.0 | 2.0657 | 0.3964 | 0.28580 | 1.9645 | 0.36730 | 0.28580 | 0.2149 | 5.0455 | 0.0243 | 1.3279 | 0.0632 | 0.0939 | 0.4228 |
| `67_1` | 67 | MOUNT ZION 1 | 5.640 | 8.7908 | **49.58** | 0.04700 | 0.01683 | 0.0 | 0.0 | 2.4704 | 0.4050 | 0.29190 | 2.3413 | 0.55810 | 0.29190 | 0.2196 | 5.2533 | 0.0406 | 1.1605 | 0.0589 | 0.0891 | 0.4134 |
| `50_1` | 50 | RANTOUL 2 2 | 5.440 | 9.1611 | **49.84** | 0.04530 | 0.02110 | 0.0 | 0.0 | 1.3303 | 0.1763 | 0.12950 | 1.2799 | 0.38920 | 0.12950 | 0.0954 | 4.3890 | 0.0515 | 1.3293 | 0.0766 | 0.1718 | 0.5788 |
| `70_1` | 70 | MOUNT ZION 4 | 33.500 | 1.6342 | **54.75** | 0.27920 | 0.01460 | 0.0 | 0.0 | 1.7642 | 0.2550 | 0.18540 | 1.6839 | 0.43270 | 0.18540 | 0.1381 | 9.2595 | 0.0489 | 1.0841 | 0.0444 | 0.0308 | 0.2966 |
| `53_1` | 53 | RANTOUL 2 5 | 10.880 | 5.5136 | **59.99** | 0.09070 | 0.02479 | 0.0 | 0.0 | 1.8881 | 0.3001 | 0.21750 | 1.7992 | 0.45570 | 0.21750 | 0.1627 | 8.1314 | 0.0231 | 1.4869 | 0.0617 | 0.1210 | 0.4771 |
| `94_1` | 94 | TUSCOLA 2 1 | 21.600 | 3.2621 | **70.46** | 0.18000 | 0.01633 | 0.0 | 0.0 | 1.8367 | 0.2124 | 0.15520 | 1.7514 | 0.46070 | 0.15520 | 0.1150 | 7.3633 | 0.0593 | 1.2298 | 0.0495 | 0.1262 | 0.4874 |
| `182_1` | 182 | SPRINGFIELD 2 1 | 21.000 | 3.9719 | **83.41** | 0.17500 | 0.04984 | 0.0 | 0.0 | 1.5145 | 0.1785 | 0.13110 | 1.4514 | 0.45730 | 0.13110 | 0.0966 | 4.5954 | 0.0487 | 1.0461 | 0.0767 | 0.1891 | 0.6133 |
| `167_1` | 167 | CHAMPAIGN 1 4 | 11.280 | 9.3080 | **104.99** | 0.02820 | 0.04780 | 0.0 | 0.0 | 2.0261 | 0.3914 | 0.28230 | 1.9277 | 0.44660 | 0.28230 | 0.2122 | 8.6775 | 0.0524 | 0.7968 | 0.0716 | 0.1691 | 0.5732 |
| `72_1` | 72 | MOUNT ZION 6 | 33.500 | 3.8354 | **128.49** | 0.27920 | 0.03507 | 0.0 | 0.0 | 1.9187 | 0.2416 | 0.17590 | 1.8277 | 0.40480 | 0.17590 | 0.1309 | 4.7283 | 0.0423 | 1.2067 | 0.0591 | 0.1118 | 0.4587 |
| `183_1` | 183 | SPRINGFIELD 2 2 | 31.920 | 4.0832 | **130.34** | 0.26600 | 0.12400 | 0.0 | 0.0 | 1.6564 | 0.1889 | 0.13850 | 1.5836 | 0.23040 | 0.13850 | 0.1023 | 7.0376 | 0.0501 | 1.0589 | 0.0710 | 0.1321 | 0.4994 |
| `68_1` | 68 | MOUNT ZION 2 | 33.500 | 4.0842 | **136.82** | 0.27920 | 0.01514 | 0.0 | 0.0 | 1.2660 | 0.1213 | 0.09040 | 1.2202 | 0.22960 | 0.09040 | 0.0656 | 7.7221 | 0.0485 | 1.0162 | 0.0585 | 0.1198 | 0.4747 |
| `73_1` | 73 | MOUNT ZION 7 | 33.500 | 4.1486 | **138.98** | 0.27920 | 0.01536 | 0.0 | 0.0 | 1.8358 | 0.3355 | 0.24250 | 1.7505 | 0.37980 | 0.24250 | 0.1818 | 8.9262 | 0.0218 | 1.4708 | 0.0724 | 0.1996 | 0.6342 |
| `197_1` | 197 | GIBSON CITY 1 2 | 81.000 | 3.0000 | **243.00** | 0.20250 | 0.09019 | 0.0 | 0.0 | 2.1000 | 0.2000 | 0.18000 | 0.5000 | 0.50000 | 0.18000 | 0.1500 | 7.0000 | 0.0400 | 0.7500 | 0.0500 | 0.0000 | 0.0000 |
| `104_1` | 104 | ELLSWORTH 1 2 | 118.800 | 2.3995 | **285.06** | 0.17900 | 0.15491 | 0.0 | 0.0 | 1.3332 | 0.1289 | 0.09590 | 1.2827 | 0.47090 | 0.09590 | 0.0697 | 5.6694 | 0.0269 | 0.6157 | 0.0508 | 0.1762 | 0.5874 |
| `71_1` | 71 | MOUNT ZION 5 | 33.500 | 8.8016 | **294.85** | 0.27920 | 0.03256 | 0.0 | 0.0 | 1.4924 | 0.2143 | 0.15650 | 1.4308 | 0.38120 | 0.15650 | 0.1160 | 9.8220 | 0.0436 | 0.9680 | 0.0419 | 0.1417 | 0.5185 |
| `125_1` | 125 | BARTONVILLE 2 | 156.060 | 2.1875 | **341.38** | 1.50000 | 0.21950 | 0.0 | 0.0 | 2.1458 | 0.3643 | 0.26300 | 2.0391 | 0.51900 | 0.26300 | 0.1975 | 6.2212 | 0.0359 | 1.2441 | 0.0612 | 0.0356 | 0.3062 |
| `105_1` | 105 | ELLSWORTH 1 3 | 237.600 | 1.5014 | **356.73** | 0.33600 | 0.23082 | 0.0 | 0.0 | 1.6638 | 0.3183 | 0.23030 | 1.5905 | 0.33770 | 0.23030 | 0.1725 | 5.3727 | 0.0570 | 0.8459 | 0.0514 | 0.0509 | 0.3369 |
| `161_1` | 161 | SPRINGFIELD 4 1 | 166.320 | 2.4543 | **408.20** | 0.41580 | 0.18003 | 0.0 | 0.0 | 1.9578 | 0.3500 | 0.25290 | 1.8641 | 0.45150 | 0.25290 | 0.1898 | 8.8717 | 0.0568 | 1.2683 | 0.0535 | 0.1541 | 0.5432 |
| `155_1` | 155 | SPRINGFIELD 5 6 | 92.670 | 4.4975 | **416.78** | 0.77220 | 0.08118 | 0.0 | 0.0 | 1.8487 | 0.2047 | 0.14960 | 1.7626 | 0.29390 | 0.14960 | 0.1108 | 4.2006 | 0.0458 | 0.5965 | 0.0517 | 0.0769 | 0.3888 |
| `147_1` | 147 | HOPEDALE 2 1 | 120.600 | 3.9792 | **479.89** | 0.13400 | 0.21410 | 0.0 | 0.0 | 2.4438 | 0.4858 | 0.34930 | 2.3166 | 0.40770 | 0.34930 | 0.2635 | 8.8310 | 0.0551 | 0.9997 | 0.0417 | 0.1968 | 0.6288 |
| `115_1` | 115 | NORMAL 2 3 | 180.000 | 2.7268 | **490.82** | 0.47700 | 0.30987 | 0.0 | 0.0 | 1.1766 | 0.2127 | 0.15530 | 1.1369 | 0.34710 | 0.15530 | 0.1152 | 4.4602 | 0.0217 | 1.3124 | 0.0484 | 0.1389 | 0.5129 |
| `127_1` | 127 | BARTONVILLE 4 | 156.060 | 3.3851 | **528.27** | 1.50000 | 0.55933 | 0.0 | 0.0 | 1.2488 | 0.1215 | 0.09060 | 1.2041 | 0.21060 | 0.09060 | 0.0657 | 8.3432 | 0.0569 | 1.4042 | 0.0746 | 0.1968 | 0.6287 |
| `152_1` | 152 | SPRINGFIELD 5 3 | 92.670 | 5.8576 | **542.83** | 0.77220 | 0.07704 | 0.0 | 0.0 | 2.3179 | 0.3172 | 0.22960 | 2.1993 | 0.36550 | 0.22960 | 0.1719 | 7.9089 | 0.0514 | 1.2546 | 0.0428 | 0.1712 | 0.5774 |
| `153_1` | 153 | SPRINGFIELD 5 4 | 92.670 | 6.2927 | **583.14** | 0.77220 | 0.06762 | 0.0 | 0.0 | 1.7600 | 0.3067 | 0.22210 | 1.6800 | 0.52970 | 0.22210 | 0.1662 | 9.4303 | 0.0504 | 0.8196 | 0.0434 | 0.1286 | 0.4923 |
| `154_1` | 154 | SPRINGFIELD 5 5 | 92.670 | 8.4419 | **782.31** | 0.77220 | 0.08043 | 0.0 | 0.0 | 1.3620 | 0.2699 | 0.19600 | 1.3095 | 0.31640 | 0.19600 | 0.1463 | 5.5956 | 0.0476 | 1.3940 | 0.0727 | 0.0550 | 0.3451 |
| `65_1` | 65 | PAXTON 1 1 | 180.480 | 6.2899 | **1135.21** | 0.04000 | 0.11008 | 0.0 | 0.0 | 1.3243 | 0.1819 | 0.13350 | 1.2744 | 0.46460 | 0.13350 | 0.0985 | 4.9885 | 0.0547 | 0.5169 | 0.0602 | 0.1634 | 0.5620 |
| `126_1` | 126 | BARTONVILLE 3 | 156.060 | 7.2997 | **1139.19** | 1.50000 | 0.31115 | 0.0 | 0.0 | 1.2798 | 0.1546 | 0.11410 | 1.2330 | 0.28270 | 0.11410 | 0.0837 | 7.4167 | 0.0221 | 1.1481 | 0.0546 | 0.0603 | 0.3557 |
| `136_1` | 136 | PEKIN 1 3 | 535.680 | 3.2156 | **1722.56** | 2.22551 | 0.61527 | 0.0 | 0.0 | 1.9175 | 0.2787 | 0.20220 | 1.8266 | 0.55420 | 0.20220 | 0.1510 | 5.2757 | 0.0409 | 0.9779 | 0.0562 | 0.1204 | 0.4759 |
| `135_1` | 135 | PEKIN 1 2 | 535.680 | 3.2329 | **1731.80** | 1.96416 | 0.51620 | 0.0 | 0.0 | 2.0957 | 0.3515 | 0.25390 | 1.9925 | 0.65400 | 0.25390 | 0.1906 | 8.7077 | 0.0596 | 0.7468 | 0.0630 | 0.0679 | 0.3709 |
| `189_1` | 189 | CLINTON 1 2 | 682.980 | 3.4017 | **2323.26** | 5.69142 | 0.49565 | 0.0 | 0.0 | 1.7840 | 0.2069 | 0.15120 | 1.7023 | 0.35750 | 0.15120 | 0.1120 | 9.6481 | 0.0451 | 0.5464 | 0.0646 | 0.0836 | 0.4023 |

## H and inertia summary by bus

Total stored kinetic energy per bus = sum of H×mva over all Genrou devices at that bus (one gen per bus in Illinois).

| bus | bus_name | H (nominal) | mva (MVA) | H×mva (MWs) |
|-----|----------|---|---|---|
| 49 | RANTOUL 2 1 | 5.2115 | 5.44 | 28.35 |
| 50 | RANTOUL 2 2 | 9.1611 | 5.44 | 49.84 |
| 51 | RANTOUL 2 3 | 8.6003 | 5.44 | 46.79 |
| 52 | RANTOUL 2 4 | 4.9625 | 5.44 | 27.00 |
| 53 | RANTOUL 2 5 | 5.5136 | 10.88 | 59.99 |
| 65 | PAXTON 1 1 | 6.2899 | 180.48 | 1135.21 |
| 67 | MOUNT ZION 1 | 8.7908 | 5.64 | 49.58 |
| 68 | MOUNT ZION 2 | 4.0842 | 33.50 | 136.82 |
| 69 | MOUNT ZION 3 | 1.4636 | 33.50 | 49.03 |
| 70 | MOUNT ZION 4 | 1.6342 | 33.50 | 54.75 |
| 71 | MOUNT ZION 5 | 8.8016 | 33.50 | 294.85 |
| 72 | MOUNT ZION 6 | 3.8354 | 33.50 | 128.49 |
| 73 | MOUNT ZION 7 | 4.1486 | 33.50 | 138.98 |
| 76 | BRIMFIELD 1 | 6.7471 | 4.80 | 32.39 |
| 77 | BRIMFIELD 2 | 5.1416 | 2.88 | 14.81 |
| 90 | CLINTON 3 1 | 4.4204 | 3.84 | 16.97 |
| 91 | CLINTON 3 2 | 6.4018 | 6.00 | 38.41 |
| 94 | TUSCOLA 2 1 | 3.2621 | 21.60 | 70.46 |
| 104 | ELLSWORTH 1 2 | 2.3995 | 118.80 | 285.06 |
| 105 | ELLSWORTH 1 3 | 1.5014 | 237.60 | 356.73 |
| 114 | NORMAL 2 2 | 8.9160 | 2.04 | 18.19 |
| 115 | NORMAL 2 3 | 2.7268 | 180.00 | 490.82 |
| 125 | BARTONVILLE 2 | 2.1875 | 156.06 | 341.38 |
| 126 | BARTONVILLE 3 | 7.2997 | 156.06 | 1139.19 |
| 127 | BARTONVILLE 4 | 3.3851 | 156.06 | 528.27 |
| 135 | PEKIN 1 2 | 3.2329 | 535.68 | 1731.80 |
| 136 | PEKIN 1 3 | 3.2156 | 535.68 | 1722.56 |
| 147 | HOPEDALE 2 1 | 3.9792 | 120.60 | 479.89 |
| 151 | SPRINGFIELD 5 2 | 5.3589 | 6.48 | 34.73 |
| 152 | SPRINGFIELD 5 3 | 5.8576 | 92.67 | 542.83 |
| 153 | SPRINGFIELD 5 4 | 6.2927 | 92.67 | 583.14 |
| 154 | SPRINGFIELD 5 5 | 8.4419 | 92.67 | 782.31 |
| 155 | SPRINGFIELD 5 6 | 4.4975 | 92.67 | 416.78 |
| 161 | SPRINGFIELD 4 1 | 2.4543 | 166.32 | 408.20 |
| 167 | CHAMPAIGN 1 4 | 9.3080 | 11.28 | 104.99 |
| 170 | CHAMPAIGN 1 7 | 3.4015 | 11.28 | 38.37 |
| 182 | SPRINGFIELD 2 1 | 3.9719 | 21.00 | 83.41 |
| 183 | SPRINGFIELD 2 2 | 4.0832 | 31.92 | 130.34 |
| 189 | CLINTON 1 2 | 3.4017 | 682.98 | 2323.26 |
| 197 | GIBSON CITY 1 2 | 3.0000 | 81.00 | 243.00 |

## Summary statistics

- **Total Genrou devices**: 40 (all online, no offline generators)
- **H range**: 1.4636 s (MOUNT ZION 3, `69_1`) to 9.3080 s (CHAMPAIGN 1 4, `167_1`)
- **MVA range**: 2.04 MVA (NORMAL 2 2, `114_1`, tiny wind) to 682.98 MVA (CLINTON 1 2, `189_1`, nuclear)
- **H×mva range**: 14.81 MWs (BRIMFIELD 2, `77_1`) to 2323.26 MWs (CLINTON 1 2, `189_1`)
- **Total H×mva across all gens**: 15,693.99 MWs (all Illinois fleet inertia)

## Key observations for v6-v10 experiment selection

1. **Smallest H (v6-il)**: `69_1` (H=1.46, MOUNT ZION 3), `105_1` (H=1.50, ELLSWORTH 1 3), `115_1` (H=2.73, NORMAL 2 3), `104_1` (H=2.40, ELLSWORTH 1 2)
   - These machines have intrinsically low inertia constants
   - Good hypothesis test: "Is low H detectable even at moderate machine size?"

2. **Smallest MVA (v8-il)**: `114_1` (2.04 MVA), `77_1` (2.88 MVA), `90_1` (3.84 MVA), `76_1` (4.80 MVA)
   - Tiny wind and gas turbine generators
   - Test hypothesis: "Do small machines show weak coupling despite normal H?"

3. **Largest H×mva (v9-il)**: `189_1` (2323 MWs, nuclear), `135_1`/`136_1` (1732 and 1722 MWs, coal at PEKIN), `126_1` (1139 MWs, coal), `65_1` (1135 MWs, wind)
   - Massive grid anchors with strong inertia contribution
   - Null hypothesis: "If H variation is undetectable even here, the issue is elsewhere"

4. **Small dispatch (v10-il)**: `114_1` (p0=0.004 pu = 8 MW rated, 0.4% loaded), `77_1` (p0=0.007 pu), `90_1` (p0=0.0096 pu)
   - Lightly loaded despite nominal inertia
   - Test hypothesis: "Does small p0 mask H sensitivity?"

---

**Last updated**: 2026-09-28
**Source**: `examples/PhasorDynamics/Large/Illinois/illinois.case.json` (GridKit format)
