# Hawaii dynamic controls

| EMT class | PowerWorld model | Connections |
| --- | --- | --- |
| `Esst1a` | ESST1A | Terminal RMS voltage `vt`, exciter-base field current `ifd`, voltage reference `vref`, stabilizer `vs` and limiters `vuel`, `voel`. Drive `Machine.efd` from `EFD`. |
| `Exst1Pti` | EXST1_PTI | Same signals. Both field-voltage limits include the commutation drop. |
| `Esst4b` | ESST4B | Also supply phase voltages `v` from the bus and net injected phase currents `i` from `Machine.i` for the compound rectifier. |
| `Ggov1` | GGOV1 | Machine speed `omega` and terminal power `p` in watts. `PM` supplies mechanical power on the machine base. |
| `Ieeeg1` | IEEEG1 | Machine speed `omega`. Odd turbine fractions drive `PM`, even fractions drive `PM_LP`, both on the configured machine base. |
| `Reeca` | REECA1 / REEC_A | Terminal RMS voltage `vt`, active and reactive powers `p`, `q`, plant references and optional turbine speed. Connect `IPCMD`, `IQCMD` to the converter. |

Supply `vt` in line-to-line RMS volts. Connect each exciter's `ifd` to
`Machine.ifd_exc`. For an unused overexcitation limiter, set the held
`voel` above the regulator range on ST1A and ST4B. EXST1_PTI uses a
summed limiter signal with zero as its neutral value.

For the governors, `S` is the machine rating in VA. Zero `Trate` selects
that rating. GGOV1 uses the 1 ms fuel-command memory `Tfsr` from gridparse
in its acceleration selector and requires zero engine transport delay
`Teng`. State files and events set the physical fuel command memory in
`MEMORY`. The governor integrates `MEMORY_DEV`, its deviation from the
initial `MEMORY`. IEEEG1 uses two six-point valve arrays, `Gv` and `Pgv`. All
zeros select a straight line. Hawaii uses zero hysteresis `Eps` and
backlash `db2`, the settings supported here.

Initialize REECA's `PREF_RATE` at the device-base active power reference.
The reference ramps between the specified rate limits and tracks its
target with a 1 ms time constant inside those limits. The power-order
lag follows the ramp limiter. Set `VTERM` to the initial per-unit
terminal voltage. Four-point `Vq`/`Iq` and `Vp`/`Ip` tables limit
current as a function of filtered voltage. All-zero tables select the
converter current circle. When voltage recovers, `Thld2` holds the
fault-end active-current limit for the specified duration. [PowerWorld current-limit equations](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Exciter%20REEC_A.htm)
