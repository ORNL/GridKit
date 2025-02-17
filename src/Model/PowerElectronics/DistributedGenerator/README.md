
The Distributed Generator Component found in references 1 and 2.

Parameters:
 + $\omega_b$ - Reference Rotating Frame
 + $\omega_c$ - Cutoff Frequency
 + $m_p$ - Drop Gain, Frequency Range
 + $V_n$ - Nominal Set Point of D-Axis Output Voltage
 + $n_q$ - Voltage Range
 + $F$ - PI Controller Parameter in 1 & 2
 + $K_{pv}$ - PI Controller Parameter in 1 & 2
 + $K_{iv}$ - PI Controller Parameter in 1 & 2
 + $K_{pc}$ - PI Controller Parameter in 1 & 2
 + $C_f$ - Shunt??
 + $r_{Lf}$ - Resistance of line f
 + $L_{f}$ - Inductance of line f
 + $r_{Lc}$ - Resistance of line c
 + $L_{c}$ - Inductance of line c

Variables (External):
 + $\omega_{ref}$   - Network reference $\omega$
 + $V_{a}$          - Incoming Bus Voltage (a)
 + $V_{b}$          - Incoming Bus Voltage (b)

Variables (Internal):
 + $\delta$     - Rotor difference from reference
 + $\omega$     - Frequency
 + $P$          - Real Power
 + $Q$          - Reactive Power
 + $\phi_{d}$   - Output Voltage Control PI Variable
 + $\phi_{q}$   - Output Voltage Control PI Variable
 + $\gamma_{d}$ - Output Current Control PI Variable
 + $\gamma_{q}$ - Output Current Control PI Variable
 + $i_{ld}$     - Current of Line l (dq-space)
 + $i_{lq}$     - Current of Line l (dq-space)
 + $v_{od}$     - Voltage of Bus o (dq-space)
 + $v_{oq}$     - Voltage of Bus o (dq-space)
 + $i_{od}$     - Current of Line o (dq-space)
 + $i_{oq}$     - Current of Line o (dq-space)

Equations (External, Residuals):
 + $\omega_{com} - \omega$ &emsp;&emsp;&emsp;(If this generator is considered the reference one, otherwise 0)
 + $\cos(\delta) i_{od} - \sin(\delta) i_{oq}$
 + $\sin(\delta) i_{od} + \cos(\delta) i_{oq}$

Equations (Internal):
 + $\omega_{com} = \omega_{b} - m_{p} P$
 + $\frac{d\delta}{dt} = \omega_{com} - \omega$
 + $\frac{dP}{dt} = \omega_c ( v_{od} i_{od} + v_{oq} i_{oq} - P)$
 + $\frac{dQ}{dt} = \omega_c ( v_{od} i_{oq} + v_{oq} i_{od} - Q)$
 + $v_{od}^* = V_{n} - n_q Q$
 + $v_{oq}^* = 0$
 + $\frac{d\phi_{d}}{dt} = v_{od}^* - v_{od}$
 + $\frac{d\phi_{q}}{dt} = v_{oq}^* - v_{oq}$
 + $i_{ld}^* = F i_{od} - \omega_{b} C_{f} v_{oq} + K_{pv} (v_{od}^* - v_{od}) + K_{iv} \phi_{d}$
 + $i_{lq}^* = F i_{oq} - \omega_{b} C_{f} v_{od} + K_{pv} (v_{oq}^* - v_{oq}) + K_{iv} \phi_{q}$
 + $\frac{d\gamma_{d}}{dt} = i_{ld}^* - i_{ld}$
 + $\frac{d\gamma_{q}}{dt} = i_{lq}^* - i_{lq}$
 + $v_{id}^* = -\omega_{b} L_{f} i_{lq} + K_{pc} ( i_{ld}^* - i_{ld}) + K_{ic} \gamma_{d}$
 + $v_{iq}^* = -\omega_{b} L_{f} i_{ld} + K_{pc} ( i_{lq}^* - i_{lq}) + K_{ic} \gamma_{q}$
 + $\frac{di_{ld}}{dt} = -(\frac{r_{Lf}}{L_{f}}) i_{ld} + \omega_{com} i_{lq} + \frac{v_{id}^* - v_{id}}{L_f}$
 + $\frac{di_{lq}}{dt} = -(\frac{r_{Lf}}{L_{f}}) i_{lq} + \omega_{com} i_{ld} + \frac{v_{iq}^* - v_{iq}}{L_f}$
 + $\frac{dv_{od}}{dt} = \omega_{com} v_{oq} + \frac{i_{ld} - i_{od}}{C_f}$
 + $\frac{dv_{oq}}{dt} = -\omega_{com} v_{od} + \frac{i_{lq} - i_{oq}}{C_f}$
 + $V_{bd,in} = \cos(\delta) V_{a} - \sin(\delta) V_{b}$
 + $V_{bq,in} = -\sin(\delta) V_{a} + \cos(\delta) V_{b}$
 + $\frac{di_{od}}{dt} = -(\frac{r_{Lc}}{L_{c}}) i_{od} + \omega_{com} i_{oq} + \frac{v_{od} - V_{bd,in}}{L_f}$
 + $\frac{di_{oq}}{dt} = -(\frac{r_{Lc}}{L_{c}}) i_{oq} + \omega_{com} i_{ld} + \frac{v_{oq} - V_{bq,in}}{L_f}$

Note all internal direct equalities are simplified into the differential equations.


1. Pogaku, Nagaraju, Milan Prodanovic, and Timothy C. Green. "Modeling, analysis and testing of autonomous operation of an inverter-based microgrid." IEEE Transactions on power electronics 22.2 (2007): 613-625.
2. Bidram, Ali, Frank L. Lewis, and Ali Davoudi. "Distributed control systems for small-scale power networks: Using multiagent cooperative control theory." IEEE Control systems magazine 34.6 (2014): 56-77.
