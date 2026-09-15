#!/usr/bin/env python3
"""Generate the PWM figure data from the continuous sigmoid pulse sum."""
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
M, fm, fc, alignment = .8, 60., 900., .5
t = np.linspace(0, 1 / fm, 4001)
duty = .5 * (1 + M * np.sin(2 * np.pi * fm * t[:, None] + [0, -2*np.pi/3, 2*np.pi/3]))
phase = np.remainder(t * fc, 1)[:, None]
on = alignment * (1 - duty)
off = alignment + (1 - alignment) * duty
columns = [t]
for mu in [200000., 1000.]:
    radius = int(np.ceil(np.log(4 / np.finfo(float).eps) * fc / mu)) + 1
    s = np.zeros_like(duty)
    for k in range(-radius, radius + 1):
        x, y = .5 * mu / fc * (phase - on - k), .5 * mu / fc * (phase - off - k)
        s += .5 * (np.tanh(x) - np.tanh(y))
    columns.extend(s.T)
np.savetxt(HERE / 'waveforms.dat', np.column_stack(columns), fmt='%.12g',
           header='t sa_top sb_top sc_top sa_bottom sb_bottom sc_bottom', comments='')
