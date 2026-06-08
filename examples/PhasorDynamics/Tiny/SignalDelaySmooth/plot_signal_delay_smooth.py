#!/usr/bin/env python3
"""Plot the SignalDelaySmooth example output."""

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

csv_path = sys.argv[1] if len(sys.argv) > 1 else "SignalDelaySmooth.out.csv"
png_path = sys.argv[2] if len(sys.argv) > 2 else "SignalDelaySmooth.png"

data = np.genfromtxt(csv_path, delimiter=",", names=True)
t = data["t"]
src = data["SampledSignalSource_SRC1_out"]
delayed = data["DelaySmooth_DS1_out"]
vss = data["Ieeest_PSS1_vss"]
tau = 0.25
ideal = np.interp(t - tau, t, src, left=src[0], right=src[-1])

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(t, src, label="source (signal 1)", lw=2)
ax.plot(t, ideal, label="ideal 0.25 s shift", lw=1.5, ls=":")
ax.plot(t, delayed, label="DelaySmooth output", lw=2, ls="--")
ax.plot(t, vss, label="IEEEST $V_{ss}$ (signal 3)", lw=1.5, alpha=0.8)

ax.set_title("SignalDelaySmooth: lag-chain approximation of a 0.25 s delay")
ax.set_xlabel("time [s]")
ax.set_ylabel("signal [-]")
ax.grid(True, alpha=0.3)
ax.legend(loc="upper right")
fig.tight_layout()
fig.savefig(png_path, dpi=150)
print(f"wrote {png_path}")
