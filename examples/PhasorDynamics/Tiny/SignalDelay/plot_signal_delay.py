#!/usr/bin/env python3
"""Plot the SignalDelay example output: a sampled source, its 0.25 s delayed
copy produced by the Delay block, and the resulting IEEEST stabilizer output.

Usage:
    python3 plot_signal_delay.py [SignalDelay.out.csv] [SignalDelay.png]

The CSV is produced by running `DynamicSimulation SignalDelay.solver.json`.
"""

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

csv_path = sys.argv[1] if len(sys.argv) > 1 else "SignalDelay.out.csv"
png_path = sys.argv[2] if len(sys.argv) > 2 else "SignalDelay.png"

data = np.genfromtxt(csv_path, delimiter=",", names=True)
t = data["t"]
src = data["SampledSignalSource_SRC1_out"]
delayed = data["Delay_DLY1_out"]
vss = data["Ieeest_PSS1_vss"]

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(t, src, label="source (signal 1)", lw=2)
ax.plot(t, delayed, label="delayed 0.25 s (signal 2)", lw=2, ls="--")
ax.plot(t, vss, label="IEEEST $V_{ss}$ (signal 3)", lw=1.5, alpha=0.8)

ax.set_title("SignalDelay: a Delay block lags the source by 0.25 s")
ax.set_xlabel("time [s]")
ax.set_ylabel("signal [-]")
ax.grid(True, alpha=0.3)
ax.legend(loc="upper right")
fig.tight_layout()
fig.savefig(png_path, dpi=150)
print(f"wrote {png_path}")
