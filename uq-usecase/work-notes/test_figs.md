# Test figures (for review, not in report)

## Figure 1: unstable vs censored, same fault onset

![Rotor-angle deviation from COI for an unstable and a censored aleatoric-1 sample](test_figs/delta_unstable_vs_censored.png)

*Figure 1. Rotor-angle deviation from COI for one unstable sample (top) and one censored sample (bottom), both aleatoric-1 baseline (no datacenter), on a shared time axis. Orange: fault on (0.25 s top, 0.5 s bottom). Top: the red trace passes $\pi$ (dotted line), a pole slip, and keeps growing, so the machine loses synchronism. Bottom: the red trace peaks at $0.47\pi$ and never reaches $\pi$, so the system recovers within the 10 s run. Grey: each machine's $\delta_i - \delta_{COI}$. Red: maximum over machines of $|\delta_i - \delta_{COI}|$ (the machine can change). Dashed: $\pm\pi$, the pole-slip threshold.*

Verification: both panels are replays of the stored `case.json` and `solver.json` with the campaign's DynamicSimulation runner. The unstable replay matches the stored 0.1 s series to 5.95e-8 rad. The censored replay has no stored series to compare against.
