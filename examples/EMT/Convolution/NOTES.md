The delay is being read from frequency samples the fit itself throws away. In MinimumPhase.cpp:27, minimumDelay takes the minimum of the phase-delay trace τ_p(ω) = ℓ·β(ω)/ω over all modes and all samples. I checked this against the actual 345 kV sweep data (LineGallery/output/345kv-horizontal/response.csv):

mode	where its τ_p minimum sits	|H| at that sample
0 (fast aerial)	334.7889 µs @ 88.6 MHz	6.0e-2 — alive
1	334.6361 µs @ 100 MHz	1.6e-18 — dead
2 (ground)	334.7232 µs @ 100 MHz	2.7e-158 — dead
The global minimum — the τ that went into the artifact — is mode 1's value at 100 MHz, where that mode's magnitude is 1.6e-18. Meanwhile the fitter's --min-mag floor (1e-4 of peak) zeroes exactly those samples before fitting. So τ is dictated by data that never participates in the fit. Physically: below its alive band, a mode's phase delay keeps sliding toward the lossless front ℓ/c = 333.56 µs, a wavefront that carries exponentially little energy. The actual aerial wavefront arrives at 334.789 µs.

The consequence is exactly what hmin_time.png shows: unwinding by 334.636 µs leaves 0.153 µs of dead time — nearly half the kernel's own width — that the rational fit must fake. A stable rational function can't represent pure delay, so it rings (the pre-peak oscillation train) and pushes 4% onto the δ(t) term. Poles are spent synthesizing delay instead of shape.

The fix
Restrict the existing rule to the samples the fit sees: per mode, only samples with |H_m| ≥ min_mag·max|H_m|, using the same floor the fitter already applies. On the real sweep data this gives:

mode 0 alive-band minimum: 334.7889 µs ← new shared τ
mode 1 (alive to 7.6 MHz): 334.8038 µs
mode 2 (alive to 172 kHz): 346.0070 µs
The Hmin peak in the current artifact sits at 334.636 + 0.153 = 334.789 µs — the corrected τ matches the observed wavefront to under a nanosecond. It's also safe by construction: since τ ≤ τ_p,m(ω) at every alive sample, ωτ ≤ ℓβ_m(ω) everywhere the fit looks — no fitted sample is ever unwound past zero phase, so no acausal content is created. It's the largest τ with that property, i.e., the tightest correct delay.

I also ran a Bode gain–phase (minimum-phase excess delay) analysis as a cross-check. It agrees for the sharp aerial mode (~334.78 µs), and it explains why you should not go further toward per-mode-optimal delays for the shared τ: the ground mode's ~78 µs of extra lateness is dispersion (the e^(−c√ω) character), which is approximately minimum-phase and absorbable by the rational part without a delay line. Min-over-modes of the alive-band phase delay is the right rule for one shared delay.

Concretely, the cleaner code is a small, contained change (these helpers are used only by the ULM app):

MinimumPhase.cpp — give minimumDelay the propagation samples and the floor: minimumDelay(tau, h, mag_floor). Per mode, find max|H_m|, then take the min of τ_p only over samples with |H_m| ≥ floor·peak. Document the invariant ("largest delay that unwinds no alive sample past zero phase") rather than the old "min over everything".
UniversalLineModel.cpp:664 — pass h and settings.h.min_mag through in fitPropagation.
Same treatment for delaySpread (UniversalLineModel.cpp:580) — this is a second real bug the same root cause hides: it currently reports an inter-mode spread of 0.153 µs (0.05%), but on the alive basis the ground mode's minimum is 346.0 µs, a spread of ~11 µs. The diagnostic that's supposed to tell you whether one shared delay is apt is understating the answer by ~70×.
After the change, regenerate the propagation artifact and rerun time_domain.py: the Hmin wavefront should stand at the origin, the pre-peak ring train should collapse, D should drop, and the order search will likely meet its target with fewer poles since it stops spending them on delay.

I've kept this as a proposal since you asked what's wrong and how to fix it — say the word and I'll make the change and regenerate the figures.
