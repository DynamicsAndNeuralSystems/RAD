"""Write the shared test series and their expected RAD values.

The reference here follows Harris, Gollo & Fulcher, PRX 14, 031021 (2024),
Eqs. (9)-(12), and is written independently of RAD.py so the tests check each
implementation against the paper, not against itself.
Rerun only to add cases: `python make_reference.py`.
"""
import numpy as np


def rad_reference(x, tau, do_centre):
    x = np.asarray(x, dtype=float)
    if do_centre:
        x = np.abs(x - np.median(x))
    lead, lag = x[tau:], x[:-tau]
    lower = lag[lag < np.median(lag)]
    upper = lag[lag >= np.median(lag)]
    if len(lower) < 2 or len(upper) < 2:
        return np.nan
    s_lower, s_upper = np.std(lower, ddof=1), np.std(upper, ddof=1)
    if s_lower == 0 or s_upper == 0:
        return np.nan
    return np.std(lead - lag, ddof=1) * (1 / s_upper - 1 / s_lower)


rng = np.random.default_rng(20241021)
n = 500
ar1 = np.zeros(n)
for t in range(1, n):
    ar1[t] = 0.9 * ar1[t - 1] + rng.standard_normal()
series = {
    "ar1": ar1,
    "sine_noise": np.sin(np.arange(n) / 10) + 0.3 * rng.standard_normal(n),
    "short": ar1[:20],
    "integer": np.round(ar1),        # tied values: degenerate lower half
    "constant": np.ones(50),
    "tiny": ar1[:4],
}

with open("series.csv", "w") as f:
    for name, x in series.items():
        f.write(name + "," + ",".join(f"{v:.17g}" for v in x) + "\n")

with open("expected.csv", "w") as f:
    f.write("series,tau,doCentre,RAD\n")
    for name, x in series.items():
        for tau in (1, 2):
            for c in (1, 0):
                f.write(f"{name},{tau},{c},{rad_reference(x, tau, c):.17g}\n")
