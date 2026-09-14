#!/usr/bin/env python
"""Per-sigma-pair check that the 60"/40" candidate pruning reproduces a full
0.12-deg cone scan. With anisotropic sigma a star far away in one coordinate
but well aligned in the other can still win, so the circular-radius guarantee
that null_ts.py relies on is not automatically valid here.
"""
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "null_resample"))
from null_multi_sigma import MultiSigmaEngine
from racs_configs import unique_configs, VARIANTS
from run_null_galactic import sample_footprint

REPO = os.path.dirname(os.path.dirname(HERE))
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300

src = pd.read_csv(os.path.join(REPO, "outputs", "racs_agn_positions.csv"))
cfg, index = unique_configs(src)
used = {v: set(index[v].tolist()) for v in VARIANTS}
d = np.load(os.path.join(REPO, "data", "gaia_null", "gaia_galactic_g20.npz"))
eng = MultiSigmaEngine(d["ra"], d["dec"], d["g_mag"])
ra, dec = sample_footprint(N, np.random.default_rng(31337))

t0 = time.time()
fast = eng.best_ts(ra, dec, cfg)
mid = time.time()
slow = eng.best_ts(ra, dec, cfg, exact=True)
print(f"pruned {mid-t0:.0f}s, exact {time.time()-mid:.0f}s ({N} events)")

diff = slow - fast
print(f"\nany pruned > exact (would be a bug): {bool((diff < -1e-9).any())}")
print(f"{'variant':18s} {'cfgs':>5s} {'sig_ra':>7s} {'sig_dec':>7s} "
      f"{'frac ev diff':>12s} {'max diff':>10s} {'q99 shift':>10s}")
for v in VARIANTS:
    js = sorted(used[v])
    fr = (diff[:, js] > 1e-9).mean()
    q_f = np.percentile(fast[:, js], 99, axis=0).mean()
    q_s = np.percentile(slow[:, js], 99, axis=0).mean()
    print(f"{v:18s} {len(js):5d} {cfg[js,0].min():7.3f} {cfg[js,1].min():7.3f} "
          f"{fr:12.4f} {diff[:, js].max():10.3g} {q_s - q_f:10.4f}")

worst = np.unravel_index(np.argmax(diff), diff.shape)
print(f"\nworst single cell: event {worst[0]} sigma {cfg[worst[1]]} "
      f"pruned {fast[worst]:.3g} exact {slow[worst]:.3g}")
