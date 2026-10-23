#!/usr/bin/env python3
"""S7 (MASCOTTE_A10_PREDICTION.md): the stopped bh_s7z.sh runs, read from their own output only (no measured value).

usage: analyse.py DIR [FIELDDIR]
  DIR holds <run>_history.csv and <run>_to<ms>_stations.csv (this directory).
  FIELDDIR, if given, holds <run>_to<ms>_field.csv (backhouse: /home/greff/crucible_s7z/runs) and the two
  figures are written into DIR.
Prints, per window, section 2's half-to-half drift at r <= 20 mm per station and its change from the previous window;
per run, the liquid ledger, protocol section 3's 1 ms drift (the RL10's settle measure) and the spray's share of the
wall time.
"""
import csv, glob, os, re, statistics as st, sys

RUNS = ['low', 'high', 'blob', 'blob_fine']
D = sys.argv[1]


def rows(p):
    return list(csv.DictReader(open(p)))


for run in RUNS:
    print(f'== {run}')
    prev = None
    wins = sorted(glob.glob(f'{D}/{run}_to*_stations.csv'), key=lambda p: int(re.search(r'_to(\d+)_', p)[1]))
    for p in wins:
        ms = int(re.search(r'_to(\d+)_', p)[1])
        S = [r for r in rows(p) if float(r['r_mm']) <= 20]
        out = []
        for z in sorted({float(r['z_mm']) for r in S}):
            d = [abs(float(r['T_first_half']) - float(r['T_second_half'])) for r in S if float(r['z_mm']) == z]
            out.append(f'z {z:g}: max {max(d):.1f} rms {st.mean(x * x for x in d) ** 0.5:.1f}')
        drift = max(abs(float(r['T_first_half']) - float(r['T_second_half'])) for r in S)
        line = f'  window {ms - 10} to {ms} ms: drift {drift:.1f} K ({"; ".join(out)})'
        T = {(r['z_mm'], r['r_mm']): float(r['T']) for r in S}
        if prev:
            c = [abs(T[k] - prev[k]) for k in T]
            line += f'; change from the previous window: max {max(c):.1f} rms {st.mean(x * x for x in c) ** 0.5:.1f} K'
        prev = T
        print(line)
    H = rows(f'{D}/{run}_history.csv')
    f = lambda r, k: float(r[k])
    t = [f(r, 't') for r in H]
    held = []
    for ms in (5, 10, 15, 20, 25):
        i = min(range(len(t)), key=lambda i: abs(t[i] - ms * 1e-3))
        held.append(f'{ms} ms {f(H[i], "liquid_mass") * 1e3:.3f} g')
    print(f'  end {t[-1] * 1e3:.2f} ms; liquid held: {", ".join(held)}, end {f(H[-1], "liquid_mass") * 1e3:.3f} g;'
          f' liquid out {f(H[-1], "liquid_out") * 1e3:.3f} g')
    late = [r for r in H if f(r, 't') > 15e-3]
    o = [f(r, 'outlet') for r in late]
    print(f'  from 15 ms: outlet {min(o) * 1e3:.1f} to {max(o) * 1e3:.1f} g/s, mean {st.mean(o) * 1e3:.2f};'
          f' inlet (gas) {f(H[-1], "inlet") * 1e3:.2f} g/s')
    ev = [(f(H[i], 'evaporated') - f(H[i - 1], 'evaporated')) / (t[i] - t[i - 1]) for i in range(1, len(H))]
    do, dp, me = [], [], []
    k = 0
    while 15e-3 + (k + 1) * 1e-3 <= t[-1] + 1e-12:
        a, b = 15e-3 + k * 1e-3, 15e-3 + (k + 1) * 1e-3
        w = [i for i in range(len(H)) if a < t[i] <= b]
        oo = [f(H[i], 'outlet') for i in w]
        pp = [f(H[i], 'p_injector') for i in w]
        do.append((max(oo) - min(oo)) / abs(oo[-1]))
        dp.append((max(pp) - min(pp)) / abs(pp[-1]))
        me.append(st.mean(ev[i - 1] for i in w))
        k += 1
    print(f'  1 ms drift (max - min) / |last| over {k} windows from 15 ms: outlet median {st.median(do) * 100:.2f}%'
          f' (min {min(do) * 100:.2f}%), injector pressure median {st.median(dp) * 100:.3f}% (min {min(dp) * 100:.3f}%);'
          f' under 0.1%: outlet {sum(x < 1e-3 for x in do)}, pressure {sum(x < 1e-3 for x in dp)}')
    print(f'  1 ms mean evaporation rate: mean {st.mean(me) * 1e3:.2f} g/s, sd {st.stdev(me) * 1e3:.2f} g/s'
          f' ({st.stdev(me) / st.mean(me) * 100:.2f}%)')
    s, w = f(H[-1], 'spray_wall_s'), f(H[-1], 'wall_s')
    print(f'  spray {s:.0f} s of {w:.0f} s wall ({s / w * 100:.1f}%; the whole step at {1 / (1 - s / w):.2f} times the gas'
          f" step's cost); parcels at the end {H[-1]['parcels']}")

if len(sys.argv) > 2:
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    latest = {r: max(glob.glob(f'{sys.argv[2]}/{r}_to*_field.csv'), key=lambda p: int(re.search(r'_to(\d+)_', p)[1]))
              for r in RUNS}
    for name, zmax in (('s7_fields.png', None), ('s7_nearfield.png', 120)):
        fig, ax = plt.subplots(len(RUNS), 2, figsize=(14, 2.2 * len(RUNS)), squeeze=False)
        for n, r in enumerate(RUNS):
            F = np.genfromtxt(latest[r], delimiter=',', names=True)
            nz, nr = int(F['i'].max()) + 1, int(F['j'].max()) + 1
            Z, R = F['z_mm'].reshape(nz, nr), F['r_mm'].reshape(nz, nr)
            ms = int(re.search(r'_to(\d+)_', latest[r])[1])
            for c, (key, cmap, lo, hi) in enumerate((('T', 'Oranges', 80, 3500), ('Y_OH', 'Blues', 0, 0.08))):
                a = ax[n, c]
                m = a.pcolormesh(Z, R, F[key].reshape(nz, nr), cmap=cmap, vmin=lo, vmax=hi, shading='nearest')
                for z in (10, 250, 410):
                    a.axvline(z, color='0.3', lw=0.8, ls=':')
                a.set_xlim(0, zmax or Z.max())
                a.set_title(f'{r}, {key}, mean over {ms - 10} to {ms} ms (model output)', fontsize=9)
                a.set_xlabel('z [mm]', fontsize=8)
                a.set_ylabel('r [mm]', fontsize=8)
                fig.colorbar(m, ax=a, label='K' if key == 'T' else 'mass fraction')
        fig.tight_layout()
        fig.savefig(os.path.join(D, name), dpi=110)
