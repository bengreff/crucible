#!/usr/bin/env python3
"""The implicit prototype's case C judgements (IMPLICIT_PROTOTYPE.md: criteria 2 to 5 with the 02:20 CDT amendment of
24 October 2026). Each history is a crucible_chamber_study <prefix>_history.csv.

  c_compare.py tau   <explicit>                        tau_C1 and the three steps
  c_compare.py et    <explicit> <dt> <dt/2>            the temporal error E_t
  c_compare.py sweep <explicit> <E_t> <E_t/1000 run> <run>...   criterion 3
  c_compare.py acc   <explicit> <run>...               criteria 2, 4 and 5 against the explicit run
  c_compare.py cost  <explicit> <run>...               M1 to M3 (each history's directory holds <prefix>.log and
                                                       <prefix>_implicit.csv; walls are the march's, from the
                                                       log's "march wall" line or else its "end ... wall" line)

The observables p_inj, the outlet mass flow and F_vac are normalised by the explicit run's mean over [6, 8] ms.
A comparison is made at the coarser run's sample times in the window, with the other run linearly interpolated there."""
import csv, gzip, math, os, re, sys

OBS = [("p_inj", "p_injector"), ("outlet flow", "outlet"), ("F_vac", "F_vac")]


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def sibling(history, suffix):
    """<prefix>_history.csv (or .csv.gz, as kept in the evidence) to <prefix><suffix>."""
    return re.sub(r"_history\.csv(\.gz)?$", suffix, history)


def load(path):
    with opener(path) as f:
        rows = list(csv.DictReader(f))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


def mean(h, col, lo, hi):
    v = [y for t, y in zip(h["t"], h[col]) if lo - 1e-12 <= t <= hi + 1e-12 and math.isfinite(y)]
    return sum(v) / len(v)


def at(h, col, t):
    ts, ys = h["t"], h[col]
    lo, hi = 0, len(ts) - 1
    if t <= ts[0]:
        return ys[0]
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ts[mid] <= t: lo = mid
        else: hi = mid
    w = (t - ts[lo]) / (ts[hi] - ts[lo])
    return ys[lo] + w * (ys[hi] - ys[lo])


def scales(explicit):
    return {c: mean(explicit, c, 6e-3, 8e-3) for _, c in OBS}


def difference(a, b, norm, lo=0.5e-3, hi=8e-3):
    """{observable: (RMS, largest |.|)} of (a - b) / norm at the coarser run's times in [lo, hi]."""
    na = sum(1 for t in a["t"] if lo - 1e-12 <= t <= hi + 1e-12)
    nb = sum(1 for t in b["t"] if lo - 1e-12 <= t <= hi + 1e-12)
    coarse, fine, sign = (a, b, 1) if na <= nb else (b, a, -1)
    times = {t for t in coarse["t"] if lo - 1e-12 <= t <= hi + 1e-12 and t <= fine["t"][-1] + 1e-12}
    out = {}
    for label, c in OBS:
        d = [sign * (y - at(fine, c, t)) / norm[c] for t, y in zip(coarse["t"], coarse[c]) if t in times]
        out[label] = (math.sqrt(sum(x * x for x in d) / len(d)), max(abs(x) for x in d), len(d))
    return out


def show(title, res):
    print(title)
    for label, (rms, big, n) in res.items():
        print(f"  {label:12s} RMS {100 * rms:.5f}%  largest {100 * big:.5f}%  ({n} times)")


cmd = sys.argv[1]
ex = load(sys.argv[2])
norm = scales(ex)
if cmd == "tau":
    settled = mean(ex, "F_vac", 6e-3, 8e-3)
    last = max(t for t, y in zip(ex["t"], ex["F_vac"]) if abs(y / settled - 1) > 1e-3)
    t01 = min(t for t in ex["t"] if t > last)
    tau = t01 / math.log(1000)
    print(f"F_vac settled (mean over [6, 8] ms) {settled:.3f} N; within 0.1% from t_0.1% = {1e3 * t01:.4f} ms")
    print(f"tau_C1 = t_0.1% / ln 1000 = {1e6 * tau:.3f} us")
    for k in (25, 63, 126, 150):
        print(f"  tau_C1/{k} = {tau / k:.6e} s")
elif cmd == "et":
    res = difference(load(sys.argv[3]), load(sys.argv[4]), norm)
    show("dt against dt/2 over [0.5, 8] ms:", res)
    et = max(r[0] for r in res.values())
    print(f"E_t = {et:.4e}; the sweep's tolerances E_t/10, /100, /1000: {et / 10:.3e} {et / 100:.3e} {et / 1000:.3e}")
elif cmd == "sweep":
    et, ref = float(sys.argv[3]), load(sys.argv[4])
    passed = None
    for path in sys.argv[5:]:
        res = difference(load(path), ref, norm)
        show(f"{path} against the E_t/1000 run over [0.5, 8] ms (bar 0.1 E_t = {100 * 0.1 * et:.5f}%):", res)
        ok = all(r[0] < 0.1 * et for r in res.values())
        print(f"  {'passes' if ok else 'fails'}")
        if ok and passed is None: passed = path
    print(f"criterion 3: {'PASS at ' + passed if passed else 'FAIL'}")
elif cmd == "acc":
    for path in sys.argv[3:]:
        h = load(path)
        print(f"== {path}")
        mass = max(abs(x) for x in h["mass_budget"])
        energy = max(abs(x) for x in h["energy_budget"])
        print(f"  the history's budgets (of the initial fill's mass and energy), largest: mass {mass:.2e}, energy {energy:.2e}")
        log = open(sibling(path, ".log")).read()
        m = re.search(r"mass \S+ kg of \S+ kg supplied \((\S+)\); energy \S+ J of \S+ J.*?\((\S+)\), and of \|E\| at "
                      r"the end \S+ J \((\S+)\)", log)
        if m:
            mt, et, ec = (float(x) for x in m.groups())
            print(f"criterion 2 (bar 1e-11 of the throughput): mass {mt:.2e} of the mass supplied, energy {et:.2e} of the"
                  f" supplied total enthalpy and the igniter's energy ({ec:.2e} of |E| at the end)"
                  f"  {'PASS' if max(mt, et) <= 1e-11 else 'FAIL'}")
        else:
            print("criterion 2: the log has no throughput line (a build before 02:38 CDT, 24 October)")
        show("  the first 0.5 ms against the explicit run (reported, not judged):", difference(h, ex, norm, 0, 0.5e-3))
        res = difference(h, ex, norm)
        show("  criterion 4 over [0.5, 8] ms against the explicit run:", res)
        ok = all(r[0] <= 5e-3 and r[1] <= 2e-2 for r in res.values())
        print(f"  criterion 4's bars (RMS 0.5%, largest 2%; judged at tau/63): {'met' if ok else 'not met'}")
        if h["t"][-1] >= 8e-3 - 1e-9:
            worst = 0
            for label, c in (("vacuum Isp", "isp_vac_supply"), ("F_vac", "F_vac")):
                e, i = mean(ex, c, 6e-3, 8e-3), mean(h, c, 6e-3, 8e-3)
                worst = max(worst, abs(i / e - 1))
                print(f"  criterion 5: {label} over [6, 8] ms {i:.6g} against {e:.6g}: {100 * (i / e - 1):+.4f}%")
            print(f"  criterion 5 (bar 0.05%): {'PASS' if worst <= 5e-4 else 'FAIL'}")
elif cmd == "cost":
    def wall(history):
        log = sibling(history, ".log")
        text = open(log).read()
        end = re.findall(r"^end \S+ ms after (\d+) steps.*wall (\d+) s", text, re.M)[-1]
        march = re.findall(r"^march wall ([\d.]+) s", text, re.M)
        return float(march[-1]) if march else float(end[1]), int(end[0])
    ew, es = wall(sys.argv[2])
    print(f"explicit: wall {ew:.1f} s to {1e3 * ex['t'][-1]:.3f} ms, {es} steps, {ew / es:.4e} s per step")
    for path in sys.argv[3:]:
        iw, isteps = wall(path)
        with open(sibling(path, "_implicit.csv")) as f:
            rows = [(float(r["t"]), int(r["iterations"])) for r in csv.DictReader(f)]
        its = sum(i for _, i in rows)
        print(f"== {path}: wall {iw:.1f} s, {len(rows)} steps, {its} inner iterations")
        print(f"  M2: {iw / its:.4e} s per inner iteration = {iw / its / (ew / es):.3f} explicit steps")
        for lo, hi, name in ((0, 0.5e-3, "start [0, 0.5] ms"), (0.5e-3, 3e-3, "fill [0.5, 3] ms"), (3e-3, 8e-3, "settled [3, 8] ms")):
            sel = [i for t, i in rows if lo < t <= hi + 1e-12]
            if sel:
                dt = (hi - lo) / len(sel)
                ratio = dt / (ex["t"][-1] / es)
                print(f"  M1 {name}: {sum(sel) / len(sel):.1f} inner iterations per step (most {max(sel)}),"
                      f" step ratio {ratio:.0f}, cost per step {sum(sel) / len(sel) * iw / its / (ew / es):.1f} explicit steps")
        if abs(rows[-1][0] - ex["t"][-1]) < 1e-9:
            print(f"  M3: explicit wall over implicit wall to {1e3 * ex['t'][-1]:.3f} ms: {ew / iw:.2f}")
