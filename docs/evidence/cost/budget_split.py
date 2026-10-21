#!/usr/bin/env python3
"""The corner's formation (COST.md, lever 3, readout (m)): a cell's specific excess split by how it changed.

usage: budget_split.py <budget.csv> <t0 s> <t1 s> <i>:<j> [<i>:<j> ...]

With beta = b + k = (rho E V + p V - h_LOX rho V - (h_H2 - h_LOX) rho f_H V) / rho V (the probe's inviscid energy
flux carries rho k, so beta is the measure its faces deliver exactly), each record interval of a cell gives, with M
its mass and beta_0 its beta at the interval's start, exactly

    M_1 (beta_1 - beta_0) = sum_faces m (beta_f - beta_0) + sum of the non-convective items - beta_0 m_spray,

where m is the face's inviscid mass gain and beta_f its b + k per unit mass. A face with m < 0 is an outflow: it adds
|m| (beta_0 - beta_f), which is positive when the face value lies below the cell's. A face with m > 0 is an inflow:
it adds m (beta_f - beta_0). Each is split into its enthalpy part (h0 + k: F_E / F_m against the cell's) and its
mixture part ((h_H2 - h_LOX) times f_H's face value against the cell's). Dividing by M_1 and summing over intervals
gives the change of beta in J/kg by mechanism. b = beta - k; k's change is printed beside it.

For the first named cell, the intervals in which its b is above that of every named neighbour across a face at both
ends (its local-maximum intervals, L) are summed apart, as registered in COST.md.
"""
import sys

from budget_read import FACES, f, load, s

MECHANISMS = ["outflow, enthalpy part", "outflow, mixture part", "inflow, enthalpy part", "inflow, mixture part",
              "molecular diffusion", "turbulent diffusion", "viscous work and unsplit", "V dp", "spray exchange",
              "sources that should be zero"]


def beta_of(r):
    return (r["rhoEV"] + r["pV"] - r["h_LOX"] * r["rhoV"] - (r["h_H2"] - r["h_LOX"]) * r["rhofHV"]) / r["rhoV"]


def interval(p, q):
    """Each mechanism's J (as M_1 times its J/kg) for one record interval, the closure residual, and M_1."""
    df = lambda k, t: f(q, k, t) - f(p, k, t)
    ds = lambda t: s(q, t) - s(p, t)
    hl, dh = 0.5 * (p["h_LOX"] + q["h_LOX"]), 0.5 * (p["h_H2"] - p["h_LOX"] + q["h_H2"] - q["h_LOX"])
    m0, m1 = p["rhoV"], q["rhoV"]
    b0, f0 = beta_of(p), p["rhofHV"] / m0
    h0 = (p["rhoEV"] + p["pV"]) / m0  # h0 + k of the cell
    out = dict.fromkeys(MECHANISMS, 0.0)
    faces = {}
    for k in range(4):
        m = df(k, "inviscid_mass")
        side = "outflow" if m < 0 else "inflow"
        enthalpy, mixture = df(k, "inviscid_energy") - m * h0, -dh * (df(k, "inviscid_tracked") - m * f0)  # m (h0_f - h0), -dh m (f_H,f - f_H)
        out[side + ", enthalpy part"] += enthalpy
        out[side + ", mixture part"] += mixture
        faces[(FACES[k], side)] = enthalpy + mixture
        parts = {t: df(k, t) for t in ("conduction", "species_enthalpy", "turbulent_conduction", "turbulent_species_enthalpy",
                                       "turbulent_k_diffusion", "tracked_diffusion", "turbulent_tracked_diffusion")}
        out["molecular diffusion"] += parts["conduction"] + parts["species_enthalpy"] - dh * parts["tracked_diffusion"]
        out["turbulent diffusion"] += (parts["turbulent_conduction"] + parts["turbulent_species_enthalpy"]
                                       + parts["turbulent_k_diffusion"] - dh * parts["turbulent_tracked_diffusion"])
        # viscous work, and whatever a wall function or an open end books beyond the parts (energy and tracked)
        out["viscous work and unsplit"] += (df(k, "viscous_energy") - sum(parts[t] for t in (
            "conduction", "species_enthalpy", "turbulent_conduction", "turbulent_species_enthalpy", "turbulent_k_diffusion"))
            - dh * (df(k, "tracked_transport") - parts["tracked_diffusion"] - parts["turbulent_tracked_diffusion"]))
    out["V dp"] = q["pV"] - p["pV"]
    out["spray exchange"] = ds("exchanged_energy") - hl * ds("exchanged_mass") - dh * ds("exchanged_tracked") - b0 * ds("exchanged_mass")
    out["sources that should be zero"] = sum(ds(t) for t in ("igniter", "heating", "body_work", "correction"))
    measured = m1 * (beta_of(q) - b0)
    return out, measured - sum(out.values()), m1, measured, faces


def split(series, t0, t1, keep=lambda p, q: True):
    inside = [r for r in series if t0 <= r["t"] <= t1]
    total = dict.fromkeys(MECHANISMS, 0.0)
    total["faces"] = {}
    change, residual, worst, n = 0.0, 0.0, 0.0, 0
    for p, q in zip(inside, inside[1:]):
        if not keep(p, q):
            continue
        out, res, m1, measured, faces = interval(p, q)
        for key in MECHANISMS:
            total[key] += out[key] / m1
        for key, v in faces.items():
            total["faces"][key] = total["faces"].get(key, 0.0) + v / m1
        change += measured / m1
        residual += res / m1
        worst = max(worst, abs(res) / max(max(abs(v) for v in out.values()), abs(measured), 1e-300))
        n += 1
    return inside, total, change, residual, worst, n


def report(label, total, change, residual, worst, n, dk=None):
    faces = total.pop("faces")
    positive = sum(v for v in total.values() if v > 0)
    print("  %s: %d intervals; beta change %+.6e J/kg (sum of mechanisms %+.6e, residual %+.2e, worst interval %.1e)%s"
          % (label, n, change, sum(total.values()), residual, worst, "" if dk is None else "; k change %+.4e J/kg" % dk))
    for key in MECHANISMS:
        v = total[key]
        print("      %-28s %+.6e J/kg%s" % (key, v, ("  (%.1f%% of the positive)" % (100 * v / positive)) if v > 0 and positive > 0 else ""))
    out_pos = sum(max(0.0, total[k]) for k in MECHANISMS if k.startswith("outflow"))
    in_pos = sum(max(0.0, total[k]) for k in MECHANISMS if k.startswith("inflow"))
    if positive > 0:
        print("      outflow share of the positive %.1f%%, inflow share %.1f%%" % (100 * out_pos / positive, 100 * in_pos / positive))
    print("      per face: " + ", ".join("%s %s %+.4e" % (k[0], k[1], v) for k, v in sorted(faces.items()) if v != 0))


def main():
    path, t0, t1 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    cells = [tuple(int(x) for x in spec.split(":")) for spec in sys.argv[4:]]
    rows = load(path)
    for n_cell, cell in enumerate(cells):
        inside, total, change, residual, worst, n = split(rows[cell], t0, t1)
        a, z = inside[0], inside[-1]
        excess = lambda r: beta_of(r) - r["rhokV"] / r["rhoV"]
        print("cell i %d j %d: t %.6f to %.6f ms, T %.0f to %.0f K, b %.4g to %.4g J/kg (probe's excess %.4g to %.4g), mass %.4g to %.4g kg"
              % (cell[0], cell[1], a["t"] * 1e3, z["t"] * 1e3, a["T"], z["T"], excess(a), excess(z), a["excess"], z["excess"],
                 a["rhoV"], z["rhoV"]))
        dk = z["rhokV"] / z["rhoV"] - a["rhokV"] / a["rhoV"]
        report("whole window", total, change, residual, worst, n, dk)
        if n_cell == 0:
            neighbours = [rows[c] for c in ((cell[0] - 1, cell[1]), (cell[0] + 1, cell[1]), (cell[0], cell[1] - 1), (cell[0], cell[1] + 1))
                          if c in rows]
            by_step = [{r["step"]: r for r in other} for other in neighbours]
            top = lambda r: all(r["step"] in b and excess(r) > excess(b[r["step"]]) for b in by_step)
            inside_l, total_l, change_l, residual_l, worst_l, n_l = split(rows[cell], t0, t1, lambda p, q: top(p) and top(q))
            first = next((r["t"] for r in inside if top(r)), None)
            print("  local maximum over its %d named face neighbours first at %s" % (
                len(neighbours), "%.6f ms" % (first * 1e3) if first is not None else "no record"))
            report("local-maximum intervals (L)", total_l, change_l, residual_l, worst_l, n_l)


if __name__ == "__main__":
    main()
