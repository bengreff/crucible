# Off-line distance of element ratios in counterflow diffusion flames with mixture-averaged
# diffusion: how far differential diffusion moves the reacting elements off the two-stream line
# z(f) = z_O + f (z_F - z_O), f from Bilger's coupling function (C 4, H 1, O -2, N 0 valences).
# Sizes the off-line range a mixture-fraction Table A must carry (docs/evidence/TABLE_A_PAIRS.md).
import sys
import numpy as np
import cantera as ct

MECH = "gri30_highT.yaml"
VALENCE = {"C": 4.0, "H": 1.0, "O": -2.0, "N": 0.0}


def element_fractions(gas, y, elements):
    out = np.zeros(len(elements))
    for k, name in enumerate(gas.species_names):
        sp = gas.species(name)
        w = gas.molecular_weights[k]
        for m, el in enumerate(elements):
            out[m] += y[k] * sp.composition.get(el, 0.0) * gas.atomic_weight(el) / w
    return out


def run(label, fuel, oxidizer, p, mdot, width=0.02):
    gas = ct.Solution(MECH)
    elements = [e for e in ("C", "H", "O", "N") if e in gas.element_names]
    gas.TPY = 300, p, fuel
    zf = element_fractions(gas, gas.Y, elements)
    gas.TPY = 300, p, oxidizer
    zo = element_fractions(gas, gas.Y, elements)
    reacting = [m for m, el in enumerate(elements) if zf[m] + zo[m] > 0]
    elements = [elements[m] for m in reacting]
    zf, zo = zf[reacting] / zf[reacting].sum(), zo[reacting] / zo[reacting].sum()
    c = np.array([VALENCE[e] / gas.atomic_weight(e) for e in elements])
    bf, bo = c @ zf, c @ zo
    f_st = -bo / (bf - bo)

    flame = ct.CounterflowDiffusionFlame(gas, width=width)
    flame.fuel_inlet.Y = fuel
    flame.fuel_inlet.T = 300
    flame.fuel_inlet.mdot = mdot
    flame.oxidizer_inlet.Y = oxidizer
    flame.oxidizer_inlet.T = 300
    flame.oxidizer_inlet.mdot = mdot
    flame.P = p
    flame.transport_model = "mixture-averaged"
    flame.set_refine_criteria(ratio=3, slope=0.08, curve=0.15, prune=0.02)
    flame.solve(loglevel=0, auto=True)
    worst, worst_at = 0.0, None
    for n in range(flame.flame.n_points):
        z = element_fractions(gas, flame.Y[:, n], [*elements])
        z = z / z.sum()
        f = (c @ z - bo) / (bf - bo)
        d = z - (zo + f * (zf - zo))
        if np.abs(d).max() > worst:
            worst, worst_at = np.abs(d).max(), (f, flame.T[n], d.copy())
    print(f"{label}: p {p/1e6:.1f} MPa, mdot {mdot} kg/m2/s, points {flame.flame.n_points}, T_max {flame.T.max():.0f} K, "
          f"f_st {f_st:.4f}; elements {elements}; max |delta| {worst:.3e} at f {worst_at[0]:.4f}, T {worst_at[1]:.0f} K, "
          f"delta {np.array2string(worst_at[2], precision=3)}")
    sys.stdout.flush()


if __name__ == "__main__":
    for mdot in (0.5, 5.0):
        run("CH4/O2", "CH4:1", "O2:1", 3e6, mdot)
        run("H2/O2", "H2:1", "O2:1", 3e6, mdot)
