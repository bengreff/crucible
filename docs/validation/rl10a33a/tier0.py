"""RL10A-3-3A tier 0 (docs/validation/RL10A33A_PROTOCOL.md, section 3 as amended 13 October 2026): the ideal
rocket with NASA CEA (RocketCEA), and the one-at-a-time brackets of section 4 that reach it.

No measured performance is read or used. Inputs: Table E1 (tm107318_table_E1.csv), P1's and P2's flows and O/F (section 2),
the plenum temperatures (tm107318_inputs_addendum.md), the radius reading of section 1's size question, Figure E1's
injector-face and exit stations.

Method, as the protocol states it:
- Wall: r = r_t sqrt(A/A_t) through Table E1's nodes, the throat (0, r_t), the face node (face, r_c) and the exit
  node (exit, r_t sqrt 61), by PCHIP (SciPy's, the same interpolant as core/contour).
- p_c solves mdot = Cd p_c A_t / c*_CEA(p_c) by fixed point (infinite-area combustor; tier 0 has no Rayleigh loss,
  so the injector-face pressure equals p_c).
- Cd: Kliegel-Levine (Johnson and Wright's coefficients, as tests/chamber_study.cpp) at the CEA chamber gamma.
- Isp_vac = lambda (Isp_CEA - P) + P with P = p_e eps c*_CEA / (p_c g0), shifting and frozen (as the engine's C1
  comparison); F_vac = Isp_vac g0 mdot; Isp on the engine-inlet flow = F_vac / (mdot_inlet g0).
- lambda: the mean of (1 + cos theta)/2 at the PCHIP exit wall angle and at the throat-to-exit cone half-angle;
  those two are the ends of its bracket.

Usage: .venv-cea/bin/python docs/validation/rl10a33a/tier0.py
"""
import csv
import importlib.metadata
import math
import os

import numpy as np
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer
from scipy.interpolate import PchipInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
IN, G0, PSIA, RANKINE = 0.0254, 9.80665, 6894.757293168, 5.0 / 9.0
EPS = 61.0
# Section 2's operating points: (name, chamber flow kg/s, engine-inlet flow kg/s, combustion O/F).
P1 = ('P1, Table 6.1.1', 16.7126, 2.7946 + 13.9525, 5.055)
P2 = ('P2, Table 2.5.1 (all burnt)', 16.9462, 16.9462, 5.26)
T_FUEL_R, T_OX_R = 359.4, 182.1  # Table 6.1.1, Fuel and Oxid. Injector Plenum (addendum part A)


def table_e1():
    rows = [line for line in open(os.path.join(HERE, 'tm107318_table_E1.csv')) if not line.startswith('#')]
    data = list(csv.reader(rows))[1:]
    return [(float(r[0]), float(r[1])) for r in data]


def o2_liquid_enthalpy(t):
    """Saturated-liquid h (kJ/kg) from the NIST table, linear in T between rows."""
    rows = [line.split('\t') for line in open(os.path.join(ROOT, 'tools', 'o2_saturation_nist.tsv'))
            if not line.startswith('#')][1:]  # after the column header; column 6 is "Enthalpy (l, kJ/kg)"
    ts = [float(r[0]) for r in rows]
    hs = [float(r[5]) for r in rows]
    return float(np.interp(t, ts, hs))


def cards(t_fuel_r=T_FUEL_R, t_ox_r=T_OX_R):
    """CEA cards at the plenum states. CEA's O2(L) is defined at 90.17 K only (h = -12.979 kJ/mol), so the LOX card
    carries its enthalpy at the plenum temperature: CEA's value plus the NIST saturated-liquid difference."""
    t_ox = t_ox_r * RANKINE
    h_ox = -12.979 + 0.0319988 * (o2_liquid_enthalpy(t_ox) - o2_liquid_enthalpy(90.17))
    return (f'fuel H2 H 2 wt%=100. t(k)={t_fuel_r * RANKINE:.3f}',
            f'oxid LOXP O 2 wt%=100. h,kj/mol={h_ox:.4f} t(k)={t_ox:.3f}')


def wall(rt, rc, face, exit_):
    """PCHIP wall (inches) through the nodes."""
    nodes = sorted([(z, rt * math.sqrt(e)) for z, e in table_e1()] + [(0.0, rt), (face, rc), (exit_, rt * math.sqrt(EPS))])
    return PchipInterpolator([p[0] for p in nodes], [p[1] for p in nodes])


def throat_curvature(p, rt):
    """R / r_t: the mean of the two one-sided radii of curvature at the throat node."""
    k = list(p.x).index(0.0)
    h = p.x[k] - p.x[k - 1]
    left = 6 * p.c[0, k - 1] * h + 2 * p.c[1, k - 1]
    right = 2 * p.c[1, k]
    return 0.5 * (1 / left + 1 / right) / rt


def kliegel_levine(r, g):
    big = 1 + r
    a2, a3 = (g + 1) / 96, (g + 1) * (8 * g - 27) / 2304
    a4 = (g + 1) * (754 * g * g - 757 * g + 3633) / 276480
    return 1 - a2 / big**2 + a3 / big**3 - a4 / big**4


def lambdas(p, rt, exit_):
    exit_angle = math.atan(float(p(exit_, 1)))
    cone = math.atan((rt * math.sqrt(EPS) - rt) / exit_)
    return (1 + math.cos(exit_angle)) / 2, (1 + math.cos(cone)) / 2, exit_angle, cone


_objects = {}


def cea(fuel, ox):
    key = (fuel, ox)
    if key not in _objects:
        name = f'P{len(_objects)}'
        add_new_fuel('F' + name, fuel)
        add_new_oxidizer('O' + name, ox)
        _objects[key] = CEA_Obj(oxName='O' + name, fuelName='F' + name)
    return _objects[key]


def predict(op, rt_in, curvature, lam, fuel, ox):
    """The tier-0 quantities at one set of inputs."""
    _, MDOT, MDOT_INLET, OF = op
    c = cea(fuel, ox)
    at = math.pi * (rt_in * IN) ** 2
    pc, cd = 3e6, 1.0
    for _ in range(100):
        psia = pc / PSIA
        gamma = c.get_Chamber_MolWt_gamma(Pc=psia, MR=OF, eps=EPS)[1]
        cd = kliegel_levine(curvature, gamma)
        cstar = c.get_Cstar(Pc=psia, MR=OF) * 0.3048
        if not (gamma > 1 and cstar > 0):  # CEA returns zeros when it rejects the dataset
            raise RuntimeError(f'CEA failed for "{fuel}" / "{ox}": gamma {gamma}, c* {cstar}')
        new = MDOT * cstar / (cd * at)
        if abs(new / pc - 1) < 1e-12:
            pc = new
            break
        pc = new
    psia = pc / PSIA
    out = {'p_c [MPa]': pc / 1e6, 'c* [m/s]': pc * at / MDOT, 'Cd': cd, 'gamma_c': gamma,
           'c*_CEA [m/s]': c.get_Cstar(Pc=psia, MR=OF) * 0.3048, 'Tc [K]': c.get_Tcomb(Pc=psia, MR=OF) / 1.8}
    for label, frozen in (('shifting', 0), ('frozen', 1)):
        isp1d = c.get_IvacCstrTc(Pc=psia, MR=OF, eps=EPS, frozen=frozen, frozenAtThroat=0)[0]
        pe = pc / c.get_PcOvPe(Pc=psia, MR=OF, eps=EPS, frozen=frozen, frozenAtThroat=0)
        part = pe * EPS * out['c*_CEA [m/s]'] / (pc * G0)
        isp = lam * (isp1d - part) + part
        out[f'Isp_vac 1-D CEA {label} [s]'] = isp1d
        out[f'Isp_vac {label}, chamber flow [s]'] = isp
        out[f'Isp_vac {label}, engine-inlet flow [s]'] = isp * MDOT / MDOT_INLET
        out[f'F_vac {label} [kN]'] = isp * G0 * MDOT / 1e3
    return out


def case(op, rt, rc, face, exit_, fuel, ox, curvature=None, lam=None):
    p = wall(rt, rc, face, exit_)
    r = throat_curvature(p, rt)
    hi, lo, _, _ = lambdas(p, rt, exit_)
    return predict(op, rt, r if curvature is None else curvature, 0.5 * (hi + lo) if lam is None else lam, fuel, ox)


def report(op):
    rt, rc, face, exit_ = 2.47, 5.13, -13.1, 46.0  # the radius reading; Figure E1's stations
    fuel, ox = cards()
    p = wall(rt, rc, face, exit_)
    r_central = throat_curvature(p, rt)
    hi, lo, exit_angle, cone = lambdas(p, rt, exit_)
    print(f'\n=== {op[0]}: O/F {op[3]}, mdot {op[1]} kg/s (engine inlet {op[2]:.4f}), eps {EPS}')
    print(f'cards: "{fuel}" / "{ox}"')
    print(f'wall: r_t {rt} in, r_c {rc} in, face {face} in, exit {exit_} in; throat R/r_t {r_central:.4f} (PCHIP); '
          f'exit wall angle {math.degrees(exit_angle):.3f} deg (lambda {hi:.5f}), cone {math.degrees(cone):.3f} deg '
          f'(lambda {lo:.5f}); lambda central {0.5 * (hi + lo):.5f}')
    # Central Cd: Kliegel-Levine at the series' edge R/r_t = 0.5, the PCHIP throat being sharper (section 4).
    central = case(op, rt, rc, face, exit_, fuel, ox, curvature=0.5)
    variations = [
        ('lambda at the exit wall angle', dict(curvature=0.5, lam=hi)),
        ('lambda at the cone half-angle', dict(curvature=0.5, lam=lo)),
        ('exit station 45.29 in', dict(exit_=45.29, curvature=0.5)),
        ('exit station 46.40 in', dict(exit_=46.40, curvature=0.5)),
        (f'Cd at the PCHIP R/r_t {r_central:.4f}', dict(curvature=r_central)),
        ('sharp throat (Cd at R/r_t 0.5)', dict(curvature=0.5)),
        ('injector face -13.4 in', dict(face=-13.4, curvature=0.5)),
        ('injector face -12.8 in', dict(face=-12.8, curvature=0.5)),
        ('fuel plenum T -0.1 R', dict(cards=cards(T_FUEL_R - 0.1, T_OX_R), curvature=0.5)),
        ('fuel plenum T +0.1 R', dict(cards=cards(T_FUEL_R + 0.1, T_OX_R), curvature=0.5)),
        ('oxidizer plenum T -0.1 R', dict(cards=cards(T_FUEL_R, T_OX_R - 0.1), curvature=0.5)),
        ('oxidizer plenum T +0.1 R', dict(cards=cards(T_FUEL_R, T_OX_R + 0.1), curvature=0.5)),
    ]
    results = []
    for label, v in variations:
        f, o = v.pop('cards', (fuel, ox))
        args = dict(rt=rt, rc=rc, face=face, exit_=exit_)
        for key in ('face', 'exit_'):
            if key in v:
                args[key] = v.pop(key)
        results.append((label, case(op, args['rt'], rc, args['face'], args['exit_'], f, o, **v)))
    keys = [k for k in central if not k.startswith(('gamma', 'Tc', 'c*_CEA', 'Isp_vac 1-D'))]
    print('\nCentral and one-at-a-time deviations')
    for k in central:
        print(f'  {k:42s} {central[k]:.6g}')
    for label, r in results:
        print(f'  {label:42s} ' + ', '.join(f'{k.split(" [")[0]} {r[k] - central[k]:+.4g}' for k in keys))
    print('\nJudged band (central, minus and plus the sums of the one-at-a-time deviations on each side; quadrature '
          'beside, not judged)')
    for k in keys:
        devs = [r[k] - central[k] for _, r in results]
        lo_sum = sum(-d for d in devs if d < 0)
        hi_sum = sum(d for d in devs if d > 0)
        quad = math.sqrt(sum(d * d for d in devs))
        print(f'  {k:42s} {central[k]:.6g}  band [{central[k] - lo_sum:.6g}, {central[k] + hi_sum:.6g}]  '
              f'(-{lo_sum:.3g} / +{hi_sum:.3g}; quadrature {quad:.3g})')
    record = case(op, 1.235, 2.565, face, exit_, fuel, ox, curvature=0.5)
    print('\nFor the record, unjudged: the diameter reading (r_t 1.235 in, r_c 2.565 in), Cd at R/r_t 0.5, lambda '
          'central at its own wall')
    for k in keys:
        print(f'  {k:42s} {record[k]:.6g}')


def main():
    print(f'RL10A-3-3A tier 0: RocketCEA {importlib.metadata.version("rocketcea")} (NASA Glenn CEA2); SciPy PCHIP')
    for op in (P1, P2):
        report(op)


if __name__ == '__main__':
    main()
