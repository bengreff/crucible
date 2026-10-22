"""Engine database batch b1_tier0 (RULES.md, committed 246d43b before any engine was run): the ideal rocket with NASA
CEA (RocketCEA) for the 27 tier0_isp engines of engines.json at cfd6f12, with the one-at-a-time members of the rules.

No measured performance is read or used: the script reads only the records' inputs (propellants, cycle, operating
point, expansion ratio, throat), never `outcome_pointers`.

Usage: <a RocketCEA venv>/bin/python docs/validation/engine_db/predictions/b1_tier0/tier0_batch.py
Writes predictions.json beside this file and prints the table.
"""
import functools
import hashlib
import importlib.metadata
import json
import math
import os
import platform

from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer

HERE = os.path.dirname(os.path.abspath(__file__))
SPEC = os.path.join(HERE, '..', '..', 'engines.json')
SPEC_SHA = 'f6cda98344ed5a230442eb8a896049c554010ee01e67d8ace5503f8860b14c47'
SPEC_COMMIT, RULES_COMMIT = 'cfd6f12', '246d43b'
G0, PSIA, PA = 9.80665, 6894.757293168, 101325.0

EXCLUDED = {'rs68a_rated', 'rd107_8d74k', 'vinci_a5me_full', 'cuip_psu_shear_coax_case1'}
GROUP_S = {'j2_sa208_504', 'lm_ascent_engine', 'lm_descent_engine', 'lr87_aj_5', 'lr91_aj_5', 'hm7b_rated',
           'aestus_eps', 'viking_4b', 'rd_253', 'h1b', 'rs27a', 'aj10_190_oms'}
THROAT_AS_IS = {'j2_sa208_504', 'le8', 'lm10mira_tca_ue2', 'r4d_apollo_rcs', 'rs25_block2a_npl', 'f1'}
THROAT_LIP = {'vulcain2_rated', 'lm_ascent_engine', 'lm_descent_engine'}
THROAT_ENVELOPE = {'j2x_100pct_pl', 'lr91_aj_5', 'aj10_190_oms', 'viking_4b', 'rd_253'}
PC_SECOND = {'j2_sa208_504': [700 * PSIA / 1e6], 'h1b': [700 * PSIA / 1e6], 'hm7b_rated': [3.5]}  # rule 8, MPa
PC_PM5 = {'rs27a'}
OPEN = {'gas_generator', 'expander_bleed'}
TX = {'dumped_overboard': (0.35, 0.2, 0.5), 'into_nozzle_extension': (0.6, 0.4, 0.8)}  # rule 4; else (0.5, 0.2, 0.8)

# Rule 2's cards (RocketCEA names, or cards added below).
OX = {'O2(L)': 'LOX', 'N2O4': 'N2O4', 'MON-3': 'MON3'}
FUEL = {'H2': 'LH2', 'CH4': 'CH4', 'LNG': 'CH4', 'RP-1': 'RP1', 'kerosene': 'RP1', 'Aerozine-50': 'A50',
        'MMH': 'MMH', 'UDMH': 'UDMH', 'UH25': 'UH25'}
add_new_fuel('UH25', 'fuel C2H8N2(L),UDMH C 2 H 8 N 2 wt%=75. h,cal=11900. t(k)=298.15\n'
                     'fuel N2H4(L) N 2 H 4 wt%=25. h,cal=12100. t(k)=298.15')
add_new_oxidizer('MON3_NO', 'oxid N2O4(L) N 2 O 4 wt%=97. h,cal=-4676. t(k)=298.15\n'
                            'oxid NO N 1 O 1 wt%=3. h,kj/mol=91.271 t(k)=298.15')

_cea = {}


def cea(ox, fuel):
    if (ox, fuel) not in _cea:
        _cea[(ox, fuel)] = CEA_Obj(oxName=ox, fuelName=fuel)
    return _cea[(ox, fuel)]


def val(field):
    return field.get('value') if isinstance(field, dict) else field


def kliegel_levine(r, g):
    big = 1 + r
    a2, a3 = (g + 1) / 96, (g + 1) * (8 * g - 27) / 2304
    a4 = (g + 1) * (754 * g * g - 757 * g + 3633) / 276480
    return 1 - a2 / big**2 + a3 / big**3 - a4 / big**4


def subsonic_mach(cr, g):
    """M at which the isentropic A/A* equals cr, on the subsonic branch."""
    def area(m):
        return ((2 / (g + 1)) * (1 + (g - 1) / 2 * m * m)) ** ((g + 1) / (2 * (g - 1))) / m
    lo, hi = 1e-6, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if area(mid) > cr else (lo, mid)
    return 0.5 * (lo + hi)


def rayleigh_ratio(cr, g):
    """Rule 3: p0 / p_inj for a constant-area chamber heated from rest to the Mach number of area ratio cr."""
    m = subsonic_mach(cr, g)
    return (1 + (g - 1) / 2 * m * m) ** (g / (g - 1)) / (1 + g * m * m)


@functools.lru_cache(maxsize=None)
def gg_ratio(c, pc_psia, t_target):
    """Rule 4: the fuel-rich O/F at which CEA's equilibrium temperature at p_c is t_target (K), by bisection from
    O/F 0.02 up to the first 0.01 step that reaches it (the target is below RULES.md's 1500 K bound, so the step
    is reached first). Returns (O/F, note)."""
    t = lambda mr: c.get_Tcomb(Pc=pc_psia, MR=mr) / 1.8
    lo = 0.02
    if t(lo) >= t_target:
        return lo, f'T({lo}) = {t(lo):.0f} K is above the target; the bisection range\'s end is taken'
    hi = lo
    while t(hi) < t_target:
        lo, hi = hi, hi + 0.01
        if hi > 10:
            raise RuntimeError(f'no O/F up to 10 reaches {t_target} K')
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if t(mid) < t_target else (lo, mid)
    return 0.5 * (lo + hi), ''


def chamber_mr(mr_e, f, r_gg):
    return (mr_e / (1 + mr_e) - f * r_gg / (1 + r_gg)) / (1 / (1 + mr_e) - f / (1 + r_gg))


def isp_vac_at(c, psia, mr, eps, frozen, lam, cstar):
    isp1d = c.get_IvacCstrTc(Pc=psia, MR=mr, eps=eps, frozen=frozen, frozenAtThroat=frozen)[0]
    pe = psia / c.get_PcOvPe(Pc=psia, MR=mr, eps=eps, frozen=frozen, frozenAtThroat=frozen)  # psia
    part = pe * PSIA * eps * cstar / (psia * PSIA * G0)
    return lam * (isp1d - part) + part, pe * PSIA


def evaluate(r, s):
    """All quantities at one set of settings s."""
    c = cea(s['ox'], s['fuel'])
    eps = val(r['geometry']['expansion_ratio'])
    pc = s['pc'] * 1e6
    mr_e = val(r['operating_point']['mixture_ratio'])
    out, diag = {}, {}
    # Rule 4: the chamber mixture ratio.
    f = s['f']
    if f > 0:
        if s['cycle'] == 'expander_bleed':
            r_gg, note = 0.0, ''
        else:
            r_gg, note = gg_ratio(c, pc / PSIA, s['t_gg'])
        mr_c = chamber_mr(mr_e, f, r_gg) if s['mr_basis'] == 'engine' else mr_e
        diag.update({'gg_O/F': r_gg, 'gg_note': note})
    else:
        mr_c = mr_e
    # Rule 3: the nozzle stagnation pressure.
    gamma = c.get_Chamber_MolWt_gamma(Pc=pc / PSIA, MR=mr_c, eps=eps)[1]
    p0 = pc * rayleigh_ratio(s['cr'], gamma) if s['station'] == 'inj' else pc
    psia = p0 / PSIA
    gamma = c.get_Chamber_MolWt_gamma(Pc=psia, MR=mr_c, eps=eps)[1]
    cstar = c.get_Cstar(Pc=psia, MR=mr_c) * 0.3048
    if not (gamma > 1 and cstar > 0):
        raise RuntimeError(f'CEA rejected {r["id"]}: {s}')
    cd = 1.0 if s['cd_r'] is None else kliegel_levine(s['cd_r'], gamma)
    lam = (1 + math.cos(math.radians(s['theta']))) / 2
    frz = s['frozen']
    isp_vac, pe = isp_vac_at(c, psia, mr_c, eps, frz, lam, cstar)
    # Rule 9: sea level, separated where the exit pressure is below s p_a.
    psep = s['s_sep'] * PA
    if pe >= psep:
        eps_sl, isp_cut = eps, isp_vac
    else:
        lo, hi = 1.0001, eps
        for _ in range(60):
            mid = math.sqrt(lo * hi)
            p_mid = isp_vac_at(c, psia, mr_c, mid, frz, lam, cstar)[1]
            lo, hi = (mid, hi) if p_mid > psep else (lo, mid)
        eps_sl = math.sqrt(lo * hi)
        isp_cut = isp_vac_at(c, psia, mr_c, eps_sl, frz, lam, cstar)[0]
    isp_sl = isp_cut - PA * eps_sl * cstar / (cd * p0 * G0)
    tx = s['tx'] * isp_vac
    out['isp_vacuum_chamber'] = isp_vac
    out['isp_sea_level_chamber'] = isp_sl
    out['isp_vacuum_engine'] = (1 - f) * isp_vac + f * tx
    out['isp_sea_level_engine'] = (1 - f) * isp_sl + f * tx
    out['c_star'] = cstar * (pc / p0) / cd
    dt = s['dt']
    if dt is not None:
        at = math.pi * dt * dt / 4
        mdot = cd * p0 * at / cstar
        out['mass_flow_chamber'] = mdot
        out['mass_flow_engine'] = mdot / (1 - f)
        for k in ('vacuum', 'sea_level'):
            fch = out[f'isp_{k}_chamber'] * G0 * mdot / 1e3
            out[f'thrust_{k}_chamber'] = fch
            out[f'thrust_{k}_engine'] = fch + f * (mdot / (1 - f)) * G0 * tx / 1e3
    diag.update({'chamber_O/F': mr_c, 'p0_MPa': p0 / 1e6, 'gamma_c': gamma, 'Cd': cd, 'lambda': lam,
                 'c*_CEA_m/s': cstar, 'Tc_K': c.get_Tcomb(Pc=psia, MR=mr_c) / 1.8, 'p_e_kPa': pe / 1e3,
                 'separated_at_sea_level': bool(pe < psep), 'eps_sea_level': float(eps_sl)})
    return out, diag


def setup(r):
    """The central settings and the members (label, overrides) that apply to this engine under the rules."""
    op, g, cyc = r['operating_point'], r['geometry'], r['cycle']
    ox, fuel = val(r['propellants']['oxidizer']), val(r['propellants']['fuel'])
    cycle = val(cyc['type'])
    station = val(op['pressure_station'])
    cr_printed = val(g.get('contraction_ratio'))
    basis = op['mixture_ratio_basis']
    dt = val(g.get('throat_diameter_m'))
    central = dict(ox=OX[ox], fuel=FUEL[fuel], pc=val(op['chamber_pressure_mpa']), frozen=0, theta=10.0, cd_r=1.0,
                   s_sep=0.35, cycle=cycle, f=0.0, t_gg=900.0, mr_basis='engine', tx=0.0, station='iac',
                   cr=cr_printed or 3.0, dt=None)
    members = [('frozen at the throat', dict(frozen=1)), ('theta 5 deg', dict(theta=5.0)),
               ('theta 15 deg', dict(theta=15.0)), ('Cd at R/r_t 0.5', dict(cd_r=0.5)), ('Cd = 1', dict(cd_r=None))]
    if fuel == 'kerosene':
        members.append(('Jet-A(L) card', dict(fuel='JetA')))
    if ox == 'MON-3':
        members.append(('MON-3 with 3% NO', dict(ox='MON3_NO')))
    if station == 'nozzle_stagnation':
        pass
    elif station == 'injector_face':
        central['station'] = 'inj'
        if not cr_printed:
            members += [('CR 2', dict(cr=2.0)), ('CR 5', dict(cr=5.0))]
    else:
        members.append(('p_c read as injector-face' + ('' if cr_printed else ' at CR 2'),
                        dict(station='inj', cr=cr_printed or 2.0)))
    if cycle in OPEN:
        exhaust = val(cyc['turbine_exhaust'])
        tx_c, tx_lo, tx_hi = TX.get(exhaust, (0.5, 0.2, 0.8))
        central.update(f=0.025, tx=tx_c)
        members += [('f 1%', dict(f=0.01)), ('f 4%', dict(f=0.04)),
                    (f'Isp_tx {tx_lo}', dict(tx=tx_lo)), (f'Isp_tx {tx_hi}', dict(tx=tx_hi))]
        if cycle == 'gas_generator':
            members += [('T_gg 800 K', dict(t_gg=800.0)), ('T_gg 1050 K', dict(t_gg=1050.0))]
        if basis.get('provenance') == 'printed' and val(basis) in ('engine', 'thrust_chamber'):
            central['mr_basis'] = 'engine' if val(basis) == 'engine' else 'chamber'
        else:
            members.append(('mixture ratio read as the chamber\'s', dict(mr_basis='chamber')))
    if dt is not None:
        rid = r['id']
        if rid in THROAT_AS_IS:
            central['dt'] = dt
        elif rid in THROAT_LIP:
            central['dt'] = dt * 0.985
            members += [('inner exit 0.97 D', dict(dt=dt * 0.97)), ('inner exit 1.00 D', dict(dt=dt))]
        elif rid in THROAT_ENVELOPE:
            central['dt'] = dt * 0.925
            members += [('inner exit 0.85 D', dict(dt=dt * 0.85)), ('inner exit 1.00 D', dict(dt=dt))]
        else:
            raise RuntimeError(f'{rid}: a throat with no rule 7 category')
    for pc2 in PC_SECOND.get(r['id'], []):
        members.append((f'second printed p_c {pc2:.4g} MPa', dict(pc=pc2)))
    if r['id'] in PC_PM5:
        members += [('p_c -5%', dict(pc=central['pc'] * 0.95)), ('p_c +5%', dict(pc=central['pc'] * 1.05))]
    return central, members


def predict(r):
    central, members = setup(r)
    c_out, c_diag = evaluate(r, central)
    sep = [evaluate(r, dict(central, s_sep=v))[1]['separated_at_sea_level'] for v in (0.25, 0.45)]
    if c_diag['separated_at_sea_level'] or any(sep):
        members += [('s 0.25', dict(s_sep=0.25)), ('s 0.45', dict(s_sep=0.45))]
    devs = []
    for label, over in members:
        m_out, _ = evaluate(r, dict(central, **over))
        devs.append((label, {k: m_out[k] - c_out[k] for k in c_out}))
    quantities = {}
    for k, v in c_out.items():
        d = {label: dv[k] for label, dv in devs}
        minus = sum(-x for x in d.values() if x < 0)
        plus = sum(x for x in d.values() if x > 0)
        quantities[k] = {'central': v, 'band': [v - minus, v + plus], 'minus': minus, 'plus': plus,
                         'quadrature': math.sqrt(sum(x * x for x in d.values())),
                         'members': {label: x for label, x in d.items() if x != 0}}
    settings = {k: v for k, v in central.items()}
    return quantities, c_diag, settings


UNITS = {'isp': 's', 'c_star': 'm/s', 'mass_flow': 'kg/s', 'thrust': 'kN'}


def main():
    raw = open(SPEC, 'rb').read()
    sha = hashlib.sha256(raw).hexdigest()
    if sha != SPEC_SHA:
        raise SystemExit(f'engines.json is {sha}, not the batch\'s {SPEC_SHA}')
    engines = [r for r in json.loads(raw)['engines'] if r['id'] not in EXCLUDED]
    assert len(engines) == 27, len(engines)
    assert GROUP_S <= {r['id'] for r in engines}
    out = {'batch': 'b1_tier0', 'rules': 'RULES.md', 'rules_commit': RULES_COMMIT,
           'spec': {'file': 'docs/validation/engine_db/engines.json', 'commit': SPEC_COMMIT, 'sha256': sha},
           'crucible_commit': None, 'table_hashes': None,
           'cea': f'RocketCEA {importlib.metadata.version("rocketcea")} (NASA Glenn CEA2, as bundled); '
                  f'Python {platform.python_version()}',
           'script_sha256': hashlib.sha256(open(__file__, 'rb').read()).hexdigest(),
           'units': {'isp_*': 's', 'c_star': 'm/s (p_c A_t / mdot)', 'mass_flow_*': 'kg/s', 'thrust_*': 'kN'},
           'band': 'central minus and plus the sums of the one-at-a-time members\' negative and positive deviations '
                   '(RULES.md, The method); quadrature beside, not used',
           'engines': []}
    print(f'b1_tier0: {out["cea"]}; engines.json {sha[:8]}...{sha[-4:]}')
    print(f'{"engine":38s} {"grp":3s} {"O/F_c":>6s} {"Isp_vac,eng [s]":>26s} {"Isp_SL,eng [s]":>26s} '
          f'{"c* [m/s]":>24s} {"F_vac,eng [kN]":>28s}')
    for r in engines:
        q, diag, settings = predict(r)
        group = 'S' if r['id'] in GROUP_S else 'P'
        out['engines'].append({'id': r['id'], 'family': r['family'], 'group': group, 'settings': settings,
                               'diagnostics': diag, 'quantities': q})
        fmt = lambda k, w: (f'{q[k]["central"]:8.1f} [{q[k]["band"][0]:7.1f},{q[k]["band"][1]:7.1f}]'
                            if k in q else '-').rjust(w)
        sep = '*' if diag['separated_at_sea_level'] else ' '
        print(f'{r["id"]:38s} {group:3s} {diag["chamber_O/F"]:6.3f} {fmt("isp_vacuum_engine", 26)} '
              f'{fmt("isp_sea_level_engine", 25)}{sep} {fmt("c_star", 24)} {fmt("thrust_vacuum_engine", 28)}'
              + (f'  gg: {diag["gg_note"]}' if diag.get('gg_note') else ''))
    with open(os.path.join(HERE, 'predictions.json'), 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=False)
        fh.write('\n')
    print('* separated at sea level (rule 9). Wrote predictions.json')


if __name__ == '__main__':
    main()
