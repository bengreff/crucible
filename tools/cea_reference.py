"""Reference ideal-rocket performance from NASA CEA (via RocketCEA, which compiles the NASA
Glenn CEA2 Fortran). Gaseous reactants at 298.15 K, defined by explicit CEA reactant cards so the
reactant state is pinned, not inherited from RocketCEA's propellant library.
Infinite-area combustor, vacuum Isp at the given area ratio, equilibrium and frozen-at-chamber.
Output: tools/cea_reference.csv (SI units). Usage: python tools/cea_reference.py
"""
import csv, os
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer

add_new_fuel('H2_G298', 'fuel H2 H 2 wt%=100. t(k)=298.15')
add_new_fuel('CH4_G298', 'fuel CH4 C 1 H 4 wt%=100. t(k)=298.15')
add_new_oxidizer('O2_G298', 'oxid O2 O 2 wt%=100. t(k)=298.15')
# Cryogenic liquids at their normal boiling points; CEA thermo.lib assigns H2(L) -9012 J/mol at
# 20.27 K and O2(L) -12979 J/mol at 90.17 K (McBride 1996). The core is given the same enthalpies.
add_new_fuel('H2_L', 'fuel H2(L) H 2 wt%=100. t(k)=20.27')
add_new_oxidizer('O2_L', 'oxid O2(L) O 2 wt%=100. t(k)=90.17')

PSIA = 6894.757293168  # Pa
CASES = [  # (pair, oxidizer card, fuel card, O/F list, Pc bar, area ratio)
    ('H2/O2', 'O2_G298', 'H2_G298', [4.0, 6.0, 8.0], 20.0, 40.0),        # Penn State GO2/GH2 class
    ('CH4/O2', 'O2_G298', 'CH4_G298', [2.6, 3.4, 4.0], 20.0, 40.0),      # TUM GOX/GCH4 class (O/F 2.6 reference point)
    ('LH2/LOX', 'O2_L', 'H2_L', [4.5, 5.5, 6.5], 32.75, 61.0),           # RL10A-3-3A class (Pc about 475 psia, area ratio 61)
]

def main():
    rows = []
    for pair, ox, fuel, ofs, PC_BAR, EPS in CASES:
        c = CEA_Obj(oxName=ox, fuelName=fuel)
        pc = PC_BAR * 1e5 / PSIA
        for of in ofs:
            for frozen in (0, 1):
                ivac, cstar, tc = c.get_IvacCstrTc(Pc=pc, MR=of, eps=EPS, frozen=frozen, frozenAtThroat=0)
                rows.append(dict(pair=pair, of=of, pc_bar=PC_BAR, eps=EPS, mode='frozen' if frozen else 'equilibrium',
                                 tc_K=float(tc) / 1.8, cstar_m_s=float(cstar) * 0.3048, isp_vac_s=float(ivac)))
    path = os.path.join(os.path.dirname(__file__), 'cea_reference.csv')
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    for r in rows: print(r)

if __name__ == '__main__':
    main()
