"""One NASA CEA ideal-rocket point (via RocketCEA, which compiles the NASA Glenn CEA2 Fortran), the
reference column of the engine's printouts. The reactants are explicit CEA cards, so their state is
pinned rather than taken from RocketCEA's propellant library; any propellant pair CEA's thermo.lib
knows can be given. Infinite-area combustor at chamber pressure pc, vacuum Isp at area ratio eps.

Usage: cea_point.py '<fuel card>' '<oxidizer card>' <O/F> <pc Pa> <eps>
  e.g. cea_point.py 'fuel H2 H 2 wt%=100. t(k)=300' 'oxid O2 O 2 wt%=100. t(k)=300' 5 3.116e6 10
Prints one line, SI: rocketcea <version> Tc <K> cstar <m/s> Ivac_shifting <s> Ivac_frozen <s> Ivac_frozen_throat <s>
(frozen at the chamber, and frozen from the throat).
Setup: python3 -m venv .venv-cea && .venv-cea/bin/pip install rocketcea
"""
import sys
import importlib.metadata
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer

PSIA = 6894.757293168  # Pa


def point(fuel_card, ox_card, of, pc, eps):
    add_new_fuel('FUEL_CARD', fuel_card)
    add_new_oxidizer('OX_CARD', ox_card)
    c = CEA_Obj(oxName='OX_CARD', fuelName='FUEL_CARD')
    psia = pc / PSIA
    ivac, cstar, tc = c.get_IvacCstrTc(Pc=psia, MR=of, eps=eps, frozen=0, frozenAtThroat=0)
    frozen = c.get_IvacCstrTc(Pc=psia, MR=of, eps=eps, frozen=1, frozenAtThroat=0)[0]
    throat = c.get_IvacCstrTc(Pc=psia, MR=of, eps=eps, frozen=1, frozenAtThroat=1)[0]
    return float(tc) / 1.8, float(cstar) * 0.3048, float(ivac), float(frozen), float(throat)


if __name__ == '__main__':
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    tc, cstar, ivac, frozen, throat = point(sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
    print(f"rocketcea {importlib.metadata.version('rocketcea')} Tc {tc:.2f} cstar {cstar:.3f} "
          f"Ivac_shifting {ivac:.3f} Ivac_frozen {frozen:.3f} Ivac_frozen_throat {throat:.3f}")
