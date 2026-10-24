"""The face inlet's numbers for the Director's 04:20 decision: dh and the same-enthalpy inlet temperature at the
plate's 199.67 K (LOW_T_THERMO.md, the Director's decision of 24 October, 04:20), using para_supply_bound.py's own functions (its source up to the performance function).
Usage: ~/venvs/crucible/bin/python docs/evidence/low_t/face_inlet_2026-10-24.py"""
import os
__file__ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../tools/mechanisms/para_supply_bound.py")
src = open(__file__).read()
exec(src[:src.index("def performance")])
ts = 199.67
for xp in (0.998, 0.95):
    target = h_h2_eq(ts) - dh(ts, xp)
    lo, hi = 20.0, ts
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if h_h2_eq(mid) > target:
            hi = mid
        else:
            lo = mid
    print("T_s %.2f K, x_p %.3f: dh %.2f kJ/kg; same-enthalpy inlet %.2f K (%.2f K colder, density %.2f%% higher)"
          % (ts, xp, dh(ts, xp) / 1e3, mid, ts - mid, 100 * (ts / mid - 1)))
zp, zo = rotation(ts, range(0, 20, 2), 1)[0], rotation(ts, range(1, 20, 2), 3)[0]
print("equilibrium para fraction at %.2f K: %.4f" % (ts, zp / (zp + zo)))
# On the engine's own continuation below 200 K (make_low_t.py: normal H2, 3:1 frozen, joined at 200 K).
def h_normal(T):
    gas.TPX = 200.0, ct.one_atm, "H2:1"
    return gas.enthalpy_mass + (2.5 * (T - 200.0) + e_rot(T, 0.25) - e_rot(200.0, 0.25)) * R / M["H2"]
for xp in (0.998, 0.95):
    target = h_h2_eq(ts) - dh(ts, xp)
    lo, hi = 20.0, ts
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if h_normal(mid) > target:
            hi = mid
        else:
            lo = mid
    print("on the normal-H2 continuation, x_p %.3f: inlet %.2f K (%.2f K colder, density %.2f%% higher)"
          % (xp, mid, ts - mid, 100 * (ts / mid - 1)))
