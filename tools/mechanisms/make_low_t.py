"""Lowest NASA-7 ranges, 20 to 200 K, for H2, O2 and H2O, joined to h2o2_highT.yaml's TM-4513 low ranges at 200 K
(docs/evidence/LOW_T_THERMO.md). Each is the least-squares fit of cp/R on a 1 K grid from 20 to 200 K to the
ideal gas's statistical mechanics, with cp(200 K) equal to the TM-4513 range's; a5 and a6 then make h and s
continuous at 200 K. Prints the fit errors and the C++ table for adapters/reaction.cpp.
  H2:  normal hydrogen, ortho and para frozen 3:1, each with its own rotational partition function (even J para,
       odd J ortho with nuclear spin weight 3) on the v = 0 levels B0 J(J+1) - D0 J^2(J+1)^2 + H0 J^3(J+1)^3,
       B0 59.3392, D0 0.04599, H0 5.2e-5 cm^-1 (Huber and Herzberg 1979); vibration (5,987 K) neglected.
  O2:  classical rotor (theta_rot 2.08 K) and a harmonic oscillator at G(1) - G(0) = 1556.385 cm^-1
       (Huber and Herzberg 1979); the ground state's fine structure (about 3 K) and the a1Delta_g state neglected.
  H2O: classical asymmetric rotor and harmonic oscillators at the fundamentals 3657.053, 1594.746 and
       3755.929 cm^-1 (Shimanouchi 1972, NSRDS-NBS 39).
Usage: python3 tools/mechanisms/make_low_t.py"""
import math, os, yaml
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
mech = yaml.safe_load(open(os.path.join(here, "h2o2_highT.yaml")))
low = {s["name"]: s["thermo"]["data"][0] for s in mech["species"]}
for s in mech["species"]:
    assert s["thermo"]["model"] == "NASA7" and s["thermo"]["temperature-ranges"][0] == 200.0
c2 = 1.438776877  # second radiation constant hc/k, cm K
E = [c2 * (59.3392 * j * (j + 1) - 0.04599 * (j * (j + 1)) ** 2 + 5.2e-5 * (j * (j + 1)) ** 3) for j in range(20)]

def cv_rotation(T, js, weight):
    z = m1 = m2 = 0.0
    for j in js:
        g = weight * (2 * j + 1) * math.exp(-E[j] / T)
        z += g; m1 += g * E[j]; m2 += g * E[j] ** 2
    m1 /= z; m2 /= z
    return (m2 - m1 * m1) / T ** 2

def cp_h2_normal(T): return 2.5 + 0.25 * cv_rotation(T, range(0, 20, 2), 1) + 0.75 * cv_rotation(T, range(1, 20, 2), 3)
def einstein(T, nu):
    x = nu * c2 / T
    return x * x * math.exp(x) / math.expm1(x) ** 2
def cp_o2(T): return 3.5 + einstein(T, 1556.385)
def cp_h2o(T): return 4.0 + sum(einstein(T, nu) for nu in (3657.053, 1594.746, 3755.929))
def cp7(a, T): return a[0] + T * (a[1] + T * (a[2] + T * (a[3] + T * a[4])))
def h7(a, T): return T * (a[0] + T * (a[1] / 2 + T * (a[2] / 3 + T * (a[3] / 4 + T * a[4] / 5)))) + a[5]
def s7(a, T): return a[0] * math.log(T) + T * (a[1] + T * (a[2] / 2 + T * (a[3] / 3 + T * a[4] / 4))) + a[6]

T0, TMIN = 200.0, 20.0
rows = []
for name, source in (("H2", cp_h2_normal), ("O2", cp_o2), ("H2O", cp_h2o)):
    T = np.arange(TMIN, T0 + 0.5, 1.0)
    y = np.array([source(t) for t in T])
    scale = np.array([T0 ** k for k in range(5)])
    A = np.vander(T, 5, increasing=True) / scale
    # Least squares with the equality cp(200 K) = TM-4513's (Lagrange multiplier), columns scaled by 200^k.
    K = np.zeros((6, 6)); K[:5, :5] = 2 * A.T @ A; K[:5, 5] = 1; K[5, :5] = 1
    rhs = np.zeros(6); rhs[:5] = 2 * A.T @ y; rhs[5] = cp7(low[name], T0)
    a = [float(v) for v in np.linalg.solve(K, rhs)[:5] / scale]
    a.append(h7(low[name], T0) - h7(a + [0.0, 0.0], T0))
    a.append(s7(low[name], T0) - s7(a + [0.0], T0))
    worst = max(abs(cp7(a, t) / source(t) - 1) for t in T)
    floor = min(cp7(a, t) for t in np.linspace(1e-3, T0, 20001))
    print(f"// {name}: largest |cp/source - 1| from 20 to 200 K {worst:.2e}; at 200 K cp {cp7(a, T0) - cp7(low[name], T0):+.1e}, "
          f"h {h7(a, T0) - h7(low[name], T0):+.1e}, s {s7(a, T0) - s7(low[name], T0):+.1e} (over R); least cp/R below 200 K {floor:.4f}")
    rows.append((name, low[name], a))
for name, joins, a in rows:
    print(f'    {{"{name}", {{{", ".join(repr(float(v)) for v in joins)}}},\n'
          f'     {{{", ".join(repr(v) for v in a)}}}}},')
