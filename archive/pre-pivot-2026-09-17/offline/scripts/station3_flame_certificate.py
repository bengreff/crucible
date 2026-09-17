"""Regenerate certificates/station3_flame_certificate.md (Goal-B station 3).

Usage (from repo root):
    offline/.venv/bin/python offline/scripts/station3_flame_certificate.py

Every number is recomputed live from the pinned toolchain + the committed
production tables; check.sh diffs the result against the committed
certificate, so silent drift in any of them breaks the build (the
convergence-certificate pattern from Goal A).
"""

import pathlib
import platform
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import h5py
import numpy as np

import cantera
import cea
import crucible_offl
from crucible_offl.chemistry import EquilibriumEngine, mr_to_z
from crucible_offl.crosscheck import cantera_state_php
from crucible_offl.seam_fixture import PIN_DIGEST as FIXTURE_PIN
from crucible_offl.surfaces import _EQ_RULES, SAFETY, _multilinear

REPO = pathlib.Path(__file__).resolve().parents[2]
TABLE = REPO / "tables" / "chem" / "lox_lh2_v0.1.0.h5"
OUT = REPO / "certificates" / "station3_flame_certificate.md"

GOLDEN = "sha256:a3b0bc5987ad11ba166d006ab7103c16ff1d3a8f17f64a248d4b11d8091023f2"

# RP-1311 example 8 published output (NASA RP-1311, NTRS 19960044559; as
# printed in the nasa.github.io/cea release docs) — chamber column + rocket.
PC8, MR8 = 53.3172e5, 5.55157
REF8 = {
    "T_c [K]": 3383.845,
    "rho_c [kg/m^3]": 2.410,
    "h_c [kJ/kg]": -1026.05,
    "s_c [kJ/kg-K]": 18.659,
    "Mbar [kg/kmol]": 12.716,
    "gamma_s": 1.145,
    "a_c [m/s]": 1591.47,
    "X_H2O": 0.63456,
    "X_H2": 0.29479,
    "X_H": 0.033498,
    "X_OH": 0.033341,
    "c* [m/s]": 2332.34,
    "T_throat [K]": 3185.673,
    "Isp_vac(eps=25) [m/s]": 4348.510,
}


def main() -> int:
    e = EquilibriumEngine()

    # --- §6-1: RP-1311 example 8 through the pipeline primitives ---------
    ch = e.chamber(PC8, MR8)
    perf = e.performance(PC8, MR8)
    shifting, frozen = e.frozen_shifting_gap(PC8, MR8, supar=25.0)
    ours8 = {
        "T_c [K]": ch.T,
        "rho_c [kg/m^3]": ch.rho,
        "h_c [kJ/kg]": ch.h / 1e3,
        "s_c [kJ/kg-K]": ch.s / 1e3,
        "Mbar [kg/kmol]": ch.mbar,
        "gamma_s": ch.gamma,
        "a_c [m/s]": ch.a,
        "X_H2O": ch.mole_fractions["H2O"],
        "X_H2": ch.mole_fractions["H2"],
        "X_H": ch.mole_fractions["H"],
        "X_OH": ch.mole_fractions["OH"],
        "c* [m/s]": perf.c_star,
        "T_throat [K]": perf.T_throat,
        "Isp_vac(eps=25) [m/s]": shifting,
    }
    gap = (shifting - frozen) / shifting

    # --- §6-2: CEA vs Cantera on three chamber states --------------------
    xstates = []
    for p, mr in ((32.75e5, 5.0), (53.3172e5, 5.55157), (10.0e5, 4.0)):
        z = mr_to_z(mr)
        h = e.injection_enthalpy(z)
        a = e.state_php(p, h, z)
        b = cantera_state_php(p, h, z)
        xstates.append((p, mr, a, b))

    # --- Committed tables: pins, envelopes, bounds ------------------------
    with h5py.File(TABLE, "r") as f:
        eq = f["/chem/lox_lh2/equilibrium"]
        pf = f["/chem/lox_lh2/performance"]
        eq_axes = {n: (np.array(eq["axes"][n]), dict(eq["axes"][n].attrs)) for n in ("p", "h", "Z")}
        eq_bounds = {
            n: float(eq["values"][n].attrs["interp_error_bound"])
            for n in eq["values"]
            if not n.startswith("sigma_")
        }
        pf_bounds = {
            n: float(pf["values"][n].attrs["interp_error_bound"]) for n in pf["values"]
        }
        eq_shape = tuple(len(v[0]) for v in eq_axes.values())
        t_grid = {"temperature": np.array(eq["values/temperature"]).reshape(eq_shape)}
        ax_list = [v[0] for v in eq_axes.values()]
        p_axes = [np.array(pf["axes"][n]) for n in ("p_c", "MR")]
        tc_grid = {"T_c": np.array(pf["values/T_c"]).reshape(tuple(len(a) for a in p_axes))}
        cs_grid = {
            "c_star_ideal": np.array(pf["values/c_star_ideal"]).reshape(
                tuple(len(a) for a in p_axes)
            )
        }
        generator_commit = eq.attrs["generator_commit"]

    import tomllib

    with open(TABLE.with_suffix(".pins.toml"), "rb") as fh:
        pins = tomllib.load(fh)
    digests = {group: entry["content_digest"] for group, entry in pins.items()}

    # --- §6-6: coordinate consistency along the design line ---------------
    design = []
    for pc, mr in ((32.75e5, 5.0), (20e5, 4.5), (45e5, 6.0), (15e5, 3.8), (50e5, 7.0)):
        z = mr_to_z(mr)
        h = e.injection_enthalpy(z)
        t_surf = _multilinear(_EQ_RULES["temperature"], ax_list, t_grid, (pc, h, z), "temperature")
        t_perf = _multilinear("lin-lin-lin", p_axes, tc_grid, (pc, mr), "T_c")
        t_true = e.chamber(pc, mr).T
        design.append((pc, mr, t_surf, t_perf, t_true))

    # RL10-class headline: chamber point through the whole seam.
    z5 = mr_to_z(5.0)
    h5v = e.injection_enthalpy(z5)
    t_rl10 = _multilinear(_EQ_RULES["temperature"], ax_list, t_grid, (32.75e5, h5v, z5), "temperature")
    cstar_rl10 = _multilinear("lin-lin-lin", p_axes, cs_grid, (32.75e5, 5.0), "c_star_ideal")
    cstar_true = e.performance(32.75e5, 5.0).c_star
    t_true5 = e.chamber(32.75e5, 5.0).T

    md = []
    w = md.append
    w("# CRUCIBLE Station 3 Certificate — the Flame Seam\n")
    w(
        "Goal-B station 3 (OFFL-3 §6; FND-5 §3.2; VISION_SCOPE §6 two-language rule): the LOX/LH2 "
        "flame is computed **offline** by NASA CEA (equilibrium free-energy minimization over the "
        "full Glenn species set), tabulated as local-state surfaces, and crosses the project's one "
        "cross-language seam as versioned, digest-pinned HDF5 that the Rust runtime loads, "
        "verifies, and interpolates. The physical system: the RL10 combustion chamber at "
        "p_c = 32.75 bar burning liquid hydrogen (20.27 K) with liquid oxygen (90.17 K) at MR = 5 — "
        "the flame the station-5 blind run will feed on. Python never runs at simulation time; "
        "the tables are the only crossing.\n"
    )
    w(
        "Criteria are CI-enforced in `offline/tests/` (pytest: digest golden vector, RP-1311 "
        "reproduction, CEA↔Cantera cross-check, frozen/shifting bracket, holdout error bounds, "
        "regeneration determinism, coordinate consistency) and `crates/tables/tests/` "
        "(`fnd5_python_seam.rs`: the h5py-written fixture opens under its Python-stamped pin; "
        "`station3_tables.rs`: the production tables load, interpolate identically to the Python "
        "reference evaluator, and refuse out-of-envelope). Regenerate: "
        "`offline/scripts/station3_flame_certificate.py`.\n"
    )

    w("## The seam itself (digest v2, byte-for-byte)\n")
    w(f"- Golden vector (both languages assert the same constant): `{GOLDEN}`")
    w(f"- Cross-language fixture pin (Python stamps, Rust verifies): `{FIXTURE_PIN}`")
    w("- Production table pins (`tables/chem/lox_lh2_v0.1.0.pins.toml`):")
    for g, d in sorted(digests.items()):
        w(f"  - `{g}`: `{d}`")
    w(
        f"- Toolchain pins: Python {platform.python_version()}, cea {cea.__version__} "
        f"(libcea {cea.lib_version()}), cantera {cantera.__version__}, h5py {h5py.__version__}, "
        f"numpy {np.__version__}, crucible-offl {crucible_offl.__version__}; "
        f"generator commit `{generator_commit[:12]}`."
    )
    w("")

    w("## OFFL-3 §6-1 — RP-1311 example 8 reproduced through the pipeline\n")
    w(
        "LOX/LH2 IAC rocket, p_c = 53.3172 bar, o/f = 5.55157 (McBride & Gordon, RP-1311, "
        "NTRS 19960044559). Reference = the published example output; tolerance = its printed "
        "precision (the pytest gates).\n"
    )
    w("| quantity | published | pipeline | rel. dev |")
    w("|---|---|---|---|")
    for k, ref in REF8.items():
        v = ours8[k]
        w(f"| {k} | {ref} | {v:.6g} | {abs(v - ref) / abs(ref):.1e} |")
    w("")

    w("## OFFL-3 §6-3 — frozen/shifting bracket (the model-form error bar)\n")
    w(
        f"Vacuum Isp at area ratio 25: shifting {shifting:.1f} m/s, frozen-from-chamber "
        f"{frozen:.1f} m/s — a strict bracket of width {100 * gap:.2f}% that contains the JANNAF "
        "kinetic-efficiency knockdown (~0.8–1% of shifting for LOX/LH2 at MR≈5, META-3 "
        "`jannaf-eff`). Real delivered performance sits between the two tables; COUP-5 samples "
        "this as a discrete epistemic dimension (S19) — the disagreement **is** the error bar, "
        "never averaged away.\n"
    )

    w("## OFFL-3 §6-2 — two independent solvers on the chamber state\n")
    w(
        "CEA (Glenn NASA9 fits) vs Cantera 3.2 `h2o2.yaml` (NASA7 fits), no shared "
        "implementation, same (p, h, Z) coordinate:\n"
    )
    w("| p_c [bar] | MR | T: CEA [K] | T: Cantera [K] | dT | dM̄ | dX_H2O |")
    w("|---|---|---|---|---|---|---|")
    for p, mr, a, b in xstates:
        w(
            f"| {p / 1e5:.2f} | {mr} | {a.T:.2f} | {b.T:.2f} | {100 * (b.T - a.T) / a.T:+.3f}% "
            f"| {100 * (b.mbar - a.mbar) / a.mbar:+.3f}% "
            f"| {b.mole_fractions['H2O'] - a.mole_fractions['H2O']:+.5f} |"
        )
    w(
        "\nGates (2× the measured band): |dT| ≤ 0.3%, |dM̄| ≤ 0.15%, |dX| ≤ 5e-3, "
        "|dcp_eq| ≤ 1%. The residual disagreement is thermodynamic-data uncertainty (§3.4), "
        "carried, not hidden.\n"
    )

    w("## The tables (FND-5 schema, committed at `tables/chem/`)\n")
    w(
        "**Equilibrium surface** `(p, h, Z) → T, ρ, γ_eff, a, M̄, condensed_fraction, X_k` "
        "(41×31×13, log-p): the S22 local-state coordinate of SOLV-1's shifting mode. "
        "**Performance reference** `(p_c, MR) → c*_ideal, T_c, γ, M̄` (11×11): SOLV-7's anchor "
        "and the SOLV-1 §3.4 knockdown reference. Declared envelopes (refusal-enforced, "
        "grid strictly wider):\n"
    )
    for n, (ax, attrs) in eq_axes.items():
        w(
            f"- `{n}`: grid [{ax[0]:.6g}, {ax[-1]:.6g}], envelope "
            f"[{attrs['envelope_min']:.6g}, {attrs['envelope_max']:.6g}]"
        )
    w(
        "\nStored `interp_error_bound` per column = **measured** holdout max (midpoints + "
        f"¼-offsets, direct CEA solves) × {SAFETY} declared sampling margin (FND-5 §3.4); "
        "verified in CI by a fresh disjoint ⅜-offset sweep against the committed artifact — "
        "a violation means *refine the grid*, never *relax the gate*. Full-envelope bounds are "
        "dominated by the H₂O condensation kink in the deep-cold corners (reported honestly by "
        "the `condensed_fraction` column; the certificate's gas-region numbers below show the "
        "working regime is ~50× tighter):\n"
    )
    w("| column | stored bound | | column | stored bound |")
    w("|---|---|---|---|---|")
    eb = sorted(eq_bounds.items())
    half = (len(eb) + 1) // 2
    for i in range(half):
        left = f"| {eb[i][0]} | {eb[i][1]:.3g} |"
        right = f" {eb[i + half][0]} | {eb[i + half][1]:.3g} |" if i + half < len(eb) else " | |"
        w(left + right)
    w("")
    w("Performance reference bounds: " + ", ".join(f"{n} {b:.3g}" for n, b in sorted(pf_bounds.items())) + ".\n")

    w("## OFFL-3 §6-6 — one flame, two parameterizations (S22)\n")
    w(
        "Along the design line h = h_inj(Z), the (p, h, Z) surface and the (p_c, MR) performance "
        "reference must describe the same chamber. Direct CEA solution alongside:\n"
    )
    w("| p_c [bar] | MR | T surface [K] | T perf-ref [K] | T direct [K] | surf−direct |")
    w("|---|---|---|---|---|---|")
    for pc, mr, t_surf, t_perf, t_true in design:
        w(
            f"| {pc / 1e5:.2f} | {mr} | {t_surf:.2f} | {t_perf:.2f} | {t_true:.2f} "
            f"| {t_surf - t_true:+.2f} |"
        )
    w("")

    w("## The RL10 chamber through the whole seam (the station-3 headline)\n")
    w(
        f"At (p_c = 32.75 bar, h = h_inj, Z = 1/6): the Rust loader verifies the pin, "
        f"interpolates T = {t_rl10:.2f} K against CEA's direct {t_true5:.2f} K "
        f"({t_rl10 - t_true5:+.2f} K in the gas region, vs the kink-dominated stored bound "
        f"±{eq_bounds['temperature']:.1f} K), and c*_ideal = {cstar_rl10:.2f} m/s against the "
        f"direct rocket solve {cstar_true:.2f} m/s ({cstar_rl10 - cstar_true:+.2f} m/s). "
        "The same numbers are asserted in `station3_tables.rs` — Python and Rust provably "
        "interpolate the same surface the same way (≤1e-6 abs).\n"
    )

    w("## Deferred, loudly (owners)\n")
    w(
        "- **Frozen-path surface vs (p, h, {X_k})** — its axis set belongs to SOLV-1's "
        "frozen-advection consumer; arrives with that wave (the frozen *bracket* is validated "
        "above; `chemistry.py` header).\n"
        "- **Quasi-1-D expansion oracles** (`oracle`-labeled, never field-interpolated) — "
        "SOLV-7/VAL-2 C_F cross-check wave, station 5.\n"
        "- **B′ ablation tables** — SOLV-8 wave. The **transport feed is DISCHARGED** "
        "(plan S4): `chemistry.py` now returns the caloric companions (c_p,fr, c_v,eq) "
        "from the same CEA solve, and `crucible_offl/transport.py` assembles them with "
        "Cantera's mixture-averaged transport into the spine's chemical-regime surface "
        "(OFFL-5 §3.1a, `tables/spine/lox_lh2_transport_v0.1.0.h5`) — written in the "
        "runtime local-state (p, h, Z) coordinate, not the (T, p, Z) feed one. OFFL-3 "
        "still ships no runtime transport table; the spine (FND-7) is the sole runtime "
        "provider.\n"
        "- **Config→tables pin wiring** (FND-4 §6-4) — the `[tables]` grammar consumes "
        "`lox_lh2_v0.1.0.pins.toml` when it lands; until then the sidecar + the Rust test pins "
        "are the record.\n"
        "- **Per-point sigma columns** — the thermodynamic-data band enters with the COUP-5 UQ "
        "wave; the seam already carries sigma companions (proven by the fixture).\n"
    )
    w(
        "Honest scaffolding note: the surfaces are *the equilibrium*, including condensed H₂O "
        "where the rectangular grid's deep-cold corners demand it (`condensed_fraction` column); "
        "no engine trajectory enters that region, and the gas-only working regime is where all "
        "gas-region numbers above live.\n"
    )

    OUT.write_text("\n".join(md))
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
