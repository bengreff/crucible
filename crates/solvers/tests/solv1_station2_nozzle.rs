//! Goal-B **Station 2** pass criteria (SOLV-1 §6-4; SOLV-7 §6-1): cold-gas
//! choked flow through the axisymmetric De Laval nozzle vs quasi-1-D
//! isentropic theory — emergent p_c, emergent choking, nothing imposed.
//! Same shared-constants contract as the other certificates
//! (`bin/station2_nozzle_certificate` renders the committed artifact from
//! the same runners and constants; check.sh gate 4 diffs it).

use crucible_solvers::station2_nozzle::{
    CD_BAND, CD_FINEST_TOL, CF_DEV_MAX, EXIT_MACH_DEV_MAX, MACH_DEV_MAX, MDOT_SPREAD_MAX,
    NOZZLE_LEVELS, PC_BAND, STEADY_RESID_FINEST, STEADY_RESID_MAX, area_ratio,
    mach_from_area_ratio, masked_uniform_fixed_point, nozzle_study, run_nozzle,
};

#[test]
fn station2_choked_flow_vs_isentropic_theory() {
    let levels = nozzle_study().expect("study runs");
    assert_eq!(levels.len(), NOZZLE_LEVELS.len());
    let finest = &levels[levels.len() - 1];

    // The choked mass flow: Cd within the band everywhere, tight at the
    // finest, |Cd−1| never increasing under refinement.
    for l in &levels {
        assert!(
            (CD_BAND.0..=CD_BAND.1).contains(&l.cd_analytic),
            "Cd = {} at ({}, {}) outside {CD_BAND:?}",
            l.cd_analytic,
            l.n_r,
            l.n_z
        );
    }
    assert!(
        (finest.cd_analytic - 1.0).abs() < CD_FINEST_TOL,
        "finest Cd = {}",
        finest.cd_analytic
    );
    for w in levels.windows(2) {
        assert!(
            (w[1].cd_analytic - 1.0).abs() <= (w[0].cd_analytic - 1.0).abs() * 1.05,
            "|Cd−1| must not grow under refinement: {} → {}",
            w[0].cd_analytic,
            w[1].cd_analytic
        );
    }

    for l in &levels {
        // Emergent chamber pressure reports the reservoir back.
        assert!(
            (PC_BAND.0..=PC_BAND.1).contains(&l.p_c_over_p0),
            "p_c/p0 = {} outside {PC_BAND:?}",
            l.p_c_over_p0
        );
        // Steady, supersonic-exit, oracle-tracking flow.
        assert!(
            l.steady_resid < STEADY_RESID_MAX,
            "resid {}",
            l.steady_resid
        );
        assert!(l.mach_dev_max < MACH_DEV_MAX, "Mach dev {}", l.mach_dev_max);
        assert!(l.exit_mach_centerline > 1.0, "exit must be supersonic");
        let m_e_1d = mach_from_area_ratio(area_ratio(6.0 - 1e-9), true);
        assert!(
            ((l.exit_mach_centerline - m_e_1d) / m_e_1d).abs() < EXIT_MACH_DEV_MAX,
            "exit M {} vs 1-D {m_e_1d}",
            l.exit_mach_centerline
        );
        // Thrust coefficient vs ideal vacuum C_F.
        assert!(
            ((l.cf - l.cf_ideal) / l.cf_ideal).abs() < CF_DEV_MAX,
            "C_F {} vs ideal {}",
            l.cf,
            l.cf_ideal
        );
        // The declared slip-wall transpiration stays bounded…
        assert!(
            l.mdot_spread < MDOT_SPREAD_MAX,
            "ṁ spread {}",
            l.mdot_spread
        );
    }
    assert!(
        finest.steady_resid < STEADY_RESID_FINEST,
        "finest resid {}",
        finest.steady_resid
    );
    // …and shrinks under refinement (it is a resolution artifact).
    for w in levels.windows(2) {
        assert!(
            w[1].mdot_spread < w[0].mdot_spread,
            "ṁ spread must shrink: {} → {}",
            w[0].mdot_spread,
            w[1].mdot_spread
        );
    }
}

#[test]
fn station2_masked_uniform_state_is_a_bitwise_fixed_point() {
    // The run-decomposed masked sweeps are exact: a uniform gas at rest in
    // the stair-stepped cavity does not move, bit for bit — grid-aligned
    // mirror and slip-ghost wall alike.
    assert!(masked_uniform_fixed_point(false), "grid-aligned mirror");
    assert!(masked_uniform_fixed_point(true), "slip-ghost wall");
}

#[test]
fn station2_rerun_is_byte_identical() {
    let bits = || -> Vec<u64> {
        let (n_r, n_z) = NOZZLE_LEVELS[0];
        let (g, f, _, _) = run_nozzle(n_r, n_z).expect("run");
        g.bricks()
            .iter()
            .flat_map(|b| {
                f.ids()
                    .into_iter()
                    .flat_map(|id| b.field(id).iter().map(|v| v.to_bits()).collect::<Vec<_>>())
            })
            .collect()
    };
    assert_eq!(bits(), bits(), "rerun must be byte-identical (S6)");
}
