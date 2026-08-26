//! SOLV-1 §3.3 / §3.6 + FND-2 §3.4 + FND-3 §3.3 (plan S9) — **cut geometry
//! at uniform N_θ > 1**: per-θ-sector κ and six face apertures consumed by
//! the Euler sweeps + State Redistribution. The gates:
//!
//! * (a) a revolved cut world (identical per-(r,z) κ/apertures through
//!   BOTH builders) marches bitwise per θ-plane against its N_θ = 1 twin —
//!   THE gate-5-in-miniature proof (the station-2/4/5 certificates march
//!   contour cut worlds at N_θ = 1; the powers-of-two metric cancellation
//!   extends to the κ/aperture products, S8's argument);
//! * (b) well-balance on a genuinely θ-varying world (κ and apertures vary
//!   by sector, all κ > 0): a uniform state at rest stays a fixed point —
//!   ρ/ρE/ρC/ρb **bitwise** (their fluxes are exactly zero), the momenta
//!   at round-off (see the in-test note for why not bitwise);
//! * (c) audited conservation through a shock transient on a θ-varying cut
//!   world with per-sector slivers (SRD active per sector);
//! * (d) SRD neighborhoods are per-sector: a sliver in ONE sector merges
//!   (through its θ-faces when r/z are sealed — the S9 θ-neighbor
//!   candidates) while the other sectors of the same (r,z) cell are
//!   untouched bitwise;
//! * (e) the six-aperture storage plane contract: on a θ-uniform cut world
//!   `aperture_cell` ≡ `aperture_rz` bitwise for every sector (and full-box
//!   worlds never allocate geometry at all);
//! * (f) a toy chamber (converging-diverging revolved profile with a
//!   θ-varying wall bump, κ > 0 per sector) marches 50+ audited cold-flow
//!   steps and conserves. NOTE: `crates/grid/src/geom3d` is absent in this
//!   worktree, so the world is built by a hand-written analytic supplier
//!   closure (the S12 `linear_wall_geom` fractions, θ-modulated); the
//!   voxelizer-output integration happens with the main S9 session.

use crucible_grid::{CellGeom, CellGeomTheta, FaceDir, Grid, GridSpec, Region};
use crucible_solvers::euler::{
    Cons, EULER_FIELDS, Euler, EulerFields, FlowBc, FlowBcs, GammaLaw, I_RHO, NCOMP, Prim,
    fill_from_prim, prim6,
};
use crucible_solvers::sdc::{FlowClass, Sdc};

const GAMMA: f64 = 1.4;
const ZERO_SRC: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];

fn closed_box_op(eos: GammaLaw) -> Euler<'static, GammaLaw> {
    Euler {
        eos,
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Reflecting,
            z_hi: FlowBc::Reflecting,
        },
        wall_normal: None,
        slip_wall_z_faces: false,
        combustion: None,
    }
}

/// Per-θ-plane bitwise snapshot: `out[bi][j]` = the plane's field bits.
fn plane_bits(g: &Grid, f: &EulerFields) -> Vec<Vec<Vec<u64>>> {
    let ids = f.ids();
    g.bricks()
        .iter()
        .map(|b| {
            (0..b.n_theta())
                .map(|j| {
                    ids.iter()
                        .flat_map(|&id| {
                            (0..64).map(move |local| b.field(id)[b.cell_index(j, local)].to_bits())
                        })
                        .collect()
                })
                .collect()
        })
        .collect()
}

// --- The revolved fixture (gates a, e): a cylinder wall at r = 0.55 --------

const H: f64 = 0.1;
const WALL_R: f64 = 0.55;

/// The ONE per-(r,z) fraction rule both builders consume: κ from the
/// cylindrical measure, r-apertures binary by face radius, z-apertures = κ.
fn cylinder_rz(i_r: usize) -> (f64, [f64; 4]) {
    let (r0, r1) = (i_r as f64 * H, (i_r + 1) as f64 * H);
    let x = WALL_R.clamp(r0, r1);
    let kappa = ((x * x - r0 * r0) / (r1 * r1 - r0 * r0)).clamp(0.0, 1.0);
    let ap_r = |rf: f64| if WALL_R > rf { 1.0 } else { 0.0 };
    (kappa, [ap_r(r0), ap_r(r1), kappa, kappa])
}

fn cylinder_spec(n_theta: u32) -> GridSpec {
    GridSpec {
        r_min: 0.0,
        dr: H,
        n_r: 8,
        z_min: 0.0,
        dz: H,
        n_z: 8,
        n_theta_max: n_theta,
        axisymmetry_assertion: n_theta == 1,
    }
}

/// The N_θ = 1 build (the certified revolved path).
fn cylinder_grid_1() -> Grid {
    Grid::build_with_geometry(cylinder_spec(1), EULER_FIELDS, |i_r, _| {
        let (kappa, aperture) = cylinder_rz(i_r);
        if kappa > 0.0 {
            CellGeom {
                region: Region::Gas,
                kappa,
                aperture,
            }
        } else {
            CellGeom {
                region: Region::Exterior,
                kappa: 0.0,
                aperture: [0.0; 4],
            }
        }
    })
    .expect("revolved cylinder builds at N_θ = 1")
}

/// The SAME per-(r,z) values through the S9 θ-builder (θ-uniform; θ-face
/// apertures = κ, matching the revolved path's declared convention).
fn cylinder_grid_theta(n_theta: u32) -> Grid {
    Grid::build_with_geometry_theta(
        cylinder_spec(n_theta),
        EULER_FIELDS,
        |i_r, _j, _i_z| {
            let (kappa, ap) = cylinder_rz(i_r);
            CellGeomTheta {
                kappa,
                aperture: [ap[0], ap[1], ap[2], ap[3], kappa, kappa],
            }
        },
        |_, _| Region::Exterior,
    )
    .expect("revolved cylinder builds through the θ-builder")
}

/// Axisymmetric Gaussian pressure pulse inside the gas region.
fn axisym_pulse(r: f64, _th: f64, z: f64) -> Prim {
    let dp = 0.4 * (-((r / 0.15).powi(2) + ((z - 0.3) / 0.15).powi(2))).exp();
    prim6(1.3, 0.0, 0.0, 0.0, 1.0e5 * (1.0 + dp), 0.4)
}

// Covers FND-3 §3.3 / SOLV-1 §3.3 (S9): THE gate-5-in-miniature proof — the
// certified N_θ = 1 revolved-cut march reproduced bitwise per θ-plane at
// N_θ = 8 through the axis machinery, the aperture-weighted sweeps, and the
// per-sector SRD (the wall ring is a κ ≈ 0.48 sliver, so SRD is active).
#[test]
fn fnd3_s33_revolved_cut_world_at_n_theta_8_matches_n_theta_1_bitwise_per_plane() {
    let eos = GammaLaw { gamma: GAMMA };
    let mut g1 = cylinder_grid_1();
    let mut g8 = cylinder_grid_theta(8);
    // The fixture must actually contain an SRD sliver (κ < 0.5).
    assert!(
        g8.kappa_at(5, 3, 0) < crucible_solvers::euler::KAPPA_SRD,
        "fixture lost its sliver ring"
    );
    let f1 = EulerFields::resolve(&g1).expect("fields");
    let f8 = EulerFields::resolve(&g8).expect("fields");
    fill_from_prim(&mut g1, &f1, &eos, axisym_pulse);
    fill_from_prim(&mut g8, &f8, &eos, axisym_pulse);
    let op = closed_box_op(eos);
    let (mut sdc1, mut sdc8) = (Sdc::new(), Sdc::new());
    let flow1 = FlowClass {
        op: &op,
        fields: &f1,
    };
    let flow8 = FlowClass {
        op: &op,
        fields: &f8,
    };
    let mut t = 0.0;
    for step in 0..30 {
        let dt = sdc8.stable_dt(&g8, &flow8, 0.4).expect("dt");
        sdc8.step_flow(&mut g8, &flow8, t, dt)
            .expect("audited 3-D step");
        sdc1.step_flow(&mut g1, &flow1, t, dt)
            .expect("audited 2-D step");
        t += dt;
        let b1 = plane_bits(&g1, &f1);
        let b8 = plane_bits(&g8, &f8);
        for (bi, planes) in b8.iter().enumerate() {
            for (j, plane) in planes.iter().enumerate() {
                assert_eq!(
                    plane, &b1[bi][0],
                    "θ-plane {j} of brick {bi} diverged from the N_θ = 1 cut march \
                     at step {step}"
                );
            }
        }
    }
}

// --- The θ-varying fixture (gate b): well-balance ---------------------------

/// A porous bumpy world: every cell gas, per-sector κ ∈ [0.55, 0.70] (no
/// SRD — this gate isolates well-balance), all apertures in (0, 1) varying
/// by sector, every shared face computed from its own indices (coherent by
/// construction).
fn bumpy_theta_grid(nt: u32) -> Grid {
    let spec = GridSpec {
        r_min: 0.1,
        dr: 0.1,
        n_r: 4,
        z_min: 0.0,
        dz: 0.1,
        n_z: 4,
        n_theta_max: nt,
        axisymmetry_assertion: false,
    };
    Grid::build_with_geometry_theta(
        spec,
        EULER_FIELDS,
        move |i_r, j, i_z| {
            let ap_r = |fi: usize| 0.8 + 0.04 * (((fi as u32 + j) % 3) as f64);
            let ap_z = |fz: usize| 0.75 + 0.05 * (((fz as u32 + j) % 3) as f64);
            let ap_th = |fi: u32| 0.7 + 0.05 * (((fi + i_r as u32 + i_z as u32) % nt) as f64 / 4.0);
            CellGeomTheta {
                kappa: 0.55 + 0.05 * (((j + i_r as u32 + i_z as u32) % 4) as f64),
                aperture: [
                    ap_r(i_r),
                    ap_r(i_r + 1),
                    ap_z(i_z),
                    ap_z(i_z + 1),
                    ap_th(j),
                    ap_th((j + 1) % nt),
                ],
            }
        },
        |_, _| Region::Exterior,
    )
    .expect("bumpy θ-varying world builds")
}

// Covers SOLV-1 §3.3 well-balance extended to the θ-limb (S9, FND-2 §3.4).
//
// WHY the momenta are round-off rather than bitwise: for a uniform state
// at rest the sweep's face contribution is `(A·ap)·p` products differenced
// then divided by κV, the geometric source is `(A_out·p − A_in·p)/V`, and
// the closure limb is `p·((a₊A₊ − a₋A₋ − κ(A_out−A_in)))/(κV)` — three
// separately-rounded product/association orders whose EXACT sum is zero
// but whose floating-point sum is ~ulp(p·A/κV). This is the identical
// round-off class the S12 cut-cell well-balance gate accepts (1e-11 over
// 25 steps there); one step here bounds it at 1e-14 of the momentum scale
// ρ·a. The zero-flux components (ρ, ρE, ρC, ρb) cancel EXACTLY (their
// face fluxes are ±0.0 products), so they are asserted bitwise.
#[test]
fn solv1_s33_uniform_state_is_a_fixed_point_on_a_theta_varying_cut_world() {
    let eos = GammaLaw { gamma: GAMMA };
    let mut g = bumpy_theta_grid(4);
    let f = EulerFields::resolve(&g).expect("fields");
    let (rho0, p0) = (1.3, 2.7e5);
    fill_from_prim(&mut g, &f, &eos, |_, _, _| {
        prim6(rho0, 0.0, 0.0, 0.0, p0, 0.4)
    });
    let before = plane_bits(&g, &f);
    let op = closed_box_op(eos);
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
    sdc.step_flow(&mut g, &flow, 0.0, dt).expect("audited step");
    let after = plane_bits(&g, &f);
    let ids = f.ids();
    // ρ, ρE, ρC, ρb: bitwise fixed point (exactly-zero fluxes).
    let n_fields = ids.len();
    for (bi, planes) in after.iter().enumerate() {
        for (j, plane) in planes.iter().enumerate() {
            let stride = plane.len() / n_fields;
            for (k, _) in ids.iter().enumerate() {
                if k == 1 || k == 2 || k == 3 {
                    continue; // the momenta are tolerance-gated below
                }
                assert_eq!(
                    &plane[k * stride..(k + 1) * stride],
                    &before[bi][j][k * stride..(k + 1) * stride],
                    "component {k} of θ-plane {j}, brick {bi} moved on a uniform state"
                );
            }
        }
    }
    // Momenta: ≤ 1e-14 of the momentum scale ρ·a (see the head note).
    let a = eos.sound_speed(rho0, p0);
    let scale = rho0 * a;
    let mut worst = 0.0f64;
    g.for_each_active_cell(|cell| {
        let b = g.brick(cell.bi);
        for k in [1usize, 2, 3] {
            worst = worst.max(b.field(ids[k])[cell.idx].abs());
        }
    });
    assert!(
        worst <= 1e-14 * scale,
        "momentum residual {worst:.3e} exceeds 1e-14·ρa = {:.3e}",
        1e-14 * scale
    );
}

// --- The per-sector-sliver fixture (gates c, d) -----------------------------

/// The cylinder world with ring 5's κ varying BY SECTOR through the SRD
/// threshold: sectors with κ ∈ {0.35, 0.45} are slivers, {0.55} are not —
/// per-sector SRD in anger. Rings 0–4 are full; the sliver ring's r− face
/// is fully open onto the κ = 1 ring (every sector's neighborhood
/// resolves in-sector through it).
fn sector_sliver_grid(nt: u32) -> Grid {
    Grid::build_with_geometry_theta(
        cylinder_spec(nt),
        EULER_FIELDS,
        move |i_r, j, _i_z| {
            if i_r < 5 {
                CellGeomTheta {
                    kappa: 1.0,
                    aperture: [1.0; 6],
                }
            } else if i_r == 5 {
                let kappa = 0.35 + 0.1 * ((j % 3) as f64);
                let ap_th = |fi: u32| 0.6 + 0.05 * ((fi % 2) as f64);
                CellGeomTheta {
                    kappa,
                    aperture: [1.0, 0.0, kappa, kappa, ap_th(j), ap_th((j + 1) % nt)],
                }
            } else {
                CellGeomTheta {
                    kappa: 0.0,
                    aperture: [0.0; 6],
                }
            }
        },
        |_, _| Region::Exterior,
    )
    .expect("sector-sliver world builds")
}

// Covers COUP-2 §3.1 (the armed audit every step) + FND-3 §3.3 (S9):
// conservation through a shock transient on a θ-varying cut world with
// per-sector SRD active.
#[test]
fn coup2_s31_transient_on_a_theta_varying_cut_world_conserves() {
    let eos = GammaLaw { gamma: GAMMA };
    let mut g = sector_sliver_grid(4);
    let f = EulerFields::resolve(&g).expect("fields");
    // A strong off-axis 3-D pressure bump drives shocks into the θ-varying
    // wall ring.
    fill_from_prim(&mut g, &f, &eos, |r, th, z| {
        let (x, y) = (r * th.cos(), r * th.sin());
        let d2 = (x - 0.2) * (x - 0.2) + y * y + (z - 0.4) * (z - 0.4);
        prim6(
            1.0,
            0.0,
            0.0,
            0.0,
            1.0e5 * (1.0 + 9.0 * (-d2 / 0.01).exp()),
            0.3,
        )
    });
    let ids = f.ids();
    let mass0 = g.reduce_kappa_volume_weighted(ids[I_RHO]);
    let en0 = g.reduce_kappa_volume_weighted(ids[4]);
    let c0 = g.reduce_kappa_volume_weighted(ids[5]);
    let op = closed_box_op(eos);
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    for _ in 0..60 {
        let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
        // The armed COUP-2 audit IS the per-step conservation gate.
        sdc.step_flow(&mut g, &flow, 0.0, dt).expect("audited step");
    }
    let mass1 = g.reduce_kappa_volume_weighted(ids[I_RHO]);
    let en1 = g.reduce_kappa_volume_weighted(ids[4]);
    let c1 = g.reduce_kappa_volume_weighted(ids[5]);
    assert!(
        ((mass1 - mass0) / mass0).abs() < 1e-12,
        "mass drift {:.3e}",
        (mass1 - mass0) / mass0
    );
    assert!(
        ((en1 - en0) / en0).abs() < 1e-12,
        "energy drift {:.3e}",
        (en1 - en0) / en0
    );
    assert!(
        ((c1 - c0) / c0).abs() < 1e-12,
        "composition drift {:.3e}",
        (c1 - c0) / c0
    );
}

// Covers SOLV-1 §3.6 (S9): SRD neighborhoods are per-(θ-sector, r, z), and
// θ-neighbors are legal merge candidates. Cell (2, 2)'s r/z faces are
// sealed (thin walls, aperture 0 — legal); its ring is open in θ. Sector 0
// is a κ = 0.3 sliver whose ONLY flow-connected neighbors are its two ring
// partners — the neighborhood must recruit through a θ-face (the fixed
// candidate order r−, r+, θ−, θ+, z−, z+ reaches θ− first at equal κ), so
// SRD moves sectors 0 and 3 and leaves sectors 1 and 2 (non-members)
// bitwise untouched, conserving Σ κV·U.
#[test]
fn solv1_s36_srd_redistributes_slivers_per_sector() {
    let nt = 4u32;
    let spec = GridSpec {
        r_min: 0.1,
        dr: 0.1,
        n_r: 4,
        z_min: 0.0,
        dz: 0.1,
        n_z: 4,
        n_theta_max: nt,
        axisymmetry_assertion: false,
    };
    let sealed = (2usize, 2usize);
    let g0 = Grid::build_with_geometry_theta(
        spec,
        EULER_FIELDS,
        move |i_r, j, i_z| {
            let touches = |a: (usize, usize), b: (usize, usize)| a == sealed || b == sealed;
            // Any r/z face adjoining the sealed cell carries aperture 0 —
            // computed from the FACE (both sides agree bitwise).
            let ap_r = |fi: usize| {
                if fi > 0 && touches((fi - 1, i_z), (fi, i_z)) {
                    0.0
                } else {
                    1.0
                }
            };
            let ap_z = |fz: usize| {
                if fz > 0 && touches((i_r, fz - 1), (i_r, fz)) {
                    0.0
                } else {
                    1.0
                }
            };
            let kappa = if (i_r, i_z) == sealed {
                if j == 0 { 0.3 } else { 0.9 }
            } else {
                1.0
            };
            CellGeomTheta {
                kappa,
                aperture: [ap_r(i_r), ap_r(i_r + 1), ap_z(i_z), ap_z(i_z + 1), 1.0, 1.0],
            }
        },
        |_, _| Region::Exterior,
    )
    .expect("sealed-ring world builds");
    let mut g = g0;
    let eos = GammaLaw { gamma: GAMMA };
    let f = EulerFields::resolve(&g).expect("fields");
    // Distinct state in the sliver sector (θ_center(0, 4) = τ/8 < τ/4).
    fill_from_prim(&mut g, &f, &eos, |r, th, z| {
        let in_sealed = (0.3..0.4).contains(&r) && (0.2..0.3).contains(&z);
        let rho = if in_sealed && th < std::f64::consts::FRAC_PI_2 {
            2.0
        } else {
            1.0
        };
        prim6(rho, 0.0, 0.0, 0.0, 1.0e5, 0.3)
    });
    let ids = f.ids();
    let op = closed_box_op(eos);
    let ws = op
        .workspace(&g)
        .expect("workspace (validate + neighborhoods)");
    let stored0 = g.reduce_kappa_volume_weighted(ids[I_RHO]);
    let before = plane_bits(&g, &f);
    op.apply_srd(&mut g, &f, &ws);
    let after = plane_bits(&g, &f);
    let bi = g.brick_index(2, 2).expect("brick");
    let b = g.brick(bi);
    let local = 2 * 8 + 2; // local_rz(2, 2) inside brick (0, 0)
    let rho_at = |j: u32| b.field(ids[I_RHO])[b.cell_index(j, local)];
    // The sliver sector merged toward its θ-neighborhood (ρ dropped from 2).
    assert!(
        rho_at(0) < 2.0 - 1e-6,
        "sliver sector did not redistribute (ρ = {})",
        rho_at(0)
    );
    // Its recruited ring partner (θ−, sector 3) moved too.
    assert_ne!(
        before[bi][3], after[bi][3],
        "the θ-recruited member (sector 3) should have moved"
    );
    // Non-member sectors of the SAME (r,z) cell: bitwise untouched.
    for j in [1usize, 2] {
        assert_eq!(
            before[bi][j], after[bi][j],
            "sector {j} is not in any neighborhood and must stay bitwise"
        );
    }
    // SRD conserves Σ κV·U exactly (round-off).
    let stored1 = g.reduce_kappa_volume_weighted(ids[I_RHO]);
    assert!(
        ((stored1 - stored0) / stored0).abs() < 1e-14,
        "SRD moved the κV-weighted total: {:.3e}",
        (stored1 - stored0) / stored0
    );
    // And the world marches (audited, finite) with the per-sector sliver.
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    for _ in 0..15 {
        let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
        sdc.step_flow(&mut g, &flow, 0.0, dt).expect("audited step");
    }
    g.for_each_active_cell(|cell| {
        let rho = g.brick(cell.bi).field(ids[I_RHO])[cell.idx];
        assert!(rho.is_finite() && rho > 0.0, "density ran away");
    });
}

// Covers FND-3 §3.3 (S9 storage contract): on a θ-uniform cut world every
// sector's `aperture_cell`/`kappa_cell` is bitwise the plane-0
// `aperture_rz`/`kappa_rz`; and no-geometry worlds never allocate
// BrickGeom (the arithmetic-identity default path).
#[test]
fn fnd3_s33_theta_uniform_cut_world_aperture_cell_matches_aperture_rz_bitwise() {
    let g = cylinder_grid_theta(8);
    let dirs = [
        FaceDir::RMinus,
        FaceDir::RPlus,
        FaceDir::ZMinus,
        FaceDir::ZPlus,
        FaceDir::ThetaMinus,
        FaceDir::ThetaPlus,
    ];
    for b in g.bricks() {
        assert!(b.has_geom());
        for local in 0..64 {
            for j in 0..b.n_theta() {
                assert_eq!(
                    b.kappa_cell(j, local).to_bits(),
                    b.kappa_rz(local).to_bits(),
                    "κ sector {j} ≠ plane 0"
                );
                for dir in dirs {
                    assert_eq!(
                        b.aperture_cell(dir, j, local).to_bits(),
                        b.aperture_rz(dir, local).to_bits(),
                        "aperture {dir:?} sector {j} ≠ plane 0"
                    );
                }
            }
        }
    }
    // Full-box worlds: no BrickGeom, ever (the S12 idiom's precondition).
    let full = Grid::build(cylinder_spec(8), EULER_FIELDS).expect("full box");
    assert!(!full.has_cut_geometry());
    for b in full.bricks() {
        assert!(!b.has_geom(), "a no-geometry world allocated BrickGeom");
    }
}

// --- Gate (f): the toy chamber (analytic supplier — geom3d absent) ----------

/// Exact fractions for a linear wall w(z) = c0 + c1·z over one cell — the
/// S12 `linear_wall_geom` oracle (solv1_cut_cells), reused verbatim as the
/// hand-written analytic supplier standing in for the absent geom3d
/// voxelizer output.
fn linear_wall_geom(c0: f64, c1: f64, r0: f64, r1: f64, z0: f64, z1: f64) -> (f64, [f64; 4]) {
    let w = |z: f64| c0 + c1 * z;
    let mut pts = vec![z0];
    for rc in [r0, r1] {
        if (w(z0) - rc) * (w(z1) - rc) < 0.0 {
            pts.push(z0 + (z1 - z0) * (rc - w(z0)) / (w(z1) - w(z0)));
        }
    }
    pts.push(z1);
    pts.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
    let mut acc = 0.0f64;
    for p in pts.windows(2) {
        let (pa, pb) = (p[0], p[1]);
        let len = pb - pa;
        if len <= 0.0 {
            continue;
        }
        let wm = w(0.5 * (pa + pb));
        if wm <= r0 {
        } else if wm >= r1 {
            acc += 0.5 * (r1 * r1 - r0 * r0) * len;
        } else {
            let (va, vb) = (w(pa), w(pb));
            acc += 0.5 * (len * (va * va + va * vb + vb * vb) / 3.0 - r0 * r0 * len);
        }
    }
    let kappa = (acc / (0.5 * (r1 * r1 - r0 * r0) * (z1 - z0))).clamp(0.0, 1.0);
    let ap_r = |rf: f64| -> f64 {
        let (wa, wb) = (w(z0), w(z1));
        if wa > rf && wb > rf {
            1.0
        } else if wa <= rf && wb <= rf {
            0.0
        } else {
            let t = (rf - wa) / (wb - wa);
            if wa > rf { t } else { 1.0 - t }
        }
    };
    let ap_z = |zf: f64| -> f64 {
        let x = w(zf).clamp(r0, r1);
        ((x * x - r0 * r0) / (r1 * r1 - r0 * r0)).clamp(0.0, 1.0)
    };
    (kappa, [ap_r(r0), ap_r(r1), ap_z(z0), ap_z(z1)])
}

// Covers FND-3 §3.1/§3.3 (S9, in-kind): a converging-diverging toy chamber
// with a θ-varying wall bump (κ > 0 in every sector) marches 50+ audited
// cold-flow steps and conserves. Built by an analytic supplier closure —
// crates/grid/src/geom3d is ABSENT in this worktree; the voxelizer-output
// integration happens with the main S9 session.
#[test]
fn fnd3_s31_toy_chamber_cold_flow_analytic_supplier() {
    let nt = 8u32;
    let (h, n_r, n_z) = (0.05, 8usize, 16usize);
    let len = h * n_z as f64;
    // The revolved profile: piecewise-linear interpolant of a sinusoidal
    // throat (radius 0.32 → 0.20 → 0.32), exact fractions per cell.
    let wall = move |z: f64| 0.32 - 0.12 * (std::f64::consts::PI * z / len).sin();
    // A mild multiplicative θ-bump on cut cells (κ stays in (0, 1]) and a
    // θ-face aperture ring profile — a genuinely 3-D wall.
    let bump = |j: u32| 0.95 + 0.00625 * f64::from(j);
    let spec = GridSpec {
        r_min: 0.0,
        dr: h,
        n_r,
        z_min: 0.0,
        dz: h,
        n_z,
        n_theta_max: nt,
        axisymmetry_assertion: false,
    };
    let g0 = Grid::build_with_geometry_theta(
        spec,
        EULER_FIELDS,
        move |i_r, j, i_z| {
            let (r0, r1) = (i_r as f64 * h, (i_r + 1) as f64 * h);
            let (z0, z1) = (i_z as f64 * h, (i_z + 1) as f64 * h);
            let (wa, wb) = (wall(z0), wall(z1));
            // The per-cell linearization runs in LOCAL z (ζ ∈ [0, h]) so
            // w(0) = wa and w(h) ≈ wb exactly at the endpoints…
            let (kappa_rz, mut ap) = linear_wall_geom(wa, (wb - wa) / h, r0, r1, 0.0, h);
            // …but a shared z-face's aperture must be the SAME BITS on both
            // sides, so it is overridden with the face-canonical fraction
            // of wall(z_face) — both adjacent cells compute the identical
            // expression (the FND-3 §3.3 "compute a face from its
            // coordinates alone" rule; consistent with the covered-face
            // rule because κ = 0 ⇔ the wall sits below r0 at both faces).
            let frac_z = |wf: f64| -> f64 {
                let x = wf.clamp(r0, r1);
                ((x * x - r0 * r0) / (r1 * r1 - r0 * r0)).clamp(0.0, 1.0)
            };
            ap[2] = frac_z(wa);
            ap[3] = frac_z(wb);
            if kappa_rz <= 0.0 {
                return CellGeomTheta {
                    kappa: 0.0,
                    aperture: [0.0; 6],
                };
            }
            let cut = kappa_rz < 1.0;
            let kappa = if cut { kappa_rz * bump(j) } else { 1.0 };
            let ap_th = |fi: u32| {
                if cut {
                    kappa_rz * (0.9 + 0.0125 * f64::from(fi % nt))
                } else {
                    1.0
                }
            };
            CellGeomTheta {
                kappa,
                aperture: [ap[0], ap[1], ap[2], ap[3], ap_th(j), ap_th((j + 1) % nt)],
            }
        },
        |_, _| Region::Exterior,
    )
    .expect("toy chamber builds from the analytic supplier");
    let mut g = g0;
    let eos = GammaLaw { gamma: GAMMA };
    let f = EulerFields::resolve(&g).expect("fields");
    // Cold flow: an axial pressure gradient drains through the throat.
    fill_from_prim(&mut g, &f, &eos, |_, _, z| {
        prim6(1.2, 0.0, 0.0, 0.0, 1.0e5 * (1.3 - 0.5 * z / 0.8), 0.3)
    });
    let ids = f.ids();
    let mass0 = g.reduce_kappa_volume_weighted(ids[I_RHO]);
    let en0 = g.reduce_kappa_volume_weighted(ids[4]);
    let op = closed_box_op(eos);
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    for step in 0..55 {
        let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt stays evaluable");
        sdc.step_flow(&mut g, &flow, 0.0, dt)
            .unwrap_or_else(|e| panic!("audited cold-flow step {step} halted: {e}"));
    }
    let mass1 = g.reduce_kappa_volume_weighted(ids[I_RHO]);
    let en1 = g.reduce_kappa_volume_weighted(ids[4]);
    assert!(
        ((mass1 - mass0) / mass0).abs() < 1e-12,
        "mass drift {:.3e}",
        (mass1 - mass0) / mass0
    );
    assert!(
        ((en1 - en0) / en0).abs() < 1e-12,
        "energy drift {:.3e}",
        (en1 - en0) / en0
    );
    g.for_each_active_cell(|cell| {
        let rho = g.brick(cell.bi).field(ids[I_RHO])[cell.idx];
        assert!(
            rho.is_finite() && rho > 1e-3,
            "cell ({}, {}, θ{}) density {rho:.3e} ran away",
            cell.i_r,
            cell.i_z,
            cell.i_theta
        );
    });
}

// --- The voxelizer-fed toy chamber (the main-session integration) -----------

/// FND-3 §3.1/§3.3 → FND-2 §3.6 → SOLV-1 §3.3, end to end (plan S9): the
/// `stl_toy_chamber` mini-sim on the REAL sampled path. A converging-
/// diverging chamber is authored as a CSG solid (a big cylinder minus the
/// revolved gas cavity, plus a shallow axial RIB protruding into the gas
/// over one θ-side — a genuinely 3-D wall), voxelized by the geom3d kernel
/// (jittered-stratified fractions, six apertures, canonical face arrays),
/// ingested through `build_with_geometry_theta` (whose per-sector
/// validation — bitwise face coherence, covered-face rule, the sector-
/// coverage scope rule — must ACCEPT the sampled output as-is), and
/// marched 50+ audited cold-flow steps.
///
/// Face order (S10 unification): the voxelizer now emits in the grid's
/// `FaceDir::index` order `[r−, r+, z−, z+, θ−, θ+]` — the single owner —
/// so the `CellCut → CellGeomTheta` ingest below is a plain identity copy
/// (the S9 seam wart is retired; no permutation to get wrong).
#[test]
fn fnd3_s33_stl_toy_chamber_cold_flow_via_the_voxelizer() {
    use crucible_grid::geom3d::{Sdf, Solid, VoxelSpec, voxelize};

    let nt = 8u32;
    let (h, n_r, n_z) = (0.05, 8usize, 16usize);
    let len = h * n_z as f64;
    // The gas cavity: a revolved polygon of the sinusoidal-throat profile
    // (the gate-f wall), sampled at cell-edge z's so the polygon is the
    // profile the fractions see. Closed via the axis (s = 0).
    let wall = |z: f64| 0.32 - 0.12 * (std::f64::consts::PI * z / len).sin();
    let mut profile: Vec<(f64, f64)> = vec![(0.0, 0.0)];
    for iz in 0..=n_z {
        let z = iz as f64 * h;
        profile.push((wall(z), z));
    }
    profile.push((0.0, len));
    let cavity = Sdf::Revolved { profile };
    // The chamber body: a big cylinder covering the whole grid; the solid
    // is body − cavity, plus the rib: a shallow box along +x protruding
    // RIB_DEPTH into the gas over a mid-z span (θ-varying by construction;
    // shallow enough that no gas ring loses a whole sector — the builder's
    // scope rule validates exactly that, loudly).
    // The rib's inner radius sits 0.019 inside the local wall (wall ≈
    // 0.200–0.209 over the rib's z-span), so it genuinely narrows the gas
    // in its sectors without severing any (≈0.4·h penetration).
    let body = Sdf::CylinderZ {
        center: [0.0, 0.0, 0.5 * len],
        radius: 1.0,
        half_height: 0.5 * len,
    };
    let rib = Sdf::BoxAxis {
        min: [0.19, -0.02, 0.30],
        max: [1.0, 0.02, 0.50],
    };
    let solid = Solid::Csg(Sdf::Union(
        Box::new(Sdf::Subtraction(Box::new(body), Box::new(cavity))),
        Box::new(rib),
    ));

    let vspec = VoxelSpec {
        r_min: 0.0,
        dr: h,
        n_r,
        z_min: 0.0,
        dz: h,
        n_z,
        n_theta: nt,
        seed: 9,
    };
    let world = voxelize(&solid, &vspec).expect("toy chamber voxelizes");
    println!(
        "stl_toy_chamber: eps_alpha = {:.3e} (declared, manifest-bound)",
        world.eps_alpha
    );

    let spec = GridSpec {
        r_min: 0.0,
        dr: h,
        n_r,
        z_min: 0.0,
        dz: h,
        n_z,
        n_theta_max: nt,
        axisymmetry_assertion: false,
    };
    let g0 = Grid::build_with_geometry_theta(
        spec,
        EULER_FIELDS,
        |i_r, j, i_z| {
            let c = &world.cells[world.cell_index(i_r, j as usize, i_z)];
            // S10 face-order unification: the voxelizer now emits in
            // FaceDir::index order, so the ingest is a plain copy.
            CellGeomTheta {
                kappa: c.kappa,
                aperture: c.aperture,
            }
        },
        |_, _| Region::Exterior,
    )
    .expect("the sampled path satisfies the per-sector ingest validation as-is");
    let mut g = g0;

    // The rib is genuinely 3-D: some (r,z) ring must carry θ-varying κ.
    let mut theta_varying = 0usize;
    for i_r in 0..n_r {
        for i_z in 0..n_z {
            let k0 = g.kappa_at(i_r, 0, i_z);
            if (0..nt).any(|j| g.kappa_at(i_r, j, i_z) != k0) && g.is_active(i_r, i_z) {
                theta_varying += 1;
            }
        }
    }
    assert!(
        theta_varying > 0,
        "the rib produced no θ-varying gas ring — the fixture is degenerate"
    );

    let eos = GammaLaw { gamma: GAMMA };
    let f = EulerFields::resolve(&g).expect("fields");
    fill_from_prim(&mut g, &f, &eos, |_, _, z| {
        prim6(1.2, 0.0, 0.0, 0.0, 1.0e5 * (1.3 - 0.5 * z / 0.8), 0.3)
    });
    let ids = f.ids();
    let mass0 = g.reduce_kappa_volume_weighted(ids[I_RHO]);
    let en0 = g.reduce_kappa_volume_weighted(ids[4]);
    let op = closed_box_op(eos);
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    for step in 0..55 {
        let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt stays evaluable");
        sdc.step_flow(&mut g, &flow, 0.0, dt)
            .unwrap_or_else(|e| panic!("audited cold-flow step {step} halted: {e}"));
    }
    let mass1 = g.reduce_kappa_volume_weighted(ids[I_RHO]);
    let en1 = g.reduce_kappa_volume_weighted(ids[4]);
    assert!(
        ((mass1 - mass0) / mass0).abs() < 1e-12,
        "mass drift {:.3e}",
        (mass1 - mass0) / mass0
    );
    assert!(
        ((en1 - en0) / en0).abs() < 1e-12,
        "energy drift {:.3e}",
        (en1 - en0) / en0
    );
    g.for_each_active_cell(|cell| {
        let rho = g.brick(cell.bi).field(ids[I_RHO])[cell.idx];
        assert!(
            rho.is_finite() && rho > 1e-3,
            "cell ({}, {}, θ{}) density {rho:.3e} ran away",
            cell.i_r,
            cell.i_z,
            cell.i_theta
        );
    });
}
