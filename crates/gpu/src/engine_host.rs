//! Host side of the resident engine (`cuda/residency_engine.cu`): the FFI
//! contract, the geometry-time device tables built from the real assembly,
//! the table marshaling, and the run SCHEDULE pieces replicated from
//! `crucible_engine::run::run` (the injector start ramp, the pump-down
//! ambient, the igniter kernel, the audit reference scales). This module is
//! the ONE owner of that replication — both the cross-check bin (which scores
//! it against the CPU `Sdc::step`) and the overnight harness use it.
#![allow(clippy::needless_range_loop)] // dense (j, rz, k) index arithmetic mirrors the device layout

use crucible_engine::assembly::EngineSpec;
use crucible_engine::eos_sel::ChemEos;
use crucible_engine::run::{IGNITER_RAMP_FRAC, open_blend_tables, open_pinned_table, plane_area};
use crucible_grid::{FaceDir, Grid};
use crucible_solvers::euler::{
    BurnBlendEos, Combustion, Cons, EosLaw, FlowLedger, I_EN, NCOMP, NPRIM,
    xcheck_srd_neighborhoods,
};
use crucible_solvers::sdc::{AuditSpec, K_AUDIT};
use crucible_tables::{ColumnMarshal, Table};

pub const GK_DOMAIN_REFLECT: i32 = 0;
pub const GK_WALL_SLIP: i32 = 3;
pub const GK_AXIS: i32 = 4;
pub const GK_INFLOW: i32 = 5;
pub const GK_OUTFLOW: i32 = 6;

/// The step's reduced audit operands (`gpu_engine_step`'s `out`): 86 scalars.
pub const N_STEP_OUT: usize = 86;

#[repr(C)]
pub struct HostCol {
    p: *const f64,
    np: i32,
    lp: i32,
    h: *const f64,
    nh: i32,
    lh: i32,
    z: *const f64,
    nz: i32,
    lz: i32,
    str: *const i32,
    data: *const f64,
    vlog: i32,
}
#[repr(C)]
pub struct HostWorld {
    n_r: i32,
    n_z: i32,
    nt: i32,
    r_min: f64,
    dr: f64,
    z_min: f64,
    dz: f64,
    has_geom: i32,
    act: *const i32,
    kappa: *const f64,
    ap: *const f64,
    rs_r: *const i32,
    rl_r: *const i32,
    klo_r: *const i32,
    khi_r: *const i32,
    nlo_r: *const f64,
    nhi_r: *const f64,
    rs_z: *const i32,
    rl_z: *const i32,
    klo_z: *const i32,
    khi_z: *const i32,
    nlo_z: *const f64,
    nhi_z: *const f64,
}
#[repr(C)]
pub struct HostSrd {
    ns: i32,
    na: i32,
    mem_off: *const i32,
    mem: *const i32,
    mem_kv: *const f64,
    cnt: *const i32,
    aff: *const i32,
    aff_owner: *const i32,
    inv_off: *const i32,
    inv: *const i32,
}
#[repr(C)]
pub struct HostTables {
    urho: HostCol,
    usnd: HostCol,
    utmp: HostCol,
    brho: HostCol,
    bsnd: HostCol,
    flame: HostCol,
    delay: HostCol,
    pu_lo: f64,
    pu_hi: f64,
    hu_floor: f64,
    hu_ceil: f64,
    pb_lo: f64,
    pb_hi: f64,
    hb_lo: f64,
    hb_hi: f64,
    h_off: f64,
    p_floor: f64,
    tu_floor: f64,
}
#[repr(C)]
pub struct HostIgniter {
    armed: i32,
    r_lo: f64,
    r_hi: f64,
    z_lo: f64,
    z_hi: f64,
    theta_gated: i32,
    nt: i32,
    jt: i32,
    t_on: f64,
    t_off: f64,
    t_ramp: f64,
    q0: f64,
}

unsafe extern "C" {
    fn gpu_engine_create(
        hw: *const HostWorld,
        hs: *const HostSrd,
        ht: *const HostTables,
        hg: *const HostIgniter,
        wrinkling: f64,
        theta: f64,
        h_total: f64,
        c_frac: f64,
    ) -> *mut std::ffi::c_void;
    fn gpu_engine_destroy(h: *mut std::ffi::c_void);
    fn gpu_engine_upload(h: *mut std::ffi::c_void, cons: *const f64, hint: *const f64, primed: i32);
    fn gpu_engine_download(h: *mut std::ffi::c_void, cons: *mut f64, hint: *mut f64) -> i32;
    fn gpu_engine_stable_dt(h: *mut std::ffi::c_void, t: f64, cfl: f64, bad: *mut i32) -> f64;
    fn gpu_engine_step(
        h: *mut std::ffi::c_void,
        t: f64,
        dt: f64,
        mdot_per_area: f64,
        p_amb_t: f64,
        p_amb_t1: f64,
        out: *mut f64,
    ) -> i32;
    fn gpu_engine_rhs(
        h: *mut std::ffi::c_void,
        t: f64,
        mdot_per_area: f64,
        p_amb: f64,
        rate: *mut f64,
        led: *mut f64,
    ) -> i32;
}

/// Owning host mirror of one marshaled column.
struct ColBuf {
    axes: [Vec<f64>; 3],
    is_log: [i32; 3],
    strides: Vec<i32>,
    data: Vec<f64>,
    vlog: i32,
}
impl ColBuf {
    fn new(c: &ColumnMarshal) -> Self {
        ColBuf {
            axes: [
                c.axis_points[0].clone(),
                c.axis_points[1].clone(),
                c.axis_points[2].clone(),
            ],
            is_log: [
                c.axis_is_log[0] as i32,
                c.axis_is_log[1] as i32,
                c.axis_is_log[2] as i32,
            ],
            strides: c.strides.iter().map(|&s| s as i32).collect(),
            data: c.data.clone(),
            vlog: c.value_is_log as i32,
        }
    }
    fn host(&self) -> HostCol {
        HostCol {
            p: self.axes[0].as_ptr(),
            np: self.axes[0].len() as i32,
            lp: self.is_log[0],
            h: self.axes[1].as_ptr(),
            nh: self.axes[1].len() as i32,
            lh: self.is_log[1],
            z: self.axes[2].as_ptr(),
            nz: self.axes[2].len() as i32,
            lz: self.is_log[2],
            str: self.strides.as_ptr(),
            data: self.data.as_ptr(),
            vlog: self.vlog,
        }
    }
}

pub fn read_rel(path: &str) -> Result<String, String> {
    std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))
}

/// The loaded run: the assembled spec + the pinned tables it needs (owned
/// here so the borrowing EOS/combustion occupants can be built by the caller).
pub struct LoadedRun {
    pub name: String,
    pub spec: EngineSpec,
    pub table: Table,
    pub unburnt: Table,
    pub ignition: Table,
}

/// Load + assemble a config exactly as the CLI does (paths repo-root
/// relative — the caller sets the working directory).
pub fn load_run(config: &str) -> Result<LoadedRun, String> {
    let author = read_rel(config)?;
    let registry = crucible_engine::registry();
    let loaded = crucible_config::load_str_with_sidecars(&author, &registry, &read_rel)
        .map_err(|d| format!("config load: {d}"))?;
    let name = if loaded.resolved.meta.name.is_empty() {
        "run".to_string()
    } else {
        loaded.resolved.meta.name.clone()
    };
    let spec = crucible_engine::assembly::assemble(&loaded, &read_rel)
        .map_err(|e| format!("assembly: {e}"))?;
    if spec.blend.is_none() || spec.igniter.is_none() {
        return Err(
            "the resident engine needs combustion_blend + spark_igniter (the ◆C3/◆C4 start core)"
                .into(),
        );
    }
    if spec.liner.is_some()
        || spec.wall_law.is_some()
        || spec.gas_diffusion
        || spec.turbopump.is_some()
    {
        return Err("the resident engine is the adiabatic flow+combustion start core (no walls, no F_visc, no turbopump)".into());
    }
    let table = open_pinned_table(&spec)?;
    let (unburnt, ignition) = open_blend_tables(&spec)?.ok_or("blend tables not selected")?;
    Ok(LoadedRun {
        name,
        spec,
        table,
        unburnt,
        ignition,
    })
}

/// The igniter kernel (run.rs `IgniterKernel`), a pure function of the spec + grid.
#[derive(Debug, Clone, Copy)]
pub struct IgniterParams {
    pub r_lo: f64,
    pub r_hi: f64,
    pub z_lo: f64,
    pub z_hi: f64,
    pub theta_gate: Option<(u32, u32)>,
    pub t_on: f64,
    pub t_off: f64,
    pub t_ramp: f64,
    pub q0: f64,
}

/// The run schedule (run.rs): the ramp, the pump-down ambient, the injector
/// face flux, the horizon, the audit reference scales, the igniter.
#[derive(Debug, Clone)]
pub struct Schedule {
    pub t_final: f64,
    pub t_ramp_window: f64,
    pub t_pump_window: f64,
    pub p_floor: f64,
    pub p_fill: f64,
    pub mdot_design: f64,
    pub a_inlet: f64,
    pub h_total: f64,
    pub c_frac: f64,
    pub a_ref: f64,
    pub span: f64,
    pub cfl: f64,
    pub t_dwell_s: Option<f64>,
    pub igniter: IgniterParams,
    pub audit: AuditSpec,
    pub n_cells: usize,
}
impl Schedule {
    /// Build from the spec AFTER the fill has been written into the grid
    /// (the audit reference scales are the fill state's stored magnitudes).
    pub fn new(spec: &EngineSpec, eos: &ChemEos<'_>, u_fill: &Cons) -> Result<Self, String> {
        let ids = spec.fields.ids();
        let a_inlet = plane_area(&spec.grid, 0, FaceDir::ZMinus);
        if a_inlet <= 0.0 {
            return Err("no active gas cells on the injector plane".into());
        }
        let span = spec.contour.z_max() - spec.contour.z_min();
        let w = eos
            .prim_checked(u_fill)
            .map_err(|e| format!("fill state sound speed: {e}"))?;
        let a_ref = eos.sound_speed_w(&w);
        let t_final = spec.flowthroughs * span / a_ref;
        let ig = spec.igniter.as_ref().ok_or("igniter")?;
        let (r_lo, r_hi) = (ig.r_m - ig.half_width_m, ig.r_m + ig.half_width_m);
        let (z_lo, z_hi) = (ig.z_m - ig.half_width_m, ig.z_m + ig.half_width_m);
        let nt_world = spec.grid.spec().n_theta_max;
        let theta_gate = if nt_world > 1 {
            let dth = std::f64::consts::TAU / f64::from(nt_world);
            Some((nt_world, ((ig.theta_rad / dth).floor() as u32) % nt_world))
        } else {
            None
        };
        let mut vol = 0.0f64;
        spec.grid.for_each_active_cell(|c| {
            let in_box = c.r > r_lo && c.r < r_hi && c.z > z_lo && c.z < z_hi;
            let in_theta = theta_gate.is_none_or(|(_, jt)| c.i_theta == jt);
            if in_box && in_theta {
                vol += spec.grid.kappa(c.i_r, c.i_z) * spec.grid.cell_volume(c.i_r, nt_world);
            }
        });
        if vol <= 0.0 {
            return Err("spark igniter kernel covers no gas cell".into());
        }
        let t_ramp = IGNITER_RAMP_FRAC * ig.window_s;
        let igniter = IgniterParams {
            r_lo,
            r_hi,
            z_lo,
            z_hi,
            theta_gate,
            t_on: ig.window_start_s,
            t_off: ig.window_start_s + ig.window_s,
            t_ramp,
            q0: ig.energy_j / (vol * (ig.window_s - 0.5 * t_ramp)),
        };
        let mass_ref = spec.grid.reduce_kappa_volume_weighted(ids[0]).abs();
        let energy_gas = spec.grid.reduce_kappa_volume_weighted(ids[4]).abs();
        let audit = AuditSpec {
            k_audit: K_AUDIT,
            ref_scale: [mass_ref, mass_ref * a_ref, energy_gas, mass_ref],
        };
        let n_cells = spec
            .grid
            .bricks()
            .iter()
            .map(|b| {
                (b.mask().count_ones() + b.solid_mask().count_ones()) as usize
                    * b.n_theta() as usize
            })
            .sum();
        Ok(Schedule {
            t_final,
            t_ramp_window: spec.injector_ramp_flowthroughs * span / a_ref,
            t_pump_window: spec.pumpdown_flowthroughs * span / a_ref,
            p_floor: spec.p_amb_floor_pa,
            p_fill: spec.fill_p_pa,
            mdot_design: spec.injector.mdot_kg_per_s,
            a_inlet,
            h_total: spec.injector.h_inj_j_per_kg,
            c_frac: spec.injector.z_frac,
            a_ref,
            span,
            cfl: spec.cfl,
            t_dwell_s: spec
                .verdict
                .as_ref()
                .map(|v| v.t_dwell_flowthroughs * span / a_ref),
            igniter,
            audit,
            n_cells,
        })
    }
    pub fn ramp_frac(&self, t: f64) -> f64 {
        if self.t_ramp_window <= 0.0 {
            1.0
        } else {
            ((t - self.t_pump_window) / self.t_ramp_window).clamp(0.0, 1.0)
        }
    }
    /// The delivered ṁ at `t` (open mode: the design ṁ under the ramp).
    pub fn mdot_delivered(&self, t: f64) -> f64 {
        self.mdot_design * self.ramp_frac(t)
    }
    pub fn mdot_per_area(&self, t: f64) -> f64 {
        self.mdot_delivered(t) / self.a_inlet
    }
    pub fn p_amb(&self, t: f64) -> f64 {
        if self.t_pump_window <= 0.0 {
            return self.p_floor;
        }
        let s = (t / self.t_pump_window).clamp(0.0, 1.0);
        self.p_fill * (self.p_floor / self.p_fill).powf(s)
    }
    /// The igniter deposit as the CPU operator's external source closure.
    pub fn source_fn(&self) -> impl Fn(f64, f64, f64, f64) -> Cons + Sync + 'static {
        let k = self.igniter;
        move |r: f64, th: f64, z: f64, t: f64| -> Cons {
            let mut src = [0.0; NCOMP];
            if t >= k.t_on && t < k.t_off && r > k.r_lo && r < k.r_hi && z > k.z_lo && z < k.z_hi {
                let in_theta = match k.theta_gate {
                    None => true,
                    Some((nt, jt)) => {
                        let dth = std::f64::consts::TAU / f64::from(nt);
                        let j = (th / dth).floor().rem_euclid(f64::from(nt)) as u32;
                        j == jt
                    }
                };
                if in_theta {
                    let ramp = if k.t_ramp > 0.0 {
                        ((t - k.t_on) / k.t_ramp).min(1.0)
                    } else {
                        1.0
                    };
                    src[I_EN] = k.q0 * ramp;
                }
            }
            src
        }
    }
}

/// The dense device layout of a grid: θ-plane-major, `c = j·NRZ + i_r·n_z + i_z`.
#[derive(Debug, Clone, Copy)]
pub struct Layout {
    pub n_r: usize,
    pub n_z: usize,
    pub nt: usize,
    pub nrz: usize,
    pub n: usize,
}
impl Layout {
    pub fn of(g: &Grid) -> Self {
        let (n_r, n_z) = (g.spec().n_r, g.spec().n_z);
        let nt = g.spec().n_theta_max as usize;
        Layout {
            n_r,
            n_z,
            nt,
            nrz: n_r * n_z,
            n: nt * n_r * n_z,
        }
    }
    #[inline]
    pub fn dense(&self, j: usize, i_r: usize, i_z: usize) -> usize {
        j * self.nrz + i_r * self.n_z + i_z
    }
    pub fn read_dense(&self, g: &Grid, ids: &[crucible_grid::FieldId; NCOMP]) -> Vec<f64> {
        let mut v = vec![0.0f64; self.n * NCOMP];
        g.for_each_active_cell(|cell| {
            let cc = self.dense(cell.i_theta as usize, cell.i_r, cell.i_z);
            let b = g.brick(cell.bi);
            for k in 0..NCOMP {
                v[cc * NCOMP + k] = b.field(ids[k])[cell.idx];
            }
        });
        v
    }
    pub fn write_dense(&self, g: &mut Grid, ids: &[crucible_grid::FieldId; NCOMP], cons: &[f64]) {
        let cells: Vec<(usize, usize, usize)> = {
            let mut v = vec![];
            g.for_each_active_cell(|c| {
                v.push((c.bi, c.idx, self.dense(c.i_theta as usize, c.i_r, c.i_z)))
            });
            v
        };
        for k in 0..NCOMP {
            for &(bi, idx, cc) in &cells {
                g.brick_field_mut(bi, ids[k])[idx] = cons[cc * NCOMP + k];
            }
        }
    }
}

/// The resident device engine: the handle + every host buffer it was built from.
pub struct DeviceEngine {
    handle: *mut std::ffi::c_void,
    pub layout: Layout,
    pub act: Vec<i32>,
    // kept alive for the FFI structs (uploaded once at create; retained for the
    // cross-check's own bookkeeping)
    _bufs: HostBufs,
    _cols: Vec<ColBuf>,
}
struct HostBufs {
    kappa: Vec<f64>,
    ap: Vec<f64>,
    rs_r: Vec<i32>,
    rl_r: Vec<i32>,
    klo_r: Vec<i32>,
    khi_r: Vec<i32>,
    nlo_r: Vec<f64>,
    nhi_r: Vec<f64>,
    rs_z: Vec<i32>,
    rl_z: Vec<i32>,
    klo_z: Vec<i32>,
    khi_z: Vec<i32>,
    nlo_z: Vec<f64>,
    nhi_z: Vec<f64>,
    mem_off: Vec<i32>,
    mem: Vec<i32>,
    mem_kv: Vec<f64>,
    cnt: Vec<i32>,
    aff: Vec<i32>,
    aff_owner: Vec<i32>,
    inv_off: Vec<i32>,
    inv: Vec<i32>,
}

impl DeviceEngine {
    /// Build the device engine for the loaded run's grid + occupants (the
    /// run's own BC/wall settings: axis, Reflecting r_outer, injector inflow
    /// at z_lo, pressure outflow at z_hi, slip walls about the contour normal).
    pub fn new(
        spec: &EngineSpec,
        blend: &BurnBlendEos<'_>,
        comb: &Combustion<'_>,
        sched: &Schedule,
    ) -> Result<Self, String> {
        let g = &spec.grid;
        let lay = Layout::of(g);
        let (n_r, n_z, nt, nrz, n) = (lay.n_r, lay.n_z, lay.nt, lay.nrz, lay.n);
        let (r_min, dr, z0, dz) = (g.spec().r_min, g.spec().dr, g.spec().z_min, g.spec().dz);
        let mut act = vec![0i32; nrz];
        let mut kappa = vec![0.0f64; n];
        let mut ap = vec![0.0f64; 6 * n];
        let dirs = [
            FaceDir::RMinus,
            FaceDir::RPlus,
            FaceDir::ZMinus,
            FaceDir::ZPlus,
            FaceDir::ThetaMinus,
            FaceDir::ThetaPlus,
        ];
        for i_r in 0..n_r {
            for i_z in 0..n_z {
                if !g.is_active(i_r, i_z) {
                    continue;
                }
                act[i_r * n_z + i_z] = 1;
                for j in 0..nt {
                    let cc = lay.dense(j, i_r, i_z);
                    kappa[cc] = g.kappa_at(i_r, j as u32, i_z);
                    for (d, &dir) in dirs.iter().enumerate() {
                        ap[d * n + cc] = g.aperture_at(i_r, j as u32, i_z, dir);
                    }
                }
            }
        }
        let contour = spec.contour.clone();
        let normal_fn = move |r: f64, z: f64| contour.wall_normal(r, z);
        let mut b = HostBufs {
            kappa,
            ap,
            rs_r: vec![0; nrz],
            rl_r: vec![0; nrz],
            klo_r: vec![0; nrz],
            khi_r: vec![0; nrz],
            nlo_r: vec![0.0; 2 * nrz],
            nhi_r: vec![0.0; 2 * nrz],
            rs_z: vec![0; nrz],
            rl_z: vec![0; nrz],
            klo_z: vec![0; nrz],
            khi_z: vec![0; nrz],
            nlo_z: vec![0.0; 2 * nrz],
            nhi_z: vec![0.0; 2 * nrz],
            mem_off: vec![0],
            mem: vec![],
            mem_kv: vec![],
            cnt: vec![1; n],
            aff: vec![],
            aff_owner: vec![],
            inv_off: vec![0],
            inv: vec![],
        };
        let on_axis = r_min == 0.0;
        for i_z in 0..n_z {
            let z = g.z_center(i_z);
            let mut i = 0usize;
            while i < n_r {
                if act[i * n_z + i_z] == 0 {
                    i += 1;
                    continue;
                }
                let start = i;
                while i < n_r && act[i * n_z + i_z] == 1 {
                    i += 1;
                }
                let len = i - start;
                let (klo, nlo) = if start == 0 && on_axis {
                    (GK_AXIS, (0.0, 0.0))
                } else if start == 0 {
                    (GK_DOMAIN_REFLECT, (0.0, 0.0))
                } else {
                    (GK_WALL_SLIP, normal_fn(r_min + start as f64 * dr, z))
                };
                let (khi, nhi) = if start + len == n_r {
                    (GK_DOMAIN_REFLECT, (0.0, 0.0))
                } else {
                    (
                        GK_WALL_SLIP,
                        normal_fn(r_min + (start + len) as f64 * dr, z),
                    )
                };
                for ii in start..start + len {
                    let rz = ii * n_z + i_z;
                    b.rs_r[rz] = start as i32;
                    b.rl_r[rz] = len as i32;
                    b.klo_r[rz] = klo;
                    b.khi_r[rz] = khi;
                    b.nlo_r[2 * rz] = nlo.0;
                    b.nlo_r[2 * rz + 1] = nlo.1;
                    b.nhi_r[2 * rz] = nhi.0;
                    b.nhi_r[2 * rz + 1] = nhi.1;
                }
            }
        }
        for i_r in 0..n_r {
            let r = g.r_center(i_r);
            let mut i = 0usize;
            while i < n_z {
                if act[i_r * n_z + i] == 0 {
                    i += 1;
                    continue;
                }
                let start = i;
                while i < n_z && act[i_r * n_z + i] == 1 {
                    i += 1;
                }
                let len = i - start;
                let (klo, nlo) = if start == 0 {
                    (GK_INFLOW, (0.0, 0.0))
                } else {
                    (GK_WALL_SLIP, normal_fn(r, z0 + start as f64 * dz))
                };
                let (khi, nhi) = if start + len == n_z {
                    (GK_OUTFLOW, (0.0, 0.0))
                } else {
                    (GK_WALL_SLIP, normal_fn(r, z0 + (start + len) as f64 * dz))
                };
                for ii in start..start + len {
                    let rz = i_r * n_z + ii;
                    b.rs_z[rz] = start as i32;
                    b.rl_z[rz] = len as i32;
                    b.klo_z[rz] = klo;
                    b.khi_z[rz] = khi;
                    b.nlo_z[2 * rz] = nlo.0;
                    b.nlo_z[2 * rz + 1] = nlo.1;
                    b.nhi_z[2 * rz] = nhi.0;
                    b.nhi_z[2 * rz + 1] = nhi.1;
                }
            }
        }
        // SRD tables in the CPU's build order.
        let small = xcheck_srd_neighborhoods(g).map_err(|e| format!("SRD neighborhoods: {e}"))?;
        let mut inv_lists: std::collections::BTreeMap<usize, (bool, Vec<i32>)> =
            std::collections::BTreeMap::new();
        for (s, (cells, kv)) in small.iter().enumerate() {
            for (m, &(rr, zz, j)) in cells.iter().enumerate() {
                let cc = lay.dense(j as usize, rr, zz);
                b.mem.push(cc as i32);
                b.mem_kv.push(kv[m]);
                if m > 0 {
                    b.cnt[cc] += 1;
                }
                let e = inv_lists.entry(cc).or_insert((false, vec![]));
                if m == 0 {
                    e.0 = true;
                }
                e.1.push(s as i32);
            }
            b.mem_off.push(b.mem.len() as i32);
        }
        for (&cc, (owner, list)) in &inv_lists {
            b.aff.push(cc as i32);
            b.aff_owner.push(i32::from(*owner));
            b.inv.extend_from_slice(list);
            b.inv_off.push(b.inv.len() as i32);
        }
        let hw = HostWorld {
            n_r: n_r as i32,
            n_z: n_z as i32,
            nt: nt as i32,
            r_min,
            dr,
            z_min: z0,
            dz,
            has_geom: i32::from(g.has_cut_geometry()),
            act: act.as_ptr(),
            kappa: b.kappa.as_ptr(),
            ap: b.ap.as_ptr(),
            rs_r: b.rs_r.as_ptr(),
            rl_r: b.rl_r.as_ptr(),
            klo_r: b.klo_r.as_ptr(),
            khi_r: b.khi_r.as_ptr(),
            nlo_r: b.nlo_r.as_ptr(),
            nhi_r: b.nhi_r.as_ptr(),
            rs_z: b.rs_z.as_ptr(),
            rl_z: b.rl_z.as_ptr(),
            klo_z: b.klo_z.as_ptr(),
            khi_z: b.khi_z.as_ptr(),
            nlo_z: b.nlo_z.as_ptr(),
            nhi_z: b.nhi_z.as_ptr(),
        };
        let hs = HostSrd {
            ns: small.len() as i32,
            na: b.aff.len() as i32,
            mem_off: b.mem_off.as_ptr(),
            mem: b.mem.as_ptr(),
            mem_kv: b.mem_kv.as_ptr(),
            cnt: b.cnt.as_ptr(),
            aff: b.aff.as_ptr(),
            aff_owner: b.aff_owner.as_ptr(),
            inv_off: b.inv_off.as_ptr(),
            inv: b.inv.as_ptr(),
        };
        let (um, bm, h_off) = blend.xcheck_branches_marshal();
        let (flame_m, p_floor_ign, tu_floor) = comb.ignition.xcheck_marshal_flame();
        let delay_m = comb.ignition.xcheck_marshal_delay();
        let cols = vec![
            ColBuf::new(&um.rho),
            ColBuf::new(&um.sound),
            ColBuf::new(&um.temperature),
            ColBuf::new(&bm.rho),
            ColBuf::new(&bm.sound),
            ColBuf::new(&flame_m),
            ColBuf::new(&delay_m),
        ];
        let ht = HostTables {
            urho: cols[0].host(),
            usnd: cols[1].host(),
            utmp: cols[2].host(),
            brho: cols[3].host(),
            bsnd: cols[4].host(),
            flame: cols[5].host(),
            delay: cols[6].host(),
            pu_lo: um.p_env.0,
            pu_hi: um.p_env.1,
            hu_floor: um.h_env.0,
            hu_ceil: um.h_env.1,
            pb_lo: bm.p_env.0,
            pb_hi: bm.p_env.1,
            hb_lo: bm.h_env.0,
            hb_hi: bm.h_env.1,
            h_off,
            p_floor: p_floor_ign,
            tu_floor,
        };
        let ig = &sched.igniter;
        let hg = HostIgniter {
            armed: 1,
            r_lo: ig.r_lo,
            r_hi: ig.r_hi,
            z_lo: ig.z_lo,
            z_hi: ig.z_hi,
            theta_gated: i32::from(ig.theta_gate.is_some()),
            nt: nt as i32,
            jt: ig.theta_gate.map_or(0, |(_, jt)| jt as i32),
            t_on: ig.t_on,
            t_off: ig.t_off,
            t_ramp: ig.t_ramp,
            q0: ig.q0,
        };
        let handle = unsafe {
            gpu_engine_create(
                &hw,
                &hs,
                &ht,
                &hg,
                comb.wrinkling,
                comb.theta,
                sched.h_total,
                sched.c_frac,
            )
        };
        if handle.is_null() {
            return Err("gpu_engine_create failed".into());
        }
        Ok(DeviceEngine {
            handle,
            layout: lay,
            act,
            _bufs: b,
            _cols: cols,
        })
    }
    pub fn n_small(&self) -> usize {
        self._bufs.mem_off.len() - 1
    }
    /// Upload a dense state; `hint` = the warm-start prim cache (None = cold).
    pub fn upload(&self, cons: &[f64], hint: Option<&[f64]>) {
        assert_eq!(cons.len(), self.layout.n * NCOMP);
        unsafe {
            match hint {
                Some(h) => {
                    assert_eq!(h.len(), self.layout.n * NPRIM);
                    gpu_engine_upload(self.handle, cons.as_ptr(), h.as_ptr(), 1)
                }
                None => gpu_engine_upload(self.handle, cons.as_ptr(), std::ptr::null(), 0),
            }
        }
    }
    /// Download the dense state and the prim cache; returns `primed`.
    pub fn download(&self, cons: &mut [f64], hint: &mut [f64]) -> bool {
        unsafe { gpu_engine_download(self.handle, cons.as_mut_ptr(), hint.as_mut_ptr()) != 0 }
    }
    pub fn stable_dt(&self, t: f64, cfl: f64) -> Result<f64, String> {
        let mut bad = 0i32;
        let dt = unsafe { gpu_engine_stable_dt(self.handle, t, cfl, &mut bad) };
        if bad != 0 || dt <= 0.0 || !dt.is_finite() {
            return Err(
                "stable_dt refused (non-physical state or the scale-separation guard)".into(),
            );
        }
        Ok(dt)
    }
    /// One audited SDC step; returns the reduced audit operands.
    pub fn step(&self, sched: &Schedule, t: f64, dt: f64) -> Result<StepOut, String> {
        let mut out = vec![0.0f64; N_STEP_OUT];
        let bad = unsafe {
            gpu_engine_step(
                self.handle,
                t,
                dt,
                sched.mdot_per_area(t),
                sched.p_amb(t),
                sched.p_amb(t + dt),
                out.as_mut_ptr(),
            )
        };
        if bad != 0 {
            return Err("device step flagged a non-physical state / projection failure".into());
        }
        Ok(StepOut::from(&out))
    }
    /// Single-shot RHS (+ ledger) at the current state.
    pub fn rhs(&self, sched: &Schedule, t: f64) -> Result<(Vec<f64>, FlowLedger), String> {
        let mut rate = vec![0.0f64; self.layout.n * NCOMP];
        let mut led = vec![0.0f64; 4 * NCOMP];
        let bad = unsafe {
            gpu_engine_rhs(
                self.handle,
                t,
                sched.mdot_per_area(t),
                sched.p_amb(t),
                rate.as_mut_ptr(),
                led.as_mut_ptr(),
            )
        };
        if bad != 0 {
            return Err("device RHS flagged a halt".into());
        }
        Ok((rate, ledger_of(&led)))
    }
}
impl Drop for DeviceEngine {
    fn drop(&mut self) {
        unsafe { gpu_engine_destroy(self.handle) }
    }
}

/// The step's reduced operands, decoded.
#[derive(Debug, Clone)]
pub struct StepOut {
    pub before: ([f64; NCOMP], [f64; NCOMP]),
    pub after: ([f64; NCOMP], [f64; NCOMP]),
    pub l0: FlowLedger,
    pub l_last: FlowLedger,
    pub burn_applied: f64,
    pub burn_gross: f64,
}
impl StepOut {
    fn from(o: &[f64]) -> Self {
        let mut before = ([0.0f64; NCOMP], [0.0f64; NCOMP]);
        let mut after = ([0.0f64; NCOMP], [0.0f64; NCOMP]);
        for k in 0..NCOMP {
            before.0[k] = o[k];
            before.1[k] = o[NCOMP + k];
            after.0[k] = o[2 * NCOMP + k];
            after.1[k] = o[3 * NCOMP + k];
        }
        StepOut {
            before,
            after,
            l0: ledger_of(&o[28..56]),
            l_last: ledger_of(&o[56..84]),
            burn_applied: o[84],
            burn_gross: o[85],
        }
    }
}
pub fn ledger_of(o: &[f64]) -> FlowLedger {
    let mut l = FlowLedger::default();
    for k in 0..NCOMP {
        l.port_net[k] = o[k];
        l.port_abs[k] = o[NCOMP + k];
        l.src_net[k] = o[2 * NCOMP + k];
        l.src_abs[k] = o[3 * NCOMP + k];
    }
    l
}
