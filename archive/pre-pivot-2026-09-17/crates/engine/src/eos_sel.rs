//! The config-selected chemical EOS occupant (S7): **one dispatch point**
//! between the shifting-equilibrium [`TableEos`] (the stations' occupant,
//! `flow_shifting` alone) and the SOLV-4 §3.6 burn-progress [`BurnBlendEos`]
//! (`combustion_blend` selected — the startup marches). Selection is config
//! data resolved at assembly; every physics method delegates to the one
//! occupant, so there is no `if(mode)` inside any operator (Rule 13 — the
//! branch lives at the seam, once).
//!
//! The engine-side helper surface (`temperature_w`, `interrogation_php`,
//! `envelopes`, `fill_cons`) mirrors the two occupants' own methods so the
//! wall law, the FND-7 spine closures, and the readout/CSV writers read one
//! type. A shifting-only config runs the `Table` arm bit-for-bit as before.

use crucible_solvers::euler::{BurnBlendEos, Cons, EosLaw, FlowError, Prim, TableEos};
use crucible_solvers::transport::MediumState;
use crucible_tables::Table;

// One instance exists per run and it is cloned once into the flow operator;
// the size skew between the arms is irrelevant next to indirecting every
// EOS call through a Box.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone)]
pub enum ChemEos<'t> {
    /// Shifting equilibrium (`c ≡ 1` corner): the pinned surface alone.
    Table(TableEos<'t>),
    /// The burn-progress blend (unburnt ⊕ burnt) — the chemical regime can
    /// ignite (S6 physics, S7 config face).
    Blend(BurnBlendEos<'t>),
}

impl<'t> ChemEos<'t> {
    /// Bind the shifting occupant (the pre-S7 path, unchanged).
    pub fn bind_table(table: &'t Table, h_offset: f64) -> Result<Self, String> {
        let mut eos = TableEos::bind(table)?;
        eos.h_offset = h_offset;
        Ok(ChemEos::Table(eos))
    }

    /// Bind the blend: unburnt + burnt branches on the same (p, h, Z)
    /// schema; the S18 `h_offset` knockdown rides the burnt branch only.
    pub fn bind_blend(unburnt: &'t Table, burnt: &'t Table, h_offset: f64) -> Result<Self, String> {
        let mut blend = BurnBlendEos::new(TableEos::bind(unburnt)?, TableEos::bind(burnt)?);
        blend.h_offset = h_offset;
        Ok(ChemEos::Blend(blend))
    }

    /// The blend occupant, when selected (the combustion operator and the
    /// class-R reaction borrow it).
    pub fn blend(&self) -> Option<&BurnBlendEos<'t>> {
        match self {
            ChemEos::Blend(b) => Some(b),
            ChemEos::Table(_) => None,
        }
    }

    /// The diagnostic temperature at a projected primitive: equilibrium T
    /// (Table) or the SOLV-4 §3.6 mass-weighted blend (Blend). The wall
    /// law's and the CSV writers' operand.
    pub fn temperature_w(&self, w: &Prim) -> Result<f64, &'static str> {
        match self {
            ChemEos::Table(t) => t
                .temperature_w(w)
                .map_err(|_| "temperature query off the pinned surface"),
            ChemEos::Blend(b) => b.temperature_w(w),
        }
    }

    /// The FND-7 §3.3 spine interrogation coordinate (one former per
    /// occupant, so the EOS and the transport surface can never be read at
    /// different states).
    pub fn interrogation_php(&self, w: &Prim) -> MediumState {
        match self {
            ChemEos::Table(t) => t.interrogation_php(w),
            ChemEos::Blend(b) => b.interrogation_php(w),
        }
    }

    /// The declared (p, h, Z) envelope box (Blend: the branches'
    /// intersection) — assembly-time sizing/refusals.
    pub fn envelopes(&self) -> [(f64, f64); 3] {
        match self {
            ChemEos::Table(t) => t.envelopes(),
            ChemEos::Blend(b) => b.envelopes(),
        }
    }

    /// The (p, h, Z) box the occupant may INTERROGATE the spine over: the
    /// pinned surface's own envelope (Table) or the UNION of the two
    /// branches' envelopes (Blend — a pure-burnt cell roams the whole burnt
    /// surface, and the cold-side partition admits states below the
    /// intersection floor). The transport-coverage gate checks THIS box,
    /// not the projection intersection (S7 review finding).
    pub fn interrogation_envelopes(&self) -> [(f64, f64); 3] {
        match self {
            ChemEos::Table(t) => t.envelopes(),
            ChemEos::Blend(b) => b.interrogation_envelopes(),
        }
    }

    /// The quiescent fill state at (p, h, Z): shifting mode fills on the
    /// equilibrium surface (the inert `b ≡ 0` slot); the blend fills PURE
    /// UNBURNT (`b = 0`) — the physical pre-start chamber contents, which is
    /// exactly what the S5 cold branch made representable.
    pub fn fill_cons(&self, p: f64, h: f64, z: f64) -> Result<Cons, String> {
        match self {
            ChemEos::Table(t) => t
                .cons_from_phz(p, h, z, [0.0, 0.0, 0.0])
                .map_err(|e| format!("fill state: {e}")),
            ChemEos::Blend(b) => b
                .cons_from_phzb(p, h, z, 0.0, [0.0, 0.0, 0.0])
                .map_err(|e| format!("fill state (unburnt branch): {e}")),
        }
    }
}

impl EosLaw for ChemEos<'_> {
    fn prim_checked(&self, u: &Cons) -> Result<Prim, &'static str> {
        match self {
            ChemEos::Table(t) => t.prim_checked(u),
            ChemEos::Blend(b) => b.prim_checked(u),
        }
    }

    fn prim_checked_hinted(&self, u: &Cons, hint: Option<f64>) -> Result<Prim, &'static str> {
        match self {
            ChemEos::Table(t) => t.prim_checked_hinted(u, hint),
            ChemEos::Blend(b) => b.prim_checked_hinted(u, hint),
        }
    }

    fn prim_to_cons(&self, w: &Prim) -> Cons {
        match self {
            ChemEos::Table(t) => t.prim_to_cons(w),
            ChemEos::Blend(b) => b.prim_to_cons(w),
        }
    }

    fn total_energy(&self, w: &Prim) -> f64 {
        match self {
            ChemEos::Table(t) => t.total_energy(w),
            ChemEos::Blend(b) => b.total_energy(w),
        }
    }

    fn sound_speed_w(&self, w: &Prim) -> f64 {
        match self {
            ChemEos::Table(t) => t.sound_speed_w(w),
            ChemEos::Blend(b) => b.sound_speed_w(w),
        }
    }

    fn roe_sound_speed(
        &self,
        wl: &Prim,
        wr: &Prim,
        h_roe: f64,
        q2_roe: f64,
        sql: f64,
        sqr: f64,
        inv: f64,
    ) -> f64 {
        match self {
            ChemEos::Table(t) => t.roe_sound_speed(wl, wr, h_roe, q2_roe, sql, sqr, inv),
            ChemEos::Blend(b) => b.roe_sound_speed(wl, wr, h_roe, q2_roe, sql, sqr, inv),
        }
    }

    fn stagnation_ghost(
        &self,
        p0: f64,
        rho0: f64,
        c_frac: f64,
        u_n: f64,
        normal: usize,
    ) -> Result<Prim, FlowError> {
        match self {
            ChemEos::Table(t) => t.stagnation_ghost(p0, rho0, c_frac, u_n, normal),
            ChemEos::Blend(b) => b.stagnation_ghost(p0, rho0, c_frac, u_n, normal),
        }
    }

    fn mass_flow_inflow_ghost(
        &self,
        mdot_per_area: f64,
        h_total: f64,
        c_frac: f64,
        p_int: f64,
        sign: f64,
        normal: usize,
    ) -> Result<Prim, FlowError> {
        match self {
            ChemEos::Table(t) => {
                t.mass_flow_inflow_ghost(mdot_per_area, h_total, c_frac, p_int, sign, normal)
            }
            ChemEos::Blend(b) => {
                b.mass_flow_inflow_ghost(mdot_per_area, h_total, c_frac, p_int, sign, normal)
            }
        }
    }
}
