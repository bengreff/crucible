//! HLLC approximate Riemann flux with **Batten wavespeed estimates**
//! (SOLV-1 §3.2; META-3 `hllc` — Toro ch. 10, Batten et al. 1997). Batten's
//! Roe-average bounds restore the contact and species-advection waves that
//! HLLE smears — mandatory with the composition block in `U` (SOLV-1 §3.2).
//!
//! The flux is written against the [`EosLaw`] seam (FND-7 constitutive
//! spine): the Riemann solve reads only (ρ, p, c, ρE) closures through the
//! trait, never a material label — the general-convex-EOS occupant supplies
//! its sound speeds (local and Roe-averaged, from the aux Γ₁ slot) without
//! touching the wave algebra (SOLV-1 §3.2, Castro/PeleC treatment).

use super::{EosLaw, I_EN, I_RB, I_RC, I_RHO, NCOMP, Prim};

/// Physical (hyperbolic) flux of the conserved vector in direction `n`
/// (component slot of the normal velocity: 1 = r, 2 = θ, 3 = z), from a
/// primitive state. Components ordered as `U` (SOLV-1 §3.1).
#[inline]
pub fn physical_flux<E: EosLaw>(w: &Prim, n: usize, eos: &E) -> [f64; NCOMP] {
    let (rho, p, c_frac) = (w[I_RHO], w[4], w[I_RC]);
    let u_n = w[n];
    let m = rho * u_n;
    let e_tot = eos.total_energy(w);
    let mut f = [0.0; NCOMP];
    f[I_RHO] = m;
    f[1] = m * w[1];
    f[2] = m * w[2];
    f[3] = m * w[3];
    f[n] += p;
    f[I_EN] = u_n * (e_tot + p);
    f[I_RC] = m * c_frac;
    // Burn progress advects with the contact, exactly like the composition
    // scalar — flux = mass flux × the upwind value (SOLV-4 §3.6; the doc's c).
    f[I_RB] = m * w[I_RB];
    f
}

/// HLLC-Batten flux across a face with normal in slot `n`, given the
/// reconstructed left/right primitive states. Deterministic: the wave-branch
/// selection is a fixed comparison ladder, no iteration.
pub fn hllc_flux<E: EosLaw>(wl: &Prim, wr: &Prim, n: usize, eos: &E) -> [f64; NCOMP] {
    let (rho_l, p_l) = (wl[I_RHO], wl[4]);
    let (rho_r, p_r) = (wr[I_RHO], wr[4]);
    let (un_l, un_r) = (wl[n], wr[n]);
    let c_l = eos.sound_speed_w(wl);
    let c_r = eos.sound_speed_w(wr);

    // Batten wavespeeds: bound each family by both the one-sided and the
    // Roe-average estimate (Batten et al., SIAM JSC 18(6), 1997).
    let sql = rho_l.sqrt();
    let sqr = rho_r.sqrt();
    let inv = 1.0 / (sql + sqr);
    let u_roe = [
        0.0, // slot 0 unused (density has no velocity)
        (sql * wl[1] + sqr * wr[1]) * inv,
        (sql * wl[2] + sqr * wr[2]) * inv,
        (sql * wl[3] + sqr * wr[3]) * inv,
    ];
    let h_l = (eos.total_energy(wl) + p_l) / rho_l;
    let h_r = (eos.total_energy(wr) + p_r) / rho_r;
    let h_roe = (sql * h_l + sqr * h_r) * inv;
    let q2_roe = u_roe[1] * u_roe[1] + u_roe[2] * u_roe[2] + u_roe[3] * u_roe[3];
    let c_roe = eos.roe_sound_speed(wl, wr, h_roe, q2_roe, sql, sqr, inv);

    let s_l = (un_l - c_l).min(u_roe[n] - c_roe);
    let s_r = (un_r + c_r).max(u_roe[n] + c_roe);

    // Contact speed (the HLLC middle wave).
    let ml = rho_l * (s_l - un_l);
    let mr = rho_r * (s_r - un_r);
    let s_m = (mr * un_r - ml * un_l + p_l - p_r) / (mr - ml);

    if s_l >= 0.0 {
        return physical_flux(wl, n, eos);
    }
    if s_r <= 0.0 {
        return physical_flux(wr, n, eos);
    }
    // Star-region flux: F_K + S_K (U*_K − U_K), K chosen by the contact.
    let (w, s_k) = if s_m >= 0.0 { (wl, s_l) } else { (wr, s_r) };
    let (rho, p, un) = (w[I_RHO], w[4], w[n]);
    let e_tot = eos.total_energy(w);
    let p_star = rho * (un - s_k) * (un - s_m) + p;

    let mut u_k = [0.0; NCOMP];
    u_k[I_RHO] = rho;
    u_k[1] = rho * w[1];
    u_k[2] = rho * w[2];
    u_k[3] = rho * w[3];
    u_k[I_EN] = e_tot;
    u_k[I_RC] = rho * w[I_RC];
    u_k[I_RB] = rho * w[I_RB];

    // Star energy factored as `fac·E + (p*S_M − p·u_n)/(S_K − S_M)` — same
    // algebra as Toro eq. 10.73, but with `fac = (S_K−u_n)/(S_K−S_M)` shared
    // with ρ*, a state at rest has fac = S_K/S_K = 1 exactly and the star
    // state IS the input state bitwise: a uniform gas at rest is a fixed
    // point of the flux to round-off-free exactness (the §3.3 well-balance
    // certificate item leans on this).
    let fac = (s_k - un) / (s_k - s_m);
    let rho_star = rho * fac;
    let mut u_star = [0.0; NCOMP];
    u_star[I_RHO] = rho_star;
    u_star[1] = rho_star * w[1];
    u_star[2] = rho_star * w[2];
    u_star[3] = rho_star * w[3];
    u_star[n] = rho_star * s_m;
    u_star[I_EN] = fac * e_tot + (p_star * s_m - p * un) / (s_k - s_m);
    u_star[I_RC] = rho_star * w[I_RC];
    u_star[I_RB] = rho_star * w[I_RB];

    let f_k = physical_flux(w, n, eos);
    let mut f = [0.0; NCOMP];
    for i in 0..NCOMP {
        f[i] = f_k[i] + s_k * (u_star[i] - u_k[i]);
    }
    f
}

// Covers SOLV-1 §3.2 (HLLC-Batten flux properties) — named per META-2 §4.
#[cfg(test)]
mod tests {
    use super::*;
    use crate::euler::{GammaLaw, prim6};

    const EOS: GammaLaw = GammaLaw { gamma: 1.4 };

    #[test]
    fn solv1_s32_identical_states_give_the_physical_flux() {
        let w: Prim = prim6(1.2, 0.3, -0.1, 0.7, 2.5, 0.4);
        for n in 1..=3 {
            let f = hllc_flux(&w, &w, n, &EOS);
            let fp = physical_flux(&w, n, &EOS);
            for i in 0..NCOMP {
                assert!(
                    (f[i] - fp[i]).abs() <= 1e-14 * fp[i].abs().max(1.0),
                    "component {i}, normal {n}: {} vs {}",
                    f[i],
                    fp[i]
                );
            }
        }
    }

    #[test]
    fn solv1_s32_supersonic_states_upwind_exactly() {
        // Both states moving right far above their sound speed: the flux is
        // exactly the left physical flux (S_L > 0).
        let wl: Prim = prim6(1.0, 0.0, 0.0, 5.0, 1.0, 1.0);
        let wr: Prim = prim6(0.5, 0.0, 0.0, 5.0, 0.5, 0.0);
        let f = hllc_flux(&wl, &wr, 3, &EOS);
        let fl = physical_flux(&wl, 3, &EOS);
        assert_eq!(f, fl);
    }

    #[test]
    fn solv1_s32_stationary_contact_is_preserved_exactly() {
        // A stationary contact (equal p, zero normal u, different ρ and C):
        // HLLC must return zero mass/species flux — the property Batten
        // wavespeeds exist to protect (HLLE smears this).
        let wl: Prim = prim6(1.0, 0.0, 0.0, 0.0, 1.0, 1.0);
        let wr: Prim = prim6(0.125, 0.0, 0.0, 0.0, 1.0, 0.0);
        let f = hllc_flux(&wl, &wr, 3, &EOS);
        assert_eq!(f[I_RHO], 0.0, "mass flux through a stationary contact");
        assert_eq!(f[I_RC], 0.0, "species flux through a stationary contact");
        assert_eq!(f[I_EN], 0.0, "energy flux through a stationary contact");
        // Momentum flux is the (equal) pressure.
        assert!((f[3] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn solv1_s32_species_ride_the_upwind_side_of_the_contact() {
        // A right-moving contact: species flux = mass flux × left C.
        let wl: Prim = prim6(1.0, 0.0, 0.0, 0.5, 1.0, 1.0);
        let wr: Prim = prim6(0.25, 0.0, 0.0, 0.5, 1.0, 0.0);
        let f = hllc_flux(&wl, &wr, 3, &EOS);
        assert!(f[I_RHO] > 0.0);
        assert!(
            (f[I_RC] - f[I_RHO]).abs() < 1e-13,
            "C = 1 upwind ⇒ species flux equals mass flux"
        );
    }
}
