//! SOLV-1 §6-2 pass criteria: the full Euler operator recovers a
//! manufactured solution at formal order — every flux direction and all
//! three geometric sources exercised with nonzero operands (which the Sod
//! run structurally cannot do: its radial velocity is identically zero).
//! Same shared-constants contract as the other certificates.

use crucible_solvers::euler::{EULER_FIELDS, NCOMP};
use crucible_solvers::euler_mms::{MMS_ORDER_MAX, MMS_ORDER_MIN, mms_axisym_swirl, mms_theta_mode};

#[test]
fn solv1_s62_mms_axisym_swirl_all_components_at_formal_order() {
    let study = mms_axisym_swirl();
    for (k, name) in EULER_FIELDS.iter().enumerate() {
        for (i, p) in study.observed_orders(k).iter().enumerate() {
            assert!(
                (MMS_ORDER_MIN..=MMS_ORDER_MAX).contains(p),
                "{}: component {k} ({name}) order {p} at refinement {i} outside \
                 [{MMS_ORDER_MIN}, {MMS_ORDER_MAX}]: {study:?}",
                study.label,
            );
        }
        // Errors must actually shrink ~4× per level, not just have noisy
        // ratios near a floor.
        for w in study.levels.windows(2) {
            assert!(w[0].l1[k] > 3.0 * w[1].l1[k]);
        }
    }
}

#[test]
fn solv1_s62_mms_theta_mode_all_components_at_formal_order() {
    let study = mms_theta_mode();
    for (k, name) in EULER_FIELDS.iter().enumerate() {
        for (i, p) in study.observed_orders(k).iter().enumerate() {
            assert!(
                (MMS_ORDER_MIN..=MMS_ORDER_MAX).contains(p),
                "{}: component {k} ({name}) order {p} at refinement {i} outside band: {study:?}",
                study.label,
            );
        }
        for w in study.levels.windows(2) {
            assert!(w[0].l1[k] > 3.0 * w[1].l1[k]);
        }
    }
}

#[test]
fn solv1_s62_mms_studies_are_reproducible() {
    let a = mms_axisym_swirl();
    let b = mms_axisym_swirl();
    for (x, y) in a.levels.iter().zip(&b.levels) {
        for k in 0..NCOMP {
            assert_eq!(x.l1[k].to_bits(), y.l1[k].to_bits());
        }
    }
}
