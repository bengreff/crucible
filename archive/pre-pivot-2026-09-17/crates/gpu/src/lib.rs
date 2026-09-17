//! crucible-gpu — the CUDA residency tier (PLAN Phase 4). The library half
//! holds the host-side ENGINE bridge shared by the cross-check bin and the
//! overnight harness: config → assembly → tables → device geometry/SRD/table
//! marshaling → the persistent device engine handle (`cuda/residency_engine.cu`),
//! plus the run schedule pieces replicated from `crucible_engine::run` (the
//! injector ramp, the pump-down ambient, the igniter kernel, the audit
//! reference scales). One owner for that replication.
pub mod engine_host;

// Every per-leg kernel entry point in the CUDA archive, referenced from a
// function every bin calls (`keep_kernels`): a static archive linked through
// an rlib surrenders only the objects the LINKED rlib code references, and a
// bin that uses none of the lib pulls nothing — so each cross-check bin calls
// this once and its own FFI prototypes then resolve. Signatures deliberately
// elided — never called through this table. (The engine entry points are
// referenced by `engine_host` with their real signatures.)
mod keep {
    unsafe extern "C" {
        fn gpu_blend_project();
        fn gpu_class_a_bench();
        fn gpu_class_a_bench_3d();
        fn gpu_class_a_march();
        fn gpu_class_a_march_3d();
        fn gpu_class_a_rhs();
        fn gpu_class_a_rhs_3d();
        fn gpu_class_d_assemble();
        fn gpu_class_d_cg();
        fn gpu_class_d_iterate();
        fn gpu_class_r_update();
        fn gpu_combustion_source();
        fn gpu_hllc();
        fn gpu_srd_3d();
        fn gpu_stable_dt();
        fn gpu_stable_dt_3d();
        fn gpu_table_project();
    }
    static KEEP_KERNELS: [unsafe extern "C" fn(); 17] = [
        gpu_blend_project,
        gpu_class_a_bench,
        gpu_class_a_bench_3d,
        gpu_class_a_march,
        gpu_class_a_march_3d,
        gpu_class_a_rhs,
        gpu_class_a_rhs_3d,
        gpu_class_d_assemble,
        gpu_class_d_cg,
        gpu_class_d_iterate,
        gpu_class_r_update,
        gpu_combustion_source,
        gpu_hllc,
        gpu_srd_3d,
        gpu_stable_dt,
        gpu_stable_dt_3d,
        gpu_table_project,
    ];
    pub fn count() -> usize {
        KEEP_KERNELS.iter().filter(|f| (**f as usize) != 0).count()
    }
}

/// Call once from any bin that reaches the per-leg kernels through its own
/// FFI prototypes (see `keep`). Returns the number of retained entry points.
pub fn keep_kernels() -> usize {
    keep::count()
}
