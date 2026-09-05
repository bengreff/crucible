//! crucible-gpu — the CUDA residency tier (PLAN Phase 4). The library half
//! holds the host-side ENGINE bridge shared by the cross-check bin and the
//! overnight harness: config → assembly → tables → device geometry/SRD/table
//! marshaling → the persistent device engine handle (`cuda/residency_engine.cu`),
//! plus the run schedule pieces replicated from `crucible_engine::run` (the
//! injector ramp, the pump-down ambient, the igniter kernel, the audit
//! reference scales). One owner for that replication.
pub mod engine_host;
