//! Will implement **FND-2 v0.6** — the natively cylindrical sparse-brick
//! world-state grid on (i_r, i_θ, i_z), adaptive N_θ with guard = 4, SoA
//! layout, CPU fixed-order reference execution (the GPU relaxed-reduction
//! path comes later, validated against this — D-H).
//!
//! Empty by design — lands test-first against FND-2 §6.
