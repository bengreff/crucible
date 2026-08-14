//! CRUCIBLE CLI. Scaffold only — run orchestration (COUP-3/4) arrives with
//! the solvers. For now it states what it is and pins its provenance.

fn main() {
    println!(
        "CRUCIBLE {} — scaffold; no engine yet. Constants: {}",
        env!("CARGO_PKG_VERSION"),
        crucible_constants::SOURCE,
    );
}
