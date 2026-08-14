//! Will implement **FND-4 v0.2** — config schema, loader pipeline
//! (parse → migrate → structural deserialize → registry dispatch → wiring
//! validation → default resolution) and the run-manifest emission.
//!
//! Empty by design: the scaffold session claims the workspace layout
//! (META-2 §4, crate-per-doc-area); the loader lands test-first against
//! FND-4 §6 in the next session. No code here may precede the doc.
