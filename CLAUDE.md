# CRUCIBLE: current instructions

Start every session with `docs/SESSION_HANDOFF.md` (current state, what's running, exact next steps).

Read `AGENTS.md` and `VISION_SCOPE.md` for what the project is and its scope, following the 17 September 2026 pivot. Read `TECHNICAL_PLAN.md` before implementation. `RESEARCH.md` holds algorithm candidates, evidence and unresolved physical models (not validated capabilities) and, in section 7, lessons from the pre-pivot attempt.

`docs/evidence/README.md` indexes verification evidence; `docs/validation/` holds real-engine cases; `docs/parked/` is parked work (the MHD spike). Everything under `archive/` is historical and does not govern new work.

Direction since 5 October 2026: the lightweight engine (TECHNICAL_PLAN, *Lightweight engine*): offline chemistry tables, wall functions, adaptive resolution, a cold start in minutes. Roadmap dates are ceilings, not a pace.

Standing rules: one project, one line. Commit directly to `main`, no PRs and no merge questions. When the Director runs two workers in parallel, the second uses a worktree and rebases onto `main` and fast-forwards it itself at every milestone; that is not a branch to ask about. Validation predictions are pre-registered blind: commit the prediction before digitising a case's measured curves, and never let a research agent read out or report measured outcome values for a case still awaiting its prediction. The backhouse engine lives at `/home/greff/crucible`; never write to `/home/greff/inquiry-project` (the obsolete pre-pivot Rust project).
