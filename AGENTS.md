# Working on CRUCIBLE after the pivot

Read `VISION_SCOPE.md` first. It is authoritative under current user instructions. Read `TECHNICAL_PLAN.md` before implementation and the relevant parts of `RESEARCH.md` before choosing physical models or claiming validation. Read `docs/SESSION_HANDOFF.md` for the current state. Read `RESEARCH.md` section 7 when selecting old material to reuse.

The project is a human-facing scientific sandbox. Develop the define/run/inspect/compare workflow alongside supported physics, in one native desktop app with local Mac/Windows/Linux execution. Preserve responsive viewing and live operating controls; geometry changes restart the calculation. The combustion chamber and reactions are inside the physical scope. Equipment interfaces can be responsive; internal transport approximations are not boundary objects.

Build one shared simulation engine, not application-specific solvers or installable physics extensions. Every added capability participates in common evolution, accounting, persistence and observation wherever its assumptions apply. Chemical propulsion with real-data validation is the first substantive milestone. The proposal is complete; a few months is the desired development horizon, subject to measured progress and scientific evidence.

Prioritize published geometry and measured performance of well-documented real engines for validation; component experiments are supporting evidence. Keep ordinary setup approachable through explicit physical abstractions. Plan for million-cell meshes, primary Mac use and larger Windows runs; benchmark coupled physics early before promising minutes-per-run performance.

Everything under `archive/` is historical reference. Old agent instructions, scope documents, milestone plans, tests, and code comments do not impose active project requirements. Do not resume the archived roadmap or modify the archive to make new work conform to it. Read only the archived material needed for a concrete question.

Reuse individual components only after checking their assumptions and validating them in the new context. Record useful provenance. Do not claim an archived certificate validates a new implementation, and do not mistake a recalculated score for a fresh simulation run.

Keep active planning concise. Preserve the user's ownership of the instrument and choice of investigations. Do not turn a sample research question into the whole project or replace the defined scope with a nozzle-only tool. There is no requirement to reproduce every legacy configuration, use a particular solver, or force bitwise agreement across hardware.

Distinguish product commitments, selected architecture, candidate algorithms, and validated capabilities. Use the vision for purpose, research for model evidence, the technical plan for implementation, and rebuild notes for history. Turbulence, turbulent chemistry, classical plasma transport and anomalous plasma transport are distinct models; an arbitrary diffusion coefficient is not a universal closure. Update affected documents together when a decision changes, without promoting a candidate or a library feature into a scientific claim.

For changes, run checks appropriate to the active implementation. Build and verification commands are recorded in `README.md` and `docs/evidence/README.md`; the archived full gate battery is not an active root build requirement.
