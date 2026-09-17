# CRUCIBLE

A human-facing technical sandbox for investigating how reacting gas and plasma, geometry, and fields interact to produce propulsion.

**Status: scope reset on 17 September 2026.** The previous implementation is archived. The new sandbox is not implemented yet.

- [Vision and scope](VISION_SCOPE.md) is the authoritative definition of the project, including the visual human workflow.
- [Rebuild notes](REBUILD_NOTES.md) preserve lessons, useful decisions, and specific archived code references.
- [Archive](archive/README.md) explains what was preserved, how to verify it, and what is local-only.

The intended workflow is **construct → operate → observe → change → compare**. The combustion chamber's reacting contents belong inside the model; equipment is represented through physical interfaces. Visualization and usability develop alongside the physical model.

There is deliberately no active solver workspace or inherited CI gate at the root. Build instructions will accompany the new implementation. Historical commands and certificates are reference material, not evidence of validation of the new project.
