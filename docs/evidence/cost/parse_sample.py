import re, sys
from collections import defaultdict
# Attribute each sample (a node's self count) to the innermost frame matching a phase rule.
RULES = [
    ("spray", r"crucible::Spray::"),
    ("chemistry: Cantera equilibrium (table fallback)", r"Cantera::|ReactionSource::equilibrate"),
    ("chemistry: Table A lookup", r"EquilibriumTable::|ReactingFlow::react"),
    ("transport properties (mu, lambda, D_k)", r"Flow::transportProperties"),
    ("gradients", r"Flow::transportGradients"),
    ("viscous and diffusive fluxes", r"Flow::transportFluxes"),
    ("SST (eddy viscosity, sources)", r"Flow::eddyViscosity|Flow::turbulence|Turbulence"),
    ("inviscid fluxes (HLLC, reconstruction)", r"Flow::axialFace|Flow::radialFace|hllc|radialProfiles|axialSlopes|Flow::faceThermal|Flow::gather"),
    ("state refresh (T from e, sound speed)", r"Flow::refresh|Flow::primitives"),
    ("time-step limit (stableDt loop)", r"Flow::stableDt"),
    ("thread waiting", r"__psynch_cvwait|__psynch_cvbroad|Pool::loop|_pthread_cond_wait"),
]
line_re = re.compile(r"^([\s+!:|]*)(\d+)\s+(.*)$")
def phase(name):
    for label, pat in RULES:
        if re.search(pat, name): return label
    return None
tree = []  # (depth, count, name)
started = False
for line in open(sys.argv[1]):
    if line.startswith("Call graph:"): started = True; continue
    if started and (line.startswith("Total number in stack") or line.strip().startswith("Sort by top")): break
    if not started: continue
    m = line_re.match(line.rstrip("\n"))
    if not m: continue
    depth = len(m.group(1)); tree.append((depth, int(m.group(2)), m.group(3)))
totals = defaultdict(int)
stack = []  # (depth, count, childsum, category)
def close(entry):
    d, c, cs, cat = entry
    totals[cat or "other"] += c - cs
for depth, count, name in tree:
    while stack and stack[-1][0] >= depth:
        close(stack.pop())
    if stack: stack[-1][2] += count
    parent = stack[-1][3] if stack else None
    cat = phase(name) or parent
    stack.append([depth, count, 0, cat])
while stack: close(stack.pop())
all_ = sum(totals.values()); busy = all_ - totals.get("thread waiting", 0)
print(f"samples {all_}, busy {busy} ({100*busy/all_:.0f}%)")
for k, v in sorted(totals.items(), key=lambda kv: -kv[1]):
    share = "" if k == "thread waiting" else f"{100*v/busy:5.1f}% of busy"
    print(f"{v:7d}  {share:>14}  {k}")
