#!/bin/sh
# AMR.md, *Fixing the tag thresholds, then criterion 4*: the rule for s*, then criterion 4 at s*, from the "c4" lines
# crucible_amr_c4 prints. Usage: sh tests/amr_c4_judge.sh LOG... (the AMR runs at s = 1, 1/2, 1/4, 1/8 and the
# uniform run, in any order).
grep -h '^c4 ' "$@" | awk '
function abs(x) { return x < 0 ? -x : x }
$2 == "uniform" { uc = $6; ui = $8; uf = $10; ul = $12; un = $18; have = 1; next }
{ n++; s[n] = $4; c[n] = $6; p[n] = $8; f[n] = $10; l[n] = $12; b[n] = $14; m[n] = $18 }
END {
  for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (s[j] > s[i]) {
    t = s[i]; s[i] = s[j]; s[j] = t; t = c[i]; c[i] = c[j]; c[j] = t; t = p[i]; p[i] = p[j]; p[j] = t
    t = f[i]; f[i] = f[j]; f[j] = t; t = l[i]; l[i] = l[j]; l[j] = t; t = b[i]; b[i] = b[j]; b[j] = t; t = m[i]; m[i] = m[j]; m[j] = t
  }
  bad = 0
  for (i = 1; i <= n; i++) {
    printf("s %-6g settled c* %.4f m/s, Isp_vac %.5f s, F_vac %.4f N; light-off %.6f ms, base step %.3e s; %d samples\n", s[i], c[i], p[i], f[i], l[i] * 1e3, b[i], m[i])
    if (m[i] != 21 || l[i] < 0) bad = 1
  }
  if (n != 4 || s[1] != 1 || s[2] != 0.5 || s[3] != 0.25 || s[4] != 0.125) bad = 1
  if (bad) { print "the runs are not the registered set (s = 1, 1/2, 1/4, 1/8, 21 samples, a light-off each): nothing judged"; exit 1 }
  star = 0
  for (i = 1; i < n; i++) {
    dc = abs(c[i] / c[i + 1] - 1); dp = abs(p[i] / p[i + 1] - 1); df = abs(f[i] / f[i + 1] - 1)
    dl = abs(l[i] - l[i + 1]); bl = b[i + 1] + 0.002 * l[i + 1]
    ok = dc <= 2e-4 && dp <= 2e-4 && df <= 2e-4 && dl <= bl
    printf("s %g against s %g: c* %.3e, Isp %.3e, F %.3e (limit 2e-4 each); light-off %.3e s (limit %.3e s): %s\n", s[i], s[i + 1], dc, dp, df, dl, bl, ok ? "meets the rule" : "does not")
    if (ok && star == 0) star = i
  }
  if (star == 0) { print "no s qualifies: the thresholds are not fixed (the Director decides); criterion 4 is not judged"; exit 1 }
  printf("s* = %g: the thresholds are fixed at f_H %g, T %g, pressure %g, vorticity %g\n", s[star], 0.02 * s[star], 0.05 * s[star], 0.10 * s[star], 0.05 * s[star])
  if (!have) { print "no uniform run: criterion 4 is not judged"; exit 1 }
  printf("uniform 64 x 12: settled c* %.4f m/s, Isp_vac %.5f s, F_vac %.4f N; light-off %.6f ms; %d samples\n", uc, ui, uf, ul * 1e3, un)
  dc = abs(c[star] / uc - 1); dp = abs(p[star] / ui - 1); df = abs(f[star] / uf - 1); dl = abs(l[star] - ul); bl = b[star] + 0.01 * ul
  fails = (dc > 1e-3) + (dp > 1e-3) + (df > 1e-3) + (dl > bl)
  printf("criterion 4 at s*: c* %+.3e (limit 1e-3) %s\n", c[star] / uc - 1, dc <= 1e-3 ? "met" : "not met")
  printf("criterion 4 at s*: Isp_vac %+.3e (limit 1e-3) %s\n", p[star] / ui - 1, dp <= 1e-3 ? "met" : "not met")
  printf("criterion 4 at s*: F_vac %+.3e (limit 1e-3) %s\n", f[star] / uf - 1, df <= 1e-3 ? "met" : "not met")
  printf("criterion 4 at s*: light-off %+.3e s (limit one base step plus 1%%, %.3e s) %s\n", l[star] - ul, bl, dl <= bl ? "met" : "not met")
  printf("criterion 4: %s\n", fails ? "not met" : "met")
}'
