import sys

for path in sys.argv[1:]:
    rows = [l.split() for l in open(path) if l.startswith("t ")]
    t = [float(r[1]) for r in rows]
    f = [float(r[10]) for r in rows]
    w = [float(r[15]) for r in rows]
    final = f[-1]
    for tol in (1e-3, 1e-4):
        s = 0
        for i in range(len(f) - 1, -1, -1):
            if abs(f[i] / final - 1) > tol:
                s = i + 1
                break
        print(f"{path}: final {final} N; within {tol:g} from t {t[s]} ms (wall {w[s]:.0f} s)")
