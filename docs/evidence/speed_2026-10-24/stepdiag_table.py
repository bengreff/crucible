import sys,re
for path in sys.argv[1:]:
    print(path)
    print(' t ms | steps | dt us | wall s/step | bound by | shares flow/ac/diff | dominant f/a/d | nu kind(eddy th) | x no-ac | x no-diff | x flow-only | LTS full | LTS no-ac | L3 cells | L3 bind cell')
    lines=open(path,errors='replace').read().splitlines()
    for k,l in enumerate(lines):
        m=re.match(r'STEPDIAG t ([\d.]+) ms: (\d+) base steps, mean dt ([\d.e+-]+) s, ([\d.]+) s wall per step \| bound by level 0..3 (\d+) (\d+) (\d+) (\d+), event (\d+), spray bound (\d+), sample/end (\d+) \| at the binding cell, rate shares flow ([\d.]+) acoustic ([\d.]+) diffusive ([\d.]+) \(dominant: flow (\d+), acoustic (\d+), diffusive (\d+)\), acoustic/flow ([\d.]+); largest nu kind \[[^\]]*\] ([\d ]+) \| speed-up of the global step: no acoustics ([\d.]+), no diffusion ([\d.]+), flow only ([\d.]+) \| work now / LTS: full ([\d.]+), no acoustics ([\d.]+), no diffusion ([\d.]+), flow only ([\d.]+)',l)
        if not m: continue
        g=m.groups(); nu=g[18].split()
        l3=[x for x in lines[k+1:k+5] if 'level 3' in x]
        c=re.search(r'(\d+) active',l3[0]).group(1) if l3 else '?'; bc=re.search(r'binding cell (\([^)]*\))',l3[0]).group(1) if l3 else '?'
        print(f' {g[0]} | {g[1]} | {float(g[2])*1e6:.3f} | {g[3]} | L0-3 {g[4]} {g[5]} {g[6]} {g[7]} ev {g[8]} spray {g[9]} smp {g[10]} | {g[11]}/{g[12]}/{g[13]} | {g[14]}/{g[15]}/{g[16]} | {nu[5]} | {g[19]} | {g[20]} | {g[21]} | {g[22]} | {g[23]} | {c} | {bc}')
