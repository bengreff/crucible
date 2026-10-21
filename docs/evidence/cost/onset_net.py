# inquiry2: net b items by class (turbulent, molecular, inviscid, spray, other) for each cell block in an onset readout file.
import re, sys
for p in sys.argv[1:]:
    blk = None; acc = {}
    for line in open(p):
        if line.startswith('## cells'):
            blk = line.strip()[3:]; continue
        m = re.match(r'\s+(above|below|low z|high z|cell)\s+(.+?)\s+([+-]\d\.\d+e[+-]\d+)\s*$', line)
        if not m or blk is None: continue
        name = m.group(2); v = float(m.group(3))
        cls = 'turbulent' if 'turbulent' in name else 'molecular' if 'molecular' in name else 'inviscid' if name == 'inviscid' else 'spray' if 'spray' in name else 'other'
        acc.setdefault(blk, {}).setdefault(cls, 0.0)
        acc[blk][cls] += v
    for b, d in acc.items():
        print(p.split('/')[-1], b, ', '.join(f'{k} {v:+.4e}' for k, v in d.items()), f'sum {sum(d.values()):+.4e}')
