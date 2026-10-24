import re,sys,collections
# Parse a macOS `sample` call graph; return per thread a list of nodes (depth, count, symbol), in order.
def parse(path):
    threads=[];cur=None
    for line in open(path,errors='replace'):
        if line.startswith('Total number in stack'): break
        m=re.match(r'^    (\d+) (Thread_\S+.*)$',line)
        if m: cur=[];threads.append((m.group(2),int(m.group(1)),cur));continue
        if cur is None: continue
        m=re.match(r'^(    [ +!:|]*?)(\d+) (.*)$',line)
        if not m: continue
        depth=(len(m.group(1))-4)//2
        sym=m.group(3)
        sym=re.sub(r'\s+\(in [^)]*\).*$','',sym)
        cur.append((depth,int(m.group(2)),sym))
    return threads
def tree(nodes):
    # children lists
    out=[];stack=[]
    for idx,(d,c,s) in enumerate(nodes):
        while stack and nodes[stack[-1]][0]>=d: stack.pop()
        parent=stack[-1] if stack else None
        out.append([d,c,s,parent,[]])
        if parent is not None: out[parent][4].append(idx)
        stack.append(idx)
    return out
if __name__=='__main__':
    th=parse(sys.argv[1]);mn=int(sys.argv[2]) if len(sys.argv)>2 else 100;maxd=int(sys.argv[3]) if len(sys.argv)>3 else 12
    name,total,nodes=th[0]
    print(name,total)
    for d,c,s in nodes:
        if c>=mn and d<=maxd: print('  '*d+f'{c} {s[:150]}')
