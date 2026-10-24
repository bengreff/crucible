import sys,re,collections
sys.path.insert(0,'/tmp'); from sample_tree import parse,tree
CATS=[('chemistry (ReactingFlow::react: CVODES, Cantera kinetics, table)',r'ReactingFlow::react\b|ReactionStep|CVode'),
      ('stable step (Flow::stableDt incl. its transport, gradients, eddy viscosity)',r'Flow::stableDt'),
      ('flow stages (Flow::firstStage/secondStage)',r'Flow::firstStage|Flow::secondStage'),
      ('flow step start (Flow::beginStep)',r'Flow::beginStep'),
      ('AMR fill',r'Amr::fill'),('AMR reflux',r'Amr::reflux'),('AMR regrid/tags',r'Amr::regrid|Amr::tag'),
      ('AMR average down',r'Amr::average|averageDown'),('AMR save/restore',r'Amr::save|Amr::restore'),
      ('spray',r'Spray::'),('checkpoint',r'[Cc]heckpoint'),
      ('driver sample (c* by Cantera equilibrate, budgets, output)',r'runCase\(.*\)::\$_'),]
SUB=[('transport properties',r'transportProperties|Medium::transport|Transport::|fits'),('gradients/eddy viscosity',r'transportGradients|eddyViscosity'),
     ('turbulence sources',r'turbulenceSource|sstSource|Sst|sourceCoefficients'),('fluxes/reconstruction',r'[Ff]lux|reconstruct|riemann|Riemann|hllc|Hllc|slope'),
     ('refresh/primitives',r'Flow::refresh|primitive'),('prescribe/boundaries',r'prescribe|boundary|Boundary|inlet|outlet|wall'),('admissibility',r'admissib|Admissib')]
def run(path,thread=0):
    th=parse(path);name,total,nodes=th[thread];t=tree(nodes)
    cat=collections.Counter();sub=collections.Counter();wait=collections.Counter()
    def walk(i,c,s,inrun):
        d,cnt,sym,parent,ch=t[i]
        for k,(n,p) in enumerate(CATS):
            if re.search(p,sym): c=n;s=None
        if c and c.startswith(('stable','flow stages')):
            for n,p in SUB:
                if re.search(p,sym): s=n
        if 'Pool::run' in sym: inrun=True
        if 'Pool::execute' in sym: inrun=False
        self_=cnt-sum(t[j][1] for j in ch)
        key=c or 'other'
        cat[key]+=self_
        if inrun: wait[key]+=self_
        if s or (c and c.startswith(('stable','flow stages'))): sub[(key,s or 'other')]+=self_
        for j in ch: walk(j,c,s,inrun)
    for i,n in enumerate(t):
        if n[3] is None: walk(i,None,None,False)
    print(f'{path}: {name.split()[0]}, {total} samples')
    for k,v in cat.most_common(): print(f'  {100*v/total:5.1f}%  {k}   (of which waiting in Pool::run for other workers {100*wait[k]/total:.1f}%)')
    print('  inside the stable step and the stages:')
    for (k,s),v in sorted(sub.items(),key=lambda x:-x[1]):
        if v/total>0.003: print(f'    {100*v/total:5.1f}%  {k.split(" (")[0]} / {s}')
    print(f'  waiting in Pool::run in total {100*sum(wait.values())/total:.1f}%')
if __name__=='__main__': run(sys.argv[1],int(sys.argv[2]) if len(sys.argv)>2 else 0)
