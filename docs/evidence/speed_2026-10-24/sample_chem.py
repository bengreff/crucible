import sys,re,collections
sys.path.insert(0,'/tmp'); from sample_tree import parse,tree
SUB=[('Cantera rates (kinetics/thermo update)',r'Cantera::|getNetProductionRates|updateROP|_update_rates|updateThermo|GasKinetics|IdealGas'),
     ('CVODES internals/linear algebra',r'cv[A-Z]|CVode|N_V|SUNLin|SUNMat|dgetr|Dense|dense'),
     ('ReactionStep rhs/jacobian (own)',r'ReactionStep::(rhs|jac|Rhs|Jac)'),
     ('ReactionStep setup/advance (own)',r'ReactionStep'),
     ('table lookup',r'Manifold|Table|manifold|table'),
     ('react own',r'ReactingFlow::react')]
th=parse(sys.argv[1]);name,total,nodes=th[0];t=tree(nodes)
cat=collections.Counter();leaf=collections.Counter()
def walk(i,inchem,s):
    d,cnt,sym,parent,ch=t[i]
    if re.search(r'ReactingFlow::react\b',sym): inchem=True
    if inchem:
        for n,p in SUB:
            if re.search(p,sym): s=n;break
    self_=cnt-sum(t[j][1] for j in ch)
    if inchem:
        cat[s or 'other']+=self_
        leaf[re.sub(r'\(.*','',sym)[:90]]+=self_
    for j in ch: walk(j,inchem,s)
for i,n in enumerate(t):
    if n[3] is None: walk(i,False,None)
c=sum(cat.values())
print(f'chemistry {c} of {total} samples ({100*c/total:.1f}%)')
for k,v in cat.most_common(): print(f'  {100*v/c:5.1f}% of chemistry  {k}')
print('top self frames in chemistry:')
for k,v in leaf.most_common(14): print(f'  {100*v/c:5.1f}%  {k}')
