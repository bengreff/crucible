import csv,sys,math
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
def cells(p): return list(csv.DictReader(open(p)))
def nodes(p):
    rows=list(csv.DictReader(open(p))); ni=max(int(r['i']) for r in rows)+1; nj=max(int(r['j']) for r in rows)+1
    Z=np.zeros((ni,nj));R=np.zeros((ni,nj));P=np.zeros((ni,nj))
    for r in rows: i,j=int(r['i']),int(r['j']);Z[i,j]=float(r['z']);R[i,j]=float(r['r']);P[i,j]=float(r['psi'])
    return Z,R,P
def panel(ax,rows,f,title,cmap,log=False):
    polys=[];vals=[]
    for r in rows:
        z0,z1=float(r['z0']),float(r['z1']);rb0,rb1,rt0,rt1=[float(r[k]) for k in ['rb0','rb1','rt0','rt1']]
        polys.append([(z0,rb0),(z1,rb1),(z1,rt1),(z0,rt0)]);vals.append(f(r))
    v=np.array(vals); 
    if log: v=np.log10(np.maximum(v,1e-6))
    pc=PolyCollection(polys,array=v,cmap=cmap,edgecolors='none')
    ax.add_collection(pc);ax.set_xlim(0,0.6);ax.set_ylim(0,0.036);ax.set_aspect('equal');ax.set_title(title,fontsize=9)
    plt.colorbar(pc,ax=ax,fraction=0.02,pad=0.01); ax.plot([0.216],[0.035],'rx')
if __name__!="__main__": pass
else:
  out=sys.argv[1];base=sys.argv[2];ref=sys.argv[3]
  rows=cells(base);Z,R,P=nodes(base+'.psi');Z0,R0,P0=nodes(ref+'.psi') if ref!='-' else (None,None,None)
  fig,axes=plt.subplots(4,1,figsize=(12,9))
  panel(axes[0],rows,lambda r:float(r['p']),base.split('/')[-1]+': pressure [Pa]','viridis')
  panel(axes[1],rows,lambda r:float(r['mach']),'Mach','magma')
  panel(axes[2],rows,lambda r:math.hypot(float(r['bz']),float(r['br'])),'log10 |B_poloidal| [T] now','cividis',True)
  panel(axes[3],rows,lambda r:math.hypot(float(r['b0z']),float(r['b0r'])),'log10 |B| applied vacuum field [T]; lines: psi applied (grey), psi now (red)','cividis',True)
  lv=np.linspace(P0.min(),P0.max(),14)[1:-1] if P0 is not None else 12
  if P0 is not None: axes[3].contour(Z0,R0,P0,levels=lv,colors='0.7',linewidths=0.8)
  for ax in axes[2:]: ax.contour(Z,R,P,levels=lv,colors='r',linewidths=0.8)
  plt.tight_layout();plt.savefig(out,dpi=100)
