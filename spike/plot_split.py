# Split-field failure state: log pressure, total |B| = |b0 + b1|, log beta, radial velocity.
import sys,math,numpy as np,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'spike'); from plot_coil import cells,panel
out,base=sys.argv[1],sys.argv[2]; rows=cells(base); mu0=4e-7*math.pi
def btot(r): return math.hypot(float(r['bz'])+float(r['b0z']),float(r['br'])+float(r['b0r']))
fig,axes=plt.subplots(4,1,figsize=(12,9))
panel(axes[0],rows,lambda r:float(r['p']),base.split('/')[-1]+': log10 pressure [Pa]','viridis',True)
panel(axes[1],rows,btot,'total |B| [T]','cividis')
panel(axes[2],rows,lambda r:2*mu0*float(r['p'])/max(btot(r)**2,1e-30),'log10 beta','coolwarm',True)
panel(axes[3],rows,lambda r:float(r['ur']),'u_r [m/s]','RdBu_r')
plt.tight_layout(); plt.savefig(out,dpi=100)
