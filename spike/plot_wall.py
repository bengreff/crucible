# 1 T split runs with dual energy, by magnetic wall condition: total |B| and axial velocity.
import sys,math,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'spike'); from plot_coil import cells,panel
out,bases=sys.argv[1],sys.argv[2:]
def btot(r): return math.hypot(float(r['bz'])+float(r['b0z']),float(r['br'])+float(r['b0r']))
fig,axes=plt.subplots(len(bases),2,figsize=(16,2.4*len(bases)),squeeze=False)
for k,base in enumerate(bases):
    rows=cells(base); name=base.split('/')[-1].replace('.csv','')
    panel(axes[k][0],rows,btot,name+': total |B| [T]','cividis')
    panel(axes[k][1],rows,lambda r:float(r['uz']),name+': u_z [m/s]','RdBu_r')
plt.tight_layout(); plt.savefig(out,dpi=90)
