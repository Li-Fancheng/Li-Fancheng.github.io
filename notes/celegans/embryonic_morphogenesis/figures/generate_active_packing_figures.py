"""Analytical and schematic figures for active matter, packing and topology.
No fitted embryo data. Equations, network construction and parameters are in note 06.
Requires NumPy and Matplotlib.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle,Polygon
OUT=Path(__file__).resolve().parent
C=['#287f83','#b2576f','#b88929','#6677aa']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.titlesize':14,
 'axes.labelsize':12,'legend.fontsize':10,'axes.spines.top':False,
 'axes.spines.right':False,'svg.fonttype':'path','savefig.facecolor':'white',
 'figure.facecolor':'white','text.color':'#243746','axes.labelcolor':'#243746'})
def save(fig,name):
 fig.savefig(OUT/(name+'.svg'),bbox_inches='tight',pad_inches=.18)
 fig.savefig(OUT/(name+'.png'),dpi=160,bbox_inches='tight',pad_inches=.18)
 plt.close(fig)

def rigidity():
 pos=np.array([[0,0],[1,0],[1,1],[0,1]],float)
 sets=[[(0,1),(1,2),(2,3),(3,0)],[(0,1),(1,2),(2,3),(3,0),(0,2)],
       [(0,1),(1,2),(2,3),(3,0),(0,2),(1,3)]]
 fig,ax=plt.subplots(1,3,figsize=(12,4.25),layout='constrained')
 result=[]
 for a,edges,title in zip(ax,sets,['A   Square','B   One diagonal','C   Two diagonals']):
  M=np.zeros((len(edges),8))
  for k,(i,j) in enumerate(edges):
   e=pos[j]-pos[i];e/=np.linalg.norm(e)
   M[k,2*i:2*i+2]=-e;M[k,2*j:2*j+2]=e
   a.plot(pos[[i,j],0],pos[[i,j],1],lw=3,c=C[0] if k<4 else C[1],zorder=1)
  rank=int(np.linalg.matrix_rank(M,tol=1e-10))
  N0=8-rank;Nss=len(edges)-rank
  a.scatter(pos[:,0],pos[:,1],s=150,c='#243746',zorder=2)
  for k,p in enumerate(pos):a.text(p[0],p[1]-.14,chr(65+k),ha='center',va='top')
  a.set(title=title,aspect='equal',xlim=(-.25,1.25),ylim=(-.45,1.25));a.axis('off')
  a.text(.5,-.35,rf'$N_b={len(edges)},\ N_0={N0},\ N_{{ss}}={Nss}$',ha='center')
  assert N0-Nss==8-len(edges)
  result.append([len(edges),N0,Nss])
 assert result==[[4,4,0],[5,3,0],[6,3,1]]
 save(fig,'20_packing_rigidity')
 return result

def cage():
 fig,ax=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
 t=np.geomspace(.01,100,501)
 for Dr,col in zip([.2,1,5],C):
  z=Dr*t;msd=2*(z+np.expm1(-z))/Dr**2
  assert np.all(msd>0)
  ax[0].loglog(t,msd,lw=2.5,c=col,label=rf'$D_r={Dr:g}$')
 ax[0].set(xlabel='Time',ylabel='Mean squared displacement',
  title='A   Free active motion, no cage')
 ax[0].legend(frameon=False)
 tau=np.geomspace(.01,100,501);exact=1/(1+1/tau);white=tau
 ax[1].loglog(tau,exact,lw=2.6,c=C[0],label='Exact harmonic cage')
 ax[1].loglog(tau,white,lw=2.2,c=C[1],ls='--',label='Fast-turning approximation')
 ax[1].axhline(1,c='#9aa7b1',lw=1,ls=':')
 ax[1].set(xlabel=r'Persistence time $\tau_p$',ylabel=r'Steady variance $\langle r^2\rangle$',
  title=r'B   Harmonic confinement, $\kappa=v_0=1$')
 ax[1].legend(frameon=False)
 assert exact[-1]<1 and np.isclose(exact[0]/white[0],1,rtol=.02)
 save(fig,'21_active_cage')

def t1():
 pos={'A':(-1,0),'B':(1,0),'C':(0,1),'D':(0,-1)}
 outer=[('A','C'),('C','B'),('B','D'),('D','A')]
 fig,ax=plt.subplots(1,2,figsize=(11.5,4.9),layout='constrained')
 for a,pair,title in zip(ax,[('A','B'),('C','D')],['A   Before: contact A-B','B   After: contact C-D']):
  for i,j in outer:
   a.plot([pos[i][0],pos[j][0]],[pos[i][1],pos[j][1]],c='#bac5cb',lw=2.5)
  i,j=pair
  a.plot([pos[i][0],pos[j][0]],[pos[i][1],pos[j][1]],lw=4,c=C[1] if pair==('A','B') else C[0])
  for name,(x,y) in pos.items():
   a.scatter([x],[y],s=430,c='#edf3f4',edgecolors='#243746',linewidths=1.5,zorder=3)
   a.text(x,y,name,ha='center',va='center',weight='bold',zorder=4)
  a.set(title=title,aspect='equal',xlim=(-1.65,1.65),ylim=(-1.5,1.5));a.axis('off')
 ax[0].text(0,-1.4,'A-B disappears',ha='center',c=C[1])
 ax[1].text(0,-1.4,r'$\Delta n_A=\Delta n_B=-1;\quad\Delta n_C=\Delta n_D=+1$',ha='center')
 save(fig,'22_t1_topology')

def defects():
 fig,ax=plt.subplots(1,2,figsize=(11.7,5.15),layout='constrained')
 vals=np.linspace(-2,2,23)
 for a,s,title in zip(ax,[.5,-.5],[r'A   $+1/2$ defect',r'B   $-1/2$ defect']):
  segs=[];colors=[]
  for x in vals:
   for y in vals:
    r=np.hypot(x,y)
    if r<.15 or r>2.35:continue
    theta=s*np.arctan2(y,x);S=np.tanh(r/.22)
    d=.073*S*np.array([np.cos(theta),np.sin(theta)])
    segs.append([np.array([x,y])-d,np.array([x,y])+d])
    colors.append((theta/np.pi)%1)
  lc=LineCollection(segs,cmap='twilight',array=np.array(colors),linewidths=2.1,clim=(0,1))
  a.add_collection(lc)
  a.add_patch(Circle((0,0),.15,color='#e5e9eb'))
  a.add_patch(Circle((0,0),1.5,fill=False,ec='#6c7c85',lw=1,ls='--'))
  a.set(title=title,aspect='equal',xlim=(-2.35,2.35),ylim=(-2.35,2.35))
  a.set_xticks([]);a.set_yticks([])
  a.text(.5,-.06,rf'$\Delta\theta={s*2:+g}\pi,\quad s={s:+g}$',ha='center',transform=a.transAxes)
 phi=np.linspace(0,2*np.pi,1001)
 for s in [.5,-.5]:
  angle=.5*np.unwrap(2*s*phi)
  assert np.isclose((angle[-1]-angle[0])/(2*np.pi),s)
 save(fig,'23_nematic_defects')

def topology():
 pos=np.array([[0,0],[1,0],[1,1],[0,1],[.5,.5]])
 edges=[(0,1,.9),(1,2,.7),(2,3,.4),(3,0,.8),(4,0,.6),(4,1,.2),(4,2,.9),(4,3,.3)]
 def stats(selected):
  parent=list(range(5))
  def find(i):
   while parent[i]!=i:i=parent[i]
   return i
  for i,j,w in selected:
   parent[find(i)]=find(j)
  b0=len({find(i) for i in range(5)})
  return b0,len(selected)-5+b0
 fig,ax=plt.subplots(1,3,figsize=(12.8,4.35),layout='constrained')
 checks=[]
 for a,h,title in zip(ax[:2],[.8,.3],[r'A   Keep $w\geq0.8$',r'B   Keep $w\geq0.3$']):
  selected=[e for e in edges if e[2]>=h-1e-10]
  for i,j,w in selected:
   a.plot(pos[[i,j],0],pos[[i,j],1],c=C[0],lw=1.7+2*w,zorder=1)
  a.scatter(pos[:,0],pos[:,1],s=120,c='#243746',zorder=2)
  for i,(x,y) in enumerate(pos):a.text(x,y-.11,str(i),ha='center',va='top',fontsize=10)
  b0,b1=stats(selected);checks.append([b0,b1])
  a.text(.5,-.34,rf'$\beta_0={b0},\quad\beta_1={b1}$',ha='center')
  a.set(title=title,aspect='equal',xlim=(-.25,1.25),ylim=(-.5,1.35));a.axis('off')
 a=ax[2]
 for offset,fill in [(0,False),(1.35,True)]:
  p=np.array([[offset,0],[offset+.95,0],[offset+.475,.82]])
  a.add_patch(Polygon(p,fill=fill,fc='#b9ded7' if fill else 'none',ec=C[0],lw=2.5))
  a.scatter(p[:,0],p[:,1],s=90,c='#243746',zorder=3)
  a.text(offset+.475,-.3,rf'$\beta_1={0 if fill else 1}$',ha='center')
 a.set(title='C   Ring versus filled face',aspect='equal',xlim=(-.2,2.5),ylim=(-.55,1.4));a.axis('off')
 assert checks==[[2,0],[1,3]],checks
 save(fig,'24_force_network_topology')
 return checks

if __name__=='__main__':
 results={'spring_networks':rigidity()}
 cage();t1();defects();results['graph_betti_numbers']=topology()
 print(json.dumps(results,indent=2))
