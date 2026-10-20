"""Calculated figures for active packing and topology; no fitted experimental data."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle,Polygon,FancyArrowPatch

ROOT=Path(__file__).resolve().parent/'active_packing'/'figures'
BLUE,ORANGE,GREEN,INK='#277da8','#d67837','#429582','#263e50'
for name in ['Arial.ttf','Arial Bold.ttf']:
 p=Path('/System/Library/Fonts/Supplemental')/name
 if p.exists():font_manager.fontManager.addfont(p)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':14,'axes.titlesize':15,'axes.labelsize':14,'legend.fontsize':12,'mathtext.fontset':'cm','axes.spines.top':False,'axes.spines.right':False,'text.color':INK,'axes.labelcolor':INK,'axes.titlecolor':INK,'lines.linewidth':2.2,'svg.fonttype':'path','svg.hashsalt':'vertex-active-packing','savefig.facecolor':'white'})

def save(fig,name):
 ROOT.mkdir(parents=True,exist_ok=True)
 fig.savefig(ROOT/(name+'.svg'),bbox_inches='tight',pad_inches=.2,metadata={'Date':None})
 fig.savefig(ROOT/(name+'.png'),bbox_inches='tight',pad_inches=.2,dpi=300)
 plt.close(fig)

def motility_instability():
 fig,axs=plt.subplots(1,3,figsize=(12.5,3.7),layout='constrained')
 x=np.linspace(0,1,401);axs[0].plot(x,1-x,color=BLUE)
 axs[0].set(xlabel=r'Density, $x=\rho/\rho_*$',ylabel=r'$v/v_0$',title='A  Crowding slows motion',ylim=(-.04,1.05))
 for delta,c in [(0,BLUE),(.06,ORANGE),(.15,GREEN)]:
  deff=delta+(1-x)*(1-2*x);axs[1].plot(x,deff,color=c,label=rf'$\delta={delta:g}$')
 axs[1].axhline(0,color=INK,lw=1);axs[1].axvspan(.5,1,color=BLUE,alpha=.07)
 axs[1].set(xlabel=r'$x=\rho/\rho_*$',ylabel=r'$2D_r D_{\rm coll}/v_0^2$',title='B  Collective diffusivity',ylim=(-.17,1.2));axs[1].legend()
 q=np.linspace(0,.34,300)
 for D,c in [(-.08,ORANGE),(.06,BLUE)]:axs[2].plot(q,-D*q*q-q**4,color=c,label=rf'$D_{{\rm coll}}={D:g}$')
 axs[2].axhline(0,color=INK,lw=1);axs[2].set(xlabel=r'Wave number, $q$',ylabel=r'Growth rate, $\lambda(q)$',title='C  A growing density mode',ylim=(-.012,.0022));axs[2].legend(loc='lower left');save(fig,'motility_instability')

def order_correlations():
 fig=plt.figure(figsize=(12,4.2),layout='constrained');grid=fig.add_gridspec(1,3)
 for j,theta in enumerate([0,np.pi/12]):
  ax=fig.add_subplot(grid[0,j]);angles=np.arange(6)*np.pi/3+theta
  points=np.c_[np.cos(angles),np.sin(angles)];ax.scatter(points[:,0],points[:,1],s=95,color=BLUE);ax.scatter(0,0,s=95,color=INK)
  for p in points:ax.plot([0,p[0]],[0,p[1]],color=BLUE,alpha=.5,lw=1.5)
  ax.set(aspect='equal',xlim=(-1.3,1.3),ylim=(-1.35,1.3),title=['A  Reference bonds','B  Rotated by 15 degrees'][j]);ax.axis('off')
  ax.text(0,-1.25,[r'$|\psi_6|=1,\quad \arg\psi_6=0$',r'$|\psi_6|=1,\quad \arg\psi_6=\pi/2$'][j],ha='center',fontsize=16)
 ax=fig.add_subplot(grid[0,2]);r=np.geomspace(1,60,300)
 ax.loglog(r,np.exp(-(r-1)/3),color=BLUE,label='Exponential');ax.loglog(r,r**(-.25),color=ORANGE,label=r'Power law: $r^{-1/4}$');ax.loglog(r,.6+.4*np.exp(-(r-1)/3),color=GREEN,label='Nonzero limit')
 ax.set(xlabel=r'Distance, $r/a$',ylabel=r'$g_6(r)$',title='C  Spatial correlations',ylim=(.008,1.15));ax.legend(loc='lower left');save(fig,'order_correlations')

def edge(ax,a,b,color=INK,lw=2,ls='-'):
 ax.plot([a[0],b[0]],[a[1],b[1]],color=color,lw=lw,ls=ls,zorder=2)

def graph_topology():
 fig=plt.figure(figsize=(12.5,7),layout='constrained');grid=fig.add_gridspec(2,6,height_ratios=[1.15,1])
 ax=fig.add_subplot(grid[0,:3]);p=np.array([[0,.9],[.9,0],[0,-.9],[-.9,0]])
 for offset,diag,title in [(np.array([-1.25,0]),(0,2),'Old edge'),(np.array([1.25,0]),(1,3),'New edge')]:
  ps=p+offset
  for i in range(4):edge(ax,ps[i],ps[(i+1)%4])
  edge(ax,ps[diag[0]],ps[diag[1]],ORANGE,3)
  ax.scatter(*ps.T,s=55,color=INK,zorder=3);ax.text(offset[0],1.18,title,ha='center')
 ax.text(0,-1.25,r'Old-edge ends: $\Delta z=-1$',ha='center',fontsize=13)
 ax.text(0,-1.65,r'New-edge ends: $\Delta z=+1$',ha='center',fontsize=13)
 ax.set(aspect='equal',xlim=(-2.5,2.5),ylim=(-1.95,1.5),title='A  Neighbour exchange');ax.axis('off')
 ax=fig.add_subplot(grid[0,3:]);tri=np.array([[0,0],[1,0],[.5,.86]])
 for offset,filled in [(np.array([-.9,0]),False),(np.array([.9,0]),True)]:
  ps=tri+offset
  if filled:ax.add_patch(Polygon(ps,facecolor='#d9ece8',edgecolor='none'))
  for i in range(3):edge(ax,ps[i],ps[(i+1)%3],GREEN if filled else BLUE)
  ax.scatter(*ps.T,color=INK,s=55,zorder=3)
  ax.text(offset[0]+.5,-.22,r'$\beta_1=0$' if filled else r'$\beta_1=1$',ha='center',fontsize=17)
  ax.text(offset[0]+.5,1.12,'Filled face' if filled else 'Edges only',ha='center')
 ax.set(aspect='equal',xlim=(-1.1,2.1),ylim=(-.45,1.45),title='B  What is counted as a hole?');ax.axis('off')
 sq=np.array([[0,0],[1,0],[1,1],[0,1]]);weights=[4,3,2,1]
 for j,(threshold,beta0,beta1) in enumerate([(3.5,3,0),(1.5,1,0),(.5,1,1)]):
  ax=fig.add_subplot(grid[1,2*j:2*j+2])
  for i,w in enumerate(weights):
   if w>=threshold:edge(ax,sq[i],sq[(i+1)%4],BLUE,lw=2.5)
   else:edge(ax,sq[i],sq[(i+1)%4],'#c9cfd5',lw=1,ls=':')
  ax.scatter(*sq.T,s=70,color=INK,zorder=3);ax.text(.5,-.3,rf'$\beta_0={beta0},\quad\beta_1={beta1}$',ha='center',fontsize=16)
  ax.set(aspect='equal',xlim=(-.28,1.28),ylim=(-.55,1.28),title=rf'C{j+1}  Keep edges with $F\geq {threshold:g}$');ax.axis('off')
 save(fig,'graph_topology')

def nematic_charges():
 fig,axs=plt.subplots(1,3,figsize=(12,4.1),layout='constrained')
 for ax,s in zip(axs,[.5,-.5,1.]):
  for x in np.linspace(-1,1,17):
   for y in np.linspace(-1,1,17):
    rr=np.hypot(x,y)
    if rr<.18 or rr>1.1:continue
    theta=s*np.arctan2(y,x);vec=.046*np.array([np.cos(theta),np.sin(theta)])
    edge(ax,np.array([x,y])-vec,np.array([x,y])+vec,INK,1.35)
  ax.add_patch(Circle((0,0),.78,facecolor='none',edgecolor=ORANGE,ls='--',lw=1.5));ax.add_patch(Circle((0,0),.07,color=BLUE))
  ax.add_patch(FancyArrowPatch((.78,0),(.765,.15),arrowstyle='-|>',mutation_scale=14,color=ORANGE,lw=1.5))
  ax.text(0,-1.25,[r'$\Delta\theta=+\pi$',r'$\Delta\theta=-\pi$',r'$\Delta\theta=+2\pi$'][[.5,-.5,1.].index(s)],ha='center',fontsize=17)
  ax.set(aspect='equal',xlim=(-1.25,1.25),ylim=(-1.4,1.2),title=rf'$s={s:+g}$');ax.axis('off')
 save(fig,'nematic_charges')

if __name__=='__main__':
 motility_instability();order_correlations();graph_topology();nematic_charges()
 print('Generated four calculated figures (SVG and PNG).')
