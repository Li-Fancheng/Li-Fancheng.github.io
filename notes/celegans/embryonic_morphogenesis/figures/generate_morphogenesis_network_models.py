"""Coarse-grained spatial examples for note 07, not fits to embryo data.

Run with NumPy, SciPy, Matplotlib. Grid nodes are numerical material points,
not actin molecules or individual cells. Three figures connect local activity
with tissue deformation, volume-constrained elongation, and retained shape.
"""
from pathlib import Path
import json
import numpy as np
from scipy.spatial import Delaunay
from scipy.linalg import solve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Ellipse

OUT=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','font.family':'DejaVu Sans','savefig.facecolor':'white'})
INK='#304857';TEAL='#27858b';ROSE='#b65376';GOLD='#bd8b24'
report={}
def save(fig,name):
 fig.savefig(OUT/(name+'.svg'),bbox_inches='tight');fig.savefig(OUT/(name+'.png'),dpi=170,bbox_inches='tight');plt.close(fig)

def tissue_patch():
 x=np.linspace(-3,3,15);y=np.linspace(-2,2,11);X,Y=np.meshgrid(x,y);pts=np.c_[X.ravel(),Y.ravel()]
 pts[1::2,0]+=.09*np.sin(2.5*pts[1::2,1]);n=len(pts)
 tris=Delaunay(pts).simplices;edges=sorted({tuple(sorted((t[j],t[(j+1)%3]))) for t in tris for j in range(3)})
 results=[]
 for coupling in [1.,.03]:
  K=.002*np.eye(2*n);f=np.zeros(2*n);rows=[]
  for i,j in edges:
   d=pts[j]-pts[i];l=np.linalg.norm(d);q=d/l;mid=(pts[i]+pts[j])/2;r=np.linalg.norm(mid)
   k=coupling if .8<r<1.5 else 1.
   ids=[2*i,2*i+1,2*j,2*j+1];B=np.r_[-q,q]
   K[np.ix_(ids,ids)]+=k*np.outer(B,B)
   dl=-.08*l if r<.8 else 0.;f[ids]+=k*dl*B;rows.append((i,j,r,k))
  boundary=(abs(pts[:,0])>2.95)|(abs(pts[:,1])>1.95)
  for i in np.where(boundary)[0]:K[2*i,2*i]+=.5;K[2*i+1,2*i+1]+=.5
  u=solve(K,f,assume_a='pos').reshape(n,2)
  error=np.linalg.norm(K@u.ravel()-f)/np.linalg.norm(f);assert error<1e-10
  assert np.linalg.norm(f.reshape(n,2).sum(axis=0))<1e-12
  rad=np.linalg.norm(pts,axis=1);local=(rad>.1)&(rad<.8);far=(rad>1.8)&(~boundary)
  radial=np.sum(pts[local]*u[local],axis=1)/rad[local]
  results.append((coupling,u,rows,{'relative_force_residual':float(error),'mean_inward_displacement':float(-radial.mean()),'far_field_mean_displacement':float(np.linalg.norm(u[far],axis=1).mean())}))
 assert results[0][3]['far_field_mean_displacement']>results[1][3]['far_field_mean_displacement']
 fig,axs=plt.subplots(1,3,figsize=(13.5,4.5),layout='constrained');mag=3
 ref=np.array([[pts[i],pts[j]] for i,j in edges]);axs[0].add_collection(LineCollection(ref,colors=[ROSE if np.linalg.norm((pts[i]+pts[j])/2)<.8 else '#a9bdc8' for i,j in edges],linewidths=1))
 axs[0].set_title('Same local active shortening')
 for ax,(c,u,rows,vals) in zip(axs[1:],results):
  new=pts+mag*u;ax.add_collection(LineCollection(ref,colors='#dce3e8',lw=.6))
  ax.add_collection(LineCollection([[new[i],new[j]] for i,j in edges],colors=[ROSE if r<.8 else TEAL if .8<r<1.5 else '#6f8794' for i,j,r,k in rows],lw=1))
  ax.quiver(pts[::3,0],pts[::3,1],mag*u[::3,0],mag*u[::3,1],angles='xy',scale_units='xy',scale=1,color=INK,width=.004)
  ax.set_title('Continuous coupling' if c==1 else 'Weak coupling around active patch')
  ax.text(.03,.03,f'Far-field motion: {vals["far_field_mean_displacement"]:.4f}',transform=ax.transAxes,color=INK)
 for ax in axs:ax.set_xlim(-3.35,3.35);ax.set_ylim(-2.35,2.35);ax.set_aspect('equal');ax.set_xlabel('Tissue coordinate x');ax.set_ylabel('Tissue coordinate y')
 fig.suptitle('Local activity reaches distant tissue through mechanical connections; displacements ×3',color=INK,fontsize=12)
 save(fig,'37_network_to_tissue');report['tissue_patch']={'nodes':n,'active_preferred_shortening':.08,'basal_anchor':.002,'boundary_anchor':.5,'display_magnification':mag,'cases':{str(c):v for c,u,rows,v in results}}

def cylinder(ax,L,R,x0=0,y0=0,color=TEAL,alpha=.18,label=None):
 # Projected flat-ended circular cylinder; volume = pi R^2 L.
 depth=.22*R
 ax.fill([x0-L/2,x0+L/2,x0+L/2,x0-L/2],[y0-R,y0-R,y0+R,y0+R],color=color,alpha=alpha)
 for xx in np.linspace(-L/2,L/2,9):ax.add_patch(Ellipse((x0+xx,y0),2*depth,2*R,fill=False,edgecolor=color,lw=.8,alpha=.8))
 ax.plot([x0-L/2,x0+L/2],[y0-R,y0-R],color=color,lw=1);ax.plot([x0-L/2,x0+L/2],[y0+R,y0+R],color=color,lw=1)
 if label:ax.text(x0,y0-R-.22,label,ha='center',va='top',fontsize=10,color=INK)

def closed_shell():
 Cz0=.4;Ctheta=2.;tau_theta=1.2;zeta=1.;L0=3.;R0=.8;times=[0,.5,2,6]
 cases=[('Hoop activity',Cz0,0),('Hoop + axial activity',Cz0,.6),('Greater axial resistance',2.,0)]
 fig,axs=plt.subplots(3,4,figsize=(13.5,7.5),layout='constrained');vals={}
 for row,(name,Cz,tau_z) in enumerate(cases):
  A=Cz+Ctheta/4;equil=(tau_theta/2-tau_z)/A;vals[name]={'Cz':Cz,'Ctheta':Ctheta,'tau_theta':tau_theta,'tau_z':tau_z,'equilibrium_log_strain':equil,'snapshots':[]}
  for ax,t in zip(axs[row],times):
   e=equil*(1-np.exp(-A*t/zeta));lam=np.exp(e);L=L0*lam;R=R0/np.sqrt(lam)
   assert abs((R/R0)**2*(L/L0)-1)<1e-12
   cylinder(ax,L,R,color=TEAL if row!=1 else GOLD,label=fr'$L/L_0={lam:.2f}$')
   ax.set_xlim(-3.3,3.3);ax.set_ylim(-1.3,1.3);ax.set_aspect('equal');ax.axis('off');ax.set_title(f't = {t:g}')
   vals[name]['snapshots'].append({'t':t,'lambda':float(lam),'radius_ratio':float(R/R0)})
  axs[row,0].text(-.08,.5,name,rotation=90,transform=axs[row,0].transAxes,va='center',ha='right',color=INK)
 fig.suptitle('A closed volume converts circumferential shortening into axial elongation',color=INK,fontsize=13)
 save(fig,'38_closed_volume_elongation');report['closed_shell']={'L0':L0,'R0':R0,'zeta':zeta,'energy':'0.5(Cz+Ctheta/4)e^2+(tau_z-tau_theta/2)e','cases':vals,'note':'Uniform cylindrical mode only; illustrative, not a heterogeneous embryo simulation.'}

def cyclic_shape():
 dt=.001;T=4.;E=2.;drag=.4;yieldstress=.3;mobility=.12;duration=20.;steps=int(duration/dt)
 time=np.arange(steps+1)*dt;load=np.where(time<16,.9*np.sin(np.pi*time/T)**2,0.)
 histories=[]
 for gp in [0.,mobility]:
  e=np.zeros_like(time);p=np.zeros_like(time)
  for j in range(steps):
   stress=E*(e[j]-p[j]);e[j+1]=e[j]+dt*(load[j]-stress)/drag;p[j+1]=p[j]+dt*gp*max(stress-yieldstress,0)
  histories.append((e,p));assert np.isfinite(e).all()
 assert histories[0][0][-1]<1e-6
 assert abs(histories[1][0][-1]-histories[1][1][-1])<1e-6
 snaps=[0,2,4,10,16,20]
 fig=plt.figure(figsize=(13.5,7.4),layout='constrained');gs=fig.add_gridspec(3,6,height_ratios=[1,1,1.2])
 for row,(e,p) in enumerate(histories):
  for col,t in enumerate(snaps):
   ax=fig.add_subplot(gs[row,col]);j=int(t/dt);lam=np.exp(e[j]);cylinder(ax,2.3*lam,.6/np.sqrt(lam),color=TEAL if row==0 else ROSE,label=f'λ = {lam:.2f}')
   ax.set_xlim(-2.7,2.7);ax.set_ylim(-1.,1.);ax.set_aspect('equal');ax.axis('off');ax.set_title(f't = {t:g}')
   if col==0:ax.text(-.08,.5,'Reversible' if row==0 else 'Reference state updates',rotation=90,transform=ax.transAxes,va='center',ha='right',color=INK)
 ax=fig.add_subplot(gs[2,:]);ax.plot(time,load,ls='--',color='#92a3ae',label='Effective periodic load')
 ax.plot(time,histories[0][0],color=TEAL,label='Reversible strain');ax.plot(time,histories[1][0],color=ROSE,label='Strain with remodeling');ax.plot(time,histories[1][1],color=GOLD,label='Reference strain p')
 ax.axvline(16,color='#b8c2cb',lw=1);ax.set(xlabel='Time t',ylabel='Dimensionless load / strain');ax.legend(ncol=2,loc='upper left',fontsize=9)
 fig.suptitle('Repeated deformations accumulate only when the load changes the reference state',color=INK,fontsize=13)
 save(fig,'39_cyclic_network_shape');report['cyclic_shape']={'dt':dt,'period':T,'E':E,'drag':drag,'yield':yieldstress,'plastic_mobility':mobility,'load_amplitude':.9,'load_stops':16,'final_elastic_strain':float(histories[0][0][-1]),'final_total_strain':float(histories[1][0][-1]),'final_reference_strain':float(histories[1][1][-1]),'note':'Load is a generalized elongation-mode load, not direct longitudinal muscle stress.'}

if __name__=='__main__':
 tissue_patch();closed_shell();cyclic_shape();(OUT/'morphogenesis_network_parameters.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
