"""Spatial morphogenesis examples for note 02. No embryo data are fitted.
Requires NumPy, SciPy and Matplotlib. Original paper images are not modified.
All coordinates, times and coefficients below are dimensionless.
"""
from pathlib import Path
import json
import numpy as np
from scipy.ndimage import laplace
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
OUT=Path(__file__).resolve().parent
C=["#287f83","#b2576f","#b88929","#6677aa"]
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
 "axes.titlesize":12,"axes.labelsize":10,"legend.fontsize":9,
 "axes.spines.top":False,"axes.spines.right":False,"svg.fonttype":"path",
 "savefig.facecolor":"white","figure.facecolor":"white","text.color":"#243746",
 "axes.labelcolor":"#243746"})
def save(fig,stem):
 fig.savefig(OUT/(stem+".svg"),bbox_inches="tight",pad_inches=.18)
 fig.savefig(OUT/(stem+".png"),dpi=160,bbox_inches="tight",pad_inches=.18)
 plt.close(fig)

def phase_layers():
 h=.05;dt=.01;eps=.14;Gamma=.08;KA=12.;rep=2.;v0=.30
 x=np.arange(-3,3.001,h);y=np.arange(-2,2.001,h);X,Y=np.meshgrid(x,y)
 phi0=.5*(1-np.tanh((np.sqrt((X+1.6)**2+Y**2)-.66)/.16))
 neighbors=np.array([.5*(1-np.tanh((np.sqrt((X-.6)**2+(Y-yy)**2)-.8)/.16)) for yy in [-1.05,1.05]])
 psi2=(neighbors**2).sum(axis=0);A0=(phi0**2).sum()*h*h
 fig,ax=plt.subplots(2,4,figsize=(13.5,7.2),layout="constrained")
 records={};max_error=0.
 for row,kab in enumerate([.1,.8]):
  a=phi0.copy();b=phi0.copy();frames=[];hist=[]
  def area(q):return (q*q).sum()*h*h
  def centroid(q):return (X*q*q).sum()*h*h/area(q)
  def chemical(q,other):
   return (-Gamma*eps*laplace(q,mode="nearest")/h**2
    +2*Gamma/eps*q*(1-q)*(1-2*q)
    +2*KA/A0*(area(q)-A0)*q+2*rep*q*psi2+kab*(q-other))
  for step in range(1201):
   if step%400==0:frames.append((a.copy(),b.copy()))
   if step%20==0:
    hist.append([step*dt,area(a),area(b),centroid(a),centroid(b)])
    max_error=max(max_error,abs(area(a)/A0-1),abs(area(b)/A0-1))
   da=-chemical(a,b)
   upwind=(b-np.pad(b,((0,0),(1,0)),mode="edge")[:,:-1])/h
   db=-chemical(b,a)-v0*upwind
   a=np.clip(a+dt*da,0,1.1);b=np.clip(b+dt*db,0,1.1)
  records[str(kab)]=hist
  for col,(a,b) in enumerate(frames):
   aa=ax[row,col]
   for mask in neighbors:aa.contourf(X,Y,mask,levels=[.5,1.01],colors=["#dce2e5"])
   aa.contourf(X,Y,b,levels=[.5,1.1],colors=["#b5dbd6"],alpha=.9)
   aa.contour(X,Y,b,levels=[.5],colors=[C[0]],linewidths=1.8)
   aa.contour(X,Y,a,levels=[.5],colors=[C[1]],linewidths=1.9)
   aa.axvline(.6,c="#a9b4bb",ls=":",lw=.8)
   aa.set(xlim=(-2.65,3.05),ylim=(-1.95,1.95),aspect="equal",
          title=rf"$k_{{ab}}={kab:g}$, $t={col*4}$")
   if row==1:aa.set_xlabel("Inward coordinate")
   if col==0:aa.set_ylabel("Weak coupling" if row==0 else "Strong coupling")
   aa.set_xticks([-2,0,2]);aa.set_yticks([-1,0,1])
 fig.legend([plt.Line2D([],[],color=C[0],lw=2),plt.Line2D([],[],color=C[1],lw=2),
  plt.Rectangle((0,0),1,1,fc="#dce2e5")],["Basal outline","Apical outline","Fixed neighbor masks"],
  loc="outside upper center",ncol=3,frameon=False)
 assert max_error<.03
 assert records["0.1"][-1][4]-records["0.1"][-1][3]>2.5
 assert records["0.8"][-1][4]-records["0.8"][-1][3]<.3
 save(fig,"29_basal_apical_fields")
 return dict(h=h,dt=dt,epsilon=eps,Gamma=Gamma,KA=KA,repulsion=rep,v_basal=v0,
  initial_area=float(A0),max_relative_area_error=float(max_error),
  final_centroid_lag={k:float(v[-1][4]-v[-1][3]) for k,v in records.items()})

def polygons(q):
 nw=np.array([-1.,1.]);ne=np.array([1.,1.]);sw=np.array([-1.,-1.]);se=np.array([1.,-1.])
 if q>=0:
  u=np.array([0,q/2]);v=np.array([0,-q/2])
  return [np.array([sw,nw,u,v]),np.array([ne,se,v,u]),np.array([nw,ne,u]),np.array([se,sw,v])]
 u=np.array([q/2,0]);v=np.array([-q/2,0])
 return [np.array([sw,nw,u]),np.array([ne,se,v]),np.array([nw,ne,v,u]),np.array([se,sw,u,v])]
def area_perimeter(poly):
 area=.5*abs(np.sum(poly[:,0]*np.roll(poly[:,1],-1)-poly[:,1]*np.roll(poly[:,0],-1)))
 per=np.linalg.norm(np.roll(poly,-1,axis=0)-poly,axis=1).sum()
 return area,per
def patch_energy(q,lam=0):
 vals=np.array([area_perimeter(v) for v in polygons(q)])
 return np.sum(.5*(vals[:,0]-1)**2+.5*.5*(vals[:,1]-3.8)**2)+lam*max(q,0)
def t1_patch():
 fig=plt.figure(figsize=(13.5,7),layout="constrained")
 gs=fig.add_gridspec(2,5,height_ratios=[1.05,1])
 for j,q in enumerate([.8,.3,0,-.3,-.8]):
  a=fig.add_subplot(gs[0,j])
  for k,poly in enumerate(polygons(q)):
   a.add_patch(Polygon(poly,facecolor=["#d1e8e3","#dce1f3","#f1ddce","#ead5de"][k],ec="#566971",lw=1.5))
   v=poly.mean(axis=0);a.text(*v,["A","B","C","D"][k],ha="center",va="center",weight="bold")
  if q>0:
   a.plot([0,0],[-q/2,q/2],c=C[1],lw=3.2)
   if j==0:
    for sign in [-1,1]:
     a.annotate("",xy=(.10,sign*q/2),xytext=(.10,sign*(q/2+.32)),
      arrowprops={"arrowstyle":"->","color":C[1],"lw":1.5})
  if q<0:a.plot([q/2,-q/2],[0,0],c=C[0],lw=3.2)
  a.set(aspect="equal",xlim=(-1.08,1.08),ylim=(-1.08,1.08),title=rf"$q={q:g}$")
  a.axis("off")
 q=np.linspace(-1.5,1.5,1001);en=np.array([patch_energy(z) for z in q]);emin=en.min()
 a=fig.add_subplot(gs[1,:3])
 for lam,c in zip([0,.25,.55],C):
  a.plot(q,[patch_energy(z,lam)-emin for z in q],c=c,lw=2,label=rf"$\Lambda_{{AB}}={lam}$")
 a.axvline(0,c="#a9b4bb",ls=":",lw=1)
 a.set(xlabel=r"Signed junction coordinate $q$",ylabel="Constrained patch energy",
  title="Geometry-derived energy along a T1 path")
 a.legend(frameon=False,fontsize=9)
 inset=a.inset_axes([.29,.46,.37,.43])
 qq=np.linspace(-.6,.6,301)
 for lam,c in zip([0,.25,.55],C):
  inset.plot(qq,[patch_energy(z,lam)-emin for z in qq],c=c,lw=1.3)
 inset.axvline(0,c="#a9b4bb",ls=":",lw=.7)
 inset.set(xlim=(-.6,.6),ylim=(-.005,.25),title="Near the exchange")
 inset.tick_params(labelsize=7);inset.title.set_fontsize(9)
 b=fig.add_subplot(gs[1,3:]);lams=np.linspace(0,.6,121);barriers=[];mins=[]
 for lam in lams:
  m=minimize_scalar(lambda z:patch_energy(z,lam),bounds=(0,1.5),method="bounded")
  barriers.append(max(patch_energy(0,lam)-m.fun,0));mins.append(m.x)
 b.plot(lams,barriers,c=C[1],lw=2.4)
 b.set(xlabel=r"Extra tension on A-B, $\Lambda_{AB}$",ylabel="Barrier to the four-way vertex",
  title="Old-contact tension lowers the barrier")
 assert max(abs(sum(area_perimeter(v)[0] for v in polygons(z))-4) for z in q)<1e-12
 assert barriers[0]>.01 and barriers[-1]<1e-5
 assert np.max(np.diff(barriers))<1e-6
 save(fig,"30_four_cell_rearrangement")
 return dict(KA=1,KP=.5,A0=1,P0=3.8,initial_barrier=float(barriers[0]),
  high_tension_barrier=float(barriers[-1]),outer_box=[-1,1,-1,1])

def enclosure():
 s=np.linspace(-2,2,161);ds=s[1]-s[0];dt=.001;w0=.65*(1-(s/2)**2)
 fig,ax=plt.subplots(2,4,figsize=(13.5,6.5),layout="constrained")
 results={}
 for row,kon in enumerate([3.,.03]):
  w=w0.copy();c=np.zeros_like(w);frames=[]
  for step in range(8001):
   t=step*dt
   if step in [0,2000,5000,8000]:frames.append((w.copy(),c.copy()))
   lap=np.r_[0,np.diff(w,2)/ds**2,0]
   fp=(.18+.35*np.exp(-((s+.35)/.35)**2)
     +(.12/(1+np.exp(-12*(s-.2))) if t>2 else 0)) if t<5 else np.zeros_like(w)
   c+=dt*(kon*(w<.035)*(1-c)-.05*c)
   w=np.maximum(0,w+dt*(.025*lap-fp+.12-.32*c*np.exp(-w/.06)))
   w[[0,-1]]=0
  for col,(w,c) in enumerate(frames):
   a=ax[row,col];outer=1.2-.04*s*s
   boundaries=np.array([-2,-1.4,-.8,-.2,.4,1.,1.5,2.])
   for k in range(len(boundaries)-1):
    z=(s>=boundaries[k])&(s<=boundaries[k+1]);ss=s[z]
    a.fill_between(ss,w[z],outer[z],fc=["#d1e8e3","#bcd7d0"][k%2],ec="#6f8a88",lw=.7)
    a.fill_between(ss,-outer[z],-w[z],fc=["#dce1f3","#c5cde6"][k%2],ec="#7c85a1",lw=.7)
   a.fill_between(s,-w,w,fc="#f5f7f8")
   a.plot(s,w,c="#506a70",lw=1.1);a.plot(s,-w,c="#506a70",lw=1.1)
   contact=(w<.035)&(abs(s)<1.9)
   # Color only interfaces which remain close, using the evolving adhesion variable.
   a.scatter(s[contact],np.zeros(contact.sum()),c=c[contact],cmap="magma",vmin=0,vmax=1,
    s=8,linewidths=0,zorder=4)
   a.axvline(-.35,c="#82949d",ls=":",lw=.7)
   a.set(aspect="equal",xlim=(-2.05,2.05),ylim=(-1.28,1.28),
    title=rf"$k_{{on}}={kon:g}$, $t={[0,2,5,8][col]}$")
   a.set_xticks([]);a.set_yticks([])
   if col==0:a.set_ylabel("Fast adhesion" if row==0 else "Slow adhesion")
   if row==1:a.set_xlabel("Anterior  \u2190        AP position        \u2192  Posterior")
  leading=abs(s+.35)<.1
  results[str(kon)]=dict(final_leading_half_gap=float(w[leading].mean()),
   final_max_half_gap=float(w.max()),final_leading_adhesion=float(c[leading].mean()),
   half_gap_at_drive_stop=float(frames[2][0][leading].mean()))
 assert results["3.0"]["final_leading_half_gap"]<.001
 assert results["0.03"]["final_leading_half_gap"]>.2
 assert results["0.03"]["half_gap_at_drive_stop"]<.001
 fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0,1),cmap="magma"),
  ax=ax,shrink=.65,pad=.02,label=r"Contact maturation $c_j$")
 fig.suptitle("Front positions and contact maturation; protrusive drive stops at t=5",fontsize=12)
 save(fig,"31_enclosure_and_adhesion")
 return dict(dt=dt,ds=float(ds),smoothing=.025,outward_force=.12,adhesion_force=.32,
  adhesion_range=.06,contact_threshold=.035,koff=.05,results=results)
if __name__=="__main__":
 result={"two_layers":phase_layers(),"four_cell_patch":t1_patch(),"enclosure":enclosure()}
 (OUT/"epithelial_spatial_model_parameters.json").write_text(json.dumps(result,indent=2)+"\n")
 print(json.dumps(result,indent=2))
