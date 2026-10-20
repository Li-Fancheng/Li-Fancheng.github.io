"""Mechanism comparisons for note 06, not fitted C. elegans data.
Run with NumPy and Matplotlib. Includes analytical and numerical checks.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
OUT=Path(__file__).resolve().parent
C=["#287f83","#b2576f","#b88929","#6677aa"]
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,
 "axes.titlesize":13,"axes.labelsize":11,"legend.fontsize":9,
 "axes.spines.top":False,"axes.spines.right":False,"svg.fonttype":"path",
 "savefig.facecolor":"white","figure.facecolor":"white","text.color":"#243746",
 "axes.labelcolor":"#243746"})
def save(fig,name):
 fig.savefig(OUT/(name+".svg"),bbox_inches="tight",pad_inches=.18)
 fig.savefig(OUT/(name+".png"),dpi=160,bbox_inches="tight",pad_inches=.18)
 plt.close(fig)

def ingression():
 fig,ax=plt.subplots(1,3,figsize=(14,4.25),layout="constrained")
 x=np.linspace(-1.6,1.8,701)
 for f,c in zip([0,.2,.5],C):
  ax[0].plot(x,.25*(x*x-1)**2-f*x,c=c,lw=2.3,label=rf"$f={f:g}$")
 fc=2/(3*np.sqrt(3))
 assert np.isclose((-1/np.sqrt(3))**3+1/np.sqrt(3),fc)
 assert len([z for z in np.roots([1,0,-1,-.2]) if abs(z.imag)<1e-9])==3
 assert len([z for z in np.roots([1,0,-1,-.5]) if abs(z.imag)<1e-9])==1
 ax[0].set(title="A   Packing landscape + drive",xlabel=r"Inward coordinate $q/a$",
  ylabel=r"Energy $U/E_0$",ylim=(-1,2))
 ax[0].legend(frameon=False)
 ax[0].text(.03,.94,rf"$f_c=2/(3\sqrt{{3}})\simeq{fc:.3f}$",transform=ax[0].transAxes,va="top")
 ky=np.geomspace(.55,6,400);kq=1.;coupling=.7;keff=kq-coupling**2/ky
 ax[1].semilogx(ky,np.ones_like(ky),c=C[1],lw=2.2,ls="--",label="Neighbors held fixed")
 ax[1].semilogx(ky,keff,c=C[0],lw=2.5,label="Neighbors allowed to relax")
 ax[1].set(title="B   Collective relaxation",xlabel=r"Neighbor stiffness $H_{yy}$",
  ylabel=r"Local stiffness along $q$",ylim=(0,1.18))
 ax[1].legend(frameon=False,loc="lower right")
 ax[1].set_xticks([.6,1,2,4,6],labels=["0.6","1","2","4","6"])
 ax[1].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
 assert np.all(keff>0) and np.all(keff<kq)
 for k in [.6,1.,5.]:
  H=np.array([[kq,coupling],[coupling,k]])
  assert np.linalg.eigvalsh(H).min()>0
  q=.12;y=-coupling*q/k
  energy=.5*np.array([q,y])@H@np.array([q,y])
  assert np.isclose(energy,.5*(kq-coupling**2/k)*q*q)
 xi=np.geomspace(.01,100,500);g=xi/(1+xi);v=(g-.2)/(1+g)
 ax[2].semilogx(xi,v,c=C[0],lw=2.5)
 ax[2].axhline(0,c="#9aa7b1",lw=1)
 ax[2].axvline(.25,c=C[1],lw=1.5,ls=":")
 ax[2].set(title="C   Motor-to-boundary coupling",xlabel=r"Clutch friction $\xi_c$",
  ylabel=r"Boundary velocity $\dot q$",ylim=(-.24,.48))
 ax[2].text(.05,.95,r"$F_m=\zeta_c=\zeta_j=1$"+"\n"+r"$F_{\rm load}=0.2$",
  transform=ax[2].transAxes,va="top")
 assert np.isclose((.25/1.25-.2)/(1+.25/1.25),0)
 assert np.isclose((1-.2)/(1+1),.4)
 save(fig,"25_worm_ingression")
 return {"fold_force":float(fc),"clutch_zero_velocity_xi":.25,
 "fixed_neighbor_stiffness":kq,"coupling":coupling}

def material():
 fig,ax=plt.subplots(1,3,figsize=(14,4.3),layout="constrained")
 angles=np.arange(6)*np.pi/3
 verts=np.c_[np.cos(angles),np.sin(angles)]
 area=3*np.sqrt(3)/2
 stretch=np.linspace(1,2.25,400);shape=[]
 for l in stretch:
  v=verts*np.array([l,1/l]);edges=np.roll(v,-1,axis=0)-v
  shape.append(np.linalg.norm(edges,axis=1).sum()/np.sqrt(area))
 shape=np.array(shape)
 assert np.isclose(shape[0],np.sqrt(24*np.tan(np.pi/6)))
 assert np.all(np.diff(shape)>0)
 ax[0].plot(stretch,shape,c=C[0],lw=2.5)
 ax[0].axhline(3.81,c=C[1],lw=1.4,ls="--",label="3.81 reference")
 ax[0].set(title="A   Fixed neighbors, changing shape",
  xlabel=r"Area-preserving stretch $\lambda$",ylabel=r"Actual shape $q=P/\sqrt{A}$")
 ax[0].legend(frameon=False,loc="upper left")
 inset=ax[0].inset_axes([.07,.43,.39,.24])
 inset.add_patch(Polygon(verts,fill=False,lw=1.4,ec="#87979e"))
 inset.add_patch(Polygon(verts*np.array([1.6,1/1.6]),fill=False,lw=1.8,ec=C[0]))
 inset.set(xlim=(-1.85,1.85),ylim=(-1.1,1.1),aspect="equal");inset.axis("off")
 E=1.;eta=.12;sigmay=.55;etap=1.;dt=.0005
 t=np.arange(0,7+dt/2,dt)
 def load(tt):
  return max(np.sin(2*np.pi*tt),0.) if tt<5 else 0.
 def rhs(tt,z):
  sig=load(tt);ep=z[2]
  return np.array([(sig-E*z[0])/eta,
   (sig-E*(z[1]-ep))/eta,max(E*(z[1]-ep)-sigmay,0)/etap])
 z=np.zeros((len(t),3))
 for k in range(len(t)-1):
  tt=t[k];old=z[k]
  k1=rhs(tt,old);k2=rhs(tt+dt/2,old+dt*k1/2)
  k3=rhs(tt+dt/2,old+dt*k2/2);k4=rhs(tt+dt,old+dt*k3)
  z[k+1]=old+dt*(k1+2*k2+2*k3+k4)/6
 sig=np.array([load(tt) for tt in t])
 elastic_stress=E*(z[:,1]-z[:,2])
 ax[1].plot(t,sig,c=C[3],lw=2,label="Applied stress")
 ax[1].plot(t,elastic_stress,c=C[0],lw=1.8,label="Elastic stress")
 ax[1].axhline(sigmay,c=C[1],lw=1.4,ls="--",label=r"Elastic yield $\sigma_y=0.55$")
 ax[1].axvspan(5,7,color="#edf3f4")
 ax[1].set(title="B   Repeated load, then release",xlabel="Time (pulse periods)",
  ylabel=r"Stress / $E$",ylim=(-.06,1.4))
 ax[1].legend(frameon=False,fontsize=8,loc="upper right")
 for ind,c,ls,label in [(0,C[3],"--","Kelvin-Voigt"),(1,C[0],"-","With reference update"),
  (2,C[1],":",r"Reference strain $\epsilon_p$")]:
  ax[2].plot(t,z[:,ind],c=c,lw=2.1,ls=ls,label=label)
 ax[2].axvspan(5,7,color="#edf3f4")
 ax[2].set(title="C   Recoverable and retained strain",xlabel="Time (pulse periods)",ylabel="Strain",ylim=(-.03,1.3))
 ax[2].legend(frameon=False,loc="upper left")
 assert abs(z[-1,0])<1e-6
 assert z[-1,2]>.1 and abs(z[-1,1]-z[-1,2])<1e-6
 assert np.min(np.diff(z[:,2]))>-1e-12
 # Passive energy balance: external work = elastic storage + dissipation.
 rate=(sig-elastic_stress)/eta
 plastic_rate=np.maximum(elastic_stress-sigmay,0)/etap
 dissipation=eta*rate**2+elastic_stress*plastic_rate
 assert np.min(dissipation)>-1e-12
 work=np.trapezoid(sig*rate,t)
 loss=np.trapezoid(dissipation,t)
 storage=.5*E*(z[-1,1]-z[-1,2])**2
 assert abs(work-loss-storage)<5e-6
 ref=np.trapezoid(plastic_rate,t)
 assert abs(ref-z[-1,2])<2e-6
 save(fig,"26_worm_material_memory")
 return {"E":E,"eta":eta,"sigma_y":sigmay,"eta_p":etap,"dt":dt,"energy_balance_error":float(work-loss-storage),
 "reference_strain_after_five_pulses":float(z[-1,2]),
 "kelvin_voigt_strain_after_release":float(z[-1,0])}
if __name__=="__main__":
 print(json.dumps({"ingression":ingression(),"material":material()},indent=2))
