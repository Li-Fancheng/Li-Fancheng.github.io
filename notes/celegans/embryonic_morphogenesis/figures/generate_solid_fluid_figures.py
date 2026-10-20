"""Analytic rheology comparisons for note 06; no fitted embryo measurements.
NumPy and Matplotlib. Solid/fluid terminology refers to the defined models.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent
C=["#287f83","#b2576f","#6677aa","#b88929"]
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,
 "axes.titlesize":13,"axes.labelsize":11,"legend.fontsize":9,
 "axes.spines.top":False,"axes.spines.right":False,"svg.fonttype":"path",
 "savefig.facecolor":"white","figure.facecolor":"white","text.color":"#243746",
 "axes.labelcolor":"#243746"})
def save(fig,stem):
 fig.savefig(OUT/(stem+".svg"),bbox_inches="tight",pad_inches=.18)
 fig.savefig(OUT/(stem+".png"),dpi=160,bbox_inches="tight",pad_inches=.18)
 plt.close(fig)

def rheology():
 t=np.geomspace(.001,1000,1000)
 # (G_infinity, relaxing modulus, stress-relaxation time)
 cases=[(0,1,1,"Maxwell fluid",C[0]),
        (0,1,100,"Slow Maxwell fluid",C[2]),
        (.3,.7,1,"Viscoelastic solid",C[1])]
 fig,ax=plt.subplots(1,3,figsize=(14.4,4.5),layout="constrained")
 for gi,gr,tau,label,c in cases:
  ax[0].semilogx(t,gi+gr*np.exp(-t/tau),c=c,lw=2.3,label=label)
 ax[0].set(title="A   Stress relaxation",xlabel="Time",ylabel=r"$G(t)$",ylim=(-.025,1.15))
 ax[0].legend(frameon=False,loc="upper right",fontsize=8.5)
 w=np.geomspace(.001,1000,1000)
 for gi,gr,tau,label,c in [cases[0],cases[2]]:
  x=w*tau;gp=gi+gr*x*x/(1+x*x);gpp=gr*x/(1+x*x)
  ax[1].loglog(w,gp,c=c,lw=2.3,label=label+r" $G'$")
  ax[1].loglog(w,gpp,c=c,lw=1.9,ls="--",label=label+r" $G''$")
 ax[1].set(title="B   Oscillatory response",xlabel=r"Angular frequency $\omega$",
  ylabel="Shear modulus",ylim=(1e-5,2))
 ax[1].legend(frameon=False,loc="lower right",fontsize=8)
 tc=np.linspace(0,15,501)
 for gi,gr,tau,label,c in cases:
  g0=gi+gr
  if gi==0:
   j=1/g0+tc/(gr*tau)
  else:
   retard=tau*g0/gi
   j=1/gi-(1/gi-1/g0)*np.exp(-tc/retard)
  ax[2].plot(tc,j,c=c,lw=2.3,label=label)
 ax[2].set(title="C   Creep under fixed stress",xlabel="Time",ylabel=r"Compliance $J(t)$")
 ax[2].legend(frameon=False,loc="upper left",fontsize=8.5)
 # Initial compliance and steady modulus/viscosity, plus creep-relaxation duality.
 gi=.3;gr=.7;tau=1.;g0=gi+gr;retard=tau*g0/gi
 assert np.isclose(gi+gr,g0) and np.isclose(1/gi-(1/gi-1/g0),1/g0)
 assert np.isclose((1/gi-(1/gi-1/g0)*np.exp(-1e4/retard)),1/gi)
 # Laplace transforms must obey G_tilde(s)*J_tilde(s)=1/s^2.
 for lap in [.01,.1,1,10]:
  gt=gi/lap+gr/(lap+1/tau)
  jt=1/(gi*lap)-(1/gi-1/g0)/(lap+1/retard)
  assert np.isclose(gt*jt,1/lap**2,rtol=1e-11)
 assert np.isclose(1*100,100)  # Slow Maxwell viscosity eta=G_r*tau.
 save(fig,"27_solid_fluid_rheology")
 return {"G0":1,"Maxwell_tau":[1,100],"solid_G_infinity":gi,
         "solid_tau_relax":tau,"solid_tau_creep":retard}

def yielding():
 fig,ax=plt.subplots(1,2,figsize=(11.6,4.5),layout="constrained")
 rate=np.geomspace(1e-4,3,700);K=.6;n=.5
 for sy,c in zip([0,.4,.8],C):
  sig=sy+K*rate**n
  ax[0].semilogx(rate,sig,c=c,lw=2.4,label=rf"$\sigma_y={sy:g}$")
  ax[0].axhline(sy,c=c,ls=":",lw=.9,alpha=.65)
 ax[0].set(title="A   Steady flow curves",xlabel=r"Shear rate $\dot\gamma$",
  ylabel=r"Stress $\sigma$",ylim=(-.04,2))
 ax[0].legend(frameon=False,loc="upper left")
 # Elastic spring in series with a Bingham element.
 t=np.linspace(0,8,401);G=1.;eta_p=1.;sy=.4
 for sig,c in zip([.2,.6],C):
  pl_rate=max(sig-sy,0)/eta_p
  strain=sig/G+pl_rate*t
  ax[1].plot(t,strain,c=c,lw=2.4,
   label=rf"$\sigma_0={sig:g}$, "+("below yield" if sig<sy else "above yield"))
  assert np.all(np.diff(strain)>=-1e-12)
 ax[1].set(title="B   One solid, two applied loads",xlabel="Time",
  ylabel=r"Strain $\gamma$",ylim=(0,2.4))
 ax[1].legend(frameon=False,loc="upper left")
 assert max(.2-sy,0)==0 and np.isclose(max(.6-sy,0)/eta_p,.2)
 save(fig,"28_yield_and_flow")
 return {"Herschel_Bulkley_K":K,"n":n,"yield_stress_values":[0,.4,.8],
         "creep_G":G,"creep_eta_p":eta_p,"creep_yield_stress":sy}
if __name__=="__main__":
 print(json.dumps({"linear_rheology":rheology(),"yielding":yielding()},indent=2))
