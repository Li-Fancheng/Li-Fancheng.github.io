"""Dimensionless model comparisons for the five morphogenesis notes.
Equations and assumptions are given beside each figure. No embryo data are fitted.
Run with NumPy, SciPy and Matplotlib. Existing paper figures are not modified.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
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

def gastrulation():
 fig,ax=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
 r=5.
 n=np.linspace(.0001,r/(1+r),3000)
 a=n*np.log(r*(1-n)/n)
 critical=minimize_scalar(lambda x:-x*np.log(r*(1-x)/x),
  bounds=(.0001,r/(1+r)),method='bounded')
 nc=critical.x; ac=-critical.fun
 ax[0].plot(a[n>=nc],n[n>=nc],c=C[0],lw=2.5,label='Stable bound fraction')
 ax[0].plot(a[n<nc],n[n<nc],c=C[1],lw=2.5,ls='--',label='Unstable branch')
 ax[0].scatter([ac],[nc],c=C[2],s=45,zorder=3)
 ax[0].axvline(ac,c=C[2],lw=1,ls=':')
 ax[0].text(ac*.97,.8,r'Capacity $\alpha_c$',ha='right',color=C[2])
 ax[0].set(xlabel=r'Total load $\alpha=F/(Nf_d)$',ylabel=r'Bound fraction $n$',
  title='A   Force sharing and bond failure',xlim=(0,ac*1.17),ylim=(0,1))
 ax[0].legend(frameon=False,loc='lower left')
 # Stability is dg/dn<0; implicit relation gives dg/dn = exp(alpha/n)*alpha'(n).
 aprime=np.log(r*(1-n)/n)-1/(1-n)
 assert np.all(aprime[n>nc+1e-3]<0)
 q=np.linspace(0,5,501)
 for aa,col in zip([1,2,4],C):
  s=aa*q*q/(1+q*q)-.25*q*q-1
  ax[1].plot(q,s,lw=2.5,c=col,label=rf'$a={aa}$')
 ax[1].axhline(0,c='#9aa7b1',lw=1)
 ax[1].set(xlabel=r'Wavenumber $q=k\ell$',ylabel=r'Growth rate $s\tau$',
  title='B   Cortical aggregation instability',ylim=(-4,1.6))
 ax[1].legend(frameon=False)
 assert (4*q*q/(1+q*q)-.25*q*q-1).max()>0
 save(fig,'15_gastrulation_mechanics')
 return {'clutch_capacity':float(ac),'critical_bound_fraction':float(nc)}

def epithelial():
 fig,ax=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
 k=np.linspace(0,2,200)
 for f,col in zip([0,.4,.8],C):
  ax[0].plot(k,k+f,lw=2.5,c=col,label=rf'$f_m={f}$')
 ax[0].set(xlabel=r'Inward curvature $\kappa$',ylabel=r'Inward speed $v_n$',
  title='A   Curvature and migration contributions')
 ax[0].legend(frameon=False)
 q=np.linspace(-1.6,1.6,501)
 for f,col in zip([0,.2,.42],C):
  ax[1].plot(q,(q*q-1)**2/4-f*q,lw=2.5,c=col,label=rf'Drive $f={f}$')
 ax[1].set(xlabel=r'Rearrangement coordinate $q$',ylabel='Effective energy',
  title='B   A driven barrier crossing example',ylim=(-.6,1.2))
 ax[1].legend(frameon=False)
 save(fig,'16_epithelial_models')

def early():
 fig,ax=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
 lam=np.linspace(.65,1.4,401);I4=lam**2
 ax[0].plot(lam,.5*(I4-1)**2,lw=2.5,c=C[1],label='Bilateral fiber')
 ax[0].plot(lam,.5*np.maximum(I4-1,0)**2,lw=2.5,c=C[0],ls='--',label='Tension-only fiber')
 ax[0].axvline(1,c='#9aa7b1',lw=1)
 ax[0].set(xlabel=r'Fiber stretch $\lambda_f$',ylabel=r'Energy $W_f/k_f$',
  title='A   Does a compressed fiber carry load?')
 ax[0].legend(frameon=False)
 t=np.linspace(0,14,501)
 for size,col in zip([1,2,4],C):
  ax[1].plot(t,np.exp(-t/size**2),lw=2.5,c=col,label=rf'Size $a={size}$')
 ax[1].set(xlabel='Time',ylabel=r'Pressure mode $p/p_0$',
  title=r'B   Diffusive relaxation: $\tau\propto a^2$',ylim=(0,1.05))
 ax[1].legend(frameon=False)
 save(fig,'17_elongation_materials')

def late():
 fig,ax=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
 w=np.geomspace(.02,50,401)
 ax[0].semilogx(w,1+w*w/(1+w*w),lw=2.5,c=C[0],label=r'Storage $E^\prime$')
 ax[0].semilogx(w,w/(1+w*w),lw=2.5,c=C[1],label=r'Loss $E^{\prime\prime}$')
 ax[0].axvline(1,c='#9aa7b1',lw=1,ls=':')
 ax[0].set(xlabel=r'Frequency $\omega\tau$',ylabel='Modulus',
  title='A   Standard linear solid')
 ax[0].legend(frameon=False)
 T=np.geomspace(.01,100,1001)
 for tr,tl,col in [(.2,1,C[0]),(1,1,C[1]),(1,3,C[2])]:
  rate=(-np.expm1(-.5*T/tr))*(-np.expm1(-.5*T/tl))/T
  ax[1].semilogx(T,rate,lw=2.5,c=col,label=rf'$\tau_r={tr:g},\ \tau_l={tl:g}$')
  assert 0<np.argmax(rate)<len(T)-1
 ax[1].set(xlabel=r'Cycle duration $T$',ylabel=r'Mean reference-strain rate $g$',
  title='B   Two steps require time to complete')
 ax[1].legend(frameon=False)
 save(fig,'18_late_model_comparison')

def geometry():
 fig,ax=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
 q=np.linspace(0,2.4,401)
 for N,col in zip([1,2,3],C):
  ax[0].plot(q,N*q*q-q**4-1,lw=2.5,c=col,label=rf'Compression $N={N}$')
 ax[0].axhline(0,c='#9aa7b1',lw=1)
 ax[0].set(xlabel=r'Wavenumber $q$',ylabel=r'Growth rate $s(q)$',
  title='A   Bending against an elastic background',ylim=(-4,1.5))
 ax[0].legend(frameon=False)
 t=np.linspace(.05,10,40);z=t/2;u=1-np.exp(-z);v=z*np.exp(-z)
 J=np.stack([u,-u+v,-v],axis=1)/.05
 ev0=np.linalg.eigvalsh(J.T@J)
 Jplus=np.vstack([J,np.array([0,1,0])/.05])
 ev1=np.linalg.eigvalsh(Jplus.T@Jplus)
 assert np.linalg.norm(J@np.ones(3))<1e-10
 assert abs(ev0[0])<1e-9 and ev1[0]>1
 vals0=np.maximum(ev0,1e-3)
 x=np.arange(3)
 ax[1].bar(x-.18,vals0,width=.36,color=C[1],label='Release curve only')
 ax[1].bar(x+.18,ev1,width=.36,color=C[0],label=r'+ independent $\ln K$')
 ax[1].set_yscale('log')
 ax[1].set(xticks=x,xticklabels=['Smallest','Middle','Largest'],
  ylabel='Information eigenvalue',title='B   Breaking a parameter symmetry',ylim=(.0003,50000))
 ax[1].annotate('Exact zero',xy=(-.18,.001),xytext=(.15,.02),
  arrowprops={'arrowstyle':'->','color':'#243746'},fontsize=10)
 ax[1].legend(frameon=False,loc='upper left')
 save(fig,'19_rod_inference')
 return {'release_information_eigenvalues':ev0.tolist(),
         'with_stiffness_eigenvalues':ev1.tolist()}

if __name__=='__main__':
 import json
 results=gastrulation();epithelial();early();late();results.update(geometry())
 print(json.dumps(results,indent=2))
