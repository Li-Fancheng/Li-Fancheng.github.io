"""Analytical model comparisons used in the C. elegans review series.

These are dimensionless calculations from equations stated in the notes,
not fitted embryo data. Run with NumPy and Matplotlib.
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np

OUT=Path(__file__).resolve().parent
C=['#287f83','#b2576f','#b88929','#6677aa']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.titlesize':14,'axes.labelsize':12,'xtick.labelsize':11,'ytick.labelsize':11,'legend.fontsize':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'path','savefig.facecolor':'white','figure.facecolor':'white','text.color':'#243746','axes.labelcolor':'#243746'})

def save(fig,name):
 fig.savefig(OUT/(name+'.svg'),bbox_inches='tight',pad_inches=.16)
 fig.savefig(OUT/(name+'.png'),dpi=160,bbox_inches='tight',pad_inches=.16)
 plt.close(fig)

def clutch():
 fig,ax=plt.subplots(1,2,figsize=(11.8,4.5),layout='constrained')
 b=np.linspace(0,8,401);vj=b/(1+1.5*b);vc=1-.5*vj
 for y,l,c in [(vc,'Cortex',C[1]),(vj,'Junction',C[0]),(vc-vj,'Relative slip',C[2])]:ax[0].plot(b,y,label=l,c=c,lw=2.5)
 ax[0].set(xlabel=r'Coupling / drag $\xi/\zeta$',ylabel=r'Speed / $v_0$',ylim=(-.02,1.08),title='A   Motor load and boundary motion');ax[0].legend(frameon=False)
 t=np.linspace(0,8,1601);v=.085*(1+.7*np.sin(2*np.pi*t));travel=np.r_[0,np.cumsum((v[1:]+v[:-1])/2*np.diff(t))]
 for beta,c in zip([.1,1,6],C):
  radius=1-beta/(1+beta)*travel
  ax[1].plot(t,radius**2,c=c,lw=2.5,label=rf'$\xi/\zeta={beta}$')
 ax[1].set(xlabel='Time',ylabel=r'Apical area / $A_0$',title='B   Same imposed cortical pulses',ylim=(0,1.05));ax[1].legend(frameon=False)
 save(fig,'08_clutch_model')

def closure():
 fig,ax=plt.subplots(1,2,figsize=(11.8,4.5),layout='constrained');t=np.linspace(0,1,301)
 ax[0].plot(t,1-t,c=C[0],lw=2.5,label='Drag distributed along edge');ax[0].plot(t,(1-t)**2,c=C[1],lw=2.5,label='Fixed total drag')
 ax[0].set(xlabel=r'Time / $T_c$',ylabel=r'Opening area / $A_0$',title='A   Same initial size and closing time');ax[0].legend(frameon=False)
 r=np.linspace(.3,2.4,400);speed=1-1/r
 ax[1].axhline(0,c='#9aa7b1',lw=1);ax[1].axvline(1,c='#9aa7b1',ls='--',lw=1);ax[1].plot(r,speed,c=C[0],lw=2.5)
 ax[1].text(.38,-1.8,'Closes',color=C[0]);ax[1].text(1.6,.32,'Expands',color=C[1]);ax[1].set(xlabel=r'Radius / $r_c$',ylabel=r'$\dot r$ (scaled)',title='B   Outward load creates an unstable radius')
 save(fig,'09_closure_model')

def anisotropy():
 fig,ax=plt.subplots(1,2,figsize=(12.2,4.9),layout='constrained');beta=.35
 k=np.linspace(1,7,241);r=np.linspace(.3,4,241);kk,rr=np.meshgrid(k,r);ez=1-beta*rr;et=rr/kk-beta
 cat=np.where((ez>0)&(et<0),0,np.where(ez>0,1,2))
 ax[0].pcolormesh(kk,rr,cat,cmap=ListedColormap(['#6bb5a0','#e9e0ba','#ca93a3']),vmin=-.5,vmax=2.5,shading='nearest',rasterized=True)
 ax[0].plot(k,beta*k,c='#243746',lw=1.5);ax[0].axhline(1/beta,c='#243746',ls='--',lw=1.3)
 ax[0].set(xlabel=r'Stiffness ratio $k=E_\theta/E_z$',ylabel=r'Stress ratio $r=\sigma_\theta/\sigma_z$',title=r'A   Orthotropic response, $\beta=0.35$')
 ax[0].legend(handles=[Patch(color='#6bb5a0',label='Longer and thinner'),Patch(color='#e9e0ba',label='Both directions extend'),Patch(color='#ca93a3',label='Axial shortening')],loc='upper left',frameon=True,fontsize=9)
 ax[1].plot(k,np.full_like(k,1-2*beta),lw=2.5,c=C[0],label=r'Axial strain $\varepsilon_z$');ax[1].plot(k,2/k-beta,lw=2.5,c=C[1],label=r'Hoop strain $\varepsilon_\theta$');ax[1].axhline(0,c='#9aa7b1',lw=1)
 ax[1].set(xlabel=r'Stiffness ratio $E_\theta/E_z$',ylabel=r'Strain / $(\sigma_z/E_z)$',title=r'B   Fixed stress ratio $r=2$');ax[1].legend(frameon=False)
 assert np.all(1/k-beta**2>0)
 save(fig,'10_anisotropy_model')

def plasticity():
 fig,ax=plt.subplots(1,3,figsize=(14,4.5),layout='constrained');t=np.linspace(0,6,6001);sigma=np.cos(2*np.pi*t);rate=np.maximum(sigma-.4,0);ep=np.r_[0,np.cumsum((rate[1:]+rate[:-1])/2*np.diff(t))]
 ax[0].plot(t,sigma,c='#9aa7b1',lw=1.5,label='Elastic response');ax[0].plot(t,sigma+ep,c=C[0],lw=2,label='With reference update');ax[0].plot(t,ep,c=C[2],lw=2,label='Reference strain')
 ax[0].set(xlabel='Time / period',ylabel='Strain',title='A   Directional memory');ax[0].legend(frameon=False,fontsize=9)
 q=np.logspace(-2,2,401);x=1/(1+q);J=q/(1+q)
 ax[1].semilogx(q,x,c=C[0],lw=2.5);ax[1].set(xlabel=r'Severing / rebundling $k_s/k_b$',ylabel=r'Bundled fraction $x_*$',title='B   Material available to carry load',ylim=(-.02,1.02))
 ax[2].semilogx(q,J,c=C[1],lw=2.5);ax[2].set(xlabel=r'Severing / rebundling $k_s/k_b$',ylabel=r'Remodeling flux $J_*/k_b$',title='C   Opportunities for renewal',ylim=(-.02,1.02))
 save(fig,'11_plasticity_model')

def frap():
 fig,ax=plt.subplots(1,2,figsize=(11.8,4.5),layout='constrained');t=np.linspace(0,8,501)
 for fm,tau,c in [(.9,1,C[0]),(.5,1,C[1]),(.9,3,C[2])]:ax[0].plot(t,fm*(1-np.exp(-t/tau)),c=c,lw=2.5,label=rf'$F_m={fm},\ \tau={tau}$')
 ax[0].set(xlabel='Time',ylabel='Recovered fluorescence',title='A   Exchange time and recoverable fraction',ylim=(0,1));ax[0].legend(frameon=False)
 for g,c in zip([0,.3,.8],C):ax[1].plot(t,(1-np.exp(-(1+g)*t))/(1+g),c=c,lw=2.5,label=rf'Growth rate $g={g}$')
 ax[1].set(xlabel='Time',ylabel='Fluorescence per unit length',title=r'B   Geometric dilution, $s=k=1$',ylim=(0,1.03));ax[1].legend(frameon=False)
 save(fig,'12_frap_model')

def rod():
 fig,ax=plt.subplots(1,3,figsize=(13.5,4.4),layout='constrained')
 for theta,c in zip([.2,1.2,2.4,4.],C):
  u=np.linspace(-theta/2,theta/2,301);R=1/theta
  ax[0].plot(R*np.sin(u),R*(np.cos(u)-np.cos(theta/2)),c=c,lw=2.5,label=rf'$\theta={theta}$')
 ax[0].set(xlabel='Position / arc length',ylabel='Position / arc length',title='A   All arcs have length 1');ax[0].set_aspect('equal');ax[0].legend(frameon=False,fontsize=9)
 theta=np.linspace(.01,5,301);ax[1].plot(theta,2*np.sin(theta/2)/theta,c=C[0],lw=2.5);ax[1].set(xlabel=r'Arc angle $\theta$ (rad)',ylabel='End distance / arc length',title='B   Chord and arc diverge',ylim=(0,1.05))
 f=np.linspace(0,2,401);pos=np.sqrt(np.maximum(f-1,0));ax[2].plot(f[f<=1],pos[f<=1],c=C[0],lw=2.5)
 ax[2].plot(f[f>=1],pos[f>=1],c=C[0],lw=2.5);ax[2].plot(f[f>=1],-pos[f>=1],c=C[0],lw=2.5);ax[2].plot([1,2],[0,0],c=C[1],ls='--',lw=1.8)
 ax[2].set(xlabel=r'Compressive load $F/F_c$',ylabel='Scaled bending amplitude',title='C   A local buckling normal form')
 save(fig,'13_rod_model')

def identifiability():
 fig,ax=plt.subplots(1,2,figsize=(11.8,4.5),layout='constrained');t=np.linspace(0,8,401)
 for c,color,style in zip([1,2,4],C,['-',(0,(6,4)),(0,(1,3))]):
  F,K,eta=c,c,2*c;u=F/K*(1-np.exp(-K*t/eta));ax[0].plot(t,u,c=color,ls=style,lw=2.8,label=rf'$(F,K,\eta)=({F},{K},{eta})$')
 ax[0].set(xlabel='Time',ylabel='Released displacement',title='A   Three parameter sets, one curve');ax[0].legend(frameon=False)
 x=np.arange(3)
 for i,c in enumerate([1,2,4]):ax[1].bar(x+(i-1)*.23,[c,c,c],width=.21,color=C[i],label=rf'Scale $c={c}$')
 ax[1].set(xticks=x,xticklabels=[r'$F/F_0$',r'$K/K_0$',r'$\eta/\eta_0$'],ylabel='Parameter / baseline',title='B   Absolute scales remain different');ax[1].legend(frameon=False)
 save(fig,'14_identifiability_model')

if __name__=='__main__':
 for f in [clutch,closure,anisotropy,plasticity,frap,rod,identifiability]:f()
 print('Generated 7 model comparisons as SVG and PNG.')
