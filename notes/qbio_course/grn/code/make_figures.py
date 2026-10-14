#!/usr/bin/env python3
"""Reproduce the five numerical figures in the GRN note. Requires numpy, scipy, matplotlib."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.stats import poisson
OUT=Path(__file__).resolve().parents[1]/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
# 1. single-gene activation and linear loss, illustrating 3 fixed points
x=np.linspace(0,4,1200); a,b,k,n,g=.05,3.5,1.,4,1.
prod=a+b*x**n/(k**n+x**n); loss=g*x
fig,ax=plt.subplots(figsize=(5.7,3.5));ax.plot(x,prod,label='production');ax.plot(x,loss,label='degradation')
f=lambda z:a+b*z**n/(1+z**n)-g*z
roots=[]
for l,r in zip(x[:-1],x[1:]):
 if f(l)*f(r)<0:roots.append(brentq(f,l,r))
for z in roots:
 slope=b*n*z**(n-1)/(1+z**n)**2-g
 ax.plot(z,z,'ko',ms=6,mfc='black' if slope<0 else 'white')
ax.set(xlabel='protein concentration x',ylabel='rate',title='Autoregulation: intersections and bistability',xlim=(0,4));ax.legend();fig.savefig(OUT/'self_activation.png',dpi=180);plt.close(fig)
# 2 toggle nullclines
alpha,n=3.0,3
xx=np.linspace(0,3.2,1200)
fig,ax=plt.subplots(figsize=(5.3,4.2));ax.plot(xx,alpha/(1+xx**n),label='y-nullcline: y=g(x)');ax.plot(alpha/(1+xx**n),xx,label='x-nullcline: x=f(y)');ax.set(xlim=(0,3.2),ylim=(0,3.2),xlabel='x',ylabel='y',title='Mutual repression: nullclines');ax.legend(loc='upper right',fontsize=8)
sym=brentq(lambda q:q+q**4-3,0,3)
sum_xy=brentq(lambda q:q-q**(-2)-3,3,4)
low=(sum_xy-np.sqrt(sum_xy**2-4/sum_xy))/2;high=sum_xy-low
toggle_points=[(sym,sym),(low,high),(high,low)]
for px,py in toggle_points:
 jac=np.array([[-1,-9*py**2/(1+py**3)**2],[-9*px**2/(1+px**3)**2,-1]])
 eig=np.linalg.eigvals(jac); stable=np.max(eig)<0
 ax.plot(px,py,'ko',ms=6,mfc='black' if stable else 'white')
 print('Toggle equilibrium / eigenvalues:',(px,py),eig)
# Flow direction and trajectories distinguish basins from the nullclines.
mesh=np.linspace(.12,3.1,15);gx,gy=np.meshgrid(mesh,mesh)
vx=alpha/(1+gy**n)-gx;vy=alpha/(1+gx**n)-gy
speed=np.hypot(vx,vy);ax.quiver(gx,gy,vx/np.maximum(speed,1e-12),vy/np.maximum(speed,1e-12),color='.75',alpha=.6,width=.002,scale=35,zorder=0)
for initial in [[.5,.8],[.8,.5]]:
 trajectory=solve_ivp(lambda t,v:[alpha/(1+v[1]**n)-v[0],alpha/(1+v[0]**n)-v[1]],(0,15),initial,rtol=1e-9,atol=1e-11)
 ax.plot(*trajectory.y,color='darkgreen',lw=1,alpha=.8)
fig.savefig(OUT/'toggle_nullclines.png',dpi=180);plt.close(fig)
# 3 three-node repressilator sample trajectories
n=3; a=20; gamma=1
fun=lambda t,y: [a/(1+y[2]**n)-y[0],a/(1+y[0]**n)-y[1],a/(1+y[1]**n)-y[2]]
t=np.linspace(0,80,1800);sol=solve_ivp(fun,(0,80),[.1,1.2,3.0],t_eval=t,rtol=1e-9,atol=1e-11)
fig,ax=plt.subplots(figsize=(6.1,3.4))
for i,label in enumerate(['X1','X2','X3']):ax.plot(sol.t,sol.y[i],label=label,lw=1.6)
ax.set(xlabel='time',ylabel='concentration',title='Three-gene negative-feedback ring');ax.legend(ncol=3);fig.savefig(OUT/'repressilator.png',dpi=180);plt.close(fig)
# 4 birth death Poisson distributions
fig,ax=plt.subplots(figsize=(5.7,3.3));ns=np.arange(0,34)
for mu in [2,8,17]:ax.plot(ns,poisson.pmf(ns,mu),'o-',ms=2.7,label=f'mean={mu}')
ax.set(xlabel='molecule count n',ylabel='stationary probability',title='Birth-death reaction: Poisson fluctuations');ax.legend();fig.savefig(OUT/'birth_death.png',dpi=180);plt.close(fig)
print('Activation fixed points:',roots)

# 5. Incoherent feed-forward step responses and steady input-output curves.
fig,axes=plt.subplots(1,2,figsize=(8.4,3.3))
time=np.linspace(0,16,1601)
for ideal,label,color in [(True,'H=u/z: exact adaptation','C0'),(False,'H=u/(1+z): shifted baseline','C1')]:
 baseline=1. if ideal else .5
 def rhs(t,state):
  output,z=state
  return [(2/(z if ideal else 1+z)-output)/.25,(2-z)/2]
 sol=solve_ivp(rhs,(0,16),[baseline,1.],t_eval=time,rtol=1e-9,atol=1e-11)
 if not sol.success:raise RuntimeError(sol.message)
 axes[0].plot(time,sol.y[0],label=label,color=color)
 axes[0].hlines(baseline,0,16,colors=color,linestyles=':',alpha=.5)
inputs=np.linspace(.1,4,200)
axes[1].plot(inputs,np.ones_like(inputs),label='H=u/z')
axes[1].plot(inputs,inputs/(1+inputs),label='H=u/(1+z)')
axes[0].set(xlabel='time after step u: 1 -> 2',ylabel='output y',title='Transient responses')
axes[1].set(xlabel='constant input u',ylabel='steady output',title='Adaptation requires a functional constraint')
for ax in axes:ax.legend(fontsize=7)
fig.tight_layout();fig.savefig(OUT/'feedforward_adaptation.png',dpi=180);plt.close(fig)
# Numerical consistency checks for the stated demonstration, not a proof of cycles.
assert len(roots)==3
for z in roots:
 slope=3.5*4*z**3/(1+z**4)**2-1
 print('Activation equilibrium / slope:',z,slope)
assert low>0 and high>low
assert np.allclose([low+low*high**3,high+high*low**3],[3,3])
print('Hill ring threshold alpha_c:',3*2**(1/3))

print('Generated figures:',sorted(p.name for p in OUT.iterdir()))
