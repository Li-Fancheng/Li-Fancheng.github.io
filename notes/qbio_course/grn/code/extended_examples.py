#!/usr/bin/env python3
"""Reproduce illustrative two-node oscillation and adaptation mechanism figures.
These models explain mechanisms; they do not reproduce Ma et al.'s full scan.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures'
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
H=lambda x:x**4/(1+x**4)
Hp=lambda x:4*x**3/(1+x**4)**2
# The same R8 topology has a stable equilibrium or sustained numerical oscillation.
a0,b0,d0=.4,4.,2.
xs=np.linspace(0,a0+b0,8000)
f=lambda x:(a0+b0*H(x))/(1+(d0*H(x))**4)-x
roots=[brentq(f,l,r) for l,r in zip(xs[:-1],xs[1:]) if f(l)*f(r)<0]
assert len(roots)==1
x0=roots[0];y0=d0*H(x0)
fig,axes=plt.subplots(1,2,figsize=(8.5,3.4))
for eps,ax in zip([1.,.05],axes):
    jac=np.array([[b0*Hp(x0)/(1+y0**4)-1,-4*y0**3*(a0+b0*H(x0))/(1+y0**4)**2],[eps*d0*Hp(x0),-eps]])
    eigen=np.linalg.eigvals(jac)
    sol=solve_ivp(lambda t,v:[(a0+b0*H(v[0]))/(1+v[1]**4)-v[0],eps*(d0*H(v[0])-v[1])],(0,500),[1.1*x0,.9*y0],t_eval=np.linspace(0,500,10000),rtol=1e-9,atol=1e-11)
    if not sol.success:raise RuntimeError(sol.message)
    ax.plot(sol.t,sol.y[0],label='X');ax.plot(sol.t,sol.y[1],label='Y')
    ax.set(xlim=(200,500),xlabel='time',ylabel='concentration',title=f'epsilon={eps:g}: '+('recovery' if max(eigen.real)<0 else 'numerical oscillation'))
    ax.legend()
    if eps==1:assert max(eigen.real)<0 and np.ptp(sol.y[0,-2000:])<1e-6
    else:
        assert min(eigen.real)>0 and np.ptp(sol.y[0,-2000:])>1
        assert abs(np.ptp(sol.y[0,-2000:])-np.ptp(sol.y[0,-4000:-2000]))<1e-3
    print('R8 example: epsilon=',eps,'trace sign=',np.sign(np.trace(jac)),'det sign=',np.sign(np.linalg.det(jac)))
fig.tight_layout();fig.savefig(OUT/'two_node_oscillation.png',dpi=180);plt.close(fig)

# Simplified activation fractions and saturating-buffer approximation.
u0,u1=.5,.8;Cref=.4;kappa=.3
A0=u0/(1+u0);B0=A0*(1-Cref)/Cref
T=np.linspace(0,200,4001)
curves={}
for label,leak in [('buffer',0.),('leaky buffer',.08)]:
    def rhs(t,v):
        A,B,C=v
        return [3*(u1*(1-A)-A),kappa*(C-Cref)-leak*(B-B0),2*(A*(1-C)-B*C)]
    sol=solve_ivp(rhs,(0,200),[A0,B0,Cref],t_eval=T,rtol=1e-9,atol=1e-11)
    if not sol.success:raise RuntimeError(sol.message)
    assert np.min(sol.y)>0 and np.max(sol.y)<1
    if leak==0:assert abs(sol.y[2,-1]-Cref)<1e-6
    else:assert sol.y[2,-1]>Cref+.01
    curves[label]=(sol.y[2]/Cref,sol.y[1]/B0)

def ff_rhs(t,v):
    A,B,C=v
    return [3*(u1*(1-A)-A),.25*(A-B),2*(A*(1-C)-B*C)]
sol=solve_ivp(ff_rhs,(0,200),[A0,A0,.5],t_eval=T,rtol=1e-9,atol=1e-11)
assert sol.success and np.all(sol.y>0) and np.all(sol.y<1)
assert abs(sol.y[2,-1]-.5)<1e-6
curves['proportioner']=(sol.y[2]/.5,sol.y[1]/A0)
fig,axes=plt.subplots(1,3,figsize=(10.5,3.3))
for label,(C,B) in curves.items():
    axes[0].plot(T,C,label=label);axes[1].plot(T,B,label=label)
axes[0].axhline(1,color='.5',ls=':',lw=1)
axes[0].set(xlim=(0,80),xlabel='time after step',ylabel='C / prestimulus C',title='Response and return')
axes[1].set(xlim=(0,80),xlabel='time after step',ylabel='B / prestimulus B',title='Regulator changes its baseline')
inputs=np.linspace(.2,1.2,160)
axes[2].plot(inputs,np.ones_like(inputs),label='buffer / proportioner')
leak_vals=[]
for u in inputs:
    A=u/(1+u)
    C=brentq(lambda c:A*(1-c)-(B0+kappa/.08*(c-Cref))*c,.1,.8)
    leak_vals.append(C/Cref)
axes[2].plot(inputs,leak_vals,label='leaky buffer')
axes[2].set(xlabel='constant input u',ylabel='C* / reference C*',title='Steady input-output relationship')
for ax in axes:ax.legend(fontsize=7)
fig.tight_layout();fig.savefig(OUT/'adaptation_mechanisms.png',dpi=180);plt.close(fig)
print('Adaptation examples: positive states verified; ideal baselines recovered; leakage shifts baseline')
