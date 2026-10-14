#!/usr/bin/env python3
"""Illustrate fast-slow frequency tuning in the note's R8 model.

This is a mechanism illustration, not a reproduction of Tsai et al. (2008).
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.signal import find_peaks
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures'
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
a0,b0,d0=.4,4.,2.
H=lambda x:x**4/(1+x**4)
Hp=lambda x:4*x**3/(1+x**4)**2
stationary=lambda x:(a0+b0*H(x))/(1+(d0*H(x))**4)-x
xeq=brentq(stationary,.5,2);yeq=d0*H(xeq)
clock=np.linspace(0,1200,24001)
results=[]
for eps in [.015,.025,.04,.065,.1]:
    sol=solve_ivp(lambda t,v:[(a0+b0*H(v[0]))/(1+v[1]**4)-v[0],eps*(d0*H(v[0])-v[1])],(0,1200),[1.1*xeq,.9*yeq],t_eval=clock,rtol=1e-9,atol=1e-11)
    if not sol.success:raise RuntimeError(sol.message)
    assert np.all(sol.y>=0)
    late=clock>600
    peaks,_=find_peaks(sol.y[0,late],prominence=.5)
    pt=clock[late][peaks]
    assert len(pt)>=4
    intervals=np.diff(pt[-4:]);assert np.std(intervals)/np.mean(intervals)<.01
    period=np.mean(intervals)
    cycle_mask=(clock>=pt[-2])&(clock<=pt[-1])
    amp=np.ptp(sol.y[0,cycle_mask]);assert amp>1
    results.append(dict(eps=eps,sol=sol,period=period,amplitude=amp,frequency=1/period))

fig,axes=plt.subplots(1,3,figsize=(10.7,3.5))
# Fast-equilibrium branches are independent of epsilon.
xx=np.linspace(.035,4.399,6000)
ratio=(a0+b0*H(xx))/xx-1
valid=ratio>=0;xx=xx[valid];yy=ratio[valid]**.25
fast_slope=b0*Hp(xx)/(1+yy**4)-1
stable=fast_slope<0
axes[0].plot(np.where(stable,yy,np.nan),np.where(stable,xx,np.nan),color='C0',lw=1.4,label='stable fast branches')
axes[0].plot(np.where(~stable,yy,np.nan),np.where(~stable,xx,np.nan),'--',color='.5',lw=1.2,label='unstable branch')
example=results[1];sol=example['sol']
mask=(clock>=1000)&(clock<=1000+example['period'])
axes[0].plot(sol.y[1,mask],sol.y[0,mask],color='C1',lw=1.3,label='late numerical orbit')
axes[0].set(xlim=(.5,1.35),ylim=(0,4.4),xlabel='slow inhibitor y',ylabel='fast variable x',title='Switch geometry and orbit')
axes[0].legend(fontsize=6)
for item in [results[1],results[-1]]:
    axes[1].plot(clock,item['sol'].y[0],label=f"epsilon={item['eps']:g}",lw=1.1)
axes[1].set(xlim=(800,1000),xlabel='time',ylabel='x',title='Different periods, similar excursions')
axes[1].legend(fontsize=7)
eps=np.array([v['eps'] for v in results]);frequencies=np.array([v['frequency'] for v in results]);amplitudes=np.array([v['amplitude'] for v in results])
axes[2].plot(eps,frequencies/frequencies[0],'o-',label='frequency / reference')
axes[2].plot(eps,amplitudes/amplitudes[0],'s-',label='amplitude / reference')
axes[2].set(xlabel='epsilon',ylabel='relative change',title='Tuning one slow time scale')
axes[2].legend(fontsize=7)
fig.tight_layout();fig.savefig(OUT/'positive_feedback_tuning.png',dpi=180);plt.close(fig)
print('Fast-slow example: positive trajectories and repeated late cycles verified')
assert np.all(np.diff(frequencies)>0)
assert np.ptp(amplitudes)/amplitudes[0]<.15
print('Frequency increases:',np.all(np.diff(frequencies)>0))
print('Relative amplitude span:',np.ptp(amplitudes)/amplitudes[0])
