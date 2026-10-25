"""Two falsifiable, dimensionless examples for the closure-to-elongation synthesis.

These are candidate constitutive models, not embryo fits or reconstructions.
Run: python generate_mechanism_synthesis.py (NumPy, SciPy, Matplotlib).
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

OUT=Path(__file__).resolve().parent
INK='#304857'; TEAL='#27858b'; ROSE='#b65376'; GOLD='#bd8b24'
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','font.family':'DejaVu Sans','savefig.facecolor':'white'})
report={'scope':'Dimensionless candidate models; no experimental data or fitted embryo parameters.'}

def save(fig,name):
    fig.savefig(OUT/(name+'.svg'),bbox_inches='tight')
    fig.savefig(OUT/(name+'.png'),dpi=170,bbox_inches='tight')
    plt.close(fig)

def mechanical_history():
    # A permanent reference-strain branch in parallel with a Maxwell branch.
    K=1.; eta=1.; Kh=.7; tau_s=2.; force=.03
    time=np.linspace(0,40,1601); cases=[]
    for q in [0.,.04]:
        def rhs(t,y):
            e,s=y; de=(force-K*(e-q)-s)/eta
            return [de,Kh*de-s/tau_s]
        sol=solve_ivp(rhs,[0,40],[0,K*q],t_eval=time,rtol=2e-10,atol=1e-12)
        assert sol.success
        e,s=sol.y
        equilibrium=q+force/K
        initial_total=-K*q+K*q
        assert abs(initial_total)<1e-14
        assert abs(rhs(0,[0,K*q])[0]-force/eta)<1e-14
        assert abs(e[-1]-equilibrium)<1e-4
        # Mechanical energy: elastic branch + Maxwell spring. Its balance is
        # P_input - eta*de^2 - s^2/(Kh*tau_s), with q held fixed.
        de=(force-K*(e-q)-s)/eta; ds=Kh*de-s/tau_s
        dW=K*(e-q)*de+(s/Kh)*ds
        rhs_energy=force*de-eta*de**2-s**2/(Kh*tau_s)
        assert np.max(np.abs(dW-rhs_energy))<1e-13
        cases.append((q,e,s,{'reference_strain':q,'initial_branch_stress':K*q,'initial_total_stress':initial_total,'initial_strain_rate':force/eta,'long_time_strain':equilibrium,'computed_final_strain':float(e[-1]),'max_energy_balance_error':float(np.max(np.abs(dW-rhs_energy)))}))
    fig,axs=plt.subplots(1,3,figsize=(13.3,4.2),layout='constrained')
    ax=axs[0];ax.axis('off');ax.set_xlim(0,1);ax.set_ylim(0,1)
    for y in [.68,.28]:
        ax.plot([.10,.22],[y,y],color=INK,lw=2);ax.plot([.78,.9],[y,y],color=INK,lw=2)
    ax.plot([.1,.1],[.16,.81],color=INK,lw=3);ax.plot([.9,.9],[.16,.81],color=INK,lw=3)
    ax.text(.5,.78,'Reference branch',ha='center',color=TEAL)
    ax.text(.5,.62,r'$K(e-q)$',ha='center',fontsize=15,color=TEAL)
    ax.text(.5,.36,'Relaxing branch',ha='center',color=ROSE)
    ax.text(.5,.22,r'$s;\quad \dot{s}=K_h\dot{e}-s/\tau_s$',ha='center',fontsize=12,color=ROSE)
    ax.text(.5,.04,'Same shape, net stress and initial rate;\ndifferent internal stresses and reference state.',ha='center',va='bottom',fontsize=10,color=INK)
    ax.set_title('A  Two load-bearing branches',loc='left',color=INK)
    for q,e,s,meta in cases:
        col=TEAL if q==0 else ROSE
        axs[1].plot(time,e,color=col,lw=2.3,label=fr'$q={q:g}$')
        axs[1].axhline(meta['long_time_strain'],color=col,ls=':',lw=1)
        axs[2].plot(time,s,color=col,lw=2.3,label=fr'$q={q:g}$')
    axs[1].set(title='B  Same subsequent drive',xlabel='Time after matched state',ylabel='Axial strain e')
    axs[1].legend(frameon=False);axs[1].set_ylim(-.002,.078)
    axs[2].set(title='C  Transient stress relaxes',xlabel='Time after matched state',ylabel='Relaxing branch stress s')
    for ax in axs[1:]:ax.grid(alpha=.16)
    save(fig,'40_closure_mechanical_history')
    report['mechanical_history']={'K':K,'eta':eta,'K_h':Kh,'tau_s':tau_s,'common_step_drive':force,'q_evolution':'held fixed during the short assay','time_interval':[0,40],'initialization':'e(0)=0, s(0)=K*q','cases':[c[3] for c in cases]}

def contact_capacity():
    # n0=0, n_infty=1, tau_m=1. Tension ramps linearly to rho.
    rho=.75; time=np.linspace(0,4,1601); n=1-np.exp(-time)
    protocols=[('Early rapid loading',0.,.2,ROSE),('Delayed gradual loading',1.5,1.,TEAL)]
    fig,axs=plt.subplots(1,2,figsize=(11.6,4.1),layout='constrained')
    axs[0].plot(time,n,color=INK,lw=2.3,label='Maturing capacity')
    cases=[]
    for label,delay,ramp,col in protocols:
        load=rho*np.clip((time-delay)/ramp,0,1)
        axs[0].plot(time,load,color=col,lw=2,ls='--' if delay==0 else '-',label=('Attempted rapid load' if delay==0 else label))
        margin_end=1-np.exp(-(delay+ramp))-rho
        during=(time>delay)&(time<=delay+ramp)
        margin=np.min(n[during]-load[during])
        assert (margin>=-1e-10)==(margin_end>=0)
        if margin_end<0:
            # n0=0 and a steeper initial load ramp cause immediate failure
            # in this threshold approximation; dashed load is the attempted load.
            assert rho/ramp>1
        cases.append({'label':label,'delay':delay,'ramp_duration':ramp,'end_of_ramp_capacity_margin':float(margin_end),'predicted_stable_contact':bool(margin_end>=0)})
    axs[0].set(title='A  Identical arrival and maturation',xlabel='Time since first contact / maturation time',ylabel='Load or capacity / maximum capacity',ylim=(0,1.05))
    axs[0].legend(frameon=False,fontsize=9,loc='lower right')
    x=np.linspace(0,4,501);frontier=1-np.exp(-x)
    axs[1].fill_between(x,0,frontier,color=TEAL,alpha=.14)
    axs[1].fill_between(x,frontier,1,color=ROSE,alpha=.1)
    axs[1].plot(x,frontier,color=INK,lw=2)
    for case,(_,delay,ramp,col) in zip(cases,protocols):axs[1].scatter(delay+ramp,rho,s=65,color=col,zorder=3)
    axs[1].text(2.05,.25,'Capacity sufficient',color=TEAL)
    axs[1].text(.25,.91,'Capacity exceeded',color=ROSE)
    axs[1].set(title='B  Testable capacity boundary',xlabel=r'$(t_d+\tau_L)/\tau_m$',ylabel=r'$F_{max}/(N f_b n_\infty)$',xlim=(0,4),ylim=(0,1))
    for ax in axs:ax.grid(alpha=.16)
    save(fig,'41_contact_maturation_loading')
    # Verify the endpoint criterion over independently sampled loading ramps.
    rng=np.random.default_rng(937);largest_error=0.
    for _ in range(150):
        delay=float(rng.uniform(0,3));ramp=float(rng.uniform(.02,3));r=float(rng.uniform(.03,.98))
        t=np.linspace(delay,delay+ramp,2001)
        m=1-np.exp(-t)-r*(t-delay)/ramp
        endpoint=min(1-np.exp(-delay),1-np.exp(-(delay+ramp))-r)
        largest_error=max(largest_error,abs(float(m.min())-endpoint))
    assert largest_error<1e-12
    report['contact_capacity']={'n_0':0,'n_infty':1,'tau_m':1,'rho':rho,'capacity':'N*f_b*n(t)','load':'F_max*clip((t-t_d)/tau_L,0,1)','cases':cases,'critical_available_time':float(-np.log(1-rho)),'random_protocols_checked':150,'max_endpoint_criterion_error':largest_error,'omitted_processes':['force-dependent unbinding','damage accumulation','force-dependent maturation','recovery after rupture']}

if __name__=='__main__':
    mechanical_history();contact_capacity()
    (OUT/'mechanism_synthesis_parameters.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
