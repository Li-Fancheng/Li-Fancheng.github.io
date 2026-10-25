#!/usr/bin/env python3
"""A single-mode invagination comparison; no cell-resolved 3D solver or fitted data.

Run after generate_figures.py and generate_depth_figures.py. CSV files contain every
plotted trajectory. A standalone copy writes figures and CSV files next to itself.
Only NumPy and Matplotlib are required. RK4 convergence, common equilibrium and hydraulic stability are
checked before output. The prescribed spatial mode is shared by construction.
"""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'morphogenesis/figures'
DATA = ROOT/'mechanism_comparison'
COLORS = ['#2563eb', '#dc6b2f', '#16856b']
DT = .002
PARAMS = {'tau_p': 4., 'kappa': 2., 'mu': 1., 'initial_x': .03,
          'pulse_force': -.25, 'pulse_on': 2., 'pulse_off': 4.,
          'mode_height_over_radius': .12, 'dt': DT}


def rhs(y, model, drive=1., external=0., mu=1.):
    x, p = y
    if model == 'active':
        return np.array([drive-x+external, 0.])
    if model == 'growth':
        return np.array([mu*x-x**3+external, 0.])
    v = p-x+external
    return np.array([v, (drive-p-PARAMS['kappa']*v)/PARAMS['tau_p']])


def integrate(model, initial, stop, dt=DT, protocol='formation'):
    n = round(stop/dt)
    times = np.arange(n+1)*dt
    ys = np.empty((n+1, 2)); ys[0] = initial
    # Discontinuous protocols use exactly aligned step boundaries. Stage values
    # use the interval midpoint so that no RK4 stage samples the next interval.
    for j in range(n):
        middle = times[j]+dt/2
        ext = PARAMS['pulse_force'] if protocol == 'probe' and 2 <= middle < 4 else 0.
        drive = 0. if protocol == 'off' else 1.
        mu = -1. if protocol == 'release' else 1.
        f = lambda y: rhs(y, model, drive=drive, external=ext, mu=mu)
        y = ys[j]; k1 = f(y); k2 = f(y+dt*k1/2); k3 = f(y+dt*k2/2); k4 = f(y+dt*k3)
        ys[j+1] = y+dt*(k1+2*k2+2*k3+k4)/6
    return times, ys


def export(name, t, series):
    with (DATA/(name+'.csv')).open('w', newline='') as f:
        writer=csv.writer(f);writer.writerow(['time']+list(series))
        writer.writerows(zip(t, *series.values()))


def main():
    DATA.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.fonttype': 'none', 'svg.hashsalt': 'cell3d-mechanism-comparison'})
    models=['active', 'growth', 'hydraulic']
    labels=['Active moment', 'Constrained growth', 'Hydraulic drive']
    formation={};probe={};off={};convergence={}
    for model in models:
        t, ys=integrate(model, [.03,.03], 40)
        formation[model]=ys[:,0]
        tp, yp=integrate(model, [1.,1.], 18, protocol='probe');probe[model]=yp[:,0]
        _, yh=integrate(model, [1.,1.], 18, dt=DT/2, protocol='probe')
        error=float(np.max(np.abs(yh[::2]-yp)));convergence[model]=error
        assert error < 1e-9
        assert np.max(np.abs(rhs(np.array([1.,1.]),model))) < 1e-13
        to, yo=integrate(model,[1.,1.],18,protocol='off');off[model]=yo[:,0]
    _, release=integrate('growth',[1.,1.],18,protocol='release')
    off['growth_load_released']=release[:,0]
    matrix=np.array([[-1,1],[.5,-.75]])
    eig=np.linalg.eigvals(matrix)
    assert np.max(eig) < 0
    # Growth arrest leaves the accumulated compressive load in this model.
    assert np.max(np.abs(off['growth']-1)) < 1e-13
    assert abs(formation['hydraulic'][-1]-1) < .003
    export('formation',t,formation);export('common_probe',tp,probe);export('targeted_interventions',to,off)
    report={'model':'prescribed axisymmetric mode, dimensionless ODEs',
            'parameters':PARAMS,'common_equilibrium_x':1.,
            'rk4_step_halving_max_error':convergence,
            'hydraulic_eigenvalues':eig.tolist(),
            'scope':'A projected mode comparison, not a cross-representation 3D tissue benchmark.'}
    (DATA/'parameters_and_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    fig,axs=plt.subplots(2,2,figsize=(11.5,8.2),layout='constrained')
    a,b,c,d=axs.flat
    r=np.linspace(-1,1,400);phi=(1-r*r)**2
    for j,(col,label) in enumerate(zip(COLORS,labels)):
        center=(j-1)*2.55
        a.plot(r+center,-.12*phi,color=col,lw=2.1)
        a.plot(r+center,.035-.12*phi,color=col,lw=1.2)
        for rr in np.linspace(-.9,.9,13):
            z=-.12*(1-rr*rr)**2;a.plot([rr+center]*2,[z,z+.035],color=col,lw=.65)
        a.text(center,-.20,label,ha='center',fontsize=9,color=col)
    a.set(xlim=(-4,4),ylim=(-.24,.12));a.set_aspect(6)
    a.axis('off');a.set_title('A  The same prescribed final cross-section ($x=1$)',loc='left',fontsize=11)
    a.text(0,.085,'Spatial mode and location are given, not predicted',ha='center',fontsize=9,color='#475569')
    for model,label,col in zip(models,labels,COLORS):
        b.plot(t,formation[model],label=label,color=col,lw=2)
        c.plot(tp,probe[model],label=label,color=col,lw=2)
    b.set(title='B  Formation after a step in the independent input',xlabel=r'Time $t/\tau_m$',ylabel='Mode amplitude $x$',xlim=(0,30),ylim=(0,1.05));b.legend(frameon=False,fontsize=9,loc='lower right')
    c.axvspan(2,4,color='#cbd5e1',alpha=.45,zorder=-1)
    c.set(title='C  The same small flattening force in all three models',xlabel=r'Time $t/\tau_m$',ylabel='Mode amplitude $x$',ylim=(.72,1.06),xlim=(0,18))
    c.text(2.2,.735,r'$f_{ext}=-0.25$',fontsize=9)
    d.plot(to,off['active'],color=COLORS[0],label='Moment switched off',lw=2)
    d.plot(to,off['growth'],color=COLORS[1],label='Growth stopped; load retained',lw=2)
    d.plot(to,off['growth_load_released'],color=COLORS[1],ls='--',label='Growth load mechanically released',lw=1.8)
    d.plot(to,off['hydraulic'],color=COLORS[2],label='Hydraulic source switched off',lw=2)
    d.set(title='D  Distinct interventions on the proposed source',xlabel=r'Time after intervention $t/\tau_m$',ylabel='Mode amplitude $x$',ylim=(-.03,1.08),xlim=(0,18));d.legend(frameon=False,fontsize=8.5,loc='upper right')
    for ax in [b,c,d]:ax.grid(alpha=.18)
    fig.savefig(OUT/'mechanism_comparison.svg',metadata={'Date':None})
    fig.savefig(OUT/'mechanism_comparison.png',dpi=300)
    plt.close(fig)
    policy_path=ROOT.parents[1]/'data/cell_simulation_3d_figures.json'
    if policy_path.is_file():
        policy=json.loads(policy_path.read_text());policy['morphogenesis/figures/mechanism_comparison']={'width':1000,'kind':'calculation','intrinsic_width':1104,'intrinsic_height':787}
        policy_path.write_text(json.dumps(policy,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__':main()
