#!/usr/bin/env python3
"""Reproduce the 20 numerical figures in these notes. Not a benchmark for large neural networks.

Run from anywhere: python experiments.py --note 9; python experiments.py --all
Keep this file beside series.json. Each experiment uses a local, fixed RNG seed.
Figures and numerical diagnostics are written beside the corresponding body.md.
"""
import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/neural-theory-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/neural-theory-cache")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_continuous_lyapunov
from scipy.special import erf, expit, gammainc, gamma

ROOT = Path(__file__).resolve().parent
SERIES = json.loads((ROOT / "series.json").read_text())
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.prop_cycle": plt.cycler(color=["#2563eb", "#e87924", "#059669", "#9333ea"]),
    "svg.hashsalt": "neural-theory-2026", "mathtext.fontset": "cm"})


def canvas():
    return plt.subplots(1, 2, figsize=(10.2, 3.7), layout="constrained")


def finish(n, fig, stats):
    out = ROOT / SERIES[n-1]["slug"] / "figures"
    out.mkdir(parents=True, exist_ok=True)
    for ax in fig.axes:
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8, frameon=False)
        ax.grid(alpha=.16)
    fig.savefig(out / "experiment.svg", metadata={"Date": None})
    fig.savefig(out / "experiment.png", dpi=220)
    plt.close(fig)
    (out / "diagnostics.json").write_text(json.dumps(stats, indent=2, allow_nan=False)+"\n")
    print(f"Note {n:02}: {SERIES[n-1]['slug']} — {stats}", flush=True)


def note01():
    t = np.geomspace(.01, 3000, 500); s = np.array([1., .1, .01])
    filt = -np.expm1(-t[:, None]*s)
    bias = np.sum(np.exp(-2*t[:, None]*s), axis=1)
    variance = np.sum(.25/100/s*filt**2, axis=1)
    fig, ax = canvas()
    for j in range(3): ax[0].semilogx(t, filt[:,j], label=f"s={s[j]}")
    ax[0].set(xlabel="Optimization time", ylabel="Signal recovery")
    for y, lab in [(bias,"Bias"),(variance,"Variance"),(bias+variance,"Excess risk")]:
        ax[1].loglog(t,y,label=lab)
    ax[1].set(xlabel="Optimization time", ylabel="Expected squared error")
    ax[1].set_ylim(1e-4,4)
    finish(1,fig,{"best_time":float(t[np.argmin(bias+variance)]), "late_variance":float(variance[-1])})


def note02():
    rng=np.random.default_rng(202602); D=80; noise=.25
    fig,ax=canvas()
    for g in [np.linspace(.1,.97,200),np.linspace(1.03,5,300)]:
        r=np.where(g<1,noise*g/(1-g),(1-1/g)+noise/(g-1))
        ax[0].plot(g,r,color="#2563eb")
    med=[]
    for d in [16,40,64,76,84,96,120,160,240,400]:
        risks=[]
        for _ in range(24):
            X=rng.normal(size=(D,d)); w=rng.normal(size=d); w/=np.linalg.norm(w)
            y=X@w+rng.normal(scale=np.sqrt(noise),size=D)
            wh=np.linalg.lstsq(X,y,rcond=None)[0]; risks.append(np.sum((wh-w)**2))
        med.append(float(np.median(risks)))
    ax[0].scatter(np.array([16,40,64,76,84,96,120,160,240,400])/D,med,s=20,color="#e87924",label="Finite-D median")
    ax[0].set(xlabel=r"$d/D$",ylabel="Excess risk",yscale="log",ylim=(.02,100))
    for g in [.3,1.,3.]:
        lo,hi=(1-np.sqrt(g))**2,(1+np.sqrt(g))**2
        s=np.linspace(max(lo,1e-4),hi,600)
        rho=np.sqrt(np.maximum((hi-s)*(s-lo),0))/(2*np.pi*g*s)
        ax[1].plot(s,rho,label=f"g={g}; zero mass={max(1-1/g,0):.2f}")
    ax[1].set(xlabel="Eigenvalue of X.T X / D",ylabel="MP continuous density",ylim=(0,2))
    finish(2,fig,{"D":D,"repeats":24,"median_risks":med})


def rbf(x,z,ell): return np.exp(-((x[:,None]-z[None,:])/ell)**2/2)


def note03():
    rng=np.random.default_rng(202603); x=np.linspace(-np.pi,np.pi,101); ell=.7
    K=rbf(x,x,ell); sizes=np.array([8,16,32,64,128,256,512,1024,2048]); means=[]; std=[]
    for n in sizes:
        errs=[]
        for _ in range(24):
            om=rng.normal(scale=1/ell,size=n); b=rng.uniform(0,2*np.pi,n)
            Z=np.sqrt(2/n)*np.cos(x[:,None]*om+b)
            errs.append(np.linalg.norm(Z@Z.T-K)/np.linalg.norm(K))
        means.append(np.mean(errs)); std.append(np.std(errs))
    xt=rng.uniform(-np.pi,np.pi,40); y=np.sin(xt)+.3*np.cos(3*xt)+rng.normal(scale=.1,size=40)
    om=rng.normal(scale=1/ell,size=256); b=rng.uniform(0,2*np.pi,256)
    Z=np.sqrt(2/256)*np.cos(xt[:,None]*om+b); Zs=np.sqrt(2/256)*np.cos(x[:,None]*om+b)
    pred=Zs@Z.T@np.linalg.solve(Z@Z.T+.04*np.eye(40),y)
    exact=rbf(x,xt,ell)@np.linalg.solve(rbf(xt,xt,ell)+.04*np.eye(40),y)
    fig,ax=canvas(); ax[0].errorbar(sizes,means,yerr=std,fmt="o-",capsize=2)
    ax[0].set(xscale="log",yscale="log",xlabel="Random features N",ylabel="Relative kernel error")
    ax[1].plot(x,np.sin(x)+.3*np.cos(3*x),label="Teacher",color="black")
    ax[1].plot(x,pred,label="N=256"); ax[1].plot(x,exact,"--",label="RBF kernel"); ax[1].scatter(xt,y,s=10,alpha=.5)
    ax[1].set(xlabel="x",ylabel="Prediction")
    finish(3,fig,{"kernel_error_N2048":float(means[-1]),"D":40,"lambda":.001})


def note04():
    rng=np.random.default_rng(202604); alphas=np.geomspace(.1,30,16); c=np.sqrt(2/np.pi); d=300
    vals=[]
    for alpha in alphas:
        rr=[]
        for _ in range(20):
            X=rng.normal(size=(max(1,round(alpha*d)),d)); y=np.sign(X[:,0]); v=y@X/len(y)
            rr.append(v[0]/np.linalg.norm(v))
        vals.append(np.mean(rr))
    fig,ax=canvas(); theory=c/np.sqrt(c*c+1/alphas)
    ax[0].semilogx(alphas,theory,label="High-d limit"); ax[0].scatter(alphas,vals,s=20,label="d=300; 20 repeats")
    ax[0].set(xlabel="D/d",ylabel="Direction overlap R")
    for p in [0,.2]:
        cp=(1-2*p)*c; r=cp/np.sqrt(cp*cp+1/alphas)
        ax[1].loglog(alphas,p+(1-2*p)*np.arccos(r)/np.pi,label=f"Label flip p={p}")
    ax[1].set(xlabel="D/d",ylabel="Test classification error")
    finish(4,fig,{"R_at_alpha1":float(c/np.sqrt(c*c+1)),"error_at_alpha1":float(np.arccos(c/np.sqrt(c*c+1))/np.pi)})


def note05():
    rng=np.random.default_rng(202605); z=rng.normal(size=(2,200000)); cs=np.linspace(-1,1,61)
    exact=lambda c:(np.sqrt(np.maximum(1-c*c,0))+(np.pi-np.arccos(c))*c)/(2*np.pi)
    mc=[np.mean(np.maximum(z[0],0)*np.maximum(c*z[0]+np.sqrt(max(1-c*c,0))*z[1],0)) for c in cs]
    fig,ax=canvas(); ax[0].plot(cs,exact(cs),label="Analytic"); ax[0].scatter(cs[::4],np.array(mc)[::4],s=16,label="Gaussian integral MC")
    ax[0].set(xlabel="Input correlation",ylabel="ReLU covariance")
    for c0 in [-.5,0,.5,.9]:
        c=c0; hist=[c]
        for _ in range(60): c=float(2*exact(np.clip(c,-1,1))); hist.append(c)
        ax[1].plot(hist,label=f"c0={c0}")
    ax[1].set(xlabel="Layer",ylabel="Correlation")
    finish(5,fig,{"max_MC_error":float(np.max(np.abs(np.array(mc)-exact(cs))))})


def note06():
    xt=np.linspace(-np.pi,np.pi,24); xs=np.linspace(-np.pi,np.pi,301); y=np.sin(xt)+.3*np.cos(3*xt)
    K=rbf(xt,xt,.8)/24; vals,U=np.linalg.eigh(K); proj=U.T@y
    fig,ax=canvas(); t=np.geomspace(.01,1e5,300)
    for j in [23,21,18,14]: ax[0].semilogx(t,np.exp(-vals[j]*t),label=f"Mode {24-j}; k={vals[j]:.2g}")
    ax[0].set(xlabel="Optimization time",ylabel="Residual fraction")
    ax[1].plot(xs,np.sin(xs)+.3*np.cos(3*xs),color="black",label="Teacher")
    for tt in [1,10,1000]:
        filt=np.divide(-np.expm1(-vals*tt),vals,out=np.full_like(vals,tt),where=vals>1e-14)
        pred=rbf(xs,xt,.8)/24@U@(filt*proj); ax[1].plot(xs,pred,label=f"t={tt}")
    ax[1].set(xlabel="x",ylabel="Prediction")
    finish(6,fig,{"largest_eigenvalue":float(vals[-1]),"D":24})


def sequence_risk(eta,a2,D,lam,sigma=.3):
    return np.sum((lam/(eta+lam))**2*a2+sigma**2/D*(eta/(eta+lam))**2)


def note07():
    k=np.arange(1,10001,dtype=float); eta=k**-2; Ds=np.geomspace(10,1e6,60); lams=np.geomspace(1e-10,1,200)
    fig,ax=canvas(); ax[0].loglog(k,eta,"--",color="black",label="Kernel eigenvalue")
    curves=[]
    for q in [2.,4.]:
        a2=k**-q; a2/=a2.sum(); ax[0].loglog(k,a2,label=f"Target q={q}")
        curve=[min(sequence_risk(eta,a2,D,l) for l in lams) for D in Ds]
        ax[1].loglog(Ds,curve,label=f"Target q={q}"); curves.append(float(curve[-1]))
    ax[0].set(xlabel="Mode i",ylabel="Spectral energy"); ax[1].set(xlabel="Sequence sample budget D",ylabel="Oracle-filter excess risk")
    finish(7,fig,{"last_risks":curves,"modes":10000,"sigma":.3})


def note08():
    rng=np.random.default_rng(202608); x=np.linspace(-1,1,201); y=np.tanh(3*x); n=128
    w0=rng.normal(scale=.4,size=n); a0=rng.normal(scale=.05,size=n); histories=[]; final=[]
    for train_w in [False,True]:
        w=w0.copy(); a=a0.copy(); hist=[]
        for step in range(5001):
            h=np.tanh(x[:,None]*w); r=h@a/n-y
            if step%20==0: hist.append([step*.05,np.mean(r*r)/2])
            da=h.T@r/len(x); dw=((1-h*h)*(r*x)[:,None]).mean(axis=0)*a
            a-=.05*da
            if train_w: w-=.05*dw
        histories.append(np.array(hist)); final.append((w,a))
    fig,ax=canvas()
    for h,lab in zip(histories,["Fixed hidden layer","Train a and w"]): ax[0].semilogy(h[:,0],h[:,1],label=lab)
    wf,af=final[1]; ax[1].plot(np.vstack([w0,wf]),np.vstack([a0,af]),color="#94a3b8",alpha=.3,lw=.6)
    ax[1].scatter(w0,a0,s=12,label="Initial"); ax[1].scatter(wf,af,s=12,label="Final")
    ax[0].set(xlabel="Mean-field time",ylabel="Population-grid loss"); ax[1].set(xlabel="Hidden weight w",ylabel="Readout a")
    finish(8,fig,{"fixed_loss":float(histories[0][-1,1]),"trained_loss":float(histories[1][-1,1]),"N":n})


def note09():
    fig,ax=canvas(); t=np.linspace(0,22,500); p0=1e-4; errs=[]
    for s in [.5,1,2]:
        p=s/(1+(s/p0-1)*np.exp(-2*s*t)); ax[0].plot(t,p/s,label=f"sigma={s}")
        sol=solve_ivp(lambda t,p:2*p*(s-p),[0,22],[p0],t_eval=t,rtol=1e-9,atol=1e-11)
        ax[0].scatter(t[::30],sol.y[0,::30]/s,s=12); errs.append(float(np.max(np.abs(sol.y[0]-p))))
    ts=np.geomspace(.001,1e4,600)
    for L in [1,2,3,4]:
        sol=solve_ivp(lambda t,p:L*np.maximum(p,0)**(2-2/L)*(1-p),[0,1e4],[1e-3],t_eval=ts,rtol=1e-8,atol=1e-10)
        ax[1].semilogx(ts,sol.y[0],label=f"L={L}")
    ax[0].set(xlabel="Time",ylabel="Fraction of target"); ax[1].set(xlabel="Time; equal effective p0",ylabel="Effective weight")
    finish(9,fig,{"max_logistic_ODE_error":max(errs),"half_time_sigma1":float(np.log(1/p0-1)/2)})


def note10():
    Z=np.array([[1.,0.],[1.,1.],[2.,-.2]]); t=np.geomspace(.01,1e7,600)
    # Integrate log time, avoiding millions of explicit Euler steps.
    sol=solve_ivp(lambda u,w:np.exp(u)*(Z.T@expit(-Z@w)),[np.log(t[0]),np.log(t[-1])],np.zeros(2),t_eval=np.log(t),rtol=1e-9,atol=1e-11)
    W=sol.y.T; losses=np.logaddexp(0,-W@Z.T).sum(axis=1); norms=np.linalg.norm(W,axis=1)
    angle=np.where(norms>1e-10,np.arctan2(np.abs(W[:,1]),W[:,0])*180/np.pi,np.nan)
    fig,ax=canvas(); ax[0].loglog(t,losses,label="Logistic loss"); ax[0].loglog(t,norms,label="Parameter norm")
    ax[1].semilogx(t,angle); ax[1].set(xlabel="Optimization time",ylabel="Angle to max margin (degrees)")
    ax[0].set(xlabel="Optimization time",ylabel="Loss / norm")
    finish(10,fig,{"last_angle_degrees":float(angle[-1]),"last_norm":float(norms[-1])})


def erf_corr(A,B,C): return 2/np.pi*np.arcsin(np.clip(C/np.sqrt((1+A)*(1+B)),-1,1))


def committee_error(Q,R,T):
    return .5*(erf_corr(np.diag(Q)[:,None],np.diag(Q)[None,:],Q).sum()
        -2*erf_corr(np.diag(Q)[:,None],np.diag(T)[None,:],R).sum()
        +erf_corr(np.diag(T)[:,None],np.diag(T)[None,:],T).sum())


def note11():
    rng=np.random.default_rng(202611); d=200; eta=.15
    B=np.linalg.qr(rng.normal(size=(d,2)))[0].T*np.sqrt(d)
    W=.25*(B[0]+B[1])[None,:]+rng.normal(scale=.08,size=(2,d)); history=[]; T=B@B.T/d
    for step in range(120001):
        if step%400==0:
            Q=W@W.T/d; R=W@B.T/d; S=(R[0,0]-R[0,1]-R[1,0]+R[1,1])/2
            history.append([step/d,*R.ravel(),committee_error(Q,R,T),abs(S)])
        x=rng.normal(size=d); h=W@x/np.sqrt(d); v=B@x/np.sqrt(d)
        residual=erf(v/np.sqrt(2)).sum()-erf(h/np.sqrt(2)).sum()
        F=residual*np.sqrt(2/np.pi)*np.exp(-h*h/2); W+=eta/np.sqrt(d)*F[:,None]*x
    h=np.array(history); fig,ax=canvas()
    for j,lab in enumerate(["R11","R12","R21","R22"]): ax[0].plot(h[:,0],h[:,j+1],label=lab)
    ax[1].semilogy(h[:,0],np.maximum(h[:,5],1e-16),label="Generalization loss")
    ax[1].plot(h[:,0],h[:,6],label="|Specialization|")
    ax[0].set(xlabel="Fresh examples / d",ylabel="Teacher overlap"); ax[1].set(xlabel="Fresh examples / d",ylabel="Loss / specialization")
    finish(11,fig,{"last_error":float(h[-1,5]),"last_specialization":float(h[-1,6]),"examples":120000})


def note12():
    from itertools import product
    X=np.array(list(product([-1,1],repeat=4))); z1=X[:,0]*X[:,1]; z2=X[:,2]*X[:,3]; y=z1*z2
    correlations=np.array([*(np.mean(y[:,None]*X,axis=0)),np.mean(y*z1),np.mean(y*z2),np.mean(y*z1*z2)])
    t=np.linspace(0,100,600); p0=.03**3
    sol=solve_ivp(lambda t,p:3*np.maximum(p,0)**(4/3)*(1-p),[0,100],[p0],t_eval=t,rtol=1e-9,atol=1e-12)
    a=1-(1-.03)*np.exp(-t)
    fig,ax=canvas(); ax[0].bar(["x1","x2","x3","x4","z1","z2","z1 z2"],correlations)
    ax[0].set(ylabel="Correlation with final label",ylim=(-.05,1.1))
    ax[1].semilogy(t,np.maximum((1-sol.y[0])**2/2,1e-15),label="Train three factors")
    ax[1].semilogy(t,np.maximum((1-a)**2/2,1e-15),label="b=c=1; train a")
    ax[1].set(xlabel="Optimization time",ylabel="Population loss")
    finish(12,fig,{"correlations":correlations.tolist(),"joint_final_loss":float((1-sol.y[0,-1])**2/2)})


def note13():
    t=np.geomspace(.01,2e4,700); q=100*np.exp(-.001*t)
    train=(np.exp(-2*t)+.001*np.exp(-.002*t))/2
    test=(np.exp(-2*t)+10000*np.exp(-.002*t))/2
    output_std=np.sqrt(2*test)
    acc=erf(.5/(np.sqrt(2)*output_std))
    fig,ax=canvas(); ax[0].loglog(t,train,label="Training risk"); ax[0].loglog(t,test,label="Test risk")
    ax[0].set(xlabel="Time",ylabel="Continuous squared risk",ylim=(1e-10,1e5))
    ax[1].semilogx(t,acc,label="Tolerance accuracy; delta=0.5")
    ax[1].set(xlabel="Time",ylabel="P(|prediction error| < delta)",ylim=(0,1.02))
    finish(13,fig,{"slow_rate":.001,"initial_slow_error":100.,"half_accuracy_time":float(np.log(100*.6744897501960817/.5)/.001)})


def note14():
    rng=np.random.default_rng(202614); H=np.diag([1.,4.]); A=np.array([[1.,.8],[.8,2.]])
    S=solve_continuous_lyapunov(H,2*A); dt=.002; n=5000; W=np.zeros((n,2)); L=np.linalg.cholesky(2*A*dt)
    for _ in range(3500): W-=dt*W@H.T; W+=rng.normal(size=W.shape)@L.T
    emp=np.cov(W.T); fig,ax=canvas()
    ax[0].scatter(W[:,0],W[:,1],s=2,alpha=.13); vals,V=np.linalg.eigh(S); theta=np.linspace(0,2*np.pi,250)
    ellipse=2*V@np.diag(np.sqrt(vals))@np.array([np.cos(theta),np.sin(theta)])
    ax[0].plot(*ellipse,color="#e87924",label="2-sigma covariance ellipse")
    ax[0].set(xlabel="Coordinate 1",ylabel="Coordinate 2",aspect="equal")
    ax[1].bar(np.arange(3)-.18,[S[0,0],S[1,1],S[0,1]],width=.36,label="Continuous OU")
    ax[1].bar(np.arange(3)+.18,[emp[0,0],emp[1,1],emp[0,1]],width=.36,label="Euler simulation")
    ax[1].set(xticks=np.arange(3),xticklabels=["S11","S22","S12"],ylabel="Stationary covariance")
    finish(14,fig,{"analytic_covariance":S.tolist(),"empirical_covariance":emp.tolist(),"relative_covariance_error":float(np.linalg.norm(emp-S)/np.linalg.norm(S))})


def note15():
    rng=np.random.default_rng(202615); n=400; P=12; xi=rng.choice([-1,1],size=(P,n))
    J=xi.T@xi/n; np.fill_diagonal(J,0); s=xi[0].copy(); flip=rng.choice(n,int(.2*n),replace=False); s[flip]*=-1
    energy=[]; overlap=[]
    for sweep in range(31):
        energy.append(float(-s@J@s/2)); overlap.append(float(s@xi[0]/n))
        for i in rng.permutation(n):
            h=J[i]@s
            if h!=0: s[i]=1 if h>0 else -1
    betas=np.linspace(.1,3,150); retrieval=[]
    for beta in betas:
        m=1.
        for _ in range(2000): m=np.tanh(beta*m)
        retrieval.append(m)
    fig,ax=canvas(); ax[0].plot(energy,label="Energy / arbitrary units")
    ax[0].set(xlabel="Asynchronous sweeps",ylabel="Energy")
    ax[1].plot(betas,retrieval,label="Low-load mean-field m")
    ax[1].set(xlabel="Inverse temperature beta",ylabel="Retrieval overlap")
    finish(15,fig,{"N":n,"P":P,"initial_overlap":overlap[0],"final_overlap":overlap[-1],"largest_energy_increase":float(np.max(np.diff(energy)))})


def note16():
    rng=np.random.default_rng(202616); n=400; dt=.05; J=rng.normal(size=(n,n))/np.sqrt(n)
    fig,ax=canvas(); summary={}
    for g in [.8,1.8]:
        x=rng.normal(scale=.1,size=n); tangent=rng.normal(size=n); tangent/=np.linalg.norm(tangent)
        states=[]; trace=[]; logs=[]
        for step in range(9000):
            ph=np.tanh(x); V=-tangent+g*J@((1-ph*ph)*tangent)
            x+=dt*(-x+g*J@ph); tangent+=dt*V
            norm=np.linalg.norm(tangent); tangent/=norm
            if step>=3000: logs.append(np.log(norm)); states.append(x.copy())
            if step%5==0: trace.append([step*dt,x[0]])
        trace=np.array(trace); A=np.array(states); ac=[]
        for lag in range(200): ac.append(np.mean(A[:len(A)-lag or None]*A[lag:]))
        scale=max(ac[0],1e-300); ax[0].plot(trace[:,0],trace[:,1],label=f"g={g}")
        ax[1].plot(np.arange(200)*dt,np.array(ac)/scale,label=f"g={g}")
        summary[str(g)]={"late_variance":float(np.mean(A*A)),"Euler_tangent_exponent":float(np.mean(logs)/dt)}
    ax[0].set(xlabel="Time",ylabel="One unit state"); ax[1].set(xlabel="Lag",ylabel="Raw correlation / C(0)")
    finish(16,fig,summary)


def note17():
    t=np.geomspace(.01,1e6,800); fig,ax=canvas()
    for nu in [.5,1.,2.]:
        # Integral from zero to one: 1/2 * int lambda^(nu-1) exp(-2 lambda t) dlambda.
        exact=.5*gamma(nu)*gammainc(nu,2*t)/(2*t)**nu
        ax[0].loglog(t,exact,label=f"nu={nu}")
    rates=np.geomspace(1e-3,1,1500); weighted=np.ones_like(rates)
    finite=np.trapezoid(weighted[None,:]*np.exp(-2*t[:,None]*rates),rates,axis=1)/2
    ideal=.5*(-np.expm1(-2*t))/(2*t)
    ax[1].loglog(t,finite,label="Lower cutoff 0.001"); ax[1].loglog(t,ideal,"--",label="No lower cutoff")
    ax[0].set(xlabel="Time",ylabel="Loss"); ax[1].set(xlabel="Time",ylabel="Loss",ylim=(1e-12,1))
    finish(17,fig,{"nu_values":[.5,1.,2.],"cutoff":.001})


def note18():
    k=np.arange(1,100001,dtype=float); p=2.; q=2.; eta=k**-p; a2=k**-q; a2/=a2.sum(); sigma=.3
    Ds=np.geomspace(10,1e7,70); lams=np.geomspace(1e-12,1,240); risk=[]; ls=[]
    for D in Ds:
        rr=np.array([sequence_risk(eta,a2,D,l,sigma) for l in lams]); idx=rr.argmin(); risk.append(rr[idx]); ls.append(lams[idx])
    risk=np.array(risk); fig,ax=canvas(); ax[0].loglog(Ds,risk,label="Exact sequence sum")
    ax[0].loglog(Ds,risk[20]*(Ds/Ds[20])**(-.5),"--",label="D^(-1/2)")
    ax[1].loglog(Ds,np.array(ls)**(-1/p),label="Effective cutoff lambda^(-1/p)")
    ax[1].loglog(Ds,np.sqrt(Ds),"--",label="D^(1/2), arbitrary prefactor")
    ax[0].set(xlabel="Sequence budget D",ylabel="Oracle excess risk"); ax[1].set(xlabel="Sequence budget D",ylabel="Effective modes")
    slope=np.polyfit(np.log(Ds[15:50]),np.log(risk[15:50]),1)[0]
    finish(18,fig,{"p":p,"q":q,"predicted_beta":.5,"measured_slope":float(slope),"modes":len(k)})


def note19():
    q=2.5; k=np.arange(1,100001,dtype=float); a2=k**-q; a2/=a2.sum(); tail=np.maximum(1-np.cumsum(a2),0)
    Ns=np.unique(np.geomspace(1,10000,200).astype(int)); sigma=.3; fig,ax=canvas(); best={}
    ax[0].loglog(Ns,tail[Ns-1],label="Representation tail")
    for D in [1e3,1e4,1e5]:
        error=tail[Ns-1]+sigma**2*Ns/D; ax[1].loglog(Ns,error,label=f"D={D:g}")
        j=error.argmin(); ax[1].scatter([Ns[j]],[error[j]],s=18); best[str(D)]=int(Ns[j])
    ax[0].set(xlabel="Retained orthogonal modes N",ylabel="Approximation error")
    ax[1].set(xlabel="Retained orthogonal modes N",ylabel="Total sequence-model risk")
    finish(19,fig,{"target_q":q,"best_ranks_on_grid":best,"sigma":sigma})


def note20():
    alpha=.4; beta=.3; A=2.; B=3.; k=6.; C=np.geomspace(1e5,1e15,100)
    N=(alpha*A/(beta*B))**(1/(alpha+beta))*(C/k)**(beta/(alpha+beta)); D=C/(k*N)
    fig,ax=canvas(); ax[0].loglog(C,N,label="Optimal N"); ax[0].loglog(C,D,label="Optimal D")
    ratios=np.geomspace(.01,100,301); excess=[]
    for cc in [1e7,1e10,1e13]:
        nn=(alpha*A/(beta*B))**(1/(alpha+beta))*(cc/k)**(beta/(alpha+beta)); dd=cc/(k*nn)
        loss=A*(nn*ratios)**(-alpha)+B*(dd/ratios)**(-beta); opt=A*nn**(-alpha)+B*dd**(-beta)
        ax[1].semilogx(ratios,loss/opt,label=f"C={cc:g}"); excess.append(float(loss.min()/opt))
    ax[0].set(xlabel="Compute C",ylabel="Resource size (model units)")
    ax[1].set(xlabel="N / N*; D adjusts to fixed C",ylabel="Excess loss / optimum")
    balance=alpha*A*N**(-alpha)/(beta*B*D**(-beta))
    finish(20,fig,{"alpha":alpha,"beta":beta,"N_exponent":beta/(alpha+beta),"D_exponent":alpha/(alpha+beta),"max_balance_error":float(np.max(np.abs(balance-1)))})


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--note",type=int,choices=range(1,21)); group.add_argument("--all",action="store_true")
    args=parser.parse_args()
    for number in (range(1,21) if args.all else [args.note]): globals()[f"note{number:02d}"]()
