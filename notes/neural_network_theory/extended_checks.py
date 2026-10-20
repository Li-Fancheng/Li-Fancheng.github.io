#!/usr/bin/env python3
"""Independent calculations for the expanded notes; run without arguments.

Exports four additional figures and extended_diagnostics.json. All examples
are explicitly defined models; no network weights are updated by this script.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/neural-theory-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/neural-theory-cache")
from pathlib import Path
import json
import shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.linalg import solve_continuous_lyapunov
from scipy.special import erf

ROOT=Path(__file__).resolve().parent
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
    "axes.spines.top":False,"axes.spines.right":False,
    "mathtext.fontset":"cm","svg.hashsalt":"neural-expanded-2026"})
report={}


def save(fig,slug,name):
    dest=ROOT/slug/"figures"; dest.mkdir(parents=True,exist_ok=True)
    for ax in fig.axes:
        ax.grid(alpha=.15)
        if ax.get_legend_handles_labels()[0]: ax.legend(frameon=False,fontsize=8)
    fig.savefig(dest/(name+".svg"),metadata={"Date":None})
    fig.savefig(dest/(name+".png"),dpi=220)
    plt.close(fig)
    merged_owner = {
        ("time_scaling", "critical_spectrum"): "../scaling_laws",
        ("implicit_bias", "initialization_bias"): "deep_optimization",
    }.get((slug, name))
    if merged_owner:
        public = ROOT/merged_owner/"figures"
        public.mkdir(parents=True,exist_ok=True)
        for extension in ("svg","png"):
            shutil.copyfile(dest/(name+"."+extension), public/(name+"."+extension))


# Scalar early stopping: numerical optimum of the full risk.
s=.07; signal=1.4; D=120; sigma=.4
t_star=np.log1p(D*s*signal**2/sigma**2)/s
risk=lambda t:np.exp(-2*s*t)*signal**2+sigma**2/(D*s)*(-np.expm1(-s*t))**2
numerical=minimize_scalar(risk,bounds=(0,10*t_star),method="bounded").x
assert abs(t_star-numerical)<1e-3
report["single_mode_stopping"]={"analytic":float(t_star),"numerical":float(numerical)}

# Unbalanced deep-linear model: Riccati solution versus the scalar product ODE.
c0=.2; target=.8; p0=-.1; ts=np.linspace(0,12,400)
omega=np.sqrt(4*target**2+c0**2)
zp=(2*target+omega)/c0; zm=(2*target-omega)/c0
z0=np.exp(np.arcsinh(2*p0/c0)); ratio=(z0-zp)/(z0-zm)*np.exp(-omega*ts)
zz=(zp-ratio*zm)/(1-ratio); closed=c0/4*(zz-1/zz)
sol=solve_ivp(lambda t,p:np.sqrt(c0**2+4*p*p)*(target-p),[0,12],[p0],t_eval=ts,rtol=1e-10,atol=1e-12)
error=float(np.max(np.abs(sol.y[0]-closed))); assert error<1e-7
report["unbalanced_deep_linear_max_error"]=error

# MP resolvent: direct quadrature, including the atom at zero.
mp_errors=[]
for g in [.4,1.,2.5]:
    lo,hi=(1-np.sqrt(g))**2,(1+np.sqrt(g))**2
    for lam in [.02,.3,2.]:
        density=lambda x:np.sqrt(max((hi-x)*(x-lo),0))/(2*np.pi*g*x)
        integral=quad(lambda x:density(x)/(x+lam),lo,hi,epsabs=1e-10)[0]+max(1-1/g,0)/lam
        exact=2/(lam+1-g+np.sqrt((lam+1-g)**2+4*g*lam))
        mp_errors.append(abs(integral-exact))
assert max(mp_errors)<1e-7
report["MP_resolvent_max_error"]=float(max(mp_errors))

# Signed diagonal model: mirror KKT solution versus original gradient flow.
x=np.array([1.,2.]); eps_values=np.geomspace(1e-4,10,90); solutions=[]; flow_errors=[]
for eps in eps_values:
    c=eps**2
    residual=lambda u:np.dot(x,2*c*np.sinh(4*u*x))-1
    u=brentq(residual,0,10,xtol=1e-13)
    beta=2*c*np.sinh(4*u*x); solutions.append(beta)
for eps in [1e-3,.1,1.]:
    c=eps**2
    sol=solve_ivp(lambda t,b:-4*np.sqrt(b*b+4*c*c)*x*(x@b-1),[0,200],np.zeros(2),rtol=1e-9,atol=1e-11)
    u=brentq(lambda u:x@(2*c*np.sinh(4*u*x))-1,0,10)
    flow_errors.append(float(np.linalg.norm(sol.y[:,-1]-2*c*np.sinh(4*u*x))))
assert max(flow_errors)<1e-6
solutions=np.array(solutions); fig,ax=plt.subplots(1,2,figsize=(10.2,3.7),layout="constrained")
ax[0].semilogx(eps_values,solutions[:,0],label="beta 1")
ax[0].semilogx(eps_values,solutions[:,1],label="beta 2")
ax[0].axhline(.2,ls="--",color="grey",lw=.8); ax[0].axhline(.4,ls="--",color="grey",lw=.8)
ax[0].set(xlabel="Factor initialization epsilon",ylabel="Selected interpolator")
ax[1].plot(solutions[:,0],solutions[:,1],label="Mirror solution path")
ax[1].scatter([0,.2],[.5,.4],s=35,color=["#e87924","#059669"],label="L1 and L2 limits")
ax[1].set(xlabel="beta 1",ylabel="beta 2",title="Constraint: beta 1 + 2 beta 2 = 1")
save(fig,"implicit_bias","initialization_bias")
report["diagonal_model_max_flow_error"]=max(flow_errors)

# Exact noise-filter covariance and probability-current production.
H=np.diag([1.,4.]); A=np.array([[1.,.8],[.8,2.]])
S=solve_continuous_lyapunov(H,2*A); B=-H+A@np.linalg.inv(S)
production=np.trace(B.T@np.linalg.solve(A,B)@S)
assert np.linalg.norm(H@S+S@H-2*A)<1e-12 and production>0
report["OU_current"]={"covariance":S.tolist(),"HS_minus_A":(H@S-A).tolist(),"production":float(production)}

# Hopfield zero-temperature RS retrieval branch, not the complete phase diagram.
def alpha_of_y(y):
    return (erf(y)/(np.sqrt(2)*y)-np.sqrt(2/np.pi)*np.exp(-y*y))**2
opt=minimize_scalar(lambda y:-alpha_of_y(y),bounds=(.1,5),method="bounded")
ys=np.linspace(.02,5,1000); alphas=alpha_of_y(ys); peak=float(-opt.fun)
assert .137<peak<.139
fig,ax=plt.subplots(1,2,figsize=(10.2,3.7),layout="constrained")
ax[0].plot(ys,alphas); ax[0].scatter([opt.x],[peak],color="#e87924",label=f"Max alpha = {peak:.4f}")
ax[0].set(xlabel="Field ratio y",ylabel="Load alpha(y)")
ax[1].plot(alphas,erf(ys)); ax[1].scatter([peak],[erf(opt.x)],color="#e87924",label="Retrieval branch endpoint")
ax[1].set(xlabel="Load alpha",ylabel="Retrieval overlap m")
save(fig,"hopfield","retrieval_capacity")
report["Hopfield_RS_branch"]={"alpha_max":peak,"field_ratio":float(opt.x),"overlap_at_endpoint":float(erf(opt.x))}

# Critical MP spectrum: quadrature after s=4 sin(theta)^2 removes the hard-edge singularity.
nodes,weights=np.polynomial.legendre.leggauss(1100)
theta=(nodes+1)*np.pi/4; ds_weight=weights*np.pi/4*4/np.pi*np.cos(theta)**2
rates=4*np.sin(theta)**2; times=np.geomspace(1,1e5,160)
decay=np.exp(-2*times[:,None]*rates)
train_signal=.5*np.sum(ds_weight*rates*decay,axis=1)
train_noise=.5*np.sum(ds_weight*decay,axis=1)
test_variance=np.sum(ds_weight/rates*(-np.expm1(-times[:,None]*rates))**2,axis=1)
slopes=[float(np.polyfit(np.log(times[-45:]),np.log(v[-45:]),1)[0]) for v in [train_signal,train_noise,test_variance]]
assert np.allclose(slopes,[-1.5,-.5,.5],atol=.025)
constant=4/np.sqrt(np.pi)*(1-1/np.sqrt(2))
assert abs(test_variance[-1]/np.sqrt(times[-1])-constant)<.01
fig,ax=plt.subplots(1,2,figsize=(10.2,3.7),layout="constrained")
ax[0].loglog(times,train_signal,label="Training signal: t^(-3/2)")
ax[0].loglog(times,train_noise,label="Training noise: t^(-1/2)")
ax[0].loglog(times,test_variance,label="Test variance: t^(1/2)")
ax[0].set(xlabel="Gradient-flow time",ylabel="Error contribution; unit amplitudes")
ax[1].semilogx(times,test_variance/np.sqrt(times),label="Test variance / sqrt(t)")
ax[1].axhline(constant,color="black",ls="--",label="Analytic large-time constant")
ax[1].set(xlabel="Gradient-flow time",ylabel="Normalized variance")
save(fig,"time_scaling","critical_spectrum")
report["critical_MP_time_slopes"]=slopes
report["critical_variance_constant"]={"analytic":float(constant),"last_quadrature":float(test_variance[-1]/np.sqrt(times[-1]))}

# Data-scaling leading constant, using one million modes and the analytic lambda.
indices=np.arange(1,1000001,dtype=float); eta=indices**-2; coeff=6/np.pi**2; a2=coeff*eta
sigma=.3; normalized=[]
for D in [1e3,1e4,1e5,1e6]:
    lam=sigma**2/(coeff*D)
    risk=np.sum((lam/(eta+lam))**2*a2+sigma**2/D*(eta/(eta+lam))**2)
    normalized.append(float(np.sqrt(D)*risk))
expected=np.sqrt(6)/2*sigma
assert abs(normalized[-1]/expected-1)<.01
report["data_scaling_constant"]={"analytic":float(expected),"sqrt_D_times_risk":normalized,"modes":len(indices)}

# Function-space symmetry of the saved modular-addition model.
weights=np.load(ROOT/"grokking/figures/modular_weights.npz"); p=19
aa,bb=np.meshgrid(np.arange(p),np.arange(p),indexing="ij"); X=np.zeros((p*p,2*p))
X[np.arange(p*p),aa.ravel()]=1; X[np.arange(p*p),p+bb.ravel()]=1
F=((X@weights["W"])**2@weights["V"]).reshape(p,p,p); F-=F.mean(axis=2,keepdims=True)
offset=(aa[:,:,None]+bb[:,:,None]-np.arange(p)[None,None,:])%p
mean=np.array([F[offset==k].mean() for k in range(p)])
std=np.array([F[offset==k].std() for k in range(p)]); projected=mean[offset]
fraction=float(np.sum(projected**2)/np.sum(F**2)); relative=float(np.linalg.norm(F-projected)/np.linalg.norm(F))
assert abs(fraction+relative**2-1)<1e-12
fig,ax=plt.subplots(1,2,figsize=(10.2,3.7),layout="constrained")
ax[0].errorbar(np.arange(p),mean,yerr=std,fmt="o-",capsize=2)
ax[0].set(xlabel="a + b - c mod 19",ylabel="Centered logit mean +/- std")
power=np.abs(np.fft.rfft(mean))**2; ax[1].bar(np.arange(1,len(power)),power[1:]/power[1:].sum())
ax[1].set(xlabel="Fourier frequency of projected logits",ylabel="Relative non-DC power")
save(fig,"grokking","logit_symmetry")
report["modular_logit_symmetry"]={"projected_energy_fraction":fraction,"relative_projection_error":relative}

(ROOT/"extended_diagnostics.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
print(json.dumps(report,indent=2))
