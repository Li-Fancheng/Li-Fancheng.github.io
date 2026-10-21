"""Independent checks for the derivations and preserved trained state."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parent
rng=np.random.default_rng(20261022);report={}
# Exact quadratic competition flow against an independent ODE integrator.
cg=.05;cm=5.;l=.02;gg=.015;gm=1.
A=np.array([[gg*(1+l*cg),gg],[gm,gm*(1+l*cm)]]);drive=np.array([gg,gm]);steady=np.linalg.solve(A,drive)
ts=np.geomspace(.01,2000,120);sol=solve_ivp(lambda t,z:drive-A@z,[0,ts[-1]],[0,0],t_eval=ts,rtol=1e-10,atol=1e-12)
exact=np.array([steady-expm(-A*t)@steady for t in ts]);err=float(np.max(abs(sol.y.T-exact)));assert err<1e-7;report['competition_ODE_max_error']=err
# Gaussian quadratic risk identity including nonzero trace.
d=5;B=rng.normal(size=(d,d));M=(B+B.T)/2;x=rng.normal(size=(500000,d));q=np.einsum('ni,ij,nj->n',x,M,x);theory=2*np.trace(M@M)+np.trace(M)**2;relative=abs(q.dot(q)/len(q)/theory-1);assert relative<.02
report['Gaussian_quadratic_risk']={'theory':float(theory),'Monte_Carlo':float(q.dot(q)/len(q)),'relative_error':float(relative),'samples':len(q)}
# Rayleigh gradient and invariants, in the centered representation subspace.
p=7;row=np.zeros((p-2,p))
for i in range(p-2):row[i,i:i+3]=[1,-2,1]
M=row.T@row;E=rng.normal(size=(p,2));E-=E.mean(0);S=np.sum(E*E);energy=np.sum(E*(M@E))/S;grad=2/S*(M@E-energy*E);v=rng.normal(size=E.shape);v-=v.mean(0);h=1e-6
f=lambda z:np.sum(z*(M@z))/np.sum(z*z)
fd=(f(E+h*v)-f(E-h*v))/(2*h);error=abs(fd-np.sum(grad*v));assert error<1e-7
assert abs(np.sum(E*grad))<1e-10 and np.max(abs(grad.mean(0)))<1e-10
assert np.linalg.matrix_rank(row)==p-2 and np.max(abs(row@np.arange(p)))==0
report['Rayleigh_directional_gradient_error']=float(error);report['Rayleigh_radial_derivative']=float(np.sum(E*grad));report['constraint_nullity']=p-int(np.linalg.matrix_rank(row))
# Scalar L1 soft threshold satisfies first-order subgradient conditions.
z=rng.normal(size=100);ll=.4;w=np.sign(z)*np.maximum(abs(z)-ll,0);active=w!=0
assert np.max(abs((w-z)[active]+ll*np.sign(w[active])))<1e-12;assert np.max(abs(z[~active]))<=ll
report['soft_threshold_optimality']='passed'
# Independent verification of saved weights, split and history endpoint.
weights=np.load(ROOT/'figures/modular_weights.npz');W,V=weights['W'],weights['V'];p=V.shape[1];aa,bb=np.meshgrid(np.arange(p),np.arange(p),indexing='ij');y=(aa+bb)%p;F=(W[aa]+W[p+bb])**2@V
ce=logsumexp(F,axis=2)-np.take_along_axis(F,y[:,:,None],axis=2)[:,:,0];correct=F.argmax(2)==y
tr,te=weights['train'],weights['test'];assert not set(tr)&set(te) and set(tr)|set(te)==set(range(p*p));assert correct.all()
meta=json.loads((ROOT/'figures/modular_run.json').read_text());assert abs(ce.ravel()[tr].mean()-meta['final'][1])<1e-11;assert abs(ce.ravel()[te].mean()-meta['final'][2])<1e-11
F-=F.mean(2,keepdims=True);offset=(aa[:,:,None]+bb[:,:,None]-np.arange(p))%p;curve=np.array([F[offset==r].mean() for r in range(p)]);PF=curve[offset];res=F-PF
risk=lambda Z:float((logsumexp(Z,axis=2)-np.take_along_axis(Z,y[:,:,None],axis=2)[:,:,0]).mean())
assert risk(PF)<=risk(F)+1e-12;orth=abs(np.sum(F*F)-np.sum(PF*PF)-np.sum(res*res))/np.sum(F*F);assert orth<1e-12
report['saved_training']={'train_ce':float(ce.ravel()[tr].mean()),'test_ce':float(ce.ravel()[te].mean()),'full_accuracy':float(correct.mean()),'weights_sha256':hashlib.sha256((ROOT/'figures/modular_weights.npz').read_bytes()).hexdigest()}
report['symmetry_projection']={'full_risk_before':risk(F),'full_risk_after':risk(PF),'Pythagoras_relative_error':float(orth)}
report['status']='passed';(ROOT/'review_diagnostics.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
