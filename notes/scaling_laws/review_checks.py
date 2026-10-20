"""Independent checks for formulas introduced by this topic."""
from pathlib import Path
import json
import numpy as np
from scipy.special import gammaincc,gamma
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parent
rng=np.random.default_rng(20261021)
report={}
# Differentiate the two-layer prediction in parameter and function coordinates.
n=7;B=rng.normal(size=(n,n));v=rng.normal(size=n);g=rng.normal(size=n)
dv=-B@g;dB=-np.outer(v,g);dtheta=dB.T@v+B.T@dv
metric=-(B.T@B+(v@v)*np.eye(n))@g
error=float(np.max(np.abs(dtheta-metric)));assert error<1e-10
report['adaptive_metric_max_error']=error
# A fixed row, positive support, Dirichlet smoothing: CE excess tends to (k-1)/(2D).
p=np.array([.1,.2,.3,.4]);records=[]
for D in [1000,10000,100000]:
 counts=rng.multinomial(D,p,size=60000);phat=(counts+.5)/(D+2)
 kl=np.sum(p*np.log(p/phat),axis=1);scaled=D*kl
 records.append({'D':D,'D_times_mean_KL':float(scaled.mean()),'standard_error':float(scaled.std(ddof=1)/np.sqrt(len(scaled)))})
assert abs(records[-1]['D_times_mean_KL']-1.5)<.025
report['multinomial_cross_entropy']={'support':len(p),'predicted_constant':1.5,'smoothing':.5,'trials':60000,'results':records}
# Poisson multiple-observation integral coefficient.
s=1.6;power=1-1/s;checks=[]
for tau in [1,5]:
 val=quad(lambda u:u**(-1/s)*gammaincc(tau,u),0,np.inf,epsabs=1e-8)[0]
 exact=sum(gamma(j+power)/gamma(j+1) for j in range(tau))
 checks.append(abs(val-exact)/exact)
assert max(checks)<1e-7
report['Poisson_integral_relative_errors']=checks
# Explicit geometric bias/variance optimum and compute optimum.
d=4;r=2;D=4000.;h=(d/(2*r*D))**(1/(2*r+d))
num=minimize_scalar(lambda z:np.exp(z)**(2*r)+1/D/np.exp(z)**d,bounds=(-10,0),method='bounded').x
assert abs(np.log(h)-num)<1e-5
report['noisy_resolution_log_optimum_error']=float(abs(np.log(h)-num))
a=.4;b=.3;A=2.;Q=3.;C=1e10;k=6
N=(a*A/(b*Q))**(1/(a+b))*(C/k)**(b/(a+b))
num=minimize_scalar(lambda z:A*np.exp(-a*z)+Q*(k*np.exp(z)/C)**b,bounds=(0,np.log(C/k)),method='bounded').x
assert abs(np.log(N)-num)<1e-4
report['compute_log_optimum_error']=float(abs(np.log(N)-num))
report['status']='passed'
(ROOT/'review_diagnostics.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
