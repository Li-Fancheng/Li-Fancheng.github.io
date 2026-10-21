"""Independent Grokking model calculations and terminal-weight diagnostics.

No new neural-network training is performed by this file.
"""
from pathlib import Path
from itertools import combinations
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.special import logsumexp
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.prop_cycle':plt.cycler(color=['#2563eb','#e87924','#059669','#9333ea']),'svg.fonttype':'none','svg.hashsalt':'grokking-series-2026','savefig.facecolor':'white'})
reports={}
def save(fig,name,stats):
 fig.tight_layout(pad=1.5);fig.savefig(OUT/(name+'.svg'),metadata={'Date':None});fig.savefig(OUT/(name+'.png'),dpi=200);plt.close(fig)
 (OUT/(name+'.json')).write_text(json.dumps(stats,indent=2,allow_nan=False)+'\n');reports[name]=stats

def random_design():
 z,w=leggauss(1600);theta=(z+1)*np.pi/2;w=w*np.pi/2;ts=np.geomspace(.01,1e5,380)
 fig,ax=plt.subplots(1,3,figsize=(13.5,4.1));stats={}
 for a,g in zip(ax,[.8,1.,1.2]):
  lam=1+g-2*np.sqrt(g)*np.cos(theta);mass=2*np.sin(theta)**2/(np.pi*lam)*w;atom=max(1-1/g,0)
  assert abs(mass.sum()+atom-1)<1e-8
  decay=np.exp(-2*ts[:,None]*lam);tr=.5*(decay@(mass*lam));te=.5*(decay@mass+atom)
  a.loglog(ts,np.maximum(tr,1e-25),label='Training risk');a.loglog(ts,np.maximum(te,1e-25),label='Population test risk')
  if g==1:
   tr_ref=tr[-1]*(ts/ts[-1])**(-1.5);te_ref=te[-1]*(ts/ts[-1])**(-.5)
   a.loglog(ts,tr_ref,'k:',lw=1);a.loglog(ts,te_ref,'k:',lw=1)
   slopes=[float(np.polyfit(np.log(ts[-70:]),np.log(v[-70:]),1)[0]) for v in [tr,te]]
   assert np.allclose(slopes,[-1.5,-.5],atol=.003)
  else:slopes=None
  if atom:a.axhline(.5*atom,color='gray',ls=':',label='Unobserved-space floor')
  a.set(xlabel='Gradient-flow time',ylabel='Half squared risk',ylim=(1e-12,1),title=f'Same distribution: d/D={g:g}');a.legend(fontsize=8)
  stats[str(g)]={'continuous_mass':float(mass.sum()),'zero_mass':atom,'lower_edge':float((1-np.sqrt(g))**2),'last_test_risk':float(te[-1]),'critical_slopes_train_test':slopes}
 save(fig,'random_design',{'quadrature_nodes':1600,'initial_error_energy':1,'kind':'MP spectral integrals, not trained nonlinear networks','cases':stats})

def competition():
 cg=.05;cm=5.;gg=.015;gm=1.;l=.02;t=np.geomspace(.01,2e4,450)
 A=np.array([[gg*(1+l*cg),gg],[gm,gm*(1+l*cm)]]);steady=np.linalg.solve(A,[gg,gm]);paths=np.array([steady-expm(-A*tt)@steady for tt in t]);g,m=paths.T
 tr=.5*(g+m-1)**2;te=.5*(g-1)**2
 fig,ax=plt.subplots(1,3,figsize=(13.5,4.1));ax[0].semilogx(t,g,label='Shared-rule amplitude g');ax[0].semilogx(t,m,label='Training-only amplitude m');ax[0].set(xlabel='Time',ylabel='Output amplitude',title='Fast fitting, slow solution selection');ax[0].legend(fontsize=8)
 ax[1].loglog(t,tr,label='Training risk');ax[1].loglog(t,te,label='Test risk');ax[1].set(xlabel='Time',ylabel='Half squared risk',title='Both curves follow a smooth linear ODE');ax[1].legend(fontsize=8)
 ls=np.geomspace(1e-4,10,120);taus=[]
 for ll in ls:
  mat=np.array([[gg*(1+ll*cg),gg],[gm,gm*(1+ll*cm)]]);taus.append(1/min(np.linalg.eigvals(mat).real))
 ax[2].loglog(ls,taus,label='Exact slow relaxation time');ax[2].set(xlabel='Regularization lambda',ylabel='Slow time',title='Stronger bias also changes the endpoint');a2=ax[2].twinx();a2.semilogx(ls,cm/(cg+cm+ls*cg*cm),color='#e87924',ls='--',label='Final rule amplitude');a2.set(ylabel='Final g',ylim=(0,1.05));ax[2].legend(fontsize=8,loc='lower left');a2.legend(fontsize=8,loc='upper right')
 analytic=np.array([cm,cg])/(cg+cm+l*cg*cm);assert np.max(abs(steady-analytic))<1e-10
 save(fig,'solution_competition',{'cg':cg,'cm':cm,'gamma_g':gg,'gamma_m':gm,'lambda':l,'steady_g_m':steady.tolist(),'exact_slow_time':float(1/min(np.linalg.eigvals(A).real)),'quasistatic_slow_time':float(1/(gg*l*(cg+cm/(1+l*cm)))),'kind':'fixed-feature convex quadratic model'})

def fourier():
 p=19;k=3;d=np.arange(p);curve=np.cos(2*np.pi*k*d/p);scales=np.geomspace(.02,300,200)
 ce=np.array([logsumexp(s*curve)-s*curve[0] for s in scales]);margin=scales*(1-np.cos(2*np.pi/p))
 fig,ax=plt.subplots(1,3,figsize=(13.5,4.1));ax[0].plot(d,curve,'o-');ax[0].scatter([0],[1],color='#e87924',s=80,label='Correct residual d=0');ax[0].set(xlabel='Residual a+b-c (mod 19)',ylabel='Unit-scale single-pair logit',title='A compact cyclic algorithm');ax[0].legend(fontsize=8)
 ax[1].semilogx(scales,np.ones_like(scales),label='Classification accuracy');ax[1].set(xlabel='Positive logit amplitude',ylabel='Accuracy',ylim=(0,1.1),title='Correct at every positive scale');ax[1].legend(fontsize=8)
 ax[2].loglog(scales,ce,label='Cross entropy');ax[2].loglog(scales,margin,label='Minimum class margin');ax[2].set(xlabel='Logit amplitude',ylabel='Cross entropy / margin',title='Probability quality still changes');ax[2].legend(fontsize=8)
 all_logits=np.cos(2*np.pi*k*((np.arange(p)[:,None,None]+np.arange(p)[None,:,None]-np.arange(p)[None,None,:])%p)/p)
 assert np.all(all_logits.argmax(2)==(np.arange(p)[:,None]+np.arange(p)[None,:])%p)
 save(fig,'fourier_algorithm',{'modulus':p,'frequency':k,'minimum_unit_scale_margin':float(1-np.cos(2*np.pi/p)),'verified_complete_domain_accuracy':1.0,'kind':'explicitly constructed Fourier function, not training'})

def constraints():
 rng=np.random.default_rng(20261021);p=12;pairs=list(combinations(range(p),2))+[(i,i) for i in range(p)];fractions=np.linspace(.15,1,18);trials=80
 def matrix(sub):
  groups={}
  for i,j in sub:groups.setdefault(i+j,[]).append((i,j))
  rows=[]
  for group in groups.values():
   for (i,j),(m,n) in combinations(group,2):
    row=np.zeros(p);row[i]+=1;row[j]+=1;row[m]-=1;row[n]-=1;rows.append(row)
  return np.array(rows).T@np.array(rows)/len(rows) if rows else np.zeros((p,p))
 probs=[];means=[];q25=[];q75=[];nullities=[];example=None
 for frac in fractions:
  gaps=[];nulls=[]
  for _ in range(trials):
   inds=rng.choice(len(pairs),int(round(frac*len(pairs))),replace=False);M=matrix([pairs[i] for i in inds]);ev=np.linalg.eigvalsh(M);null=int(np.sum(ev<1e-10));nulls.append(null);gaps.append(max(float(ev[2]),0) if null==2 else 0.)
   if example is None and null==2 and frac>=.6:example=M
  probs.append(np.mean(np.array(nulls)==2));means.append(np.mean(gaps));q25.append(np.quantile(gaps,.25));q75.append(np.quantile(gaps,.75));nullities.append(float(np.mean(nulls)))
 assert probs[-1]==1 and example is not None
 vals,U=np.linalg.eigh(example);E0=rng.normal(size=(p,2));E0-=E0.mean(0);E0/=np.linalg.norm(E0);ts=np.geomspace(.01,3000,200);history=[];energies=[]
 for t in ts:
  raw=U@(np.exp(-2*vals*t)[:,None]*(U.T@E0));E=raw/np.linalg.norm(raw);history.append(E);energies.append(float(np.sum(E*(example@E))))
 assert np.max(np.diff(energies))<1e-9
 assert max(abs(np.linalg.norm(E)-1) for E in history)<1e-12
 assert np.max(abs(example@np.ones(p)))<1e-12 and np.max(abs(example@np.arange(p)))<1e-12
 fig,ax=plt.subplots(1,3,figsize=(13.5,4.1));ax[0].plot(fractions,probs,'o-');ax[0].set(xlabel='Fraction of unordered training pairs',ylabel='Probability nullity = 2',ylim=(-.03,1.03),title=f'Structure identifiability, p={p}, {trials} draws')
 ax[1].plot(fractions,means,'o-',label='Mean structural gap');ax[1].fill_between(fractions,q25,q75,alpha=.18,label='25-75% across datasets');ax[1].set(xlabel='Training fraction',ylabel='Third eigenvalue (zero if underconstrained)',title='Data also control formation speed');ax[1].legend(fontsize=8)
 ax[2].scatter(E0[:,0],E0[:,1],color='gray',marker='x',label='Initial');last=history[-1];sc=ax[2].scatter(last[:,0],last[:,1],c=np.arange(p),cmap='viridis',s=35,label='Normalized late-time shape');ax[2].plot(last[:,0],last[:,1],color='#2563eb',lw=.8);ax[2].set(xlabel='Embedding coordinate 1',ylabel='Embedding coordinate 2',title='Shape becomes affine without collapse');ax[2].legend(fontsize=7);fig.colorbar(sc,ax=ax[2],label='Symbol index',fraction=.045)
 save(fig,'representation_constraints',{'seed':20261021,'symbols':p,'unordered_pairs':len(pairs),'trials_per_fraction':trials,'fractions':fractions.tolist(),'probability_unique_affine_shape':probs,'mean_nullity':nullities,'mean_structural_gap':means,'flow_initial_energy':energies[0],'flow_final_energy':energies[-1],'flow_norm_conserved':True,'kind':'constraint matrices and normalized Rayleigh flow, no decoder training'})

def interventions():
 data=np.load(OUT/'modular_weights.npz');W,V=data['W'],data['V'];p=V.shape[1];aa,bb=np.meshgrid(np.arange(p),np.arange(p),indexing='ij');y=(aa+bb)%p
 F=((W[aa]+W[p+bb])**2)@V;F-=F.mean(2,keepdims=True);offset=(aa[:,:,None]+bb[:,:,None]-np.arange(p))%p
 curve=np.array([F[offset==r].mean() for r in range(p)]);P=curve[offset];Q=F-P;freq=np.fft.fft(curve);k=int(np.argmax(abs(freq[1:(p+1)//2])**2)+1);topfreq=np.zeros(p,complex);topfreq[k]=freq[k];topfreq[-k]=freq[-k];top=np.fft.ifft(topfreq).real[offset]
 def metric(Z):
  ce=logsumexp(Z,axis=2)-np.take_along_axis(Z,y[:,:,None],axis=2)[:,:,0];acc=Z.argmax(2)==y;bad=Z.copy();np.put_along_axis(bad,y[:,:,None],-np.inf,axis=2);margin=np.take_along_axis(Z,y[:,:,None],axis=2)[:,:,0]-bad.max(2)
  return {'full_ce':float(ce.mean()),'full_accuracy':float(acc.mean()),'train_ce':float(ce.ravel()[data['train']].mean()),'test_ce':float(ce.ravel()[data['test']].mean()),'train_accuracy':float(acc.ravel()[data['train']].mean()),'test_accuracy':float(acc.ravel()[data['test']].mean()),'min_margin':float(margin.min())}
 models={s:metric(Z) for s,Z in [('Original',F),('Invariant',P),('Residual',Q),('Top pair',top),('Without pair',F-top)]};ss=np.linspace(0,4,161);sweep=[metric(P+s*Q) for s in ss]
 assert models['Original']['full_accuracy']==1 and models['Invariant']['full_ce']<=models['Original']['full_ce']
 orth=abs(np.sum(F*F)-np.sum(P*P)-np.sum(Q*Q))/np.sum(F*F);assert orth<1e-12
 fig,ax=plt.subplots(1,2,figsize=(11.5,4.5));names=list(models);ax[0].bar(np.arange(5),[models[n]['full_ce'] for n in names],color=['#2563eb','#059669','#9333ea','#e87924','#64748b']);ax[0].set_yscale('log');ax[0].set(xticks=np.arange(5),xticklabels=names,ylabel='Full-domain cross entropy',title='Different functions of the same saved weights');ax[0].tick_params(axis='x',labelrotation=18);a2=ax[0].twinx();a2.plot(np.arange(5),[models[n]['full_accuracy'] for n in names],'ko--',label='Accuracy');a2.set(ylabel='Full-domain accuracy',ylim=(-.05,1.1));a2.legend(fontsize=8,loc='upper right')
 ax[1].plot(ss,[v['train_accuracy'] for v in sweep],label='Train accuracy');ax[1].plot(ss,[v['test_accuracy'] for v in sweep],label='Test accuracy');ax[1].axvline(1,color='gray',ls=':',label='Original function');ax[1].set(xlabel='Residual scale s in P(F)+s*(I-P)(F)',ylabel='Accuracy',ylim=(-.03,1.05),title='Residual amplification harms test before train');ax[1].legend(fontsize=8)
 h=np.genfromtxt(OUT/'modular_history.csv',delimiter=',',names=True);times={}
 for a in [.5,.9,.95,.99,1.]:
  valid=np.flatnonzero(h['test_accuracy']>=a);below=np.flatnonzero(h['test_accuracy']<a);times[str(a)]={'first_test':int(h['step'][valid[0]]) if len(valid) else None,'sustained_test':int(h['step'][below[-1]+1]) if len(below) and below[-1]+1<len(h) else (0 if not len(below) else None)}
 report={'projection_energy':float(np.sum(P*P)/np.sum(F*F)),'projection_Pythagoras_relative_error':float(orth),'dominant_frequency':k,'models':models,'residual_sweep_at_selected_scales':{str(s):metric(P+s*Q) for s in [-1,0,.25,.5,1,2,4]},'threshold_times':times,'protocol':'Known-task, full-grid terminal function analysis. Original logit scale; no temperature fitting, new gradients, intermediate checkpoints or hyperparameter selection.'}
 save(fig,'function_interventions',report)

if __name__=='__main__':
 random_design();competition();fourier();constraints();interventions()
 (ROOT/'numerical_checks.json').write_text(json.dumps(reports,indent=2,allow_nan=False)+'\n')
 print(json.dumps({name:'passed' for name in reports}))
