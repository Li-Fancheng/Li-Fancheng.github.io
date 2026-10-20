"""Independent mechanism calculations. No transformer training is performed."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import gammaincc,gamma
from scipy.optimize import brentq
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.prop_cycle':plt.cycler(color=['#2563eb','#e87924','#059669','#9333ea']),'svg.fonttype':'none','savefig.facecolor':'white'})
report={}
def save(fig,name,stats):
    fig.tight_layout(pad=1.6)
    fig.savefig(OUT/(name+'.svg'))
    fig.savefig(OUT/(name+'.png'),dpi=200)
    plt.close(fig)
    (OUT/(name+'.json')).write_text(json.dumps(stats,indent=2,ensure_ascii=False)+'\n')
    report[name]=stats

def statistical():
    fig,ax=plt.subplots(1,3,figsize=(13.5,4.1))
    D=np.geomspace(10,1e6,100)
    ax[0].loglog(D,1/D,label='Mean squared error: $1/D$')
    ax[0].loglog(D,np.sqrt(2/np.pi/D),label='Mean absolute error: $D^{-1/2}$')
    ax[0].set(xlabel='Independent samples D',ylabel='Estimation error',title='Same estimator, different loss');ax[0].legend(fontsize=8)
    xx=(np.arange(200003)+.5)/200003; yy=np.sin(2*np.pi*xx)
    sizes=np.array([8,12,16,24,32,48,64,96,128,192,256]); er0=[];er1=[]
    for m in sizes:
        lo=np.floor(xx*m)/m
        er0.append(np.mean((yy-np.sin(2*np.pi*(lo+.5/m)))**2))
        nodes=np.linspace(0,1,m+1)
        er1.append(np.mean((yy-np.interp(xx,nodes,np.sin(2*np.pi*nodes)))**2))
    slopes=[float(np.polyfit(np.log(sizes[-5:]),np.log(e[-5:]),1)[0]) for e in [er0,er1]]
    assert abs(slopes[0]+2)<.01 and abs(slopes[1]+4)<.01
    ax[1].loglog(sizes,er0,'o-',label='Piecewise constant')
    ax[1].loglog(sizes,er1,'s-',label='Piecewise linear')
    ax[1].set(xlabel='Uniform intervals M',ylabel='Integrated squared error',title='Noiseless resolution on a circle');ax[1].legend(fontsize=8)
    r=2;d=4;A=B=1
    h=(d*B/(2*r*A*D))**(1/(2*r+d)); risk=A*h**(2*r)+B/D/h**d
    ax[2].loglog(D,risk,label='Optimized noisy risk')
    ax[2].loglog(D,1/D,ls='--',label='Noiseless $D^{-4/d}$')
    ax[2].set(xlabel='Independent samples D',ylabel='Model risk',title='Noise changes the exponent (d=4)');ax[2].legend(fontsize=8)
    noisy_slope=float(np.polyfit(np.log(D),np.log(risk),1)[0]); assert abs(noisy_slope+.5)<1e-10
    save(fig,'statistical_geometry',{'interpolation_slopes':slopes,'noisy_slope':noisy_slope,'dimension':d,'local_approximation_order':r,'integration_points':len(xx)})

def skills():
    fig,ax=plt.subplots(1,3,figsize=(13.5,4.1));a=.6;s=1+a;K=100000
    ks=np.arange(1,K+1,dtype=float);p=ks**(-s);p/=p.sum();Ds=np.geomspace(10,1e7,130)
    curves=[]
    for tau in [1,5]:
        risk=np.array([np.sum(p*gammaincc(tau,D*p)) for D in Ds]);curves.append(risk)
        ax[0].loglog(Ds,risk,label=f'Poisson threshold {tau}')
    expected=a/(1+a);ref=curves[0][30]*(Ds/Ds[30])**(-expected)
    ax[0].loglog(Ds,ref,'k:',label=f'Infinite-tail slope -{expected:.3f}')
    ax[0].set(xlabel='Independent examples D',ylabel='Unmastered test mass',title='Random coverage of a long tail');ax[0].legend(fontsize=8)
    scales=np.geomspace(1,1e7,300);inds=np.geomspace(1,100,10).astype(int);sel=sorted(set(inds));pk=p[:200]
    mastering=-np.expm1(-scales[:,None]*pk[None,:])
    for i in sel:ax[1].semilogx(scales,mastering[:,i-1],lw=.8,alpha=.55)
    ax[1].semilogx(scales,mastering@pk/pk.sum(),color='black',lw=2.5,label='Frequency-weighted aggregate')
    ax[1].set(xlabel='Resource scale',ylabel='Mastered fraction',title='Many local transitions, smooth average');ax[1].legend(fontsize=8)
    finite={}
    for cap in [100,1000,100000]:
        probs=p[:cap]/p[:cap].sum();v=np.array([np.sum(probs*np.exp(-D*probs)) for D in Ds]);finite[str(cap)]=v[-1]
        ax[2].loglog(Ds,np.maximum(v,1e-15),label=f'K={cap:,}')
    ax[2].set(xlabel='Independent examples D',ylabel='Unseen test mass',ylim=(1e-8,1),title='Finite vocabulary cuts off the law');ax[2].legend(fontsize=8)
    # Infinite-support integral and finite sum agree before truncation dominates.
    c=1/np.sum(np.arange(1,2000001,dtype=float)**(-s));constant=c**(1/s)/s*gamma(1-1/s)
    sampleD=10000.;finite_sum=np.sum(p*np.exp(-sampleD*p));integral=constant*sampleD**(-expected)
    assert abs(finite_sum/integral-1)<.1
    save(fig,'discrete_skills',{'tail_a':a,'skills':K,'predicted_data_exponent':expected,'example_D':sampleD,'finite_sum':float(finite_sum),'infinite_tail_integral_approx':float(integral),'finite_table_last_risk':finite})

def features():
    fig,ax=plt.subplots(1,3,figsize=(13.5,4.1));bs=np.linspace(.05,2,200);chi=np.where(bs<1,2*bs/(1+bs),bs)
    ax[0].plot(bs,bs,'--',label='Frozen kernel');ax[0].plot(bs,chi,label='Adaptive asymptotic branch');ax[0].axvline(1,color='gray',ls=':')
    ax[0].set(xlabel='Source exponent beta',ylabel='Time exponent chi',title='Improvement depends on the task');ax[0].legend(fontsize=8)
    beta=.5;ts=np.geomspace(1,1e9,220);clock_stats={}
    for g in [0.,.01,1.]:
        losses=[]
        for t in ts:
            # Positive root of R=(t+g*t^2*R)^(-beta), a scaling illustration.
            logR=brentq(lambda z:z+beta*np.log(t+g*t*t*np.exp(z)),-100,10)
            losses.append(np.exp(logR))
        losses=np.array(losses);sl=-np.gradient(np.log(losses),np.log(ts));clock_stats[str(g)]=float(sl[-1])
        ax[1].semilogx(ts,sl,label=f'gamma={g:g}')
    assert abs(clock_stats['0.0']-.5)<1e-6 and abs(clock_stats['1.0']-2/3)<.001
    ax[1].set(xlabel='Training time',ylabel='Local effective exponent',title='Self-consistent clock illustration');ax[1].legend(fontsize=8)
    alpha=2;nu=alpha*bs
    ax[2].plot(bs,nu*bs/(nu+bs),'--',label='Frozen: chi=beta')
    ax[2].plot(bs,nu*chi/(nu+chi),label='Adaptive chi')
    ax[2].set(xlabel='Source exponent beta',ylabel='Conditional compute exponent',title='Two-term cost model, alpha=2');ax[2].legend(fontsize=8)
    save(fig,'feature_learning',{'kind':'asymptotic formulas and self-consistent clock illustration, not full DMFT','clock_beta':beta,'last_local_slopes':clock_stats,'compute_assumptions':'R=t^-chi+N^-nu, C=N*t, nu=alpha*beta; hard-task branch is the relevant projected-model case','alpha':alpha})

def graphs():
    rng=np.random.default_rng(20261020);n=1600;m=6
    BA=np.zeros((n,n),dtype=bool);BA[:m+1,:m+1]=True;np.fill_diagonal(BA,False);degrees=BA.sum(0).astype(float)
    for v in range(m+1,n):
        targets=rng.choice(v,m,replace=False,p=degrees[:v]/degrees[:v].sum());BA[v,targets]=True;BA[targets,v]=True;degrees[targets]+=1;degrees[v]=m
    edges=int(BA.sum()//2);ii,jj=np.triu_indices(n,1);which=rng.choice(len(ii),edges,replace=False)
    ER=np.zeros_like(BA);ER[ii[which],jj[which]]=True;ER|=ER.T.copy()
    fig,ax=plt.subplots(1,3,figsize=(13.5,4.1));stats={}
    for name,mat in [('ER',ER),('BA',BA)]:
        k=mat.sum(1);assert np.all(k>0);assert connected_components(csr_matrix(mat),directed=False)[0]==1
        pi=k/k.sum();entropy=float(pi@np.log(k));unique=np.unique(k);ccdf=np.array([np.mean(k>=v) for v in unique]);rank=np.sort(pi)[::-1]
        ax[0].loglog(unique,ccdf,'o-',markersize=3,label=name);ax[1].loglog(np.arange(1,n+1),rank,label=name)
        order=np.argsort(k)[::-1];cum=np.cumsum(pi[order]*np.log(k[order]))/entropy
        ax[2].semilogx(np.arange(1,n+1)/n,cum,label=name)
        # Detailed balance and row-normalization checks directly from graph.
        assert np.allclose((pi[:,None]*mat/k[:,None]),(pi[:,None]*mat/k[:,None]).T)
        stats[name]={'nodes':n,'edges':edges,'mean_degree':float(k.mean()),'max_degree':int(k.max()),'entropy_rate_nats':entropy,'known_support_dimensions':2*edges-n}
    ax[0].set(xlabel='Degree k',ylabel='Fraction with degree >= k',title=f'Degree tails: V={n}, mean degree~12');ax[0].legend()
    ax[1].set(xlabel='Token rank',ylabel='Stationary visitation probability',title='Token rank is not degree-density rank');ax[1].legend()
    ax[2].set(xlabel='Fraction of highest-degree nodes',ylabel='Cumulative share of entropy rate',title='Weighted contributions to uncertainty');ax[2].legend()
    save(fig,'random_graphs',{'seed':20261020,'kind':'independent graph statistics, no transformer training','graphs':stats})

def fitting():
    fig,ax=plt.subplots(1,2,figsize=(10.8,4.1));x=np.geomspace(1,1e7,240);E=1.5;A=3.;a=.4;ex=A*x**(-a)
    ax[0].loglog(x,E+ex,label='Total loss E+A*x^-a');ax[0].loglog(x,ex,label='True excess loss');ax[0].axhline(E,color='gray',ls=':',label='True floor E')
    ax[0].set(xlabel='Resource x',ylabel='Loss',title='A floor flattens the total-loss slope');ax[0].legend(fontsize=8)
    for delta in [-.02,0.,.02]:
        den=ex-delta;valid=den>0;ax[1].semilogx(x[valid],a*ex[valid]/den[valid],label=f'Estimated floor E{delta:+.2f}')
    ax[1].set(xlabel='Resource x',ylabel='Local fitted exponent',ylim=(0,1.5),title='Floor error can mimic exponent changes');ax[1].legend(fontsize=8)
    save(fig,'fitting_floor',{'true_floor':E,'amplitude':A,'true_exponent':a,'floor_offsets':[-.02,0,.02],'kind':'analytic example, not language model data'})

if __name__=='__main__':
    statistical();skills();features();graphs();fitting()
    (ROOT/'numerical_checks.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k: 'passed' for k in report}))
