"""Fixed-dictionary vector BP versus the scalar annealed closure in Note 03.

This experiment uses the same planted soft fields as the scalar closure. It does
not simulate a Gaussian diffusion trajectory or estimate a critical exponent.
Run: python notes/hierarchical_data/bp_mean_field_comparison.py
"""
from pathlib import Path
from itertools import product
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent


def draw_rules(rng, v, m, depth, s=2):
    """One distinct ordered tuple per rule, with a shared dictionary per level."""
    result = []
    for _ in range(depth):
        code = rng.choice(v**s, size=v*m, replace=False).reshape(v, m)
        result.append(np.stack([(code//v**(s-1-j)) % v for j in range(s)], -1))
    return result


def draw_tree(rng, rules, samples):
    v, m, s = rules[0].shape
    root = rng.integers(v, size=samples)
    levels = [root[:, None]]
    for dictionary in rules[::-1]:
        state = levels[-1]
        chosen = rng.integers(m, size=state.shape)
        levels.append(dictionary[state, chosen].reshape(samples, -1))
    return root, levels[-1]


def normalize(x):
    mass = x.sum(-1, keepdims=True)
    if np.any(mass <= 0):
        raise ValueError('Zero-evidence observation; no posterior exists.')
    return x/mass


def vector_bp(rules, likelihood, downward=False):
    """likelihood shape: (field settings, samples, leaves, states).

    Upward messages are exact likelihood vectors up to a state-independent
    scale. A uniform root prior makes their normalized root the posterior.
    The optional downward pass returns exact node marginals for validation.
    """
    v, m, s = rules[0].shape
    upward = [normalize(likelihood)]
    for dictionary in rules:
        child = upward[-1]
        e, n, leaves, _ = child.shape
        child = child.reshape(e, n, leaves//s, s, v)
        weights = np.ones((e, n, leaves//s, v, m))
        for j in range(s):
            weights *= np.take(child[..., j, :], dictionary[..., j], axis=-1)
        upward.append(normalize(weights.mean(-1)))
    posterior = upward[-1][..., 0, :]
    if not downward:
        return posterior
    outside = [None]*len(upward)
    outside[-1] = np.ones_like(upward[-1])/v
    eye = np.eye(v)
    for ell in range(len(rules)-1, -1, -1):
        dictionary = rules[ell]
        child = upward[ell]
        e, n, leaves, _ = child.shape
        child = child.reshape(e, n, leaves//s, s, v)
        down = np.empty_like(child)
        for j in range(s):
            weights = np.broadcast_to(outside[ell+1][..., None], (e,n,leaves//s,v,m)).copy()/m
            for k in range(s):
                if k != j:
                    weights *= np.take(child[..., k, :], dictionary[..., k], axis=-1)
            down[..., j, :] = np.einsum('...ar,arb->...b', weights, eye[dictionary[..., j]])
        outside[ell] = normalize(down.reshape(e,n,leaves,v))
    return posterior, [normalize(u*d) for u,d in zip(upward,outside)]


def annealed(p, v, m, s=2):
    power = p**s
    return (power+(m-1)/(v**s-1)*(1-power))/(power+(m*v-1)/(v**s-1)*(1-power))


def check_bp():
    rng = np.random.default_rng(19063)
    v, m, depth = 4, 2, 2
    rules = draw_rules(rng,v,m,depth)
    likelihood = rng.uniform(.1,1,(1,1,4,v))
    root_bp, marg = vector_bp(rules,likelihood,True)
    root_exact = np.zeros(v)
    leaf_exact = np.zeros((4,v))
    # Independently enumerate all roots and all three internal rule choices.
    for y, r0, r1, r2 in product(range(v),range(m),range(m),range(m)):
        h = rules[1][y,r0]
        leaves = np.concatenate([rules[0][h[0],r1],rules[0][h[1],r2]])
        weight = np.prod(likelihood[0,0,np.arange(4),leaves])/(v*m**3)
        root_exact[y] += weight
        leaf_exact[np.arange(4),leaves] += weight
    z = root_exact.sum()
    errors = {'root_vs_enumeration':float(np.max(np.abs(root_bp[0,0]-root_exact/z))),
              'leaf_vs_enumeration':float(np.max(np.abs(marg[0][0,0]-leaf_exact/z)))}
    assert max(errors.values()) < 1e-12
    return errors


def crossing(x,y,target):
    indices = np.flatnonzero(y <= target)
    if len(indices)==0 or indices[0]==0:
        return None
    j = indices[0]
    return float(x[j-1]+(target-y[j-1])*(x[j]-x[j-1])/(y[j]-y[j-1]))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grammars',type=int,default=24)
    parser.add_argument('--samples',type=int,default=16)
    parser.add_argument('--seed',type=int,default=20261025)
    args=parser.parse_args()
    v,m,s=16,4,2
    depths=[2,4,6,8]
    epsilon=np.linspace(0,1,41)
    rng=np.random.default_rng(args.seed)
    checks=check_bp()
    observations={str(d):[] for d in depths}
    for g in range(args.grammars):
        rules=draw_rules(rng,v,m,max(depths))
        for depth in depths:
            y,x=draw_tree(rng,rules[:depth],args.samples)
            evidence=(epsilon[:,None,None,None]/v
                      +(1-epsilon[:,None,None,None])*np.eye(v)[x][None,...])
            root=vector_bp(rules[:depth],evidence)
            mass=root[:,np.arange(args.samples),y]
            observations[str(depth)].append(mass.T)
        if (g+1)%6==0:
            print(f'Finished {g+1}/{args.grammars} fixed grammars',flush=True)
    result={'seed':args.seed,'v':v,'m':m,'s':s,'depths':depths,
            'grammars':args.grammars,'reference_trees_per_grammar_and_depth':args.samples,
            'epsilon':epsilon.tolist(),
            'observable':'Posterior mass assigned to the planted root under deterministic reference-centered soft likelihoods; not a forward-backward channel overlap.',
            'band':'One sample standard deviation across grammar means (not a confidence interval).',
            'bp_checks':checks,'curves':{}}
    for depth in depths:
        raw=np.array(observations[str(depth)]) # grammar, tree, epsilon
        per_grammar=raw.mean(1)
        mean=per_grammar.mean(0)
        std=per_grammar.std(0,ddof=1)
        sem=std/np.sqrt(args.grammars)
        mf=1-epsilon+epsilon/v
        for _ in range(depth):mf=annealed(mf,v,m,s)
        q=(mean-1/v)/(1-1/v)
        qmf=(mf-1/v)/(1-1/v)
        error=mean-mf
        summary={'posterior_mass_mean':mean.tolist(),'grammar_mean_std':std.tolist(),
                 'grammar_mean_sem':sem.tolist(),'per_grammar_means':per_grammar.tolist(),
                 'annealed':mf.tolist(),'maximum_absolute_gap':float(np.max(np.abs(error))),
                 'rmse':float(np.sqrt(np.mean(error**2))),
                 'bp_half_height_epsilon':crossing(epsilon,q,.5),
                 'mf_half_height_epsilon':crossing(epsilon,qmf,.5),
                 'bp_10_90_width':crossing(epsilon,q,.1)-crossing(epsilon,q,.9),
                 'mf_10_90_width':crossing(epsilon,qmf,.1)-crossing(epsilon,qmf,.9),
                 'largest_grammar_std':float(std.max())}
        result['curves'][str(depth)]=summary
        assert np.allclose(raw[...,0],1,atol=1e-12)
        assert np.allclose(raw[...,-1],1/v,atol=1e-12)
        assert np.all((raw>=-1e-12)&(raw<=1+1e-12))
    result['bp_checks']['all_perfect_and_uniform_endpoints']='passed'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none',
                         'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,3,figsize=(13.2,4),constrained_layout=True)
    colors=['#4978a8','#309186','#ad782c','#8c5fa1']
    for depth,color in zip(depths,colors):
        row=result['curves'][str(depth)]
        mean=np.array(row['posterior_mass_mean']);std=np.array(row['grammar_mean_std'])
        mf=np.array(row['annealed'])
        axs[0].plot(epsilon,mean,color=color,label=f'L={depth}')
        axs[0].fill_between(epsilon,np.maximum(0,mean-std),np.minimum(1,mean+std),color=color,alpha=.12)
        axs[0].plot(epsilon,mf,'--',color=color,alpha=.8)
        axs[1].plot(epsilon,mean-mf,color=color)
    axs[0].set(xlabel='Soft-field parameter epsilon',ylabel='Planted-root posterior mass',title='A  Fixed dictionaries (solid) / closure (dashed)')
    axs[0].axhline(1/v,color='#888888',ls=':',lw=.8)
    axs[0].legend(frameon=False,ncol=2)
    axs[1].axhline(0,color='#888888',lw=.8)
    axs[1].set(xlabel='Soft-field parameter epsilon',ylabel='Vector BP mean - annealed closure',title='B  Closure error')
    axs[2].plot(depths,[result['curves'][str(d)]['bp_10_90_width'] for d in depths],'o-',color='#2563a6',label='Vector BP')
    axs[2].plot(depths,[result['curves'][str(d)]['mf_10_90_width'] for d in depths],'s--',color='#b56c36',label='Annealed closure')
    axs[2].set(xlabel='Tree depth L',ylabel='10%-90% crossover width',title='C  Finite-depth sharpening')
    axs[2].legend(frameon=False)
    for ax in axs:ax.grid(alpha=.16)
    out=ROOT/'figures'
    out.mkdir(exist_ok=True)
    fig.savefig(out/'bp_mean_field.svg',bbox_inches='tight')
    fig.savefig(out/'bp_mean_field.png',bbox_inches='tight',dpi=170)
    plt.close(fig)
    (out/'bp_mean_field.json').write_text(json.dumps(result,indent=2)+'\n')
    from PIL import Image
    w,h=Image.open(out/'bp_mean_field.png').size
    policy_path=ROOT.parents[1]/'data/hierarchical_data_figures.json'
    policy=json.loads(policy_path.read_text())
    policy['figures/bp_mean_field']={'width':1100,'kind':'original-bp-comparison','intrinsic_width':w,'intrinsic_height':h}
    policy_path.write_text(json.dumps(policy,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({str(d):{key:row[key] for key in ['maximum_absolute_gap','bp_half_height_epsilon','mf_half_height_epsilon','bp_10_90_width','mf_10_90_width','largest_grammar_std']} for d in depths for row in [result['curves'][str(d)]]},indent=2),flush=True)


if __name__=='__main__':main()
