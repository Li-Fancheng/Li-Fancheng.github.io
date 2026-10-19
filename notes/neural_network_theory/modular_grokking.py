#!/usr/bin/env python3
"""Full-batch MLP on modular addition, using NumPy AdamW.

This is an experiment, not a guarantee that all settings grok. Keep the same
split and report losses as well as accuracy. Results include all hyperparameters.
Run: python modular_grokking.py --steps 40000 --seed 7
"""
import argparse
import csv
import json
import os
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR", "/tmp/neural-theory-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/neural-theory-cache")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def run(args):
    rng=np.random.default_rng(args.seed); p=args.modulus; n=p*p
    a,b=np.meshgrid(np.arange(p),np.arange(p),indexing="ij"); y=((a+b)%p).ravel()
    X=np.zeros((n,2*p)); X[np.arange(n),a.ravel()]=1; X[np.arange(n),p+b.ravel()]=1
    perm=rng.permutation(n); train=perm[:int(n*args.fraction)]; test=perm[len(train):]
    W=rng.normal(scale=1/np.sqrt(2*p),size=(2*p,args.width))
    V=rng.normal(scale=1/np.sqrt(args.width),size=(args.width,p))
    params=[W,V]; m=[np.zeros_like(z) for z in params]; v=[np.zeros_like(z) for z in params]
    activate=(lambda h:h*h) if args.activation=="square" else (lambda h:np.maximum(h,0))
    derivative=(lambda h:2*h) if args.activation=="square" else (lambda h:(h>0).astype(float))
    history=[]; first_train=None; first_test=None
    for step in range(args.steps+1):
        if step%100==0 or step==args.steps:
            H=X@W; logits=activate(H)@V; logits-=logits.max(axis=1,keepdims=True)
            logZ=np.log(np.exp(logits).sum(axis=1)); ce=logZ-logits[np.arange(n),y]
            correct=logits.argmax(axis=1)==y
            # DFT power in first-layer weights; DC removed. Diagnostic only.
            power=np.abs(np.fft.rfft(W[:p],axis=0))**2+np.abs(np.fft.rfft(W[p:],axis=0))**2
            power=power[1:].sum(axis=1); concentration=float(np.max(power)/max(power.sum(),1e-30))
            row=[step,float(ce[train].mean()),float(ce[test].mean()),float(correct[train].mean()),float(correct[test].mean()),float(sum(np.sum(z*z) for z in params)),concentration]
            history.append(row)
            if first_train is None and row[3]>=.99: first_train=step
            if first_test is None and row[4]>=.99: first_test=step
            if step%5000==0: print(f"step={step} train={row[3]:.3f} test={row[4]:.3f} CE={row[1]:.4g}/{row[2]:.4g}",flush=True)
        if step==args.steps: break
        H=X[train]@W; logits=activate(H)@V; logits-=logits.max(axis=1,keepdims=True)
        prob=np.exp(logits); prob/=prob.sum(axis=1,keepdims=True); prob[np.arange(len(train)),y[train]]-=1; prob/=len(train)
        grads=[X[train].T@((prob@V.T)*derivative(H)),activate(H).T@prob]
        for j,(z,g) in enumerate(zip(params,grads)):
            m[j]=.9*m[j]+.1*g; v[j]=.999*v[j]+.001*g*g
            z*=1-args.lr*args.decay
            z-=args.lr*(m[j]/(1-.9**(step+1)))/(np.sqrt(v[j]/(1-.999**(step+1)))+1e-8)
    out=Path(args.output) if args.output else Path(__file__).resolve().parent/"grokking"/"figures"
    out.mkdir(parents=True,exist_ok=True)
    with (out/"modular_history.csv").open("w") as f:
        writer=csv.writer(f); writer.writerow(["step","train_ce","test_ce","train_accuracy","test_accuracy","norm_squared","fourier_concentration"]); writer.writerows(history)
    below=[row[0] for row in history if row[4]<.99]
    sustained=(max(below)+100 if below else 0) if history[-1][4]>=.99 else None
    meta={**vars(args),"train_examples":len(train),"test_examples":len(test),"first_train_99":first_train,"first_test_99":first_test,"sustained_test_99_from":sustained,"final":history[-1]}
    (out/"modular_run.json").write_text(json.dumps(meta,indent=2)+"\n")
    np.savez_compressed(out/"modular_weights.npz",W=W,V=V,train=train,test=test)
    data=np.array(history); fig,ax=plt.subplots(1,2,figsize=(10.2,3.7),layout="constrained")
    for col,lab in [(3,"Train"),(4,"Test")]: ax[0].plot(data[:,0],data[:,col],label=lab)
    for col,lab in [(1,"Train"),(2,"Test")]: ax[1].semilogy(data[:,0],np.maximum(data[:,col],1e-12),label=lab)
    ax[0].set(xlabel="AdamW updates",ylabel="Accuracy",ylim=(0,1.02)); ax[1].set(xlabel="AdamW updates",ylabel="Cross entropy")
    for aa in ax: aa.legend(frameon=False); aa.grid(alpha=.2); aa.spines[["top","right"]].set_visible(False)
    fig.savefig(out/"modular_arithmetic.svg",metadata={"Date":None}); fig.savefig(out/"modular_arithmetic.png",dpi=220); plt.close(fig)
    print(json.dumps(meta,indent=2),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modulus",type=int,default=19); parser.add_argument("--width",type=int,default=128)
    parser.add_argument("--fraction",type=float,default=.7); parser.add_argument("--steps",type=int,default=40000)
    parser.add_argument("--lr",type=float,default=.003); parser.add_argument("--decay",type=float,default=1.)
    parser.add_argument("--seed",type=int,default=7); parser.add_argument("--output")
    parser.add_argument("--activation",choices=["square","relu"],default="square")
    args=parser.parse_args()
    if args.modulus<2 or not 0<args.fraction<1 or int(args.modulus**2*args.fraction)<1 or args.steps<1 or args.width<1 or args.lr<=0 or args.decay<0 or args.lr*args.decay>=1:
        parser.error("Need p>=2, nonempty split, positive width/steps/lr and 0<=lr*decay<1.")
    run(args)
