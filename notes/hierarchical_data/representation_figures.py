"""Original explanatory calculations for the expanded five-note series."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
from PIL import Image

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figures"
POLICY = ROOT.parents[1] / "data/hierarchical_data_figures.json"
policy = json.loads(POLICY.read_text())
checks = {}
plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":11,
    "axes.spines.top":False, "axes.spines.right":False,
    "axes.labelcolor":"#24364b", "text.color":"#24364b",
    "axes.edgecolor":"#a5afba", "svg.fonttype":"none",
    "figure.facecolor":"white", "axes.facecolor":"white"})
blue, orange, green, purple = "#2563a6", "#c16436", "#268575", "#8a5fa7"

def save(fig, name, data):
    fig.savefig(OUT/f"{name}.svg",bbox_inches="tight")
    fig.savefig(OUT/f"{name}.png",bbox_inches="tight",dpi=180)
    w,h=Image.open(OUT/f"{name}.png").size
    policy[f"figures/{name}"]={"width":1000,"kind":"original-explanatory",
        "intrinsic_width":w,"intrinsic_height":h}
    (OUT/f"{name}.json").write_text(json.dumps(data,indent=2)+"\n")
    plt.close(fig)

# Same observed sequence, different legal boundaries and root categories.
fig,axs=plt.subplots(1,2,figsize=(12,3.7),constrained_layout=True)
for ax,groups,names,label in zip(axs,[[[0,1],[2,3,4]],[[0,1,2],[3,4]]],[["A","B"],["C","D"]],[0,1]):
    ax.text(2,1.9,f"Y = {label}",ha="center",va="center",color="white",
            bbox={"boxstyle":"round,pad=.4","fc":purple,"ec":"none"})
    for inds,name,color in zip(groups,names,[blue,orange]):
        center=np.mean(inds)
        ax.plot([2,center],[1.8,1.1],color="#b6c4cf",zorder=0)
        ax.text(center,1,name,ha="center",va="center",color="white",
                bbox={"boxstyle":"round,pad=.4","fc":color,"ec":"none"})
        for i in inds: ax.plot([center,i],[.9,.15],color=color,alpha=.6,zorder=0)
    for i,char in enumerate("abcde"):
        ax.text(i,0,char,ha="center",va="center",fontsize=15,
                bbox={"boxstyle":"round,pad=.25","fc":"#edf2f7","ec":"none"})
    ax.set(xlim=(-.5,4.5),ylim=(-.4,2.25),title="A   Two + three" if label==0 else "B   Three + two")
    ax.axis("off")
save(fig,"parsing",{"observed":"abcde","partitions":[[2,3],[3,2]],"posterior_root":[.5,.5]})

# Geometry stays identical under a global rotation; fixed-coordinate transfer does not.
fig,axs=plt.subplots(1,2,figsize=(12,4.3),constrained_layout=True)
t=np.linspace(0,2*np.pi,300)
axs[0].plot(np.cos(t),np.sin(t),color="#b6c4cf")
points=np.array([[1,0],[0,1]])
axs[0].scatter(points[:,0],points[:,1],c=[blue,orange],s=90,zorder=3)
axs[0].annotate("same item in A",xy=(1,0),xytext=(.4,-.4),arrowprops={"arrowstyle":"->","color":blue},color=blue)
axs[0].annotate("same item in B",xy=(0,1),xytext=(-1.1,1.18),arrowprops={"arrowstyle":"->","color":orange},color=orange)
axs[0].axvline(0,color="#9ba9b5",ls="--")
axs[0].arrow(0,0,.75,0,width=.012,head_width=.07,color=green,length_includes_head=True)
axs[0].text(.25,.13,"fixed readout w",color=green)
axs[0].set(aspect="equal",xlim=(-1.4,1.4),ylim=(-1.35,1.4),title="A   Coordinates rotate; the circle does not")
axs[0].axis("off")
beta=np.linspace(0,np.pi,301)
accuracy=1-beta/np.pi
axs[1].plot(beta/np.pi,np.ones_like(beta),color=blue,label="within-space neighborhood agreement")
axs[1].plot(beta/np.pi,accuracy,color=orange,label="fixed-readout transfer accuracy")
axs[1].axhline(.5,color="#9ba9b5",ls=":")
axs[1].set(xlabel=r"relative rotation $\beta/\pi$",ylabel="score",ylim=(-.04,1.1),title="B   Similarity and alignment separate")
axs[1].legend(frameon=False,fontsize=9,loc="lower left")
axs[1].grid(alpha=.18)
# Independent quadrature checks the analytical transfer-accuracy curve.
angles=(np.arange(20000)+.371)*2*np.pi/20000
err=max(abs(np.mean((np.cos(angles)>=0)==(np.cos(angles+b)>=0))-(1-b/np.pi)) for b in beta[::10])
assert err<2e-4
checks["rotation_accuracy_quadrature_max_error"]=float(err)
save(fig,"geometry_alignment",{"rotation_over_pi":(beta/np.pi).tolist(),"transfer_accuracy":accuracy.tolist(),"neighborhood_agreement":np.ones_like(beta).tolist()})

# Nonlinear semantic computation with language-specific surface signs.
rng=np.random.default_rng(20261021)
mix=np.linalg.qr(rng.normal(size=(6,6)))[0]
def paired(n):
    z=rng.normal(size=(n,2))
    us=[]
    for _ in range(2):
        signs=rng.choice([-1.,1.],size=(n,2))
        us.append(np.column_stack([signs[:,0],z[:,0]*signs[:,0],signs[:,1],z[:,1]*signs[:,1]]))
    return z,us
ztrain,utrain=paired(3000)
ztest,utest=paired(600)
def hidden(u,z): return np.column_stack([3*u,z])@mix
x=np.vstack([np.column_stack([u,np.ones(len(u))]) for u in utrain])
y=np.vstack([hidden(u,ztrain) for u in utrain])
reg=np.diag([1e-5]*4+[0.])
coef=np.linalg.solve(x.T@x+reg,x.T@y)
raw=[hidden(u,ztest) for u in utest]
res=[h-np.column_stack([u,np.ones(len(u))])@coef for u,h in zip(utest,raw)]
delta=[h-np.pad(u,((0,0),(0,2))) for u,h in zip(utest,raw)]
def metrics(a,b,k=10):
    da=cdist(a,a,"sqeuclidean");db=cdist(b,b,"sqeuclidean")
    np.fill_diagonal(da,np.inf);np.fill_diagonal(db,np.inf)
    ia=np.argsort(da,axis=1);ib=np.argsort(db,axis=1)
    knn=np.mean([len(set(x[:k])&set(y[:k]))/k for x,y in zip(ia,ib)])
    ra=np.argsort(ia,axis=1)+1;rb=np.argsort(ib,axis=1)+1
    n=len(a);ii=(np.sum(rb[np.arange(n),ia[:,0]])+np.sum(ra[np.arange(n),ib[:,0]]))/n**2
    return float(knn),float(ii)
names=["Raw hidden state","Direct difference","Affine residual"]
values=np.array([metrics(*pair) for pair in [raw,delta,res]])
fig,axs=plt.subplots(1,2,figsize=(12,4.3),constrained_layout=True)
labels=["Raw","Direct\ndifference","Affine\nresidual"]
for j,(ax,title,ylabel) in enumerate(zip(axs,["A   Shared neighborhoods (higher is better)","B   Information imbalance (lower is better)"],["mKNN (k = 10)","symmetrized II"])):
    ax.bar(np.arange(3),values[:,j],color=[blue,orange,green],width=.6)
    ax.set(xticks=np.arange(3),xticklabels=labels,ylabel=ylabel,title=title)
    ax.grid(axis="y",alpha=.18)
    if j==0:ax.axhline(10/599,color="#9ba9b5",ls=":",label="independent-neighbor baseline")
    else:ax.axhline(1,color="#9ba9b5",ls=":",label="independent-rank baseline")
    ax.legend(frameon=True,facecolor="white",edgecolor="none",framealpha=1,fontsize=9)
checks["novelty_fit_sentences"]=6000
checks["novelty_test_pairs"]=600
checks["novelty_mKNN_raw"]=float(values[0,0])
checks["novelty_mKNN_residual"]=float(values[2,0])
checks["novelty_II_raw"]=float(values[0,1])
checks["novelty_II_residual"]=float(values[2,1])
checks["semantic_product_max_error"]=float(max(np.max(np.abs(np.column_stack([u[:,0]*u[:,1],u[:,2]*u[:,3]])-ztest)) for u in utest))
assert values[2,0]>values[0,0] and values[2,1]<values[0,1]
save(fig,"novelty",{"metric":"Euclidean neighborhoods on held-out pairs; fixed generative parameters, not LLM results","names":names,"k":10,"N":600,"mKNN":values[:,0].tolist(),"II":values[:,1].tolist()})

# A positive hierarchical block kernel, whose eigenspaces are exactly known.
idx=np.arange(8)
a0,a1,a2,a3=.1,.4,.3,.2
K=a0*np.ones((8,8))+a1*(idx[:,None]//4==idx[None,:]//4)+a2*(idx[:,None]//2==idx[None,:]//2)+a3*np.eye(8)
basis=[np.ones(8)/np.sqrt(8),np.r_[np.ones(4),-np.ones(4)]/np.sqrt(8)]
for start in [0,4]:
    row=np.zeros(8);row[start:start+2]=.5;row[start+2:start+4]=-.5;basis.append(row)
for start in range(0,8,2):
    row=np.zeros(8);row[start]=1/np.sqrt(2);row[start+1]=-1/np.sqrt(2);basis.append(row)
H=np.array(basis)
eigs=np.array([8*a0+4*a1+2*a2+a3,4*a1+2*a2+a3]+[2*a2+a3]*2+[a3]*4)
error=float(np.max(np.abs(H@K@H.T-np.diag(eigs))))
assert error<1e-12
checks["hierarchical_kernel_diagonalization_error"]=error
checks["hierarchical_kernel_min_eigenvalue"]=float(np.linalg.eigvalsh(K).min())
fig,axs=plt.subplots(1,2,figsize=(12,4.4),constrained_layout=True)
im=axs[0].imshow(K,cmap="Blues",vmin=0,vmax=1)
axs[0].set(title="A   Shared ancestors induce a block kernel",xlabel="leaf concept",ylabel="leaf concept",xticks=idx,yticks=idx)
fig.colorbar(im,ax=axs[0],shrink=.8)
im=axs[1].imshow(H,cmap="RdBu_r",vmin=-.71,vmax=.71,aspect="auto")
axs[1].set(title="B   Coarse-to-fine eigendirections",xlabel="leaf concept",xticks=idx,yticks=idx,
           yticklabels=[f"common  ({eigs[0]:.1f})",f"coarse  ({eigs[1]:.1f})"]+[f"middle  ({eigs[2]:.1f})"]*2+[f"fine  ({eigs[4]:.1f})"]*4)
fig.colorbar(im,ax=axs[1],shrink=.8)
save(fig,"hierarchy_spectrum",{"coefficients":[a0,a1,a2,a3],"K":K.tolist(),"basis":H.tolist(),"eigenvalues":eigs.tolist()})

# Numerical checks for the hand-worked BP example and the small grammar.
weights=np.array([.6*.1,.3*.2,.3*.7,.1*.2]);posterior=weights/weights.sum()
checks["worked_BP_joint_error"]=float(np.max(np.abs(posterior-np.array([6,6,21,2])/35)))
assert checks["worked_BP_joint_error"]<1e-12
POLICY.write_text(json.dumps(policy,indent=2)+"\n")
(ROOT/"representation_checks.json").write_text(json.dumps(checks,indent=2)+"\n")
print(json.dumps(checks,indent=2))
