"""Worked-example figures and independent checks for the completed notes."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import quad, solve_ivp
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent
FIG = ROOT / 'figures'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.titleweight': 'bold', 'axes.titlesize': 12, 'svg.fonttype': 'path'})
BLUE, GREEN, GRAY = '#2563eb', '#0f766e', '#64748b'
report = {}

def save(fig, name):
    fig.savefig(FIG / (name + '.svg'), metadata={'Date': None}, facecolor='white')
    fig.savefig(FIG / (name + '.png'), dpi=200, facecolor='white')
    plt.close(fig)

# An independently fitted additive model recovers the entire table.
a = np.array([1., 2., 4.]); b = np.array([0., 3., 5.])
truth = a[:, None] + b[None, :]
observed = [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0)]
X = np.zeros((5, 6))
for row, (i, j) in enumerate(observed):
    X[row, i] = 1; X[row, 3+j] = 1
coeff = np.linalg.lstsq(X, np.array([truth[i, j] for i, j in observed]), rcond=None)[0]
fitted = coeff[:3, None] + coeff[None, 3:]
error = float(np.max(np.abs(fitted - truth)))
assert np.linalg.matrix_rank(X) == 5 and error < 1e-12
report['additive_composition'] = {'observations': len(observed), 'identifiable_rank': 5,
                                  'full_table': truth.tolist(), 'recovery_max_error': error}
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5), constrained_layout=True)
for ax, title, partial in [(axes[0], 'Observed entries: 5 of 9', True),
                            (axes[1], 'Recovered additive table', False)]:
    values = truth.copy()
    if partial:
        for i in range(3):
            for j in range(3):
                if (i, j) not in observed: values[i, j] = np.nan
    cmap = plt.get_cmap('Blues').copy(); cmap.set_bad('#e2e8f0')
    ax.imshow(values, vmin=0, vmax=10, cmap=cmap)
    for i in range(3):
        for j in range(3):
            text = '?' if np.isnan(values[i, j]) else str(int(values[i, j]))
            ax.text(j, i, text, ha='center', va='center', fontsize=18,
                    color='white' if not np.isnan(values[i,j]) and values[i,j] >= 6 else '#0f172a')
    ax.set(xticks=range(3), xticklabels=['relation 1', 'relation 2', 'relation 3'],
           yticks=range(3), yticklabels=['object 1', 'object 2', 'object 3'], title=title)
    ax.set_xticks(np.arange(-.5, 3, 1), minor=True); ax.set_yticks(np.arange(-.5, 3, 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=3); ax.tick_params(which='both', length=0)
K = np.arange(1, 11)
axes[2].plot(K, K*K, 'o-', color=BLUE, label='Independent: K²')
axes[2].plot(K, 2*K-1, 's-', color=GREEN, label='Additive: 2K − 1')
axes[2].set(xlabel='Objects / relations per type, K', ylabel='Independent quantities',
            title='Degrees of freedom')
axes[2].legend(fontsize=9); axes[2].grid(alpha=.2)
save(fig, 'skills_composition')

# Balanced factor dynamics: closed form against a direct two-variable ODE.
theta0 = .01
T = np.linspace(0, 100, 2001)
theta = 1 / (1 + (1/theta0-1)*np.exp(-2*T))
frozen = 1 - (1-theta0)*np.exp(-2*theta0*T)
sol = solve_ivp(lambda t,w: [w[1]*(1-w[0]*w[1]), w[0]*(1-w[0]*w[1])],
                (0., 100.), [np.sqrt(theta0)]*2, t_eval=T, rtol=1e-10, atol=1e-12)
error = float(np.max(np.abs(sol.y[0]*sol.y[1]-theta)))
assert sol.success and error < 1e-8
report['factor_dynamics'] = {'theta0': theta0, 'closed_form_vs_ode_max_error': error,
    'adaptive_half_target_time': float(np.log(1/theta0-1)/2),
    'frozen_half_target_time': float(np.log(2*(1-theta0))/(2*theta0)),
    'scope': 'One mode; late convergence is exponential, not a neural power law.'}
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5), constrained_layout=True)
axes[0].plot(T, frozen, color=BLUE, label='Frozen initial metric')
axes[0].plot(T, theta, color=GREEN, label='Balanced factor model')
axes[0].axhline(1, color=GRAY, ls=':', lw=1)
axes[0].set(xlabel='Training time', ylabel='Prediction coefficient θ', title='Same initial predictor')
axes[0].legend(fontsize=9)
short = np.linspace(0, 12, 700)
s_adapt = 1/(1+(1/theta0-1)*np.exp(-2*short))
s_frozen = 1-(1-theta0)*np.exp(-2*theta0*short)
axes[1].semilogy(short, (1-s_frozen)**2, color=BLUE)
axes[1].semilogy(short, (1-s_adapt)**2, color=GREEN)
axes[1].set(xlabel='Training time', ylabel='Excess squared risk', title='Single-mode squared risk')
axes[2].plot(short, np.full_like(short, 2*theta0), color=BLUE)
axes[2].plot(short, 2*s_adapt, color=GREEN)
axes[2].set(xlabel='Training time', ylabel='Effective metric', title='Effective update metric')
for ax in axes: ax.grid(alpha=.2)
save(fig, 'factor_dynamics')

# A non-bipartite four-state walk with exactly known entropy and support.
adj = np.zeros((4, 4))
edges = [(0,1), (1,2), (2,0), (0,3)]
for i,j in edges: adj[i,j] = adj[j,i] = 1
degree = adj.sum(axis=1); P = adj/degree[:,None]; pi = degree/degree.sum()
assert np.max(np.abs(pi @ P-pi)) < 1e-14
joint = pi[:,None]*P
assert np.max(np.abs(joint[adj>0]-1/8)) < 1e-14
entropy = float(sum(pi[i]*np.log(degree[i]) for i in range(4)))
q = np.array([[.1,.2,.3,.4],[.2,.1,.6,.1],[.4,.4,.1,.1],[.6,.1,.2,.1]])
errors = []
for i in range(4):
    keep = adj[i] > 0; legal_mass = q[i,keep].sum(); pv = P[i,keep]
    direct = np.sum(pv*np.log(pv/q[i,keep]))
    decomposed = -np.log(legal_mass)+np.sum(pv*np.log(pv/(q[i,keep]/legal_mass)))
    errors.append(float(abs(direct-decomposed)))
assert max(errors) < 1e-12
report['four_state_markov'] = {'P': P.tolist(), 'pi': pi.tolist(), 'entropy_nats': entropy,
    'uniform_predictor_excess': float(np.log(4)-entropy), 'known_support_dimension': int(sum(degree-1)),
    'cross_entropy_decomposition_max_error': max(errors)}
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5), constrained_layout=True)
pos = np.array([[-.15,.35],[-.85,-.65],[.45,-.65],[.95,.7]])
for i,j in edges: axes[0].plot(pos[[i,j],0],pos[[i,j],1],color='#94a3b8',lw=2,zorder=1)
axes[0].scatter(pos[:,0],pos[:,1],s=900,color=BLUE,zorder=2)
for i,(x,y) in enumerate(pos): axes[0].text(x,y,str(i+1),ha='center',va='center',color='white',fontsize=16,zorder=3)
axes[0].set(title='Triangle plus one leaf', xlim=(-1.2,1.25), ylim=(-1,.98)); axes[0].set_aspect('equal');axes[0].axis('off')
axes[1].imshow(P,cmap='Blues',vmin=0,vmax=1)
for i in range(4):
    for j in range(4): axes[1].text(j,i,f'{P[i,j]:.2f}' if P[i,j] else '0',ha='center',va='center',color='white' if P[i,j]>.6 else '#0f172a')
axes[1].set(xticks=range(4),xticklabels=range(1,5),yticks=range(4),yticklabels=range(1,5),
            xlabel='Next state',ylabel='Current state',title='Conditional transition matrix')
axes[2].bar(range(1,5),pi,color=[BLUE,GRAY,GRAY,GREEN])
axes[2].set(xticks=range(1,5),xlabel='State',ylabel='Stationary visit probability',ylim=(0,.45), title='Stationary visit probabilities')
for i,v in enumerate(pi): axes[2].text(i+1,v+.01,f'{v:.3f}',ha='center')
axes[2].grid(axis='y',alpha=.2)
save(fig, 'markov_example')

# Check the local interpolation constants by direct function-error integration.
M = 64; h = 1/M
const = sum(quad(lambda x,c=(i+.5)*h: (np.sin(2*np.pi*x)-np.sin(2*np.pi*c))**2,
                 i*h,(i+1)*h,epsabs=1e-14)[0] for i in range(M))
linear = 0.
for i in range(M):
    lo=i*h; hi=(i+1)*h; f0=np.sin(2*np.pi*lo); f1=np.sin(2*np.pi*hi)
    linear += quad(lambda x: (np.sin(2*np.pi*x)-(f0+(f1-f0)*(x-lo)/h))**2,lo,hi,epsabs=1e-14)[0]
relative = [abs(const/(np.pi**2*h*h/6)-1),abs(linear/(np.pi**4*h**4/15)-1)]
assert max(relative) < .002
report['interpolation_leading_constants'] = {'intervals': M,'relative_errors': [float(v) for v in relative]}

# Lifecycle cost optimum: derivative condition checked against direct minimization.
A=2.; B=3.; alpha=.4; beta=.3; C=1e6; k=6.; qcost=1.; records=[]
for Q in [0.,1e3,1e4]:
    upper=np.log(C/(qcost*Q))-.001 if Q else np.log(C/k)
    def risk(u):
        N=np.exp(u); D=C/(k*N)-qcost*Q/k
        return A*N**(-alpha)+B*D**(-beta)
    optimum=minimize_scalar(risk,bounds=(-10,upper),method='bounded',options={'xatol':1e-10})
    N=np.exp(optimum.x); D=C/(k*N)-qcost*Q/k
    lhs=alpha*A*N**(-alpha);rhs=beta*B*D**(-beta)*(1+qcost*Q/(k*D))
    rel=float(abs(lhs-rhs)/max(lhs,rhs));assert rel<2e-6
    records.append({'Q':Q,'N':float(N),'D':float(D),'risk':float(optimum.fun),'marginal_condition_relative_error':rel})
report['lifecycle_optima']=records
report['status']='passed'
(ROOT/'worked_examples_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
