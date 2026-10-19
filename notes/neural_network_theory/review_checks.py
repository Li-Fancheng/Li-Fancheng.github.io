"""Independent finite-dimensional checks for the expanded review notes.
Run from any working directory; no training or downloaded data is needed.
"""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import solve_continuous_lyapunov
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parent

def main():
    rng = np.random.default_rng(20261019)
    report = {}
    # Check deleted-sample normal equations directly, holding D*lambda fixed.
    X = rng.normal(size=(19, 7))
    y = rng.normal(size=19)
    penalty = 0.7
    H = X.T @ X + penalty * np.eye(7)
    S = X @ np.linalg.solve(H, X.T)
    predicted = (y - S @ y) / (1 - np.diag(S))
    direct = []
    for i in range(len(y)):
        keep = np.arange(len(y)) != i
        w = np.linalg.solve(X[keep].T @ X[keep] + penalty * np.eye(7), X[keep].T @ y[keep])
        direct.append(y[i] - X[i] @ w)
    report['ridge_leave_one_out_max_error'] = float(np.max(np.abs(predicted-direct)))
    assert report['ridge_leave_one_out_max_error'] < 1e-11

    # Marginal likelihood derivative: analytic trace vs finite differences.
    Q = rng.normal(size=(9, 5)); K = Q @ Q.T
    target = rng.normal(size=9)
    def evidence(a):
        A = a*K + 0.8*np.eye(9)
        return -0.5*(target @ np.linalg.solve(A,target) + np.linalg.slogdet(A)[1])
    a = 1.3; A = a*K + 0.8*np.eye(9)
    alpha = np.linalg.solve(A,target)
    analytic = 0.5*np.trace((np.outer(alpha,alpha)-np.linalg.inv(A)) @ K)
    eps = 1e-5
    fd = (evidence(a+eps)-evidence(a-eps))/(2*eps)
    report['GP_evidence_derivative_error'] = float(abs(analytic-fd))
    assert abs(analytic-fd) < 1e-8

    # Kernel mismatch initialization covariance: expand the affine map directly.
    Q = rng.normal(size=(8, 8)); C = Q @ Q.T
    Q = rng.normal(size=(8, 8)); theta = Q @ Q.T + np.eye(8)
    B = theta[5:,:5] @ np.linalg.inv(theta[:5,:5])
    L = np.concatenate((-B,np.eye(3)),axis=1)
    expanded = C[5:,5:] - B@C[:5,5:] - C[5:,:5]@B.T + B@C[:5,:5]@B.T
    err = float(np.max(np.abs(L@C@L.T-expanded)))
    report['NTK_initialization_covariance_error'] = err
    assert err < 1e-11

    # Bayesian posterior rank-one update vs a fresh inverse.
    Q = rng.normal(size=(6,6)); C = Q@Q.T+np.eye(6)
    z = rng.normal(size=6); noise = 0.4
    update = C - np.outer(C@z,C@z)/(noise+z@C@z)
    direct = np.linalg.inv(np.linalg.inv(C)+np.outer(z,z)/noise)
    report['posterior_rank_one_error'] = float(np.max(np.abs(update-direct)))
    assert report['posterior_rank_one_error'] < 1e-10

    # Matrix factorization balance: differentiate the two Gram matrices.
    W1 = rng.normal(size=(4,7)); W2 = rng.normal(size=(5,4))
    G = rng.normal(size=(5,7))
    d1 = -W2.T@G; d2 = -G@W1.T
    derivative = d2.T@W2 + W2.T@d2 - d1@W1.T - W1@d1.T
    report['deep_linear_Gram_balance_derivative'] = float(np.max(np.abs(derivative)))
    assert report['deep_linear_Gram_balance_derivative'] < 1e-11

    # OU: solve the matrix equation, not the component formula.
    H = np.diag([1.,4.]); A = np.array([[1.,.8],[.8,2.]])
    covariance = solve_continuous_lyapunov(H,2*A)
    report['OU_component_covariance_error'] = float(np.max(np.abs(covariance-np.array([[1.,.32],[.32,.5]]))))
    report['OU_mean_loss_minus_half_trace_diffusion'] = float(.5*np.trace(H@covariance)-.5*np.trace(A))
    assert report['OU_component_covariance_error'] < 1e-12
    assert abs(report['OU_mean_loss_minus_half_trace_diffusion']) < 1e-12

    # Group projection on the complete modular grid and convex CE inequality.
    p = 7; a,b,c = np.indices((p,p,p)); delta = (a+b-c)%p
    F = rng.normal(size=(p,p,p)); F -= F.mean(axis=2,keepdims=True)
    def project(f):
        curve = np.bincount(delta.ravel(),weights=f.ravel(),minlength=p)/(p*p)
        return curve[delta]
    PF = project(F)
    orthogonal_error = abs(np.sum(F*F)-np.sum(PF*PF)-np.sum((F-PF)**2))
    labels = (a[:,:,0]+b[:,:,0])%p
    def ce(f):
        return float(np.mean(logsumexp(f,axis=2)-np.take_along_axis(f,labels[:,:,None],axis=2)[:,:,0]))
    report['grokking_projection_Pythagoras_error'] = float(orthogonal_error)
    report['grokking_full_grid_cross_entropy'] = {'original':ce(F),'group_average':ce(PF)}
    assert np.max(np.abs(project(PF)-PF)) < 1e-12
    assert orthogonal_error < 1e-10
    assert ce(PF) <= ce(F) + 1e-12

    # Noisy perceptron formula vs a large independent sample.
    R = .7; delta_noise = .8; n = 600000
    v = rng.normal(size=n)
    u = R*v + np.sqrt(1-R*R)*rng.normal(size=n)
    label = np.sign(v+delta_noise*rng.normal(size=n))
    empirical = float(np.mean(np.sign(u)!=label))
    theory = float(np.arccos(R/np.sqrt(1+delta_noise**2))/np.pi)
    report['noisy_teacher_classification'] = {'theory':theory,'Monte_Carlo':empirical,'samples':n}
    assert abs(empirical-theory) < .002

    # Third-resource optimum: optimize independently and check marginal balance.
    alpha,beta,G,kappa,K,budget = .6,.4,2.,.8,1.7,1e8
    r = alpha*beta/(alpha+beta)
    def objective(T):
        return K*budget**(-r)*T**r+G*np.exp(-kappa*T)
    opt = minimize_scalar(objective,bounds=(1.,80.),method='bounded',options={'xatol':1e-10})
    T = opt.x
    lhs = r*K*budget**(-r)*T**(r-1)
    rhs = kappa*G*np.exp(-kappa*T)
    report['compute_training_pass_optimum'] = {'T':float(T),'relative_marginal_error':float(abs(lhs-rhs)/lhs)}
    assert abs(lhs-rhs)/lhs < 1e-6
    # Joint resource risk: direct coefficient ODE and independent noisy observations.
    indices = np.arange(1,26,dtype=float)
    rates = indices**-1.4
    teacher = (-1.)**indices * indices**-1.1
    N, D, sigma, time = 7, 90, .3, 8.
    observations = teacher[:N] + sigma/np.sqrt(D)*rng.normal(size=N)
    flow = solve_ivp(lambda t,w: rates[:N]*(observations-w), (0,time), np.zeros(N),
                     rtol=1e-11, atol=1e-12)
    filt = -np.expm1(-rates[:N]*time)
    ode_error = float(np.max(np.abs(flow.y[:,-1]-filt*observations)))
    expected = (np.sum(teacher[:N]**2*np.exp(-2*rates[:N]*time))
                + np.sum(teacher[N:]**2) + sigma**2/D*np.sum(filt**2))
    draws = teacher[:N] + sigma/np.sqrt(D)*rng.normal(size=(60000,N))
    empirical = float(np.mean(np.sum((draws*filt-teacher[:N])**2,axis=1))
                      + np.sum(teacher[N:]**2))
    report['joint_time_data_model_risk'] = {
        'coefficient_ODE_max_error':ode_error,'analytic_risk':float(expected),
        'Monte_Carlo_risk':empirical,'samples':60000}
    assert ode_error < 1e-9
    assert abs(empirical-expected)/expected < .01
    report['status'] = 'passed'
    (ROOT/'review_diagnostics.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__ == '__main__':
    main()
