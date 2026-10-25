"""Independent checks for the empirical entry and cross-mechanism discussion.

These are low-dimensional formula checks; no neural network is trained.
"""
from pathlib import Path
from math import comb
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma, zeta

report = {}
# A finite continuous rate measure, integrated independently of its asymptote.
rows = []
for nu in [0.3, 0.5, 2.0]:
    time = 50.0
    risk = quad(lambda rate: np.exp(-time * rate), 0.0, 1.0,
                weight='alg', wvar=(nu - 1, 0), epsabs=1e-13)[0] * nu
    predicted = gamma(nu + 1) * time ** (-nu)
    error = abs(risk / predicted - 1)
    assert error < 1e-9
    rows.append({'nu': nu, 'resource': time, 'risk': risk,
                 'asymptote': predicted, 'relative_error': error})
report['cumulative_measure_constant'] = rows

# Explicit finite Zipf sums and rigorous bounds on the omitted test mass.
# The asymptotic law refers to the infinite table; the finite cutoff is recorded.
rows = []
ranks = np.arange(1, 1000001, dtype=float)
for s in [1.6, 2.0, 3.0]:
    c = 1.0 / zeta(s)
    probabilities = c * ranks ** (-s)
    D = 10000.0
    risk = float(probabilities @ np.exp(-D * probabilities))
    predicted = c ** (1 / s) / s * gamma(1 - 1 / s) * D ** (-(s - 1) / s)
    omitted_mass_bound = c / (s - 1) * len(ranks) ** (1 - s)
    # The exact infinite sum differs from this finite sum by at most this mass.
    relative_error = abs(risk / predicted - 1)
    assert relative_error < 0.02
    assert abs(risk - predicted) < omitted_mass_bound + 1e-10
    rows.append({'s': s, 'D': D, 'skills': len(ranks), 'finite_sum': risk,
                 'infinite_asymptote': predicted,
                 'omitted_mass_upper_bound': omitted_mass_bound,
                 'relative_error': relative_error})
report['zipf_poisson_constant'] = rows

# Fixed-D binomial coverage is exactly a mixture at rate -log(1-p).
p = np.array([0.6, 0.3, 0.1])
coverage_errors = []
for D in [1, 3, 10, 100]:
    direct = np.sum(p * (1 - p) ** D)
    transformed = np.sum(p * np.exp(D * np.log1p(-p)))
    coverage_errors.append(float(abs(direct - transformed)))
assert max(coverage_errors) < 1e-14
report['fixed_count_rate_change_max_error'] = max(coverage_errors)
# Two-evidence survival has negative second derivative near zero, precluding
# a fixed positive-measure exponential mixture in general.
rate, resource = 0.2, 1.0
curvature = rate ** 2 * (rate * resource - 1) * np.exp(-rate * resource)
assert curvature < 0
report['two_evidence_negative_curvature'] = float(curvature)

# Rosenfeld's complex envelope and the real rewrite agree exactly numerically.
x = np.geomspace(1e-6, 1e6, 1000)
eta, epsilon0, floor_parameter = 0.7, 0.9, 0.03
complex_form = epsilon0 * np.abs(x / (x - 1j * eta))
real_form = epsilon0 * x / np.sqrt(x * x + eta * eta)
envelope_error = float(np.max(np.abs(complex_form - real_form)))
assert envelope_error < 1e-14 and np.all(np.diff(real_form) > 0)
report['rosenfeld_envelope'] = {'complex_vs_real_max_error': envelope_error,
    'internal_offset': floor_parameter,
    'actual_floor': float(epsilon0 * floor_parameter / np.hypot(floor_parameter, eta)),
    'small_scale_limit': epsilon0}

# Exact enumeration of a two-point gradient distribution checks the expected
# quadratic loss and the optimized single-step improvement, without sampling.
H = np.diag([1.0, 4.0]); theta = np.array([1.0, 0.5])
g = H @ theta; perturbation = np.array([0.4, 1.2])
Sigma = np.outer(perturbation, perturbation)
loss = lambda point: float(point @ H @ point / 2)
curvature_signal = float(g @ H @ g)
curvature_noise = float(np.trace(H @ Sigma))
max_improvement = float((g @ g) ** 2 / (2 * curvature_signal))
rows = []
for B in [1, 2, 4, 8]:
    step = float((g @ g) / (curvature_signal + curvature_noise / B))
    exact = sum(comb(B, k) / 2 ** B * loss(theta - step * (g + (2*k/B - 1) * perturbation))
                for k in range(B + 1))
    expected = loss(theta) - step * (g @ g) + step ** 2 / 2 * (curvature_signal + curvature_noise / B)
    predicted_improvement = max_improvement / (1 + curvature_noise / curvature_signal / B)
    error = max(abs(exact - expected), abs(loss(theta) - exact - predicted_improvement))
    assert error < 1e-13
    rows.append({'batch': B, 'exact_expected_loss': exact, 'formula_loss': float(expected),
                 'predicted_improvement': predicted_improvement, 'max_error': float(error)})
report['local_batch_quadratic'] = {'noise_scale': curvature_noise / curvature_signal,
                                  'enumerated_checks': rows}
report['scope'] = 'Independent formula checks, not empirical-paper or neural-network training reproduction.'
report['status'] = 'passed'
Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
