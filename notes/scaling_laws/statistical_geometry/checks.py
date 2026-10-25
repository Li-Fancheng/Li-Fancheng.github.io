"""Independent checks of the worked teacher-student examples.

Run from any directory with the parent topic's NumPy dependency installed.
No network training or paper experiment reproduction is performed.
"""

from pathlib import Path
import json

import numpy as np


ROOT = Path(__file__).resolve().parent
SEED = 20261025
rng = np.random.default_rng(SEED)
report = {"seed": SEED, "squared_loss_factor": 1.0}

# Compare the Gaussian risk prediction with repeated actual least-squares fits.
# Independent omitted coordinates act as extra noise for the visible student.
P, S, D, trials = 4, 7, 60, 12000
sigma2, tau2 = 0.3, 0.7
teacher_visible = np.array([0.4, -0.3, 0.2, 0.1])
teacher_omitted = np.sqrt(tau2 / (S - P)) * np.ones(S - P)
risks = []
for _ in range(trials):
    features = rng.normal(size=(D, S))
    labels = (
        features[:, :P] @ teacher_visible
        + features[:, P:] @ teacher_omitted
        + np.sqrt(sigma2) * rng.normal(size=D)
    )
    fitted = np.linalg.lstsq(features[:, :P], labels, rcond=None)[0]
    risks.append(tau2 + np.sum((fitted - teacher_visible) ** 2))
predicted = tau2 + (sigma2 + tau2) * P / (D - P - 1)
mean = float(np.mean(risks))
se = float(np.std(risks, ddof=1) / np.sqrt(trials))
z_score = abs(mean - predicted) / se
assert z_score < 5
report["misspecified_gaussian_risk"] = {
    "P": P, "S": S, "D": D, "trials": trials,
    "sigma_squared": sigma2, "omitted_signal_energy": tau2,
    "theory": predicted, "measured_mean": mean,
    "standard_error": se, "absolute_z_score": z_score,
    "test_metric": "Excess squared risk, evaluated with exact test covariance I_S",
}

# Full-rank noiseless recovery is a useful counterexample to a nonzero 1/D tail.
design = rng.normal(size=(D, P))
fit = np.linalg.lstsq(design, design @ teacher_visible, rcond=None)[0]
recovery_error = float(np.linalg.norm(fit - teacher_visible))
assert recovery_error < 1e-12
report["noiseless_matched_parameter_error"] = recovery_error

# Check the population projection against direct residuals. The finite set of
# test features below has exactly the requested second moment C.
S, P = 9, 4
matrix = rng.normal(size=(S, S))
C = matrix @ matrix.T / S + np.eye(S)
A = rng.normal(size=(P, S))
omega = rng.normal(size=S) / np.sqrt(S)
theta = np.linalg.solve(A @ C @ A.T, A @ C @ omega)
residual_matrix = C - C @ A.T @ np.linalg.solve(A @ C @ A.T, A @ C)
test_features = np.sqrt(S) * np.linalg.cholesky(C).T
direct = np.mean((test_features @ (omega - A.T @ theta)) ** 2)
quadratic = omega @ residual_matrix @ omega
projection_error = float(abs(direct - quadratic))
assert projection_error < 1e-12

# Use a deterministic set of teachers with exactly E[omega omega.T] = I_S/S
# to check the teacher average without Monte Carlo uncertainty.
teachers = np.eye(S)
coefficients = np.linalg.solve(A @ C @ A.T, A @ C @ teachers)
residuals = teachers - A.T @ coefficients
teacher_mean = np.mean(np.sum(residuals * (C @ residuals), axis=0))
trace_error = float(abs(teacher_mean - np.trace(residual_matrix) / S))
assert trace_error < 1e-12
report["population_projection"] = {
    "fixed_teacher_residual_error": projection_error,
    "teacher_average_trace_error": trace_error,
}

# For A = I_S and zero initialization, compare the kernel residual variance
# against direct errors over the same exact-covariance tests and teachers.
train_features = rng.normal(size=(4, S)) @ np.linalg.cholesky(C).T
K = train_features @ train_features.T / S
cross_kernel = train_features @ test_features.T / S
predictions = cross_kernel.T @ np.linalg.solve(K, train_features @ teachers)
targets = test_features @ teachers
direct_risk = np.mean((predictions - targets) ** 2)
kernel_risk = np.mean(
    np.sum(test_features ** 2, axis=1) / S
    - np.sum(cross_kernel * np.linalg.solve(K, cross_kernel), axis=0)
)
kernel_error = float(abs(direct_risk - kernel_risk))
interpolation_error = float(np.max(np.abs(
    K @ np.linalg.solve(K, train_features @ teachers) - train_features @ teachers
)))
assert kernel_error < 1e-12 and interpolation_error < 1e-12
report["full_feature_kernel"] = {
    "direct_average_risk": float(direct_risk),
    "kernel_residual_risk": float(kernel_risk),
    "risk_identity_error": kernel_error,
    "training_interpolation_error": interpolation_error,
}

report["status"] = "passed"
(ROOT / "checks.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
