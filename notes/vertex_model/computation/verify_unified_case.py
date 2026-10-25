"""Independent directional derivatives, embedding checks and saved-data audits."""

from pathlib import Path
import json
import numpy as np
from unified_vertex import Mesh, prepare, self_intersections, stress

HERE = Path(__file__).resolve().parent


def main():
    data = json.loads((HERE / "unified_vertex_data.json").read_text())
    checks = json.loads((HERE / "verification.json").read_text())
    assert checks["max_force_norm"] < 1e-8
    assert checks["max_hessian_decomposition_error"] < 1e-9
    assert checks["max_shear_modulus_difference"] < 2e-6
    assert checks["max_path_gradient"] < 1e-5
    assert checks["max_collapsed_branch_mismatch"] < 1e-10
    assert checks["min_path_nonselected_edge"] > 0.01
    assert checks["dynamic_intersection_samples"] == 0
    assert checks["min_dynamic_area"] > 0.2
    matrices = np.load(HERE / "unified_vertex_matrices.npz")
    mesh, _, _ = prepare()
    x, H = matrices["vertices"], matrices["H"]
    rng = np.random.default_rng(251025)
    errors = []
    for _ in range(5):
        direction = rng.normal(size=x.shape)
        direction /= np.linalg.norm(direction)
        eps = 1e-5
        gp = mesh.energy_grad(x + eps * direction, 3.75)[1].ravel()
        gm = mesh.energy_grad(x - eps * direction, 3.75)[1].ravel()
        exact = H @ direction.ravel()
        errors.append(
            np.linalg.norm((gp - gm) / (2 * eps) - exact) / np.linalg.norm(exact)
        )
    assert max(errors) < 1e-7
    eps, gamma = 1e-6, 0.07
    fd = (
        mesh.energy_grad(x, 3.75, gamma + eps)[0]
        - mesh.energy_grad(x, 3.75, gamma - eps)[0]
    ) / (2 * eps * 16)
    stress_error = abs(fd - stress(mesh, x, 3.75, gamma))
    assert stress_error < 1e-8
    area_error = 0.0
    intersections = 0
    for path in data["paths"]:
        for i, points in enumerate(path["coordinates"]):
            points = np.asarray(points)
            cells = data["preparation"]["cells"] if i < 35 else path["new_cells"]
            m = Mesh(cells, points)
            area_error = max(area_error, abs(m.geom(points)[-2].sum() - 16))
            intersections += int(self_intersections(m, points, 0.0))
    assert area_error < 1e-10
    assert intersections == 0
    for row in data["static"]:
        assert row["G_affine"] + 1e-8 >= row["G_relaxed"]
    for row in data["dynamics"]["active"]:
        assert np.isfinite(row["response"]).all()
        assert np.all(np.diff(row["events"]) >= 0)
    checks["independent_audit"] = dict(
        hessian_directional_max_relative_error=float(max(errors)),
        stress_virtual_work_absolute_error=float(stress_error),
        path_total_area_max_error=float(area_error),
        path_self_intersections=intersections,
        static_modulus_ordering=True,
        dynamic_arrays_finite=True,
    )
    (HERE / "verification.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
