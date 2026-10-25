"""One periodic Vertex network: static response, a constrained T1 path and dynamics.

All lengths use sqrt(A0), energy is 1/2 sum[(A-1)^2+(P-p0)^2].
Run this script for the saved data; plot_unified_case.py renders notes 03-05.
"""

from pathlib import Path
import json, argparse
from scipy.linalg import null_space
import numpy as np
from scipy.spatial import Voronoi
from scipy.optimize import minimize
import autograd.numpy as anp
from autograd import jacobian, hessian


class Mesh:
    def __init__(self, cells, x, L=4.0):
        self.cells = [list(c) for c in cells]
        self.L = L
        self.V = len(x)
        ids = []
        imgs = []
        nxt = []
        prev = []
        ci = []
        pairs = {}
        for c, vs in enumerate(cells):
            off = len(ids)
            points = x[vs]
            delta = np.roll(points, -1, axis=0) - points
            jumps = -np.rint(delta / L).astype(int)
            im = np.vstack([np.zeros(2, dtype=int), np.cumsum(jumps[:-1], axis=0)])
            assert np.all(jumps.sum(0) == 0)
            n = len(vs)
            ids.extend(vs)
            imgs.extend(im)
            nxt.extend(off + (np.arange(n) + 1) % n)
            prev.extend(off + (np.arange(n) - 1) % n)
            ci.extend([c] * n)
            for k, u in enumerate(vs):
                pairs.setdefault(tuple(sorted((u, vs[(k + 1) % n]))), []).append(c)
        assert all(len(v) == 2 for v in pairs.values())
        self.ids = np.array(ids)
        self.imgs = np.array(imgs)
        self.nxt = np.array(nxt)
        self.prev = np.array(prev)
        self.ci = np.array(ci)
        self.N = len(cells)
        self.edges = [(u, v, *adj) for (u, v), adj in sorted(pairs.items())]
        self.pairs = np.array([tuple(sorted((a, b))) for u, v, a, b in self.edges])
        lookup = {tuple(sorted((u, v))): k for k, (u, v, a, b) in enumerate(self.edges)}
        self.eids = np.array(
            [
                lookup[tuple(sorted((u, self.ids[j])))]
                for u, j in zip(self.ids, self.nxt)
            ]
        )

    def geom(self, x, gamma=0):
        F = np.array([[1, gamma], [0, 1.0]])
        p = (x[self.ids] + self.imgs * self.L) @ F.T
        pn = p[self.nxt]
        pp = p[self.prev]
        d = p - pn
        length = np.linalg.norm(d, axis=1)
        area = np.bincount(
            self.ci,
            weights=0.5 * (p[:, 0] * pn[:, 1] - p[:, 1] * pn[:, 0]),
            minlength=self.N,
        )
        perim = np.bincount(self.ci, weights=length, minlength=self.N)
        return F, p, pn, pp, d, length, area, perim

    def energy_grad(self, x, p0, gamma=0, noise=None):
        F, p, pn, pp, d, l, A, P = self.geom(x, gamma)
        ra = A - 1
        rp = P - p0
        E = 0.5 * np.sum(ra * ra + rp * rp)
        unit = d / np.maximum(l[:, None], 1e-12)
        gp = unit - unit[self.prev]
        ga = 0.5 * np.stack([pn[:, 1] - pp[:, 1], pp[:, 0] - pn[:, 0]], axis=1)
        gc = ra[self.ci, None] * ga + rp[self.ci, None] * gp
        if noise is not None:
            lam = noise[self.eids]
            E += 0.5 * np.sum(lam * l)
            gc += 0.5 * (lam[:, None] * unit - lam[self.prev, None] * unit[self.prev])
        result = np.zeros_like(x)
        np.add.at(result, self.ids, gc @ F)
        return E, result

    def minedge(self, x):
        return self.geom(x)[5].min()

    def relax(self, x, p0, gamma=0):
        fixed = x[0].copy()

        def fun(z):
            xx = np.vstack([fixed, z.reshape(-1, 2)])
            e, g = self.energy_grad(xx, p0, gamma)
            return e, g[1:].ravel()

        r = minimize(
            fun,
            x[1:].ravel(),
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": 4000, "gtol": 2e-10, "ftol": 1e-15, "maxls": 50},
        )
        xx = np.vstack([fixed, r.x.reshape(-1, 2)])
        return xx, r

    def derivatives(self, x, p0, gamma=0.0):
        ids, img, nxt, ci = self.ids, self.imgs, self.nxt, self.ci
        B = np.zeros((self.N, len(ids)))
        B[ci, np.arange(len(ids))] = 1

        def residual(z):
            p = anp.dot(
                z.reshape(-1, 2)[ids] + img * self.L,
                anp.array([[1.0, gamma], [0.0, 1.0]]).T,
            )
            pn = p[nxt]
            A = anp.dot(B, 0.5 * (p[:, 0] * pn[:, 1] - p[:, 1] * pn[:, 0]))
            P = anp.dot(B, anp.sqrt(anp.sum((p - pn) ** 2, axis=1)))
            return anp.concatenate([A - 1, P - p0])

        def energy(y):
            p = y[:-1].reshape(-1, 2)
            gamma = y[-1]
            F = anp.array([[1.0, gamma], [0.0, 1.0]])
            points = anp.dot(p[ids] + img * self.L, F.T)
            pn = points[nxt]
            A = anp.dot(B, 0.5 * (points[:, 0] * pn[:, 1] - points[:, 1] * pn[:, 0]))
            P = anp.dot(B, anp.sqrt(anp.sum((points - pn) ** 2, axis=1)))
            return 0.5 * anp.sum((A - 1) ** 2 + (P - p0) ** 2)

        f = np.asarray(residual(x.ravel()))
        J = np.asarray(jacobian(residual)(x.ravel()))
        HH = np.asarray(hessian(energy)(np.r_[x.ravel(), gamma]))
        H = HH[:-1, :-1]
        S = np.asarray(hessian(lambda z: anp.dot(f, residual(z)))(x.ravel()))
        return f, J, H, S, HH


def mesh_from_centers(centers, L=4.0):
    N = len(centers)
    shifts = np.array([(i, j) for i in [-1, 0, 1] for j in [-1, 0, 1]]) * L
    vor = Voronoi(np.concatenate([centers + d for d in shifts]))
    key_to_id = {}
    x = []
    cells = []
    for c in range(N):
        poly = vor.vertices[vor.regions[vor.point_region[4 * N + c]]]
        mean = poly.mean(0)
        poly = poly[np.argsort(np.arctan2(poly[:, 1] - mean[1], poly[:, 0] - mean[0]))]
        vs = []
        for r in poly:
            rr = r % L
            key = tuple(np.round(rr, 7))
            if key not in key_to_id:
                key_to_id[key] = len(x)
                x.append(rr)
            vs.append(key_to_id[key])
        cells.append(vs)
    x = np.array(x)
    return Mesh(cells, x, L), x


def flip(mesh, x, uv, newlength=0.05, gamma=0.0):
    u, v = uv
    cells = [list(c) for c in mesh.cells]
    shared = [i for i, c in enumerate(cells) if u in c and v in c]
    if len(shared) != 2:
        raise ValueError("shared")
    A, B = shared
    if cells[A][(cells[A].index(u) + 1) % len(cells[A])] != v:
        A, B = B, A
    assert cells[A][(cells[A].index(u) + 1) % len(cells[A])] == v
    C = next(i for i, c in enumerate(cells) if u in c and v not in c)
    D = next(i for i, c in enumerate(cells) if v in c and u not in c)
    if len(cells[A]) <= 3 or len(cells[B]) <= 3:
        raise ValueError("triangle")
    if any(C in adj and D in adj for adj in [e[2:] for e in mesh.edges]):
        raise ValueError("already neighbours")
    cells[A].remove(v)
    cells[B].remove(u)
    i = cells[C].index(u)
    cells[C].insert(i, v)
    i = cells[D].index(v)
    cells[D].insert(i, u)
    d = x[v] - x[u]
    shift = -np.rint(d / mesh.L) * mesh.L
    d += shift
    F = np.array([[1.0, gamma], [0.0, 1.0]])
    physical = F @ d
    n = np.linalg.solve(
        F, np.array([-physical[1], physical[0]]) / np.linalg.norm(physical)
    )
    mid = x[u] + d / 2
    xx = x.copy()
    xx[u] = mid + n * newlength / 2
    xx[v] = mid - n * newlength / 2 - shift
    return Mesh(cells, xx, mesh.L), xx, (A, B, C, D)


def constrain(mesh, x, p0, u, v, vector):
    free = [i for i in range(mesh.V) if i not in [0, v]]
    fixed = x.copy()

    def unpack(z):
        xx = fixed.copy()
        xx[free] = z.reshape(-1, 2)
        xx[v] = xx[u] + vector
        return xx

    def fun(z):
        xx = unpack(z)
        e, g = mesh.energy_grad(xx, p0)
        g[u] += g[v]
        return e, g[free].ravel()

    rr = minimize(
        fun,
        x[free].ravel(),
        jac=True,
        method="L-BFGS-B",
        options={"gtol": 1e-8, "ftol": 1e-14, "maxiter": 1500, "maxls": 40},
    )
    xx = unpack(rr.x)
    return xx, rr


def prepare():
    rng = np.random.default_rng(3)
    centers = np.array([(i + 0.5, j + 0.5) for i in range(4) for j in range(4)])
    centers += rng.uniform(-0.25, 0.25, (16, 2))
    mesh, x = mesh_from_centers(centers)
    return mesh, x, centers


def refine(mesh, x, p0, gamma=0.0):
    """Improve a fixed-topology minimum without moving its pinned vertex."""
    x, result = mesh.relax(x, p0, gamma)
    for iteration in range(12):
        energy, gradient = mesh.energy_grad(x, p0, gamma)
        if np.linalg.norm(gradient) < 2e-9:
            break
        f, J, H, S, HH = mesh.derivatives(x, p0, gamma)
        if p0 >= 3.82:
            step = -np.linalg.lstsq(J[:, 2:], f, rcond=1e-10)[0]
        else:
            w, U = np.linalg.eigh(H[2:, 2:])
            step = -U @ ((U.T @ gradient[1:].ravel()) / np.maximum(w, 1e-7))
        step *= min(1.0, 0.2 / max(np.linalg.norm(step), 1e-15))
        for a in 2.0 ** (-np.arange(18)):
            trial = x.copy()
            trial[1:] += a * step.reshape(-1, 2)
            if (
                mesh.energy_grad(trial, p0, gamma)[0] < energy
                and mesh.geom(trial, gamma)[-2].min() > 0.1
            ):
                x = trial
                break
        else:
            break
    return x


def projection(V):
    translations = np.zeros((2 * V, 2))
    translations[::2, 0] = 1.0
    translations[1::2, 1] = 1.0
    return null_space(translations.T)


def static_case(mesh, x0):
    Q = projection(mesh.V)
    out = []
    states = {}
    for p0 in [3.65, 3.70, 3.75, 3.79, 3.81, 3.84, 3.90]:
        x = refine(mesh, x0, p0)
        f, J, H, S, HH = mesh.derivatives(x, p0)
        hp = Q.T @ H @ Q
        w, U = np.linalg.eigh(hp)
        xi = Q.T @ HH[:-1, -1]
        mask = w > 1e-7
        born = HH[-1, -1] / 16
        relaxed = born - np.sum((U[:, mask].T @ xi) ** 2 / w[mask]) / 16
        direction = Q @ U[:, 0]
        states[str(p0)] = x.tolist()
        h = 1e-4
        xp = refine_at_shear(mesh, x, p0, h)
        xm = refine_at_shear(mesh, x, p0, -h)
        fd = (
            (
                mesh.energy_grad(xp, p0, h)[0]
                + mesh.energy_grad(xm, p0, -h)[0]
                - 2 * mesh.energy_grad(x, p0)[0]
            )
            / h**2
            / 16
        )
        record = dict(
            p0=p0,
            energy=mesh.energy_grad(x, p0)[0],
            force_norm=float(np.linalg.norm(mesh.energy_grad(x, p0)[1])),
            min_edge=float(mesh.minedge(x)),
            area_rms=float(np.sqrt(np.mean(f[:16] ** 2))),
            perimeter_rms=float(np.sqrt(np.mean(f[16:] ** 2))),
            h_min=float(w[0]),
            gram_min=float(np.linalg.eigvalsh(Q.T @ J.T @ J @ Q)[0]),
            s_min=float(np.linalg.eigvalsh(Q.T @ S @ Q)[0]),
            zero_modes=int((abs(w) < 1e-6).sum()),
            gram_rank=int(np.linalg.matrix_rank(J, tol=1e-7)),
            mode_gram=float(np.linalg.norm(J @ direction) ** 2),
            mode_stress=float(direction @ S @ direction),
            G_affine=float(born),
            G_relaxed=float(relaxed),
            G_finite_difference=float(fd),
            hessian_identity_max_error=float(np.max(abs(H - J.T @ J - S))),
            softest_mode=direction.reshape(-1, 2).tolist(),
        )
        out.append(record)
        print("static", p0, record["force_norm"], relaxed, fd, flush=True)
    return out, states


def refine_at_shear(mesh, x, p0, gamma):
    return refine(mesh, x, p0, gamma)


def relaxed_t1_path(mesh, x0, p0, edge=(2, 3), count=35):
    """A reproducible upper-bound path: fix edge vector, relax all remaining coordinates.

    This does not claim an unconstrained minimum-energy path or a universal T1 barrier.
    Both sides converge to the same collapsed-edge configuration at s=0.
    """
    x = refine(mesh, x0, p0)
    u, v = edge
    d = x[v] - x[u]
    shift = -np.rint(d / mesh.L) * mesh.L
    d += shift
    unit = d / np.linalg.norm(d)
    old_length = np.linalg.norm(d)
    coords = []
    energies = []
    lengths = []
    gradnorm = []
    xmin = x.copy()
    for ell in np.linspace(old_length, 0.0, count):
        xmin, r = constrain(mesh, xmin, p0, u, v, ell * unit - shift)
        coords.append(xmin.tolist())
        energies.append(float(r.fun))
        lengths.append(float(ell))
        gradnorm.append(float(np.linalg.norm(r.jac)))
    collapsed = xmin.copy()
    near = collapsed.copy()
    near[v] += unit * 0.002
    newmesh, near, adj = flip(mesh, near, edge, 0.02)
    # At zero edge length these two cell lists describe identical polygons.
    collapsed_old = mesh.energy_grad(collapsed, p0)[0]
    collapsed_new = newmesh.energy_grad(collapsed, p0)[0]
    endpoint = refine(newmesh, near, p0)
    dn = endpoint[v] - endpoint[u] + shift
    new_length = np.linalg.norm(dn)
    new_unit = dn / new_length
    xmin = collapsed.copy()
    for ell in np.linspace(0.0, new_length, count)[1:]:
        xmin, r = constrain(newmesh, xmin, p0, u, v, ell * new_unit - shift)
        coords.append(xmin.tolist())
        energies.append(float(r.fun))
        lengths.append(-float(ell))
        gradnorm.append(float(np.linalg.norm(r.jac)))
    A0 = mesh.energy_grad(x, p0)[0]
    # Audit polygons on both branches, allowing the selected edge to vanish.
    min_area = 1.0
    min_other = 1.0
    for i, xx in enumerate(coords):
        mm = mesh if i < count else newmesh
        geo = mm.geom(np.array(xx))
        min_area = min(min_area, float(geo[-2].min()))
        for j, k in enumerate(mm.nxt):
            if set((mm.ids[j], mm.ids[k])) != set(edge):
                min_other = min(min_other, float(geo[5][j]))
    return dict(
        p0=p0,
        edge=list(edge),
        adjacent_cells=list(adj),
        s=lengths,
        energy=energies,
        coordinates=coords,
        initial_energy=float(A0),
        path_barrier=float(max(energies) - A0),
        constrained_gradient_max=float(max(gradnorm)),
        collapsed_energy_mismatch=float(abs(collapsed_old - collapsed_new)),
        min_area=min_area,
        min_other_edge=min_other,
        new_cells=newmesh.cells,
    )


def stress(mesh, x, p0, gamma, noise=None):
    # The derivative is with respect to an affine physical shear, not a force on a cell center.
    F, p, pn, pp, d, l, A, P = mesh.geom(x, gamma)
    unit = d / np.maximum(l[:, None], 1e-12)
    gc = (A - 1)[mesh.ci, None] * 0.5 * np.stack(
        [pn[:, 1] - pp[:, 1], pp[:, 0] - pn[:, 0]], axis=1
    ) + (P - p0)[mesh.ci, None] * (unit - unit[mesh.prev])
    if noise is not None:
        lam = noise[mesh.eids]
        gc += 0.5 * (lam[:, None] * unit - lam[mesh.prev, None] * unit[mesh.prev])
    return float(np.sum(gc[:, 0] * p[:, 1]) / 16)


def self_intersections(mesh, x, gamma):
    p = mesh.geom(x, gamma)[1]
    for cell in range(mesh.N):
        r = p[mesh.ci == cell]
        n = len(r)

        def cross(a, b):
            return a[0] * b[1] - a[1] * b[0]

        for i in range(n):
            for j in range(i + 2, n):
                if i == 0 and j == n - 1:
                    continue
                a, b = r[i], r[(i + 1) % n]
                c, d = r[j], r[(j + 1) % n]
                if (
                    cross(b - a, c - a) * cross(b - a, d - a) < -1e-14
                    and cross(d - c, a - c) * cross(d - c, b - c) < -1e-14
                ):
                    return True
    return False


def dynamic_pair(args):
    cells, xstart, seed, amplitude, T, dt, threshold, allow_t1, *rest = args
    p0 = 3.75
    gamma_abs = rest[0] if rest else 0.15
    tau = 1.0
    rng = np.random.default_rng(seed)
    xstart = np.array(xstart)
    meshes = [Mesh(cells, xstart), Mesh(cells, xstart)]
    xs = [xstart.copy(), xstart.copy()]
    gammas = [gamma_abs, -gamma_abs]
    metrics = [
        np.linalg.inv(np.array([[1.0, g], [0.0, 1.0]]))
        @ np.linalg.inv(np.array([[1.0, 0.0], [g, 1.0]]))
        for g in gammas
    ]
    nsteps = round(T / dt)
    stride = 20
    noise = np.zeros((16, 16))
    a = np.exp(-dt / tau)
    b = amplitude * np.sqrt(1 - a * a)
    times = []
    responses = []
    events = []
    counts = [0, 0]
    maxjump = 0.0
    minarea = 1.0
    minedge = 1.0
    crossings = 0
    for step in range(nsteps + 1):
        if step % stride == 0:
            sig = []
            for j in range(2):
                m = meshes[j]
                lam = noise[m.pairs[:, 0], m.pairs[:, 1]]
                sig.append(stress(m, xs[j], p0, gammas[j], lam))
                geo = m.geom(xs[j], gammas[j])
                minarea = min(minarea, float(geo[-2].min()))
                minedge = min(minedge, float(geo[5].min()))
                if step % 200 == 0:
                    crossings += int(self_intersections(m, xs[j], gammas[j]))
            times.append(step * dt)
            responses.append((sig[0] - sig[1]) / (2 * gamma_abs))
            events.append(sum(counts) / 2)
        if step == nsteps:
            break
        z = np.triu(rng.normal(size=(16, 16)), 1)
        z += z.T
        noise = a * noise + b * z
        for j in range(2):
            m = meshes[j]
            x = xs[j]
            gamma = gammas[j]
            lam = noise[m.pairs[:, 0], m.pairs[:, 1]]
            e, g = m.energy_grad(x, p0, gamma, lam)
            x = x - dt * (g @ metrics[j])
            if allow_t1:
                ls = m.geom(x, gamma)[5]
                short = {
                    tuple(sorted((m.ids[k], m.ids[n])))
                    for k, n in enumerate(m.nxt)
                    if ls[k] < threshold
                }
                for u, v in sorted(short):
                    try:
                        mn, xn, adj = flip(m, x, (u, v), 1.5 * threshold, gamma)
                    except (ValueError, StopIteration, AssertionError):
                        continue
                    if mn.geom(xn, gamma)[-2].min() < 0.2 or self_intersections(
                        mn, xn, gamma
                    ):
                        continue
                    before = m.energy_grad(
                        x, p0, gamma, noise[m.pairs[:, 0], m.pairs[:, 1]]
                    )[0]
                    after = mn.energy_grad(
                        xn, p0, gamma, noise[mn.pairs[:, 0], mn.pairs[:, 1]]
                    )[0]
                    maxjump = max(maxjump, abs(after - before))
                    m, x = mn, xn
                    counts[j] += 1
            meshes[j], xs[j] = m, x
    return dict(
        seed=seed,
        time=times,
        response=responses,
        events=events,
        min_area=minarea,
        min_edge=minedge,
        intersection_samples=crossings,
        max_event_energy_jump=maxjump,
    )


def run_dynamics(mesh, x, repeats=24):
    from concurrent.futures import ProcessPoolExecutor

    common = (mesh.cells, x.tolist())
    passive = dynamic_pair((*common, 0, 0.0, 40.0, 0.005, 0.025, False))
    passive_t1 = dynamic_pair((*common, 0, 0.0, 40.0, 0.005, 0.025, True))
    args = [(*common, 1000 + i, 0.1, 40.0, 0.005, 0.025, True) for i in range(repeats)]
    with ProcessPoolExecutor(max_workers=4) as ex:
        active = list(ex.map(dynamic_pair, args))
    print(
        "dynamics",
        repeats,
        "active mean T1",
        np.mean([r["events"][-1] for r in active]),
        flush=True,
    )
    return dict(
        amplitude_checks=[
            dict(
                gamma=g,
                **dynamic_pair((*common, 0, 0.0, 200.0, 0.005, 0.025, False, g))
            )
            for g in [0.075, 0.0375]
        ],
        protocol=dict(
            p0=3.75,
            gamma=0.15,
            dt=0.005,
            T=40.0,
            tau=1.0,
            noise_sd=0.1,
            T1_threshold=0.025,
            new_edge_length=0.0375,
            repeats=repeats,
            noise_initialization="zero; driving and step start together",
            noise_labels="independent OU fields on unordered cell pairs, evolving also while non-contacting",
            T1_rule="geometric threshold, perpendicular new edge; reject invalid polygons; do not impose monotone total energy",
        ),
        passive_fixed=passive,
        passive_t1=passive_t1,
        active=active,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=24)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    mesh, x0, centers = prepare()
    stat, states = static_case(mesh, x0)
    paths = [relaxed_t1_path(mesh, x0, p) for p in [3.65, 3.75, 3.79, 3.84, 3.90]]
    x = refine(mesh, x0, 3.75)
    dynamic = run_dynamics(mesh, x, args.repeats)
    data = dict(
        preparation=dict(
            seed=3,
            cell_number=16,
            vertex_number=mesh.V,
            L=4.0,
            initial_centers=centers.tolist(),
            initial_vertices=x0.tolist(),
            cells=mesh.cells,
            initial_method="periodic Voronoi of a 4x4 square grid with uniform [-.25,.25] jitter; identical initial coordinates and topology for each p0; fixed-topology minimization",
        ),
        static=stat,
        states=states,
        paths=paths,
        dynamics=dynamic,
    )
    f, J, H, S, HH = mesh.derivatives(x, 3.75)
    np.savez_compressed(
        root / "unified_vertex_matrices.npz",
        p0=3.75,
        vertices=x,
        J=J,
        H=H,
        S=S,
        affine_extended_hessian=HH,
        translation_projector=projection(mesh.V),
    )
    (root / "unified_vertex_data.json").write_text(json.dumps(data, indent=2) + "\n")
    verification = dict(
        max_hessian_decomposition_error=max(
            r["hessian_identity_max_error"] for r in stat
        ),
        max_force_norm=max(r["force_norm"] for r in stat),
        max_shear_modulus_difference=max(
            abs(r["G_relaxed"] - r["G_finite_difference"]) for r in stat
        ),
        max_path_gradient=max(r["constrained_gradient_max"] for r in paths),
        max_collapsed_branch_mismatch=max(
            r["collapsed_energy_mismatch"] for r in paths
        ),
        min_path_nonselected_edge=min(r["min_other_edge"] for r in paths),
        min_path_area=min(r["min_area"] for r in paths),
        dynamic_intersection_samples=sum(
            r["intersection_samples"] for r in dynamic["active"]
        ),
        min_dynamic_area=min(r["min_area"] for r in dynamic["active"]),
        max_dynamic_event_energy_jump=max(
            r["max_event_energy_jump"] for r in dynamic["active"]
        ),
        passive_amplitude_final_values={
            str(r["gamma"]): r["response"][-1] for r in dynamic["amplitude_checks"]
        },
    )
    assert verification["max_force_norm"] < 1e-8
    assert verification["max_hessian_decomposition_error"] < 1e-9
    assert verification["max_shear_modulus_difference"] < 2e-6
    assert verification["max_path_gradient"] < 1e-5
    assert verification["max_collapsed_branch_mismatch"] < 1e-10
    assert verification["min_path_nonselected_edge"] > 0.01
    assert verification["min_dynamic_area"] > 0.2
    assert verification["dynamic_intersection_samples"] == 0
    (root / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    print(json.dumps(verification, indent=2), flush=True)


if __name__ == "__main__":
    main()
