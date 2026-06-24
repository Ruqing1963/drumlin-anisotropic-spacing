"""
stats_common.py
Shared functions for the direction-resolved drumlin spacing analysis.
Pure NumPy/SciPy; no GIS dependencies.
"""
import numpy as np
from scipy.spatial import cKDTree, ConvexHull
from matplotlib.path import Path as MplPath


def circ_axis(az_deg):
    """Axial (period-180) circular mean and mean resultant length of a set of
    long-axis azimuths. Returns (mean_axis_deg, concentration R in [0,1])."""
    a = np.radians(np.asarray(az_deg) * 2.0)
    C, S = np.cos(a).mean(), np.sin(a).mean()
    return (np.degrees(np.arctan2(S, C)) / 2.0) % 180.0, np.hypot(C, S)


def rotate_to_flow(P, axis_deg):
    """Rotate points so the flow azimuth lies along +y: x=transverse, y=along-flow."""
    th = np.radians(axis_deg)
    R = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])
    return (P - P.mean(0)) @ R.T


def clark_evans(P, rng, S=99):
    """Edge-corrected Clark-Evans R = mean NND / MC-CSR mean NND. <1 clustered."""
    if len(P) < 50:
        return np.nan
    nn = cKDTree(P).query(P, k=2)[0][:, 1]
    obs = nn.mean()
    V = P[ConvexHull(P).vertices]; pa = MplPath(V); lo, hi = V.min(0), V.max(0)
    nulls = np.empty(S)
    for i in range(S):
        o = []
        while len(o) < len(P):
            q = rng.uniform(lo, hi, (len(P) * 3, 2)); q = q[pa.contains_points(q)]
            o.extend(q.tolist())
        Q = np.array(o[:len(P)]); nulls[i] = cKDTree(Q).query(Q, k=2)[0][:, 1].mean()
    return obs / nulls.mean()


def isotropic_g(P, redges):
    """Isotropic pair correlation g(r) (no edge correction; for field interiors)."""
    D = np.sqrt(((P[:, None] - P[None]) ** 2).sum(2)); np.fill_diagonal(D, np.inf)
    rc = 0.5 * (redges[1:] + redges[:-1])
    area = (P[:, 0].max() - P[:, 0].min()) * (P[:, 1].max() - P[:, 1].min())
    rho = len(P) / area
    g = np.array([((D > redges[i]) & (D <= redges[i + 1])).sum() /
                  (len(P) * rho * np.pi * (redges[i + 1] ** 2 - redges[i] ** 2))
                  for i in range(len(rc))])
    return rc, g


def slab_g(Q, sep_axis, slab_axis, slab=0.5, dx=0.1, xmax=4.0):
    """Directional pair correlation along `sep_axis`, using pairs within a thin
    slab (height `slab`) in `slab_axis`. Isolates 1-D spacing in one direction."""
    s = Q[:, slab_axis]; x = Q[:, sep_axis]
    order = np.argsort(s); s = s[order]; x = x[order]; n = len(Q)
    dxs = []
    for i in range(n):
        j = i + 1
        while j < n and s[j] - s[i] < slab:
            dxs.append(abs(x[j] - x[i])); j += 1
    dxs = np.array(dxs); edges = np.arange(0, xmax, dx); xc = 0.5 * (edges[1:] + edges[:-1])
    cnt, _ = np.histogram(dxs, bins=edges); Lx = x.max() - x.min()
    tri = np.clip(1 - xc / Lx, 0, None)
    g = cnt / (cnt.sum() * (dx / Lx) * tri + 1e-9)
    norm = np.nanmedian(g[(xc > 2.0) & (xc < 3.5)])
    return xc, g / (norm if norm > 0 else 1)


def tonks_null_g(Q, width, rng, slab=0.5, dx=0.1, xmax=4.0, S=200):
    """Transverse g_perp expected from a 1-D equilibrium hard-rod (Tonks) gas:
    each along-flow slab is repopulated with the same number of rods of diameter
    `width`, placed at random subject to non-overlap. Returns (xc, lo95, hi95, med)."""
    s = Q[:, 1]; x = Q[:, 0]; order = np.argsort(s); s = s[order]; x = x[order]; n = len(Q)
    slabs = []; i = 0
    while i < n:
        j = i
        while j < n and s[j] - s[i] < slab:
            j += 1
        idx = np.arange(i, j)
        if len(idx) >= 2:
            slabs.append((x[idx].min(), x[idx].max(), len(idx)))
        i = j
    edges = np.arange(0, xmax, dx); xc = 0.5 * (edges[1:] + edges[:-1]); acc = np.zeros((S, len(xc)))
    meanL = np.mean([sl[1] - sl[0] for sl in slabs]) if slabs else 1.0
    for r in range(S):
        dd = []
        for (x0, x1, m) in slabs:
            L = x1 - x0; avail = L - m * width
            if avail <= 0:
                pos = np.linspace(0, L, m)
            else:
                u = np.sort(rng.uniform(0, avail, m)); pos = u + width * (0.5 + np.arange(m))
            if m >= 2:
                D = np.abs(pos[:, None] - pos[None]); dd.extend(D[np.triu_indices(m, 1)].tolist())
        cnt, _ = np.histogram(np.array(dd), bins=edges)
        tri = np.clip(1 - xc / meanL, 0, None)
        g = cnt / (cnt.sum() * (dx / meanL) * tri + 1e-9)
        norm = np.nanmedian(g[(xc > 2.0) & (xc < 3.5)]); acc[r] = g / (norm if norm > 0 else 1)
    return xc, np.nanpercentile(acc, 2.5, 0), np.nanpercentile(acc, 97.5, 0), np.nanmedian(acc, 0)


def structure_factor(P, kmax, mmax=80, nbin=30):
    """Radially averaged structure factor S(k) at the bounding-window wavevectors."""
    q = P - P.mean(0); Lx, Ly = (P.max(0) - P.min(0))
    ms = np.arange(-mmax, mmax + 1)
    KX, KY = np.meshgrid(2 * np.pi * ms / Lx, 2 * np.pi * ms / Ly)
    K = np.c_[KX.ravel(), KY.ravel()]; kk = np.hypot(K[:, 0], K[:, 1])
    sel = (kk > 0) & (kk <= kmax); K = K[sel]; kk = kk[sel]
    S = (np.abs(np.exp(-1j * (K @ q.T)).sum(1)) ** 2) / len(q)
    bins = np.linspace(0, kmax, nbin); idx = np.digitize(kk, bins); kc = 0.5 * (bins[1:] + bins[:-1])
    return kc, np.array([S[idx == i].mean() if np.any(idx == i) else np.nan
                         for i in range(1, len(bins))])
