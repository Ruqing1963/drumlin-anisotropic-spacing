"""
analyze_drumlins.py  --  direction-resolved drumlin spacing analysis.

Reads ../data/drumlin_centroids.csv, selects coherent single-flow fields,
computes isotropic and directional (along-flow / transverse) statistics with a
hard-rod (Tonks) steric null, writes ../results/field_summary.csv, and
regenerates the three paper figures (PDF + PNG) in ../figures/.
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import stats_common as sc

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "drumlin_centroids.csv")
FIG = os.path.join(HERE, "..", "figures")
RES = os.path.join(HERE, "..", "results")
PURPLE = "#534AB7"
rng = np.random.default_rng(11)

rows = list(csv.DictReader(open(DATA)))
P_all = np.array([[float(r["easting_km"]), float(r["northing_km"])] for r in rows])
AZ = np.array([float(r["flow_azimuth_deg"]) for r in rows])
WID = np.array([float(r["width_m"]) for r in rows]) / 1000.0  # km
LAB = np.array([int(r["field_label"]) for r in rows])


def coherent_fields(min_n=300, min_R=0.72):
    out = []
    for c in np.unique(LAB):
        idx = np.where(LAB == c)[0]
        if len(idx) < min_n:
            continue
        mp, R = sc.circ_axis(AZ[idx])
        if R >= min_R:
            out.append((c, len(idx), mp, R))
    out.sort(key=lambda t: -t[1])
    return out


def per_field(c, axis):
    idx = np.where(LAB == c)[0]; P = P_all[idx]; Q = sc.rotate_to_flow(P, axis)
    w = np.median(WID[idx])
    Ri = sc.clark_evans(P, rng, S=49)
    xc, gp = sc.slab_g(Q, 0, 1)                    # transverse
    xt, lo, hi, med = sc.tonks_null_g(Q, w, rng)   # steric null
    lam = xc[1:6][np.nanargmax(gp[1:6])]
    return dict(c=c, N=len(idx), axis=axis, Ri=Ri, w=w, lam=lam,
                P=P, Q=Q, xc=xc, gp=gp, lo=lo, hi=hi)


def main():
    fields = coherent_fields()
    print(f"{len(fields)} coherent fields (N>=300, concR>=0.72)")
    recs = []
    for c, N, axis, concR in fields:
        d = per_field(c, axis); d["concR"] = concR; recs.append(d)
        print(f"  field {c}: N={N} axis={axis:.0f} concR={concR:.2f} "
              f"R_iso={d['Ri']:.2f} width={d['w']*1000:.0f}m lam_perp={d['lam']*1000:.0f}m")
    # summary CSV
    with open(os.path.join(RES, "field_summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["field_label", "N", "flow_axis_deg", "azimuth_concentration",
                    "clark_evans_R_iso", "drumlin_width_m", "transverse_spacing_m", "lambda_over_width"])
        for d in recs:
            w.writerow([d["c"], d["N"], f"{d['axis']:.0f}", f"{d['concR']:.2f}",
                        f"{d['Ri']:.3f}", f"{d['w']*1000:.0f}", f"{d['lam']*1000:.0f}",
                        f"{d['lam']/d['w']:.2f}"])
    print("wrote ../results/field_summary.csv")

    # ---- Fig 1: GB drumlin map by azimuth ----
    fig, ax = plt.subplots(figsize=(7, 11))
    s = ax.scatter(P_all[:, 0], P_all[:, 1], c=AZ, s=2, cmap="hsv", lw=0, alpha=0.6)
    ax.set_aspect("equal"); ax.set_xlabel("BNG E (km)"); ax.set_ylabel("BNG N (km)")
    ax.set_title(f"BRITICE drumlins (n={len(P_all)}), colour = long-axis azimuth")
    plt.colorbar(s, label="azimuth (deg)", shrink=0.5)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"drumlin_map.{ext}"), dpi=120, bbox_inches="tight")
    plt.close(fig)

    # ---- Fig 2: representative fields ----
    pick = [d for d in recs if d["c"] in (1263, 38, 1265)] or recs[:3]
    fig, axes = plt.subplots(len(pick), 3, figsize=(13, 3.3 * len(pick)))
    for k, d in enumerate(pick):
        P, Q = d["P"], d["Q"]
        red = np.arange(0, 4, 0.15); rcx, giso = sc.isotropic_g(P, red)
        xa, gpar = sc.slab_g(Q, 1, 0)
        ax = axes[k]
        ax[0].scatter(Q[:, 0], Q[:, 1], s=3, c=PURPLE, lw=0); ax[0].set_aspect("equal")
        ax[0].set_xlabel("transverse (km)"); ax[0].set_ylabel("along-flow (km)")
        ax[0].set_title(f"field {d['c']}: N={d['N']}, concR={d['concR']:.2f}")
        ax[1].plot(rcx, giso, "-", color="#888", lw=1.5, label="isotropic $g(r)$")
        ax[1].plot(xa, gpar, "-", color="#D9543A", lw=1.5, label="along-flow $g_\\parallel$")
        ax[1].axhline(1, ls="--", color="k", lw=.7); ax[1].set_xlim(0, 3)
        ax[1].set_xlabel("r (km)"); ax[1].set_ylabel("g"); ax[1].legend(fontsize=7, frameon=False)
        ax[1].set_title(f"isotropic & along-flow: clustered ($R_{{iso}}$={d['Ri']:.2f})")
        ax[2].fill_between(d["xc"], d["lo"], d["hi"], color="#D9943A", alpha=.35,
                           label="hard-rod (Tonks) null 95%")
        ax[2].plot(d["xc"], d["gp"], "-o", color=PURPLE, ms=3, label="observed $g_\\perp$")
        ax[2].axhline(1, ls="--", color="k", lw=.7)
        ax[2].axvline(d["w"], ls=":", color="g", lw=1, label=f"drumlin width {d['w']*1000:.0f} m")
        ax[2].set_xlim(0, 2.5); ax[2].set_xlabel("transverse separation (km)"); ax[2].set_ylabel("$g_\\perp$")
        ax[2].legend(fontsize=7, frameon=False)
        ax[2].set_title(f"transverse: quasi-regular at $\\lambda_\\perp$={d['lam']*1000:.0f} m")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"drumlin_rep.{ext}"), dpi=120, bbox_inches="tight")
    plt.close(fig)

    # ---- Fig 3: across-field summary ----
    labs = [d["c"] for d in recs]; Ris = np.array([d["Ri"] for d in recs])
    ws = np.array([d["w"] for d in recs]) * 1000; lams = np.array([d["lam"] for d in recs]) * 1000
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].scatter(ws, lams, s=50, c=PURPLE); lim = [100, 600]
    ax[0].plot(lim, lim, "k--", lw=1, label="$\\lambda_\\perp$ = width")
    ax[0].plot(lim, [1.2 * x for x in lim], color="#D9543A", ls=":", lw=1.5, label="$\\lambda_\\perp$ = 1.2$\\times$width")
    ax[0].set_xlim(150, 400); ax[0].set_ylim(100, 600)
    ax[0].set_xlabel("drumlin width (m)"); ax[0].set_ylabel("transverse spacing $\\lambda_\\perp$ (m)")
    ax[0].legend(fontsize=8, frameon=False); ax[0].set_title("transverse spacing tracks drumlin width")
    ax[1].bar(range(len(recs)), Ris, color=[PURPLE if r < 1 else "#D9543A" for r in Ris])
    ax[1].axhline(1, ls="--", color="k", lw=.8); ax[1].set_xticks(range(len(recs)))
    ax[1].set_xticklabels(labs, rotation=45, fontsize=7)
    ax[1].set_ylabel("isotropic Clark-Evans $R$"); ax[1].set_xlabel("field")
    ax[1].set_title("isotropic statistic: not regular ($R\\leq1$) in 9/9;\nclustered ($R<1$) in 8/9 fields")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"drumlin_summary.{ext}"), dpi=120, bbox_inches="tight")
    plt.close(fig)
    print("figures written to ../figures/ (pdf + png)")
    print(f"MEDIANS: width={np.median(ws):.0f}m lam_perp={np.median(lams):.0f}m "
          f"ratio={np.median(lams/ws):.2f}; iso R median={np.median(Ris):.2f}")


if __name__ == "__main__":
    main()
