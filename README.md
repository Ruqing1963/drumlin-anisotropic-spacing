# Anisotropic spacing of drumlins

Direction-resolved point-pattern analysis of drumlin spacing in the **BRITICE v2**
record. Drumlin "spacing regularity", examined as a 2-D point pattern, turns out
to be strongly **anisotropic**: drumlins **cluster along ice flow** (bedforms in
trains) while being **quasi-regular across flow** at a scale set by the bedform
width. The transverse spacing is largely accounted for by steric (non-overlap)
packing, and the pattern is **not hyperuniform**.

This repository contains the code, the derived data and the figures for the paper:

> R. Chen (2026). *Anisotropic spacing of drumlins: transverse quasi-regularity at
> the bedform scale versus along-flow clustering in the BRITICE record.*
> ([`paper/drumlin_paper.pdf`](paper/drumlin_paper.pdf))

## Key results

| Quantity | Value |
|---|---|
| Drumlins analysed (Great Britain) | 41,702 |
| Coherent single-flow fields (N≥300, azimuth concentration ≥0.72) | 9 (N = 504–1547) |
| Isotropic Clark–Evans *R* | ≤1 in 9/9 fields (clustered, *R*<1, in 8/9; median 0.89) |
| Transverse spacing λ⊥ | ~250 m (median), ≈1.2× the drumlin width |
| First-neighbour transverse order vs width-matched hard-rod (Tonks) null | largely steric |
| Hyperuniform? | No (structure factor patchiness-dominated, ≈1 at the spacing scale) |

**Interpretation.** Isotropic spacing regularity is a weak diagnostic of
bed-instability theories for drumlins, because a transverse spacing close to the
bedform width is observationally degenerate with simple steric packing. Transverse,
direction-resolved statistics and size–spacing covariation are more diagnostic.

## Repository layout

```
code/      analysis pipeline (pure NumPy/SciPy/Matplotlib)
  stats_common.py      shared statistics: circular axis, flow rotation,
                       Clark-Evans, isotropic g(r), directional slab g,
                       Tonks hard-rod null, structure factor
  parse_britice.py     raw BRITICE shapefile -> centroids CSV (no GIS deps)
  segment_fields.py    density-connectivity field labelling
  analyze_drumlins.py  field selection, stats, results CSV, the 3 figures
data/
  drumlin_centroids.csv   41,702 derived drumlin centroids (see data/README.md)
figures/   the three paper figures (PDF + PNG)
results/   field_summary.csv (per-field statistics)
paper/     manuscript (PDF + LaTeX source + figures)
```

## Reproduce

```bash
pip install -r requirements.txt

# (A) From the shipped derived data — reproduces all figures and the summary table:
cd code
python analyze_drumlins.py

# (B) From the raw BRITICE shapefile (optional, full pipeline):
#   download the BRITICE v2 filtered GIS (see data/README.md), then
python parse_britice.py /path/to/Subglacial_lineations_polygons.shp
python segment_fields.py
python analyze_drumlins.py
```

Outputs are written to `../figures/` (PDF + PNG) and `../results/field_summary.csv`.

## Data source and licence

Drumlin outlines are from the **BRITICE Glacial Map, version 2** (Clark et al.,
2018, *Boreas*, [doi:10.1111/bor.12273](https://doi.org/10.1111/bor.12273)),
distributed under **CC-BY-4.0**. This repository redistributes only **derived**
centroids; see [`DATA_LICENSE.md`](DATA_LICENSE.md) for the required attribution
and the link to the raw GIS download. Code is released under the MIT License
([`LICENSE`](LICENSE)).

## Companion study

A contrasting geological spacing system — salt structures, where a vertical
buoyancy (Rayleigh–Taylor) instability with no imposed flow direction produces an
isotropic, hard-core, liquid-like regular spacing — is analysed in the companion
repository [salt-structure-spacing](https://github.com/Ruqing1963/salt-structure-spacing)
(Zenodo [10.5281/zenodo.20834957](https://doi.org/10.5281/zenodo.20834957)).
