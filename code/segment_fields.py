"""
segment_fields.py  --  assign coherent-field labels to drumlin centroids.

Links drumlins whose centroids lie within three median nearest-neighbour
distances (density connectivity) and labels the connected components as fields.
Reads ../data/drumlin_centroids_raw.csv (from parse_britice.py) and writes the
labelled ../data/drumlin_centroids.csv used by analyze_drumlins.py.
"""
import os, csv
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    src = os.path.join(HERE, "..", "data", "drumlin_centroids_raw.csv")
    rows = list(csv.DictReader(open(src)))
    P = np.array([[float(r["easting_km"]), float(r["northing_km"])] for r in rows])
    nn = cKDTree(P).query(P, k=2)[0][:, 1]; d0 = 3.0 * np.median(nn)
    pairs = cKDTree(P).query_pairs(d0, output_type="ndarray"); n = len(P)
    M = csr_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, lab = connected_components(M + M.T, directed=False)
    out = os.path.join(HERE, "..", "data", "drumlin_centroids.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["easting_km", "northing_km", "flow_azimuth_deg", "width_m", "length_m", "field_label"])
        for r, l in zip(rows, lab):
            w.writerow([r["easting_km"], r["northing_km"], r["flow_azimuth_deg"],
                        r["width_m"], r["length_m"], int(l)])
    print(f"wrote {out} ({n} drumlins, {lab.max()+1} components)")


if __name__ == "__main__":
    main()
