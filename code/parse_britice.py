"""
parse_britice.py  --  regenerate drumlin centroids from the raw BRITICE shapefile.

Pure-Python (no GIS libraries). Reads the BRITICE v2 filtered
`Subglacial_lineations_polygons.shp` (drumlin outlines, British National Grid),
reduces each polygon to an area-weighted centroid with a long-axis (flow)
azimuth, width and length from a principal-axis decomposition, and writes
`../data/drumlin_centroids_raw.csv` (without field labels).

Run `segment_fields.py` afterwards to add the coherent-field labels and produce
`../data/drumlin_centroids.csv`.

Usage:
  python parse_britice.py /path/to/Subglacial_lineations_polygons.shp
"""
import sys, struct, csv, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def parse_polys(shp):
    d = open(shp, "rb").read(); pos = 100; out = []
    while pos < len(d):
        _, clen = struct.unpack(">ii", d[pos:pos + 8]); pos += 8
        rec = d[pos:pos + clen * 2]; pos += clen * 2
        if len(rec) < 44 or struct.unpack("<i", rec[0:4])[0] != 5:
            continue
        nP, nPt = struct.unpack("<ii", rec[36:44])
        if nP < 1 or nPt < 3:
            continue
        parts = struct.unpack("<%di" % nP, rec[44:44 + 4 * nP])
        xy = np.frombuffer(rec[44 + 4 * nP:44 + 4 * nP + 16 * nPt], dtype="<f8").reshape(-1, 2)
        end = parts[1] if nP > 1 else nPt
        ring = xy[parts[0]:end]; x, y = ring[:, 0], ring[:, 1]
        A = 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
        if abs(A) > 1:
            cx = np.sum((x[:-1] + x[1:]) * (x[:-1] * y[1:] - x[1:] * y[:-1])) / (6 * A)
            cy = np.sum((y[:-1] + y[1:]) * (x[:-1] * y[1:] - x[1:] * y[:-1])) / (6 * A)
        else:
            cx, cy = x.mean(), y.mean()
        Pc = ring - ring.mean(0); Cmat = Pc.T @ Pc / len(Pc); ev, evec = np.linalg.eigh(Cmat)
        major = evec[:, 1]
        az = (np.degrees(np.arctan2(major[0], major[1])) % 180)
        elong = np.sqrt(ev[1] / ev[0]) if ev[0] > 0 else 1.0
        length = np.sqrt(ev[1]) * 4.0
        width = length / elong if elong > 0 else length
        out.append((cx / 1000.0, cy / 1000.0, az, width, length))  # E_km, N_km, az, width_m, length_m
    return out


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    rows = parse_polys(sys.argv[1])
    out = os.path.join(HERE, "..", "data", "drumlin_centroids_raw.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["easting_km", "northing_km", "flow_azimuth_deg", "width_m", "length_m"])
        for cx, cy, az, wd, ln in rows:
            w.writerow([f"{cx:.4f}", f"{cy:.4f}", f"{az:.1f}", f"{wd:.0f}", f"{ln:.0f}"])
    print(f"wrote {out} ({len(rows)} drumlins)")
