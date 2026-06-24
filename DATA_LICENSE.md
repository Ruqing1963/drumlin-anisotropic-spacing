# Data provenance and licence

This repository redistributes only **derived** data: drumlin centroids with a
flow-parallel long-axis azimuth, width and length, computed from the original
landform outlines. The raw shapefiles are **not** redistributed here.

## Source
- **Dataset:** BRITICE Glacial Map, version 2 (filtered GIS database),
  *Subglacial lineations (polygons)* layer.
- **Reference:** Clark, C.D., Ely, J.C., Greenwood, S.L., et al. (2018).
  BRITICE Glacial Map, version 2: a map and GIS database of glacial landforms
  of the last British-Irish Ice Sheet. *Boreas* 47, 11-e8.
  https://doi.org/10.1111/bor.12273
- **Download (raw GIS, 230 MB):** University of Sheffield BRITICE project page
  (https://www.sheffield.ac.uk/geography-planning/research/geography/projects/britice).
- **Licence:** Creative Commons Attribution 4.0 International (CC-BY-4.0).
  **Required attribution:** "Drumlin outlines from BRITICE Glacial Map v2
  (Clark et al., 2018), CC-BY-4.0."

## Derived fields (`data/drumlin_centroids.csv`)
`easting_km, northing_km` (British National Grid), `flow_azimuth_deg`
(long-axis), `width_m`, `length_m`, `field_label` (coherent-field id from
`code/segment_fields.py`). 41,702 drumlins.

Please retain the attribution above when using the derived data.
