# Tetra-Code
Tetra code is a novel notation system for scalp positions based on an existing widely-used geonavigation notation system
## Useful tools
`tet_code2MNI_lookup_extended.xlsx` is a lookup table that allows conversion between any given point in the 3D XYZ coordinates of the MNI152 atlas and the Euclidean nearest location on the scalp, expressed in Tetra code notation. It also contains an interactive search tool to convert between MNI coordinates on the scalp and gray matter surface, Tetra code notation, CPC coordinates, Beam-style X/Y percentages, and 10–10 EEG electrode location

`tetraConverter.m` and `tetraConverter.py` allows the above conversion in function form

The `fit_tetra.m` script allows Tetra code mapping on any SimNIBS-generated individual head mesh and stores the relevant information in a `Tetra_code` folder under the same directory as the head mesh.
