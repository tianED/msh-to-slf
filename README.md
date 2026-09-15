# Usage

1. Convert a Gmsh mesh `mesh.msh` to a SELAFIN mesh `mesh_no_bathy.slf` and a boundary condition file `boundary.cli`.

   By default, Gmsh Physical ID `1` is treated as the liquid/open boundary with prescribed water level.

```bash
python gmsh_to_slf_cli.py mesh.msh mesh_no_bathy.slf boundary.cli
```

2. Interpolate `Bathymetry.xyz` onto the input SELAFIN mesh `mesh_no_bathy.slf` and generate the final mesh `mesh.slf` with bathymetry information.

   The bathymetry file should contain `x`, `y`, and `depth` columns, for example:

```text
x y depth
345000.0 1020000.0 -35.2
345050.0 1020000.0 -34.8
345100.0 1020000.0 -34.5
```

The bathymetry and mesh coordinates should use the same coordinate system and units.

```bash
python interpolate_bathymetry.py Bathymetry.xyz mesh_no_bathy.slf mesh.slf
```
