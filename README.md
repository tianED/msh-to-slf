# Usage

1. Convert a Gmsh mesh `mesh.msh` to a SELAFIN mesh `mesh_no_bathy.slf` and a boundary condition file `boundary.cli`.

   By default, Gmsh Physical ID `1` is treated as the liquid/open boundary with prescribed water level `5 4 4`.

```bash
python gmsh_to_slf_cli.py mesh.msh mesh_no_bathy.slf boundary.cli
```

2. Interpolate `Bathymetry.xyz` onto the input SELAFIN mesh `mesh_no_bathy.slf` and generate the final mesh `mesh.slf` with bathymetry information.

```bash
python interpolate_bathymetry.py Bathymetry.xyz mesh_no_bathy.slf mesh.slf
```
