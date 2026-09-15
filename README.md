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

## Physical IDs and boundary conditions

A Gmsh physical ID is an integer label assigned to a group of 1-D line elements (Physical Curves), such as an open boundary or coastline. It is not a mesh node number or a SELAFIN `IPOBO` number. Nodes inherit the labels of the lines they belong to.

By default, physical ID `1` uses `5 4 4 2`, and all other IDs use `2 2 2 2` (solid wall). The four codes are `LIHBOR LIUBOR LIVBOR LITBOR`: depth, U velocity, V velocity, and tracer boundary conditions.

Override an ID using `--bc ID LIHBOR LIUBOR LIVBOR LITBOR`, repeating the option for additional IDs:

```bash
python gmsh_to_slf_cli.py mesh.msh mesh_no_bathy.slf boundary.cli --bc 1 5 4 4 4 --bc 7 4 6 6 5
```

Here, ID `1` has prescribed depth with free velocity and tracer; ID `7` has free depth with prescribed velocity and tracer. Codes must be integers from `0` to `6`; prescribed values are configured separately. At a node shared by several groups, the last matching `--bc` rule takes precedence over earlier rules and defaults. Without overrides, ID `1` takes precedence over walls.
