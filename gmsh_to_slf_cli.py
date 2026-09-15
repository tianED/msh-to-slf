"""Generate matching TELEMAC .cli/.slf files from Gmsh 2.2 ASCII triangles.

Usage:
  python gmsh_to_slf_cli.py mesh.msh mesh.slf boundary.cli

Requires NumPy only; all conversion functions are included in this file.
After conversion, use interpolate_bathymetry.py to add bathymetry to the .slf.
The output mesh.slf has IPOBO numbered in CLI row order.
Geometry and BOTTOM are written from Gmsh XYZ at time zero.

ID 1 -> 5 4 4; ID 1000 -> 2 2 2. Liquid wins at shared junction nodes.
Outer loops are counterclockwise, islands clockwise. N is the one-based mesh
node index (not necessarily its original Gmsh tag); K is the CLI row number.
Defaults match Fine.cli: LITBOR=2 and all real coefficients zero. --hbor sets
HBOR at liquid nodes only. Actual time-varying forcing is configured separately.
"""

import argparse
from pathlib import Path
import shlex
import struct

import numpy as np

def read_gmsh(filename):
    """Read linear triangles, mapping arbitrary Gmsh node tags to row indices."""
    node_ids, xyz, triangles = [], [], []
    seen = set()
    with open(filename, encoding="utf-8-sig") as stream:
        if stream.readline().strip() != "$MeshFormat":
            raise ValueError("Missing $MeshFormat header.")
        if stream.readline().split() != ["2.2", "0", "8"]:
            raise ValueError("Export the mesh as Gmsh 2.2 ASCII (2.2 0 8).")
        if stream.readline().strip() != "$EndMeshFormat":
            raise ValueError("Missing $EndMeshFormat.")
        for line in stream:
            section = line.strip()
            if not section:
                continue
            if not section.startswith("$") or section.startswith("$End"):
                raise ValueError(f"Unexpected section: {section}")
            end = "$End" + section[1:]
            if section in ("$Nodes", "$Elements"):
                if section in seen:
                    raise ValueError(f"Repeated section: {section}")
                seen.add(section)
                count = int(stream.readline())
                if count <= 0:
                    raise ValueError(f"Empty section: {section}")
                for _ in range(count):
                    fields = stream.readline().split()
                    if section == "$Nodes":
                        if len(fields) != 4:
                            raise ValueError("Invalid node record.")
                        node_ids.append(int(fields[0]))
                        xyz.append([float(value) for value in fields[1:]])
                    else:
                        values = [int(value) for value in fields]
                        if len(values) < 3 or values[2] < 0:
                            raise ValueError("Invalid element record.")
                        element_type = values[1]
                        # Points and lines are unnecessary for edge counting.
                        if element_type in (15, 1):
                            continue
                        if element_type != 2:
                            raise ValueError(
                                f"Unsupported element type {element_type}; "
                                "use a 2-D mesh of linear triangles."
                            )
                        nodes = values[3 + values[2]:]
                        if len(nodes) != 3:
                            raise ValueError("Invalid triangle record.")
                        triangles.append(nodes)
                if stream.readline().strip() != end:
                    raise ValueError(f"Missing {end}.")
            else:
                # Ignore optional metadata, including physical group names.
                for extra in stream:
                    if extra.strip() == end:
                        break
                else:
                    raise ValueError(f"Missing {end}.")
    if not node_ids or not triangles:
        raise ValueError("Mesh must contain nodes and triangles.")
    if min(node_ids) <= 0 or len(set(node_ids)) != len(node_ids):
        raise ValueError("Node tags must be positive and unique.")
    mapping = {tag: index for index, tag in enumerate(node_ids)}
    try:
        ikle = np.array([[mapping[tag] for tag in tri] for tri in triangles],
                        dtype=np.int64)
    except KeyError as exc:
        raise ValueError(f"Triangle references missing node {exc.args[0]}.") from exc
    xyz = np.asarray(xyz, dtype=np.float64)
    if not np.isfinite(xyz).all():
        raise ValueError("Coordinates must be finite.")
    if np.unique(ikle).size != len(xyz):
        raise ValueError("Mesh contains nodes unused by triangles; remove them first.")
    if len(np.unique(np.sort(ikle, axis=1), axis=0)) != len(ikle):
        raise ValueError("Mesh contains duplicate triangles.")
    # Ensure counterclockwise connectivity instead of blindly swapping columns.
    a, b, c = xyz[ikle, :2].transpose(1, 0, 2)
    area2 = ((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
             - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0]))
    if np.any(area2 == 0):
        raise ValueError("Mesh contains zero-area triangles.")
    clockwise = area2 < 0
    ikle[clockwise] = ikle[clockwise][:, [0, 2, 1]]
    return xyz, ikle


def write_selafin(filename, xyz, ikle, ipobo, overwrite=False):
    """Write Fortran records in the same layout as parser_gmsh/Selafin."""
    def record(stream, values, dtype=None):
        payload = values if dtype is None else np.asarray(values, dtype=dtype).tobytes()
        tag = struct.pack(">i", len(payload))
        stream.write(tag)
        stream.write(payload)
        stream.write(tag)

    with open(filename, "wb" if overwrite else "xb") as stream:
        record(stream, b"newSelafin".ljust(80))
        record(stream, [1, 0], ">i4")
        record(stream, b"BOTTOM".ljust(16) + b"M".ljust(16))
        record(stream, [1, 0, 0, 0, 0, 0, 0, 0, 0, 0], ">i4")
        record(stream, [len(ikle), len(xyz), 3, 1], ">i4")
        record(stream, ikle + 1, ">i4")  # SELAFIN uses one-based node indices.
        record(stream, ipobo, ">i4")
        record(stream, xyz[:, 0], ">f4")
        record(stream, xyz[:, 1], ">f4")
        record(stream, [0.0], ">f4")
        record(stream, xyz[:, 2], ">f4")


def read_labels(filename):
    """Return (node mapping, line records, dimension-1 physical names)."""
    nodes, lines, names = {}, [], {}
    seen = set()
    with open(filename, encoding="utf-8-sig") as stream:
        if stream.readline().strip() != "$MeshFormat":
            raise ValueError("Missing $MeshFormat.")
        if stream.readline().split() != ["2.2", "0", "8"]:
            raise ValueError("Use Gmsh 2.2 ASCII format (2.2 0 8).")
        if stream.readline().strip() != "$EndMeshFormat":
            raise ValueError("Missing $EndMeshFormat.")
        for line in stream:
            section = line.strip()
            if not section:
                continue
            if not section.startswith("$") or section.startswith("$End"):
                raise ValueError(f"Unexpected section {section}.")
            end = "$End" + section[1:]
            if section in ("$Nodes", "$Elements", "$PhysicalNames"):
                if section in seen:
                    raise ValueError(f"Repeated section {section}.")
                seen.add(section)
                count = int(stream.readline())
                if count < 0:
                    raise ValueError("Negative section count.")
                for index in range(count):
                    fields = shlex.split(stream.readline())
                    if section == "$Nodes":
                        if len(fields) != 4:
                            raise ValueError("Invalid node record.")
                        tag = int(fields[0])
                        if tag <= 0 or tag in nodes:
                            raise ValueError("Node IDs must be positive and unique.")
                        nodes[tag] = (index + 1, *map(float, fields[1:]))
                    elif section == "$PhysicalNames":
                        if len(fields) != 3:
                            raise ValueError("Invalid physical name record.")
                        if int(fields[0]) == 1:
                            names[int(fields[1])] = fields[2]
                    else:
                        values = list(map(int, fields))
                        if len(values) < 3:
                            raise ValueError("Invalid element record.")
                        eid, kind, ntags = values[:3]
                        if ntags < 0 or len(values) < 3 + ntags:
                            raise ValueError("Invalid element tags.")
                        if kind in (8, 26, 27, 28):
                            raise ValueError("Higher-order lines are unsupported; export linear elements.")
                        if kind != 1:
                            continue
                        tags = values[3:3 + ntags]
                        endpoints = values[3 + ntags:]
                        if len(endpoints) != 2:
                            raise ValueError("Invalid linear element record.")
                        # MSH2: first tag = physical ID, second = geometric ID.
                        lines.append((eid, tags[0] if tags else 0,
                                      tags[1] if len(tags) > 1 else 0, *endpoints))
                if stream.readline().strip() != end:
                    raise ValueError(f"Missing {end}.")
            else:
                for extra in stream:
                    if extra.strip() == end:
                        break
                else:
                    raise ValueError(f"Missing {end}.")
    if not nodes or not lines:
        raise ValueError("No nodes or linear 1-D elements found in the file.")
    for _, _, _, first, second in lines:
        if first not in nodes or second not in nodes:
            raise ValueError("Line references an unknown node ID.")
    if len({row[0] for row in lines}) != len(lines):
        raise ValueError("Duplicate line element IDs.")
    return nodes, lines, names


def ordered_boundary(xyz, ikle):
    """Follow directed single-use edges of counterclockwise triangles."""
    edges = np.concatenate([ikle[:, [0, 1]], ikle[:, [1, 2]], ikle[:, [2, 0]]])
    _, groups, counts = np.unique(np.sort(edges, axis=1), axis=0,
                                  return_inverse=True, return_counts=True)
    if np.any(counts > 2):
        raise ValueError('Nonmanifold mesh edge.')
    boundary = edges[counts[groups] == 1]
    successor = {int(a): int(b) for a, b in boundary}
    if (not len(boundary) or len(successor) != len(boundary)
            or len(set(successor.values())) != len(boundary)
            or set(successor) != set(successor.values())):
        raise ValueError('Boundary must consist of disjoint closed loops.')
    remaining = set(successor)
    loops = []
    while remaining:
        start = min(remaining)
        loop, node = [], start
        while node in remaining:
            remaining.remove(node)
            loop.append(node)
            node = successor[node]
        if node != start:
            raise ValueError('Boundary traversal did not close at its start.')
        xy = xyz[loop, :2] - xyz[loop[0], :2]
        area = .5 * np.sum(xy[:, 0] * np.roll(xy[:, 1], -1)
                            - xy[:, 1] * np.roll(xy[:, 0], -1))
        if area == 0:
            raise ValueError('Boundary loop has zero area.')
        loops.append((area, loop))
    # Directed CCW triangle edges keep water on the left: islands are CW.
    loops.sort(key=lambda item: (item[0] < 0, -abs(item[0]), item[1][0]))
    return loops, boundary


def make_cli(xyz, ikle, nodes, lines, hbor=0.0, liquid_litbor=2):
    loops, edges = ordered_boundary(xyz, ikle)
    labels = {}
    for _, physical, _, first, second in lines:
        edge = tuple(sorted((nodes[first][0] - 1, nodes[second][0] - 1)))
        labels.setdefault(edge, set()).add(physical)
    memberships = {}
    for first, second in edges:
        edge = tuple(sorted((int(first), int(second))))
        tags = labels.get(edge, set())
        if len(tags) != 1 or not tags <= {1, 1000}:
            raise ValueError(f'Boundary edge {edge} needs one physical ID, 1 or 1000; found {tags}.')
        for node in edge:
            memberships.setdefault(node, set()).update(tags)
    ipobo = np.zeros(len(xyz), dtype=np.int32)
    rows, liquid = [], 0
    for _, loop in loops:
        for node in loop:
            is_liquid = 1 in memberships[node]
            liquid += is_liquid
            code = '5 4 4' if is_liquid else '2 2 2'
            tracer = liquid_litbor if is_liquid else 2
            k = len(rows) + 1
            ipobo[node] = k
            rows.append(f'{code}  {hbor if is_liquid else 0.0:.9g} 0.000 0.000 0.000'
                        f'  {tracer}  0.000 0.000 0.000  {node + 1:10d} {k:10d}\n')
    return ''.join(rows), ipobo, loops, liquid


def generate(mesh, slf, cli, hbor=0.0, liquid_litbor=2, overwrite=False):
    cli, slf = Path(cli), Path(slf)
    if cli.suffix.lower() != '.cli' or slf.suffix.lower() != '.slf':
        raise ValueError('Output filenames must end in .slf and .cli respectively.')
    for output in (cli, slf):
        if output.resolve() == Path(mesh).resolve():
            raise ValueError('Output must not replace the input mesh.')
        if output.exists() and not overwrite:
            raise FileExistsError(f'{output} exists; use --overwrite to replace it.')
    if not np.isfinite(hbor):
        raise ValueError('HBOR must be finite.')
    xyz, ikle = read_gmsh(mesh)
    nodes, lines, _ = read_labels(mesh)
    content, ipobo, loops, liquid = make_cli(xyz, ikle, nodes, lines, hbor, liquid_litbor)
    write_selafin(slf, xyz, ikle, ipobo, overwrite)
    with cli.open('w' if overwrite else 'x', encoding='ascii', newline='') as dst:
        dst.write(content)
    print(f'Wrote {cli} and {slf}')
    print(f'{sum(area > 0 for area, _ in loops)} outer loop(s), '
          f'{sum(area < 0 for area, _ in loops)} island loop(s); '
          f'{liquid} liquid nodes, {np.count_nonzero(ipobo) - liquid} wall nodes.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mesh', type=Path)
    parser.add_argument('slf', type=Path, help='Output SELAFIN .slf filename')
    parser.add_argument('cli', type=Path, help='Output boundary .cli filename')
    parser.add_argument('--hbor', type=float, default=0.0)
    parser.add_argument('--liquid-litbor', type=int, choices=(2, 4, 5), default=2)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    try:
        generate(args.mesh, args.slf, args.cli, args.hbor, args.liquid_litbor, args.overwrite)
    except (OSError, ValueError) as exc:
        parser.exit(1, f'Error: {exc}\n')


if __name__ == '__main__':
    main()
