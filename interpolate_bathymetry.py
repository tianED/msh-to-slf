"""Interpolate x,y,depth XYZ data onto an existing 2-D SELAFIN mesh.

Requires numpy and scipy. Coordinates must use the same coordinate system
and units. Depth signs are preserved (use --scale -1 to reverse them).
Updates BOTTOM at every existing timestep; other records are preserved.
MATLAB telstepr exposes this variable as a column of RESULT.

Example:
    python interpolate_bathymetry.py bath_points_Fine.xyz Convert_Fine.slf Fine_bathy.slf
"""

import argparse
from pathlib import Path
import struct

import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.spatial import QhullError


def load_xyz(filename):
    """Read comma-separated or whitespace-separated data with named columns."""
    with open(filename, encoding="utf-8-sig") as stream:
        header = stream.readline().strip()
    delimiter = "," if "," in header else None
    names = [name.strip().lower() for name in
             (header.split(",") if delimiter else header.split())]
    if any(names.count(name) != 1 for name in ("x", "y", "depth")):
        raise ValueError("XYZ header must contain unique x, y, depth columns.")
    data = np.loadtxt(filename, delimiter=delimiter, skiprows=1, ndmin=2,
                      encoding="utf-8-sig")
    data = data[:, [names.index(name) for name in ("x", "y", "depth")]]
    if not len(data) or not np.isfinite(data).all():
        raise ValueError("XYZ data must be nonempty and contain only finite values.")
    # Average repeated measurements at the same XY location.
    points, groups = np.unique(data[:, :2], axis=0, return_inverse=True)
    depths = np.bincount(groups, weights=data[:, 2]) / np.bincount(groups)
    return points, depths


def interpolate(points, depths, nodes, outside="nearest"):
    if len(points) < 3:
        raise ValueError("Linear interpolation needs at least 3 distinct XY points.")
    try:
        values = np.asarray(LinearNDInterpolator(points, depths)(nodes))
    except QhullError as exc:
        raise ValueError("XYZ points must span a two-dimensional area.") from exc
    missing = ~np.isfinite(values)
    count = int(missing.sum())
    if count and outside == "error":
        raise ValueError(f"{count} mesh nodes lie outside the XYZ convex hull.")
    if count:
        values[missing] = NearestNDInterpolator(points, depths)(nodes[missing])
    return values, count


def read_record(stream, endian):
    tag = stream.read(4)
    if not tag:
        return None
    if len(tag) != 4:
        raise ValueError("Truncated SELAFIN record marker.")
    size = struct.unpack(endian + "i", tag)[0]
    if size < 0:
        raise ValueError("Invalid SELAFIN record length.")
    payload = stream.read(size)
    if len(payload) != size or stream.read(4) != tag:
        raise ValueError("Invalid or truncated SELAFIN record.")
    return payload


def write_record(stream, payload, endian):
    tag = struct.pack(endian + "i", len(payload))
    stream.write(tag + payload + tag)


def convert(xyz_path, source, target, variable="BOTTOM", outside="nearest",
            scale=1.0, overwrite=False):
    if Path(target).resolve() in (Path(source).resolve(), Path(xyz_path).resolve()):
        raise ValueError("Output must differ from both input files.")
    points, depths = load_xyz(xyz_path)
    if not np.isfinite(scale):
        raise ValueError("Scale must be finite.")
    with open(source, "rb") as src:
        tag = src.read(4)
        endian = ">" if tag == struct.pack(">i", 80) else "<"
        if tag != struct.pack(endian + "i", 80):
            raise ValueError("Expected a SELAFIN file with an 80-byte title.")
        src.seek(0)
        header = [read_record(src, endian), read_record(src, endian)]
        nbv1, nbv2 = struct.unpack(endian + "2i", header[1])
        nvar = nbv1 + nbv2
        names = []
        for _ in range(nvar):
            payload = read_record(src, endian)
            header.append(payload)
            names.append(payload[:16].decode("ascii").strip())
        matches = [i for i, name in enumerate(names) if name.upper() == variable.upper()]
        if len(matches) != 1:
            raise ValueError(f"Select one existing variable with --variable; found {names}.")
        column = matches[0]
        header.append(read_record(src, endian))
        iparam = struct.unpack(endian + "10i", header[-1])
        if iparam[9] == 1:
            header.append(read_record(src, endian))
        header.append(read_record(src, endian))
        nelem, npoin, ndp, _ = struct.unpack(endian + "4i", header[-1])
        if ndp != 3 or iparam[6] > 1 or npoin <= 0 or nelem <= 0:
            raise ValueError("Only 2-D triangular SELAFIN meshes are supported.")
        # Preserve IKLE and IPOBO byte-for-byte.
        header.extend([read_record(src, endian), read_record(src, endian)])
        xrecord, yrecord = read_record(src, endian), read_record(src, endian)
        header.extend([xrecord, yrecord])
        precision = len(xrecord) // npoin
        if precision not in (4, 8) or len(xrecord) != npoin * precision or len(yrecord) != len(xrecord):
            raise ValueError("Unsupported coordinate precision or record size.")
        dtype = np.dtype(endian + f"f{precision}")
        nodes = np.column_stack([np.frombuffer(xrecord, dtype), np.frombuffer(yrecord, dtype)])
        if not np.isfinite(nodes).all():
            raise ValueError("Mesh coordinates must be finite.")
        values, fallback = interpolate(points, depths, nodes, outside)
        values = (values * scale).astype(dtype)
        if not np.isfinite(values).all():
            raise ValueError("Interpolated depths exceed the output precision.")
        # Read the first time record before creating an output file.
        time_record = read_record(src, endian)
        if time_record is None:
            raise ValueError("Input has no timestep data to update.")
        steps = 0
        with open(target, "wb" if overwrite else "xb") as dst:
            for payload in header:
                write_record(dst, payload, endian)
            while time_record is not None:
                if len(time_record) != precision:
                    raise ValueError("Invalid timestep record size.")
                write_record(dst, time_record, endian)
                for index in range(nvar):
                    payload = read_record(src, endian)
                    if payload is None or len(payload) != npoin * precision:
                        raise ValueError("Invalid variable record size.")
                    write_record(dst, values.tobytes() if index == column else payload, endian)
                steps += 1
                time_record = read_record(src, endian)
    return npoin, fallback, steps, float(values.min()), float(values.max())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xyz", type=Path)
    parser.add_argument("input", type=Path, help="Existing SELAFIN file")
    parser.add_argument("output", type=Path, help="New SELAFIN file")
    parser.add_argument("--variable", default="BOTTOM", help="Variable name to update")
    parser.add_argument("--outside", choices=("nearest", "error"), default="nearest")
    parser.add_argument("--scale", type=float, default=1.0, help="Depth multiplier; -1 reverses sign")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    try:
        nodes, fallback, steps, low, high = convert(
            args.xyz, args.input, args.output, args.variable, args.outside,
            args.scale, args.overwrite)
    except (OSError, ValueError, struct.error) as exc:
        parser.exit(1, f"Error: {exc}\n")
    print(f"Wrote {args.output}: {nodes} nodes, {steps} timestep(s), "
          f"{fallback} nearest-point fallbacks, depth range [{low:g}, {high:g}].")


if __name__ == "__main__":
    main()
