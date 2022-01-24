import glob
import os
import sys

import numpy as np
from affine import Affine
from osgeo import gdal
from osgeo import osr


def read_index(index_path):
    int_fields = [
        "category_min",
        "category_max",
        "tile_x",
        "tile_y",
        "tile_z",
        "wordsize",
        "filename_digits",
    ]
    float_fields = ["dx", "dy", "known_x", "known_y", "known_lat", "known_lon"]
    index = {"filename_digits": 5, "tile_z": 1}

    with open(index_path) as fp:
        for line in fp:
            line = line.strip()
            if not line or line[0] == "#":
                continue
            k, v = line.split("=", maxsplit=1)
            k = k.lower()

            if k in int_fields:
                v = int(v)
            elif k in float_fields:
                v = float(v)
            elif v[0] == '"' and v[-1] == '"':
                v = v[1:-1]
            index[k] = v

    return index


def main(geog_dir, output_path):
    index = read_index(os.path.join(geog_dir, "index"))

    assert index["type"] == "categorical"
    assert index["wordsize"] == 1
    assert index["tile_z"] == 1

    if index["filename_digits"] == 5:
        tile_paths = glob.glob(os.path.join(geog_dir, "?????-?????.?????-?????"))
    else:
        tile_paths = glob.glob(os.path.join(geog_dir, "??????-??????.??????-??????"))
    tile_bounds = []
    for tile_path in tile_paths:
        basename = os.path.basename(tile_path)
        tile_bounds.append(
            [
                int(
                    basename[
                        x
                        * (index["filename_digits"] + 1) : (x + 1)
                        * (index["filename_digits"] + 1)
                        - 1
                    ]
                )
                for x in range(4)
            ]
        )
    tile_bounds = np.array(tile_bounds)

    x_size = int(np.max(tile_bounds[:, 1]))
    y_size = int(np.max(tile_bounds[:, 3]))

    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(
        output_path,
        x_size,
        y_size,
        1,
        gdal.GetDataTypeByName("Byte"),
        options=["COMPRESS=LZW"],
    )

    sr = osr.SpatialReference()
    sr.ImportFromEPSG(4326)
    ds.SetSpatialRef(sr)

    affine = (
        Affine.translation(index["known_lon"], index["known_lat"])
        * Affine.scale(index["dx"], index["dy"])
        * Affine.translation(-index["known_x"], -index["known_y"])
    )
    ulx, uly = affine * (1, y_size)
    ds.SetGeoTransform([ulx, index["dx"], 0, uly, 0, -index["dy"]])

    band = ds.GetRasterBand(1)

    for tile_path, (x_start, _, _, y_end) in zip(tile_paths, tile_bounds):
        with open(tile_path, "rb") as fp:
            tile = np.frombuffer(fp.read(), dtype="u1")
        tile = tile.reshape((index["tile_y"], index["tile_x"]))[::-1]
        band.WriteRaster(
            int(x_start - 1),
            int(y_size - y_end),
            index["tile_x"],
            index["tile_y"],
            tile,
        )


if __name__ == "__main__":
    geog_dir = sys.argv[1]
    output_path = sys.argv[2]
    main(geog_dir, output_path)
