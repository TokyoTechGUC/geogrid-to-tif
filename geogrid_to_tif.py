import argparse
import glob
import os

import numpy as np
from affine import Affine
from numba import njit
from osgeo import gdal
from osgeo import osr
from tqdm import tqdm


def read_index(index_path):
    int_fields = {
        "category_min",
        "category_max",
        "tile_x",
        "tile_y",
        "tile_z",
        "tile_z_start",
        "tile_z_end",
        "wordsize",
        "filename_digits",
        "tile_bdr",
    }
    float_fields = {
        "dx",
        "dy",
        "known_x",
        "known_y",
        "known_lat",
        "known_lon",
        "scale_factor",
        "missing_value",
    }
    boolean_fields = {"signed"}
    index = {
        "filename_digits": 5,
        "signed": False,
        "known_x": 1,
        "known_y": 1,
        "scale_factor": 1,
        "endian": "big",
        "tile_bdr": 0,
    }

    with open(index_path) as fp:
        for line in fp:
            line = line.strip()
            if not line or line[0] == "#":
                continue
            k, v = line.split("=", maxsplit=1)
            k = k.lower().strip()
            v = v.strip()

            if k in int_fields:
                v = int(v)
            elif k in float_fields:
                v = float(v)
            elif k in boolean_fields:
                v = v == "yes"
            elif v[0] == '"' and v[-1] == '"':
                v = v[1:-1]
            index[k] = v

    if "tile_z_start" in index or "tile_z_end" in index:
        assert (
            "tile_z_start" in index and "tile_z_end" in index
        ), "Both tile_z_start and tile_z_end are required."
        index["tile_z"] = index["tile_z_end"] - index["tile_z_start"] + 1

    return index


def get_gdal_type(wordsize, signed):
    if signed:
        return [gdal.GDT_Int16, gdal.GDT_Int16, gdal.GDT_Int32, gdal.GDT_Int32][
            wordsize - 1
        ]
    return [gdal.GDT_Byte, gdal.GDT_UInt16, gdal.GDT_UInt32, gdal.GDT_UInt32][
        wordsize - 1
    ]


@njit
def decode_geogrid_binary(binary_str, wordsize, signed):
    arr = []
    for start in range(0, len(binary_str), wordsize):
        num = 0
        for i in range(0, wordsize):
            num = (num << 8) | binary_str[start + i]
        if signed and num >= (1 << (8 * wordsize - 1)):
            num -= 1 << (8 * wordsize)
        arr.append(num)
    return np.array(arr)


def main(geog_dir, output_path):
    index = read_index(os.path.join(geog_dir, "index"))

    assert index["projection"] == "regular_ll"
    assert index["dx"] > 0

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
        index["tile_z"],
        get_gdal_type(index["wordsize"], index["signed"]),
        options=["COMPRESS=LZW", "TILED=YES"],
    )

    sr = osr.SpatialReference()
    sr.ImportFromEPSG(4326)
    ds.SetSpatialRef(sr)

    affine = (
        Affine.translation(index["known_lon"], index["known_lat"])
        * Affine.scale(index["dx"], index["dy"])
        * Affine.translation(-index["known_x"], -index["known_y"])
    )

    is_top_to_bottom = index["dy"] < 0

    if is_top_to_bottom:
        ulx, uly = affine * (1, 1)
        ds.SetGeoTransform([ulx, index["dx"], 0, uly, 0, index["dy"]])
    else:
        ulx, uly = affine * (1, y_size)
        ds.SetGeoTransform([ulx, index["dx"], 0, uly, 0, -index["dy"]])

    bands = [ds.GetRasterBand(i + 1) for i in range(index["tile_z"])]
    for band in bands:
        band.SetScale(index["scale_factor"])
        if "missing_value" in index:
            band.SetNoDataValue(index["missing_value"])

    for tile_path, (x_start, _, y_start, y_end) in tqdm(
        list(zip(tile_paths, tile_bounds))
    ):
        with open(tile_path, "rb") as fp:
            tile = decode_geogrid_binary(fp.read(), index["wordsize"], index["signed"])
        tile = tile.reshape(
            (
                index["tile_z"],
                index["tile_y"] + 2 * index["tile_bdr"],
                index["tile_x"] + 2 * index["tile_bdr"],
            )
        )
        if index["tile_bdr"] > 0:
            tile = tile[
                :,
                index["tile_bdr"] : -index["tile_bdr"],
                index["tile_bdr"] : -index["tile_bdr"],
            ]
        if is_top_to_bottom:
            x_off, y_off = x_start - 1, y_start - 1
        else:
            x_off, y_off = x_start - 1, y_size - y_end
            tile = tile[:, ::-1]
        for band, z_slice in zip(bands, tile):
            band.WriteRaster(
                int(x_off), int(y_off), index["tile_x"], index["tile_y"], z_slice
            )

    bands = None
    ds = None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("geog_dir", help="Path to geogrid directory")
    parser.add_argument("tif", help="Path to save converted tif")
    args = parser.parse_args()
    main(args.geog_dir, args.tif)
