# geogrid_to_tif

This is a tool for converting WRF/WPS geogrid binary format to GeoTIFF.

## Example usage

```
python geogrid_to_tif.py WPS_GEOG/modis_20class_30s_with_lakes/ modis_20class_30s_with_lakes.tif
```

Currently only categorical, 1 byte, two-dimensional dataset is supported.
Support for other format might be done in the future.
