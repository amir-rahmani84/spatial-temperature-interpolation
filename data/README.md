# Data acquisition and preparation

Evaluation needs only `data.nc` and `gadm41_IRN_1.json`. Put authorized copies here or set `SPATIAL_DATA_DIR` to their directory. Full IDW maps additionally need the **replacement** `gilan_dem_masked.csv` and `gilan_dem_masked.provenance.json`. `SPATIAL_DEM_FILE` may point to that CSV elsewhere. Never include credentials.

## ERA5

Dataset: [ERA5 monthly averaged data on single levels](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means), DOI [10.24381/cds.f17050d7](https://doi.org/10.24381/cds.f17050d7). The verified original has 120 monthly fields (2015–2024), nine latitudes 38.5–36.5°N and nine longitudes 48.5–50.5°E at 0.25° spacing. `t2m` is in K (parameter 167); surface `z` is in m² s⁻² (parameter 129). Metadata records GRIB-to-NetCDF conversion on 2025-08-05, not an original download receipt.

Follow [official CDS setup](https://cds.climate.copernicus.eu/how-to-api), accept the dataset terms, and configure your account outside this folder. Install the optional acquisition client with `python -m pip install "cdsapi>=0.7.7"`. Exact request implemented in `download_era5.py`:

```python
import cdsapi
cdsapi.Client().retrieve(
    "reanalysis-era5-single-levels-monthly-means",
    {
        "product_type": ["monthly_averaged_reanalysis"],
        "variable": ["2m_temperature", "geopotential"],
        "year": [str(y) for y in range(2015, 2025)],
        "month": [f"{m:02d}" for m in range(1, 13)],
        "time": ["00:00"],
        "area": [38.5, 48.5, 36.5, 50.5],  # north, west, south, east
        "data_format": "netcdf",
        "download_format": "unarchived",
    },
    "data/data.nc",
)
```

The default monthly reanalysis output grid is 0.25°; no custom regridding is requested. This requests monthly means across all hours, not monthly means for a selected hour of day; 00:00 is the product's nominal time selection. `python download_era5.py --request-only` prints the request without authentication. `python download_era5.py` downloads, refusing to overwrite an existing input. The current official collection form was inspected on 2026-11-27 and contains these product, variable and format values; [the collection API](https://cds.climate.copernicus.eu/api/catalogue/v1/collections/reanalysis-era5-single-levels-monthly-means) links to the live form. An authenticated download was not performed. This is a replacement request, not a recovered historical receipt.

To verify a new file despite a changed checksum:

1. Check `valid_time`, `latitude`, `longitude`, variable names/units, 120 unique month starts, 9 × 9 coordinates and no missing values. If the service splits files or changes dimension names, merge/rename explicitly and record that preparation.
2. Sort both old and new datasets by time/latitude/longitude; transpose each variable to the same dimension order; compare coordinate and `t2m`/`z` arrays numerically. Different encoding/compression/history can change the file hash without changing physical values. Scalar `number` and string `expver` differences do not alone establish a physical difference.
3. Prefer exact physical-array equality. If values differ, record maximum absolute differences in K and m² s⁻², compare climatologies and rerun evaluation before accepting the new input. Do not declare equivalence from filenames or a similar map.
4. The Gilan mask must retain 23 locations and 2,760 records. Original numerical input hashes/metadata are in `results/data_audit.json`; verified CSV hashes are in `results/verified_numerical_baseline.json`. Celsius conversion is applied exactly once.

The current official catalogue declares CC-BY-4.0. Preserve attribution, citation, source and modification notice; original acquisition terms/receipt are unavailable. Acknowledge: **Contains modified Copernicus Climate Change Service information (2015–2024). Neither the European Commission nor ECMWF is responsible for use of the information.** Raw NetCDF remains excluded by default.

## GADM

Acquire **GADM 4.1, Iran, level 1** through [GADM](https://gadm.org/download_country.html). Export an authorized local copy to EPSG:4326 GeoJSON, retaining `NAME_1`, as `gadm41_IRN_1.json`. Select the single `NAME_1 == 'Gilan'` feature. Its supplied west/south/east/north bounds are **48.4351, 36.5605, 50.5765, 38.4455**. Check the recorded hash before substituting another boundary revision.

[GADM's license](https://gadm.org/license.html) permits academic/noncommercial use and academic publication of maps; raw redistribution requires permission. Raw GADM stays excluded. The academic derived figures fall within its stated map-publication allowance, not a general permission to redistribute boundaries. The masked elevation CSV remains local by default too.

## Completed SRTM replacement

The six actual rasters were acquired on **27 November 2026** from OpenTopography's documented public SRTM GL1 mirror. This route was used after USGS registration proved unsuitable for the author. No account, API key or personal profile was supplied. The replacement is **NASA SRTM GL1 version 3**, 1-arcsecond posting (approximately 30 m), WGS84 geographic coordinates and EGM96 orthometric elevations in metres.

Actual distribution citation: **NASA Shuttle Radar Topography Mission (SRTM) (2013). Shuttle Radar Topography Mission (SRTM) Global. Distributed by OpenTopography. https://doi.org/10.5069/G9445JDF. Accessed 2026-11-27.** Select GL1, not GL3 or the GL1 ellipsoidal variant. [Official dataset/access page](https://portal.opentopography.org/raster?opentopoID=OTSRTM.082015.4326.1) documents unsigned bulk access; [citation and licensing policy](https://opentopography.org/citations) requires source/OpenTopography acknowledgment. [USGS SRTM 1 Arc-Second Global](https://doi.org/10.5066/F7PR7TFT) supplies the product overview/public-domain context; the files used here are not described as EarthExplorer downloads or claimed byte-identical to another distributor's exports.

All tiles have **3601 × 3601** point postings, EPSG:4326, spacing 1/3600 degree in both axes, NoData -32768 and `AREA_OR_POINT=Point`. Their bounds include the shared outer edge postings. Coverage was derived from the Gilan boundary and masked interpolation grid: 57,536 centres from 48.440–50.575°E and 36.565–38.445°N. Full metadata, actual source URLs, checksums and verification timestamps are in [../results/dem_provenance.json](../results/dem_provenance.json).

| Tile | SHA-256 |
|---|---|
| N36E049 | `5fdea7f218f3bdc05fe42d360d32a6c50b2b68a2c126c99a0e7da24690a70d81` |
| N36E050 | `08582dda55623edfead9b4d9b7daa468c89470969ad3b866ec4cb771303e7fa1` |
| N37E048 | `d735880f7c86ac6641510c04ca9ee71285a6034ef75b266e832eef420075cebc` |
| N37E049 | `a6871edc445b11990382fdcc00cfc3e93458eec43bb968a31c94ce36e35b6476` |
| N37E050 | `1d3a940c48d18433eff04779b3a6e5942aaead5efaa211781359dfb20c7db16e` |
| N38E048 | `03deb9d540d16b27894eb7565bbf61a2dce07de39a9a57ebd7520f0a3d66fc85` |

### Reproduce without registration

From the project directory, with the authorized local GADM/ERA5 inputs available:

```sh
python download_dem.py
python prepare_dem.py --tile-dir data/srtm --source-record data/srtm/source.json
python analysis.py
python prepare_documents.py
```

`download_dem.py` uses Python's standard library and public HTTPS URLs, checks response length, records hashes/date and writes the acquisition record automatically. It refuses to overwrite existing files; preserve their source record or use a new empty directory with `--tile-dir`. If local TLS configuration prevents Python downloads, the same direct links can be downloaded in a browser, retaining actual acquisition details:

- [N36E049.tif](https://opentopography.s3.sdsc.edu/raster/SRTM_GL1/SRTM_GL1_srtm/N36E049.tif)
- [N36E050.tif](https://opentopography.s3.sdsc.edu/raster/SRTM_GL1/SRTM_GL1_srtm/N36E050.tif)
- [N37E048.tif](https://opentopography.s3.sdsc.edu/raster/SRTM_GL1/SRTM_GL1_srtm/N37E048.tif)
- [N37E049.tif](https://opentopography.s3.sdsc.edu/raster/SRTM_GL1/SRTM_GL1_srtm/N37E049.tif)
- [N37E050.tif](https://opentopography.s3.sdsc.edu/raster/SRTM_GL1/SRTM_GL1_srtm/N37E050.tif)
- [N38E048.tif](https://opentopography.s3.sdsc.edu/raster/SRTM_GL1/SRTM_GL1_srtm/N38E048.tif)

The acquisition record's product DOI/platform pair is validated. Existing EarthExplorer acquisitions can still be used with their actual source record and USGS DOI; do not label public-mirror data as an EarthExplorer acquisition. Raw tiles, acquisition records and masked elevation CSV remain local/ignored by default. The public provenance record contains no credentials.

### Preprocessing and height interpretation

`prepare_dem.py` checks single-band north-up CRS, native resolution, aligned pixel lattices, declared tile bounds and acquisition hashes. It merges aligned tiles in recorded order, clips to a two-native-pixel halo around the target extent and bilinearly samples the existing point centres in original mesh row order. The output is `lon,lat,elevation` plus a `.provenance.json` sidecar. The processed CSV SHA-256 is:

`557f228540f3ebc2b066abfc86fb1dafb337ef091216378630d762e592db449f`


All 57,536 targets have complete valid support. Valid negative coastal elevations are retained. No void filling, extrapolation, zero replacement, artificial 12.5 m resampling or vertical transformation is performed. The 0.005° output spacing is roughly 440 × 555 m here but represents point predictions. Bilinear sampling preserves that estimand; averaging raster cells would instead define areal-mean elevations. See [Rasterio resampling guidance](https://rasterio.readthedocs.io/en/stable/topics/resampling.html).

SRTM's EGM96 orthometric elevations and ERA5 `z/9.80665` geopotential heights share a geoid reference over land but are distinct height quantities. Mean gravity and spatial support differ. Their use in a common fixed lapse correction remains an explicit approximation, not an exact vertical-height conversion. The validated numerical experiment still uses ERA5-derived heights only. No unverified ellipsoidal/geoid adjustment was applied.

The historical CSV is private; the author recalls ASF Vertex/ALOS PALSAR, but its precise source remains unverified. ASF RTC ancillary heights can be transformed and resampled; they are not automatically equivalent to native SRTM. A URL or similar numerical values does not establish that historical provenance.

### Comparison with the private historical CSV

Use private paths for the optional comparison:

```sh
python prepare_dem.py --tile-dir data/srtm --source-record data/srtm/source.json --historical-dem ../private-dem/original_gilan_dem_masked.csv --comparison-dir ../private-dem/comparison
```

The executed comparison has matching coordinates at all 57,536 centres. Replacement minus historical elevation: mean -2.0230 m, mean absolute difference 6.4319 m, RMSE 8.8067 m, maximum absolute difference 84 m, Pearson correlation 0.9999465. Historical-data point tables and spatial/histogram/scatter figures are private. Public aggregate statistics are in [../results/dem_comparison_summary.json](../results/dem_comparison_summary.json).

For unchanged training inputs, `new_IDW - old_IDW = -0.0065 * (new_height - old_height)` at each target in every period. Regenerated maps verified this relation: mean change +0.01315 °C, mean absolute change 0.04181 °C, RMSE 0.05724 °C, maximum absolute change 0.546 °C. Larger differences are localized; the unknown old datum/sampling prevents assigning them to one processing cause. The four non-IDW surfaces are exactly unchanged, as are all 575 predictions, 25 metrics and 75 rankings.

### Attribution and redistribution

Retain NASA/NGA SRTM source attribution, the OpenTopography DOI/access citation above, and the USGS product overview where relevant. Public SRTM data and OpenTopography's stated data policy permit reuse subject to source/license compliance; no third-party rights are replaced by this project's MIT software license. The current figures are regenerated from documented sources. GADM explicitly permits academic map publication; its raw boundary and the masked CSV remain excluded by default. Historical uncertain data and detailed comparison artifacts are not public.
