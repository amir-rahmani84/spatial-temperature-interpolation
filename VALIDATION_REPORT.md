# Final validation — November 2026

**SRTM replacement, all five-method maps, numerical validation and LaTeX compilation are complete locally.** No remote repository was created, changed or pushed. Original software is MIT-licensed, copyright 2026 Seyedamir Rahmani. This software license does not override third-party data rights.

## Scope and preserved inputs

This remains supervised undergraduate numerical analysis research by Seyedamir Rahmani at IUST, supervised by Mahboubeh Molavi-Arabshahi, adapting Yang and Xing (2021). The author's supplied concise README, standard MIT license and single-column A4 LaTeX report were retained as the foundation. The prior DOCX and obsolete document generator are private historical copies; the public report has one authoritative source, `report/research_report.tex`.

The original project and pre-replacement corrected project are preserved outside the clean release folder. The historical processed DEM, detailed comparison points/figures and original reports remain private. The original raster and acquisition metadata are unavailable; the author's recollection of ASF Vertex/ALOS PALSAR is recorded as unverified. It is not used to establish the replacement's provenance.

## Actual replacement acquisition

After the author requested a free alternative to the USGS registration form, the six required SRTM GL1 v3 tiles were obtained through **OpenTopography's documented unsigned public mirror on 27 November 2026**. No account, personal profile, API key or token was supplied. The actual distribution is cited as [NASA SRTM Global, distributed by OpenTopography, DOI 10.5069/G9445JDF](https://doi.org/10.5069/G9445JDF). The [USGS 1 Arc-Second Global overview](https://doi.org/10.5066/F7PR7TFT) remains product context; these are not represented as EarthExplorer downloads or asserted byte-identical to a separate USGS GeoTIFF export.

All six files are actual GeoTIFFs with 3601 × 3601 point postings, EPSG:4326 horizontal coordinates, 1/3600-degree spacing, source EGM96 orthometric elevations in metres, and NoData -32768. Headers, dimensions, spacing, shared-edge bounds, lattice alignment, complete file lengths and SHA-256 hashes were inspected. Vertical reference comes from the documented GL1 source, rather than inference from a filename or an EPSG:4326 horizontal tag alone. GL3 and GL1 Ellipsoidal were not substituted.

The required tiles are N36E049, N36E050, N37E048, N37E049, N37E050 and N38E048. They were derived from the actual Gilan mask and 57,536 target centres, not from guessed elevation values. The public `results/dem_provenance.json` contains actual acquisition date, direct URLs, file hashes, CRS, transforms, tags, shapes, NoData, Rasterio/GDAL versions and the processed CSV hash. `results/dem_acquisition_plan.json` records completed acquisition and coverage. No date, hash or elevation was invented.

## Preprocessing and physical assumptions

The existing `prepare_dem.py` pipeline was used. It merges aligned tiles on the native lattice, clips to a two-native-pixel halo around the required extent and bilinearly samples the existing grid centres in their original mesh row order. All **57,536 targets** have complete valid raster support; duplicate targets and nonfinite elevations were rejected. No void filling, extrapolation, zero substitution, vertical transformation or artificial 12.5 m upsampling was introduced. Valid negative Caspian-coast elevations were retained.

Bilinear point sampling preserves the existing point-prediction estimand. The 0.005° output spacing does not define an areal mean or provide new observational resolution. SRTM EGM96 orthometric elevations and ERA5 EGM96-referenced geopotential heights `z/9.80665` are distinct height quantities. Their common fixed lapse-rate use remains an approximation; no exact height equivalence or local-gravity conversion is claimed.

The loader accepts only declared USGS or documented OpenTopography GL1 source/DOI pairs with matching CSV hashes. This input-provenance accommodation does not change any interpolation equation, parameter, seed, training fold or criterion configuration.

## Measured replacement effects

Differences use replacement minus historical values at matching locations:

| Statistic | Elevation difference (m) | IDW temperature change (°C) |
|---|---:|---:|
| Mean signed difference | -2.023029 | +0.013150 |
| Mean absolute difference | 6.431851 | 0.041807 |
| RMSE | 8.806733 | 0.057244 |
| Maximum absolute difference | 84.000000 | 0.546000 |

Elevation Pearson correlation is **0.999946535**. The 5th–95th percentile elevation differences are -16 to +11 m. Spatial/scatter/histogram plots were inspected privately: the largest discrepancies are localized rather than a uniform offset. Unknown historical raster sampling/datum prevents attributing them to one processing cause or identifying the old source from similarity.

All five annual/seasonal IDW surfaces were regenerated. The measured changes agree with the exact unchanged-model relation `delta_temperature = -0.0065 * delta_height` at each target, identically across periods. The four non-IDW surfaces are exactly unchanged. Aggregate differences are in `results/dem_comparison_summary.json` and map comparisons in `results/map_comparison_summary.json`; historical point data and figures are not redistributed.

## Numerical experiment unchanged

| Export | Records | Maximum numerical difference | File bytes versus verified baseline |
|---|---:|---:|---|
| Cross-validation predictions | 575 | 0.0 | Identical |
| Aggregate evaluation metrics | 25 | 0.0 | Identical |
| TOPSIS rankings | 75 | 0.0 | Identical |
| TOPSIS criterion weights | 50 | 0.0 | Identical |

`results/verified_numerical_baseline.json` records numerical CSV hashes/comparison evidence. `config.json` remains byte-identical to the independently verified configuration. The interpolation and statistical evaluation definitions are unchanged. The external DEM affects map targets only; cross-validation still uses ERA5 geopotential-derived reference/withheld heights.

The verified annual Multiquadric RBF MAE is 0.322690 °C, RMSE 0.425185 °C and NSE 0.970996. Winter IDW/Multiquadric MAEs remain nearly tied. Historical versus primary/reduced TOPSIS criteria still change the winter winner. There is no new independent ground-station validation, tuning experiment or claim of general method superiority.

## Tests and execution

All **16 automated tests pass**, including:

- IDW exact matches, lapse adjustment and constant-temperature elevation behavior.
- Kriging recovery and unit-sum constraint; duplicate-coordinate rejection.
- Hand-calculated metrics, undefined constant cases and nonfinite/mismatched inputs.
- Tied normal-score inversion and negative/constant/tied TOPSIS criteria.
- Repeatable bounded sampling, independent deterministic worker streams and withheld-temperature exclusion.
- Synthetic bilinear sampling, native-lattice seams, row order, invalid CRS/NoData rejection and comparison statistics.
- The actual `analysis.py --no-maps` subprocess with only ERA5/GADM fixture inputs and no DEM file.
- Rejection of an undocumented historical CSV and safe optional-backend version reporting. Missing optional distribution metadata is recorded as JSON null, rather than failing after successful analysis.

The final notebook executed sequentially from a fresh kernel in `MAP_MODE='full'`, regenerating the complete experiment and maps. Evaluation-only mode remains available by selecting `none` or using `analysis.py --no-maps`. The test output is in `results/numerical_validation.txt`. Jupyter's local Windows event-loop/transport notices did not invalidate results. `pip check` reports no broken requirements.

Python 3.14.0 on Windows and the pinned direct dependencies were tested. Other operating systems/Python versions are not claimed compatible. The optional authenticated CDS downloader was not run; its request values were checked against the official current collection form. Actual SRTM acquisition used public HTTPS downloads with length/hash checks; the credential-free downloader and direct links reproduce that route.

## Figures and report

All **25 method-period surfaces** are represented in five temperature figures, with SRTM-adjusted IDW and four unchanged configurations. Coordinates are matched to elevations explicitly. Each period uses a common scale spanning all five surfaces. A previous display-only defect clipped cold mountain IDW values; that was corrected without changing predictions. Color limits and image hashes are in `results/map_manifest.json`.

Temperature and correlation figures use two rows so their labels remain readable in single-column A4. The supplied report's writing, equations, tables, authorship and supervisor acknowledgment were preserved. Updates concern actual replacement provenance, measured effects, figure captions, parameter-unit clarification and pagination/legibility. Its revision note explicitly distinguishes the 2026 replacement from the original 2025 study.

The authoritative LaTeX source was compiled twice with **pdfLaTeX (pdfTeX 1.40.29, TeX Live 2026)**. Because the system MiKTeX installation had missing runtime DLLs, a private workspace TinyTeX installation was used; no system TeX installation was modified. All rendered pages were inspected, including larger views of map, metric, reference and seasonal pages. The final PDF retains single-column A4, readable labels, intact figures/tables, no blank pages, no overfull boxes and no unresolved citations/references. A compact bibliography avoids an otherwise nearly empty reference spill page. Hyperlinks are retained.

`prepare_documents.py` now only compiles the supplied `.tex` source; it cannot overwrite README, notebook, source text or a DOCX. Previous DOCX/ReportLab dependencies were removed. TeX build intermediates and private rendering/runtime files are excluded from the public payload.

## Data, licensing and remaining limitations

- Original software is under the supplied **MIT LICENSE**, copyright 2026 Seyedamir Rahmani. No competing software license is introduced.
- ERA5 attribution and modification notice are retained. The [official current catalogue](https://cds.climate.copernicus.eu/api/catalogue/v1/collections/reanalysis-era5-single-levels-monthly-means) declares CC-BY-4.0. The original acquisition receipt remains unavailable; none was fabricated. Raw NetCDF remains local by default.
- NASA/NGA SRTM and OpenTopography are acknowledged with the actual distribution DOI/date. OpenTopography's published data/citation policy and upstream SRTM public-domain context were checked. Raw tiles and the GADM-masked CSV remain excluded by default; provenance and reproducible acquisition instructions are public.
- [GADM](https://gadm.org/license.html) permits academic map publication while requiring permission for raw redistribution. Raw GADM is absent from the release. The MIT software license does not grant rights over third-party data, report or figures.
- Historical DEM provenance remains unknown, but the historical data/derived comparison plots are private and are not a prerequisite for reconstructing the replacement maps. Prior unverified bibliography entries are not used as authority; their earlier audit is preserved privately.

Scientific limitations remain: 23 spatial locations; smoothed ERA5 reconstruction rather than measured ground-truth validation; spatially dependent LOOCV; equal-month pooled DJF; fixed untuned variograms/lapse rate; exclusive elevation use by IDW; distinct height quantities/spatial supports; raw-kernel RBF sensitivity; short unconverged pseudo-target chains and inverse-score clipping; redundant metrics and criterion-dependent rankings. No uncertainty model or broader research contribution was added.

The public payload is checked for correct relative paths, local links, report/figure consistency, hashes, private files, credentials and restricted raw data. Its exact contents/hashes are listed in `public_file_manifest.csv`. The parent workspace's old Git history is not certified for publication; use the clean release folder. No data-acquisition or compilation blocker remains for that folder. Remote publication still requires the author's approval and was not performed.
