"""Prepare documented USGS SRTM elevations at the existing map centres.

No download, inferred source identification, void filling or vertical transformation
is performed. EarthExplorer acquisitions must include the source record described
in data/README.md. Use --plan before downloading.
"""
from pathlib import Path
from contextlib import ExitStack
from datetime import date
import argparse
import json
import hashlib
import math
import os
import numpy as np
import pandas as pd
import rasterio
from rasterio.merge import merge
from rasterio.windows import from_bounds, bounds as window_bounds, Window
import analysis as a

DOI='10.5066/F7PR7TFT'
PUBLIC_DOI='10.5069/G9445JDF'
PUBLIC_PLATFORM='https://portal.opentopography.org/raster?opentopoID=OTSRTM.082015.4326.1'
SOURCES={
    ('https://earthexplorer.usgs.gov/',DOI): 'USGS SRTM 1 Arc-Second Global',
    (PUBLIC_PLATFORM,PUBLIC_DOI): 'NASA SRTM GL1 v3, distributed by OpenTopography',
}
PIXEL=1/3600

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def required_tiles(points):
    # SRTM tiles have shared edge postings; a point on a degree boundary is
    # covered by the adjacent named tile's first posting. No cell-area averaging.
    return sorted({f'N{math.floor(lat):02d}E{math.floor(lon):03d}' for lon,lat in np.round(points,10)})

def plan(cfg,boundary):
    _,_,_,points=a.target_grid(cfg,boundary)
    return dict(product='NASA SRTM GL1 v3, distributed by OpenTopography',product_doi=PUBLIC_DOI,
                platform=PUBLIC_PLATFORM,
                boundary_bounds_west_south_east_north=list(boundary.bounds),
                target_bounds_west_south_east_north=[float(points[:,0].min()),float(points[:,1].min()),float(points[:,0].max()),float(points[:,1].max())],
                target_count=len(points),tile_ids=required_tiles(points),
                sampling='bilinear at existing target centres; no area averaging',
                horizontal_crs='EPSG:4326',vertical_datum='EGM96 orthometric',units='metres',
                acquisition_status='coverage plan only; --plan performs no download',input_sha256={})

def bilinear(array,transform,points,nodata):
    inv=~transform
    col=inv.a*points[:,0]+inv.b*points[:,1]+inv.c-.5
    row=inv.d*points[:,0]+inv.e*points[:,1]+inv.f-.5
    # Snap negligible floating error at native posting centres before indexing.
    col=np.where(abs(col-np.round(col))<1e-7,np.round(col),col)
    row=np.where(abs(row-np.round(row))<1e-7,np.round(row),row)
    c0=np.floor(col).astype(int); r0=np.floor(row).astype(int)
    dc=col-c0; dr=row-r0; values=np.zeros(len(points)); invalid=np.zeros(len(points),bool)
    for dy,dx,w in [(0,0,(1-dr)*(1-dc)),(0,1,(1-dr)*dc),(1,0,dr*(1-dc)),(1,1,dr*dc)]:
        rr=r0+dy; cc=c0+dx; used=w>1e-12
        inside=(rr>=0)&(cc>=0)&(rr<array.shape[0])&(cc<array.shape[1])
        sample=np.full(len(points),np.nan); sample[inside]=array[rr[inside],cc[inside]]
        bad=~np.isfinite(sample) | (sample==nodata)
        invalid |= used & bad
        values += np.where(used & ~bad,w*sample,0)
    if invalid.any():
        example=np.round(points[invalid][:5],8).tolist()
        raise ValueError(f'{invalid.sum()} targets lack complete valid SRTM support (including voids); examples {example}. No silent filling.')
    return values

def sample_rasters(paths,points):
    """Merge on the original 1-arcsecond lattice, clipping to a padded region.

    The two-pixel halo supports bilinear interpolation at the outer target centres.
    Nearest resampling during the merge only copies aligned native postings.
    """
    with ExitStack() as stack:
        sources=[stack.enter_context(rasterio.open(p)) for p in paths]
        if not sources: raise FileNotFoundError('No SRTM GeoTIFF tiles supplied')
        first=sources[0]; metadata=[]
        for p,s in zip(paths,sources):
            if s.crs != rasterio.crs.CRS.from_epsg(4326): raise ValueError(f'{Path(p).name}: require WGS84 EPSG:4326')
            if s.count!=1 or s.transform.b!=0 or s.transform.d!=0 or s.transform.e>=0:
                raise ValueError(f'{Path(p).name}: require single-band north-up raster')
            if not np.allclose(s.res,(PIXEL,PIXEL),rtol=0,atol=1e-12):
                raise ValueError(f'{Path(p).name}: require 1 arc-second spacing; no substitute products')
            offset=(np.array([s.transform.c-first.transform.c,s.transform.f-first.transform.f])/PIXEL)
            if not np.allclose(offset,np.round(offset),atol=1e-6): raise ValueError('Tiles have misaligned pixel lattices')
            metadata.append(dict(filename=Path(p).name,crs=s.crs.to_string(),resolution=list(s.res),shape=[s.height,s.width],
                                 transform=list(s.transform),nodata=s.nodata,tags=s.tags(),band_tags=s.tags(1),sha256=sha256(p)))
        box=(points[:,0].min()-2*PIXEL,points[:,1].min()-2*PIXEL,points[:,0].max()+2*PIXEL,points[:,1].max()+2*PIXEL)
        w=from_bounds(*box,first.transform)
        c0=math.floor(w.col_off); r0=math.floor(w.row_off)
        w=Window(c0,r0,math.ceil(w.col_off+w.width)-c0,math.ceil(w.row_off+w.height)-r0)
        clip=window_bounds(w,first.transform)
        mosaic,transform=merge(sources,bounds=clip,res=(PIXEL,PIXEL),nodata=-32768,dtype='float32',method='first')
        values=bilinear(mosaic[0],transform,points,-32768)
    return values,metadata

def compare(historical,replacement,directory,lapse_rate):
    """Private detailed differences; public summary contains no historical values."""
    old=pd.read_csv(historical); new=replacement.copy()
    for frame in (old,new):
        frame['lon']=frame.lon.round(8); frame['lat']=frame.lat.round(8)
        if frame.duplicated(['lon','lat']).any(): raise ValueError('Duplicate comparison coordinates')
    joined=new.merge(old,on=['lon','lat'],suffixes=('_replacement','_historical'),validate='one_to_one')
    if len(joined)!=len(new) or len(old)!=len(new): raise ValueError('DEM grid locations differ')
    d=joined.elevation_replacement.to_numpy()-joined.elevation_historical.to_numpy()
    def summary(v): return dict(mean_signed_difference=float(v.mean()),mean_absolute_difference=float(abs(v).mean()),rmse=float(np.sqrt(np.mean(v*v))),max_absolute_difference=float(abs(v).max()),min=float(v.min()),max=float(v.max()),median=float(np.median(v)),p05=float(np.quantile(v,.05)),p95=float(np.quantile(v,.95)))
    stats=dict(n_locations=len(joined),difference_definition='replacement minus historical, metres',elevation=summary(d),
               elevation_pearson_correlation=float(np.corrcoef(joined.elevation_replacement,joined.elevation_historical)[0,1]) if joined.elevation_replacement.std()>0 and joined.elevation_historical.std()>0 else None,
               idw_temperature_change=summary(lapse_rate*d),
               temperature_change_rationale='For unchanged training values, coordinates and lapse rate, IDW new minus old = lapse_rate*(new height - old height), identically for every period; other four methods use no DEM.',
               cause_assessment='Inspect spatial plot and metadata. Unknown historical datum/source prevents causal attribution; an offset alone does not establish ellipsoidal/geoid correction.')
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    joined['elevation_difference']=d; joined['idw_temperature_difference']=lapse_rate*d
    joined.to_csv(directory/'dem_comparison_points.csv',index=False)
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    im=axes[0].scatter(joined.lon,joined.lat,c=d,s=1,cmap='coolwarm'); fig.colorbar(im,ax=axes[0],label='Replacement minus historical (m)')
    axes[0].set(xlabel='Longitude',ylabel='Latitude',title='Elevation differences')
    axes[1].hist(d,bins=60); axes[1].set(xlabel='Difference (m)',ylabel='Count')
    axes[2].scatter(joined.elevation_historical,joined.elevation_replacement,s=1); axes[2].set(xlabel='Historical elevation (m)',ylabel='SRTM elevation (m)')
    fig.savefig(directory/'dem_comparison.png',dpi=160); plt.close(fig)
    (directory/'summary.json').write_text(json.dumps(stats,indent=2),encoding='utf-8')
    return stats

def prepare(tile_dir,source_record,output,cfg,boundary,historical=None,comparison_dir=None):
    if historical and not comparison_dir:
        raise ValueError('Specify a private comparison directory for historical DEM comparisons')
    record=json.loads(Path(source_record).read_text(encoding='utf-8-sig'))
    source_key=(record.get('platform'),record.get('product_doi'))
    if source_key not in SOURCES:
        raise ValueError('Source record must identify documented SRTM 1-arcsecond data from USGS or OpenTopography GL1 (not GL3 or ellipsoidal GL1)')
    date.fromisoformat(record['acquired_on'])
    required=plan(cfg,boundary); entries=record['tiles']
    if sorted(e['tile_id'] for e in entries)!=required['tile_ids']: raise ValueError('Source record tile IDs do not match required coverage')
    paths=[]
    for e in entries:
        name=e['filename']
        if Path(name).name!=name: raise ValueError('Use a filename, not an external path in the source record')
        p=Path(tile_dir)/name
        if not p.is_file(): raise FileNotFoundError(p)
        if e.get('sha256') and e['sha256']!=sha256(p): raise ValueError(f'Acquisition hash mismatch: {name}')
        paths.append(p)
    _,_,_,points=a.target_grid(cfg,boundary)
    values,metadata=sample_rasters(paths,points)
    # Check each file's declared tile against its actual raster bounds.
    for entry,meta in zip(entries,metadata):
        lat=int(entry['tile_id'][1:3]); lon=int(entry['tile_id'][4:7])
        with rasterio.open(Path(tile_dir)/entry['filename']) as src:
            if not np.allclose([src.bounds.left,src.bounds.bottom,src.bounds.right,src.bounds.top],[lon,lat,lon+1,lat+1],atol=2*PIXEL,rtol=0):
                raise ValueError(f"Raster bounds inconsistent with {entry['tile_id']}")
        meta['tile_id']=entry['tile_id']
        for key in ('url','verified_at_utc'):
            if key in entry: meta[key]=entry[key]
    replacement=pd.DataFrame(dict(lon=points[:,0],lat=points[:,1],elevation=values))
    output=Path(output); output.parent.mkdir(parents=True,exist_ok=True)
    if historical is not None and output.resolve()==Path(historical).resolve(): raise ValueError('Never overwrite historical DEM')
    replacement.to_csv(output,index=False)
    provenance=dict(product=SOURCES[source_key],product_doi=record['product_doi'],platform=record['platform'],acquired_on=record['acquired_on'],
                    horizontal_crs='EPSG:4326',vertical_datum='EGM96 orthometric',units='metres',native_spacing_arcseconds=1,
                    sampling='bilinear at geographic target centres on original aligned raster lattice; no area aggregation or geoid conversion',
                    tiles=metadata,output_sha256=sha256(output),n_targets=len(points),
                    boundary_sha256=record.get('boundary_sha256'),
                    vertical_compatibility='ERA5 z/9.80665 is EGM96-referenced geopotential height over land, not orthometric height. Approximation retained only for fixed lapse correction; no ellipsoidal transformation or change to LOOCV.',
                    rights='USGS public-domain elevation data; GADM raw geometry excluded. CSV is sampled elevations at an existing study grid, not a redistributable GADM polygon.',
                    rasterio_version=rasterio.__version__,gdal_version=rasterio.__gdal_version__,
                    attribution=('NASA Shuttle Radar Topography Mission (SRTM) (2013). Shuttle Radar Topography Mission (SRTM) Global. Distributed by OpenTopography. https://doi.org/10.5069/G9445JDF' if source_key[1]==PUBLIC_DOI else 'U.S. Geological Survey, Earth Resources Observation and Science (EROS) Center. SRTM 1 Arc-Second Global. https://doi.org/10.5066/F7PR7TFT'))
    output.with_suffix('.provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    if historical:
        stats=compare(historical,replacement,comparison_dir,cfg['lapse_rate_c_per_m'])
        (a.ROOT/cfg['results_dir']).mkdir(parents=True,exist_ok=True)
        (a.ROOT/cfg['results_dir']/'dem_comparison_summary.json').write_text(json.dumps(stats,indent=2),encoding='utf-8')
    result_dir=a.ROOT/cfg['results_dir']; result_dir.mkdir(parents=True,exist_ok=True)
    required['acquisition_status']='documented rasters prepared; see output provenance'
    required['product']=provenance['product']; required['product_doi']=provenance['product_doi']; required['platform']=provenance['platform']
    required['input_sha256']={m['filename']:m['sha256'] for m in metadata}
    (result_dir/'dem_acquisition_plan.json').write_text(json.dumps(required,indent=2),encoding='utf-8')
    return provenance

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',action='store_true'); parser.add_argument('--tile-dir',type=Path)
    parser.add_argument('--source-record',type=Path); parser.add_argument('--output',type=Path)
    parser.add_argument('--historical-dem',type=Path); parser.add_argument('--comparison-dir',type=Path)
    args=parser.parse_args(); cfg=a.configuration()
    folder=Path(os.environ.get('SPATIAL_DATA_DIR',a.ROOT/cfg['data_dir']))
    path=folder/cfg['boundary']; boundary=a.load_boundary(path)
    if args.plan:
        info=plan(cfg,boundary); info['boundary_sha256']=sha256(path)
        print(json.dumps(info,indent=2)); return
    if args.tile_dir is None or args.source_record is None: parser.error('--tile-dir and --source-record are required (or use --plan)')
    # Source record is left unchanged; the generated sidecar records actual masking input.
    info=prepare(args.tile_dir,args.source_record,args.output or a.ROOT/cfg['data_dir']/cfg['dem'],cfg,boundary,args.historical_dem,args.comparison_dir)
    sidecar=(args.output or a.ROOT/cfg['data_dir']/cfg['dem']).with_suffix('.provenance.json')
    info['boundary_sha256']=sha256(path); sidecar.write_text(json.dumps(info,indent=2),encoding='utf-8')
    print(f"Prepared {info['n_targets']} documented target elevations.")

if __name__=='__main__': main()
