"""Notebook helpers for the fixed comparative experiment; distances are kilometres."""
from pathlib import Path
import json
import os
import warnings
import hashlib
import platform
import importlib.metadata
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import shape
from shapely import contains_xy
from pyproj import Transformer
from scipy.spatial.distance import cdist, pdist
from scipy.spatial import KDTree
from scipy.interpolate import Rbf
from scipy.linalg import solve, LinAlgWarning
from scipy.stats import norm, rankdata, spearmanr, pearsonr
from joblib import Parallel, delayed

ROOT = Path(__file__).resolve().parent

def configuration():
    return json.loads((ROOT / 'config.json').read_text(encoding='utf-8'))

def project(lonlat, cfg):
    tr = Transformer.from_crs('EPSG:4326', cfg['projection'], always_xy=True)
    return np.column_stack(tr.transform(*np.asarray(lonlat).T)) / 1000

def load_boundary(path):
    geo = json.loads(Path(path).read_text(encoding='utf-8'))
    feature = [f for f in geo['features'] if f['properties']['NAME_1'] == 'Gilan']
    if len(feature) != 1:
        raise ValueError('Expected one Gilan feature')
    boundary = shape(feature[0]['geometry'])
    if not boundary.is_valid:
        raise ValueError('Invalid province geometry')
    return boundary

def target_grid(cfg, boundary):
    lo,hi,la,ha=cfg['bounds']; step=cfg['grid_step_degrees']
    lon,lat=np.meshgrid(np.arange(lo,hi+step/2,step),np.arange(la,ha+step/2,step))
    mask=contains_xy(boundary,lon,lat)
    return lon,lat,mask,np.column_stack([lon[mask],lat[mask]])

def load_dem(path):
    path=Path(path); provenance=path.with_suffix('.provenance.json')
    if not path.is_file() or not provenance.is_file():
        raise FileNotFoundError('Documented replacement DEM missing. Run prepare_dem.py; see data/README.md. Evaluation uses --no-maps.')
    info=json.loads(provenance.read_text(encoding='utf-8'))
    source=(info.get('platform'),info.get('product_doi'))
    accepted_sources={('https://earthexplorer.usgs.gov/','10.5066/F7PR7TFT'),
                      ('https://portal.opentopography.org/raster?opentopoID=OTSRTM.082015.4326.1','10.5069/G9445JDF')}
    if source not in accepted_sources or info.get('output_sha256')!=hashlib.sha256(path.read_bytes()).hexdigest():
        raise ValueError('DEM provenance/hash mismatch; do not use the historical unverified CSV')
    dem=pd.read_csv(path)
    if dem.duplicated(['lon','lat']).any() or not np.isfinite(dem[['lon','lat','elevation']]).all().all():
        raise ValueError('Invalid DEM rows')
    dem.attrs['provenance']=info
    return dem,info

def load_data(cfg, require_dem=False):
    folder = Path(os.environ.get('SPATIAL_DATA_DIR', ROOT / cfg['data_dir']))
    files = [folder / cfg[k] for k in ('netcdf', 'boundary')]
    for f in files:
        if not f.is_file():
            raise FileNotFoundError(f'{f}: see data/README.md')
    boundary=load_boundary(files[1])
    with xr.open_dataset(files[0]) as source:
        ds = source.sel(valid_time=slice(cfg['start'], cfg['end'])).load()
    if ds.t2m.attrs.get('units') != 'K' or ds.z.attrs.get('units') not in ('m**2 s**-2', 'm2 s-2'):
        raise ValueError('Unexpected temperature/geopotential units')
    dates = pd.DatetimeIndex(ds.valid_time.values)
    expected = pd.date_range(cfg['start'], cfg['end'], freq='MS')
    if not dates.equals(expected):
        raise ValueError('Missing, duplicate or unexpected monthly timestamps')
    frame = ds[['t2m', 'z']].to_dataframe().reset_index()
    frame = frame[contains_xy(boundary, frame.longitude.to_numpy(), frame.latitude.to_numpy())].copy()
    frame['temperature'] = frame.t2m.astype(float) - 273.15
    frame['elevation'] = frame.z.astype(float) / 9.80665
    if not np.isfinite(frame[['temperature', 'elevation']]).all().all():
        raise ValueError('Missing/nonfinite data; no silent imputation')
    if frame.duplicated(['valid_time', 'longitude', 'latitude']).any():
        raise ValueError('Duplicate spatial-month records')
    counts = frame.groupby(['longitude', 'latitude']).size()
    if not (counts == len(expected)).all():
        raise ValueError('Unequal temporal coverage')
    periods = {}
    for name, months in {'annual':range(1,13), 'spring':[3,4,5], 'summer':[6,7,8], 'fall':[9,10,11], 'winter':[12,1,2]}.items():
        periods[name] = frame[frame.valid_time.dt.month.isin(months)].groupby(['longitude','latitude'])[['temperature','elevation']].mean().reset_index()
    dem=None; dem_info=None
    if require_dem:
        dem_path=Path(os.environ.get('SPATIAL_DEM_FILE', ROOT / cfg['data_dir'] / cfg['dem']))
        dem,dem_info=load_dem(dem_path)
    audit = dict(months=len(dates), locations=len(counts), records=len(frame),
                 first=str(dates.min()), last=str(dates.max()), missing=0,
                 temperature_units='K converted using T - 273.15',
                 netcdf_attributes=dict(ds.attrs), variable_attributes={v:dict(ds[v].attrs) for v in ('t2m','z')},
                 dem_required=require_dem,dem_rows=len(dem) if dem is not None else None,
                 dem_provenance=dem_info,
                 input_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in files})
    if dem_info is not None: audit['input_sha256'][cfg['dem']]=dem_info['output_sha256']
    return periods, boundary, dem, audit

def validate_xy(x, y, q):
    x, y, q = np.asarray(x,float), np.asarray(y,float), np.asarray(q,float)
    if x.ndim != 2 or x.shape[1] != 2 or q.ndim != 2 or q.shape[1] != 2 or y.shape != (len(x),):
        raise ValueError('Expected (n,2), (n,), (m,2)')
    if len(x)<3 or not all(np.isfinite(a).all() for a in (x,y,q)):
        raise ValueError('Insufficient or nonfinite inputs')
    if len(np.unique(x,axis=0)) != len(x):
        raise ValueError('Duplicate spatial coordinates: resolve provenance before fitting')
    return x,y,q

def kriging_weights(x, q, nugget, sill, distance_range, spherical=False):
    if nugget<0 or sill<=0 or distance_range<=0:
        raise ValueError('Invalid covariance parameters')
    def cov(d):
        r=d/distance_range
        return sill * (np.where(r<1, 1-1.5*r+0.5*r**3, 0) if spherical else np.exp(-r))
    n=len(x); k=np.zeros((n+1,n+1)); k[:n,:n]=cov(cdist(x,x))
    k[np.arange(n),np.arange(n)] += nugget
    k[n,:n]=1; k[:n,n]=1
    v=np.vstack([cov(cdist(x,q)),np.ones(len(q))])
    # Symmetric indefinite augmented covariance system: no Cholesky or diagonal jitter.
    with warnings.catch_warnings():
        warnings.simplefilter('error', LinAlgWarning)
        try:
            w=solve(k,v,assume_a='sym')
        except (np.linalg.LinAlgError,LinAlgWarning):
            w,_,rank,_=np.linalg.lstsq(k,v,rcond=1e-12)
            if rank < n+1 or not np.allclose(k@w,v,rtol=1e-8,atol=1e-10):
                raise ValueError('Singular or inconsistent kriging system')
    if not np.allclose(w[:-1].sum(axis=0),1,atol=1e-8):
        raise ValueError('Kriging weights violate unbiasedness')
    return w[:-1]

def normal_scores(y):
    return norm.ppf(rankdata(y,method='average')/(len(y)+1))

def inverse_scores(z,y):
    u=np.unique(y)
    scores=normal_scores(y)
    return np.interp(z,[scores[y==v][0] for v in u],u)

def empirical_variogram(x,z,n_lags):
    d=pdist(x); g=pdist(z[:,None])**2/2
    edges=np.linspace(0,d.max(),n_lags+1)
    bins=np.minimum(np.searchsorted(edges,d,side='right')-1,n_lags-1)
    return np.array([[d[bins==i].mean(),g[bins==i].mean()] for i in range(n_lags) if (bins==i).any()]).T

def sample_parameters(h,g,cfg,rng):
    bounds=np.array(cfg['parameter_bounds'])
    def target(p):
        if np.any(p<bounds[:,0]) or np.any(p>bounds[:,1]): return -np.inf
        model=p[0]+p[1]*(-np.expm1(-h/p[2]))
        return -np.sum(((g-model)/cfg['pseudo_residual_scale'])**2)
    current=np.array(cfg['initial_parameters'],float); score=target(current)
    samples=[]; accepted=0
    for i in range(cfg['mh_iterations']):
        proposal=current+rng.normal(0,cfg['proposal_std'],3)
        new_score=target(proposal)
        if np.log(rng.random())<new_score-score:
            current=proposal; score=new_score; accepted+=1
        if i>=cfg['mh_burn_in']: samples.append(current.copy())
    if not samples: raise ValueError('Burn-in consumes entire chain')
    samples=np.array(samples)
    return samples[rng.choice(len(samples),cfg['sample_models'],replace=True)],accepted/cfg['mh_iterations']

def local_group(idx, rows, x, y, q, z, cfg, seed):
    rng=np.random.default_rng(np.random.SeedSequence([seed,*idx]))
    h,g=empirical_variogram(x[idx],z[idx],cfg['n_lags'])
    if len(h)<3: raise ValueError('Insufficient occupied variogram bins')
    samples,acceptance=sample_parameters(h,g,cfg,rng)
    # Model-wise inverse transforms precede averaging; no predictive SD is claimed.
    predictions=[inverse_scores(z[idx]@kriging_weights(x[idx],q[rows],*p),y) for p in samples]
    return rows,np.mean(predictions,axis=0),acceptance

def predict(method,x,y,q,elev,qelev,cfg,seed=0):
    x,y,q=validate_xy(x,y,q)
    if np.ptp(y)==0 and method in ('Ordinary Kriging','EBK-inspired local kriging'):
        return np.full(len(q),y[0])
    if method=='IDW':
        residual=y-cfg['lapse_rate_c_per_m']*np.asarray(elev)
        values=[]
        for start in range(0,len(q),2000):
            d=cdist(q[start:start+2000],x); exact=d==0
            d=np.where(exact,1,d); w=(d.min(axis=1)[:,None]/d)**cfg['idw_power']
            hits=exact.any(axis=1); w[hits]=exact[hits]
            values.extend(w@residual/w.sum(axis=1))
        result=np.array(values)+cfg['lapse_rate_c_per_m']*np.asarray(qelev)
    elif method in cfg['rbf_kernels']:
        with warnings.catch_warnings():
            warnings.simplefilter('error',LinAlgWarning)
            f=Rbf(*x.T,y,function=cfg['rbf_kernels'][method],smooth=cfg['rbf_smooth'])
        result=f(*q.T)
    elif method=='Ordinary Kriging':
        result=y@kriging_weights(x,q,0,float(np.var(y)),cfg['ok_spherical_range_km'],True)
    elif method=='EBK-inspired local kriging':
        z=normal_scores(y); k=min(cfg['neighbors'],len(x))
        if k<4: raise ValueError('Insufficient local neighborhood')
        indices=KDTree(x).query(q,k=k)[1]
        groups={}
        for i,idx in enumerate(indices): groups.setdefault(tuple(sorted(idx)),[]).append(i)
        jobs=[(np.array(idx),np.array(rows)) for idx,rows in groups.items()]
        ans=Parallel(n_jobs=cfg['workers'],prefer='threads')(delayed(local_group)(idx,rows,x,y,q,z,cfg,seed) for idx,rows in jobs)
        result=np.empty(len(q))
        for rows,vals,_ in ans: result[rows]=vals
    else: raise ValueError(method)
    result=np.asarray(result).reshape(-1)
    if result.shape!=(len(q),) or not np.isfinite(result).all(): raise ValueError('Invalid prediction shape or values')
    return result

def metrics(reference,predicted):
    a,b=np.asarray(reference,float),np.asarray(predicted,float)
    if a.ndim!=1 or a.shape!=b.shape or len(a)<2: raise ValueError('Incompatible metric arrays')
    if not np.isfinite(a).all() or not np.isfinite(b).all(): raise ValueError('Nonfinite predictions; aggregate invalid')
    err=b-a; mse=np.mean(err**2); den=np.sum((a-a.mean())**2)
    return dict(MAE=np.mean(abs(err)),RMSE=np.sqrt(mse),MSE=mse,NSE=1-np.sum(err**2)/den if den>0 else np.nan,
                Spearman=spearmanr(a,b).statistic if np.ptp(a)>0 and np.ptp(b)>0 else np.nan,
                Pearson=pearsonr(a,b).statistic if np.ptp(a)>0 and np.ptp(b)>0 else np.nan,
                MAPE=100*np.mean(abs(err/a)) if np.all(a!=0) else np.nan,
                SMAPE=200*np.mean(abs(err)/(abs(a)+abs(b))) if np.all(abs(a)+abs(b)>0) else np.nan,
                min_absolute_error=np.min(abs(err)),max_absolute_error=np.max(abs(err)))

def cross_validate(periods,cfg):
    rows=[]
    for pi,(period,df) in enumerate(periods.items()):
        x=project(df[['longitude','latitude']],cfg); y=df.temperature.to_numpy(); elev=df.elevation.to_numpy()
        for i in range(len(df)):
            train=np.arange(len(df))!=i
            for mi,method in enumerate(cfg['methods']):
                failure=''; value=np.nan
                try:
                    value=float(predict(method,x[train],y[train],x[[i]],elev[train],elev[[i]],cfg,seed=cfg['seed']+1000*pi+10*i+mi)[0])
                except (ValueError,np.linalg.LinAlgError,LinAlgWarning) as exc: failure=f'{type(exc).__name__}: {exc}'
                rows.append(dict(period=period,method=method,location_id=i,longitude=df.longitude[i],latitude=df.latitude[i],reference=y[i],prediction=value,error=value-y[i],failure=failure))
    pred=pd.DataFrame(rows); aggregates=[]
    for (period,method),g in pred.groupby(['period','method'],sort=False):
        nfailed=int((g.failure!='').sum()); values=metrics(g.reference,g.prediction) if nfailed==0 else {}
        aggregates.append(dict(period=period,method=method,n_expected=len(g),n_failed=nfailed,status='complete' if nfailed==0 else 'invalid',**values))
    return pred,pd.DataFrame(aggregates)

def topsis(df,criteria):
    a=df[list(criteria)].to_numpy(float)
    if not np.isfinite(a).all(): raise ValueError('TOPSIS requires all criteria for every alternative')
    span=np.ptp(a,axis=0); active=span>0
    normalized=np.zeros_like(a)
    for j,direction in enumerate(criteria.values()):
        if active[j]: normalized[:,j]=(a[:,j]-a[:,j].min())/span[j] if direction==1 else (a[:,j].max()-a[:,j])/span[j]
    sums=normalized.sum(axis=0); p=np.divide(normalized,sums,out=np.zeros_like(a),where=sums>0)
    plog=np.zeros_like(p); mask=p>0; plog[mask]=p[mask]*np.log(p[mask])
    divergence=np.where(active,1+plog.sum(axis=0)/np.log(len(a)),0)
    weights=divergence/divergence.sum() if divergence.sum()>0 else np.zeros(a.shape[1])
    v=normalized*weights; plus=np.linalg.norm(v-v.max(axis=0),axis=1); minus=np.linalg.norm(v-v.min(axis=0),axis=1)
    c=np.divide(minus,plus+minus,out=np.full(len(a),0.5),where=(plus+minus)>0)
    return pd.DataFrame(dict(method=df.method.to_numpy(),closeness=c,rank=pd.Series(c).rank(ascending=False,method='min').to_numpy(),distance_positive=plus,distance_negative=minus)),dict(zip(criteria,weights))

def rankings(aggregate,cfg):
    rows=[]; weights=[]
    for period,df in aggregate.groupby('period',sort=False):
        if (df.n_failed>0).any(): raise ValueError('Failed methods preclude comparable ranking')
        for name,criteria in cfg['topsis_criteria'].items():
            r,w=topsis(df,criteria); r['period']=period; r['configuration']=name; rows.append(r)
            weights.extend(dict(period=period,configuration=name,criterion=k,weight=v) for k,v in w.items())
    return pd.concat(rows,ignore_index=True),pd.DataFrame(weights)

def output_dirs(cfg):
    for p in [ROOT/cfg['results_dir'],*[ROOT/cfg['figures_dir']/k for k in ('temperature_maps','error_analysis','correlation','rankings')]]: p.mkdir(parents=True,exist_ok=True)

def make_maps(periods,boundary,dem,cfg):
    import matplotlib.pyplot as plt
    if dem is not None and not dem.attrs.get('provenance'):
        raise ValueError('Map elevations require the documented SRTM provenance sidecar')
    lon,lat,mask,ll=target_grid(cfg,boundary); q=project(ll,cfg)
    methods=cfg['methods'] if dem is not None else [m for m in cfg['methods'] if m!='IDW']
    qelev=np.zeros(len(q))
    if dem is not None:
        d,idx=KDTree(dem[['lon','lat']]).query(ll)
        if d.max()>cfg['grid_step_degrees']/10: raise ValueError(f'DEM does not align with map grid: maximum offset {d.max()} degrees')
        qelev=dem.elevation.to_numpy()[idx]
    color_limits={}
    for pi,(period,df) in enumerate(periods.items()):
        fig,panel_axes=plt.subplots(2,3,figsize=(11,8),layout='constrained')
        axes=panel_axes.ravel()[:len(methods)]
        for unused in panel_axes.ravel()[len(methods):]: unused.set_visible(False)
        predictions={}
        for method in methods:
            mi=cfg['methods'].index(method)
            val=predict(method,project(df[['longitude','latitude']],cfg),df.temperature.to_numpy(),q,df.elevation.to_numpy(),qelev,cfg,cfg['seed']+1000*pi+mi)
            predictions[method]=val
        low=float(np.floor(min(v.min() for v in predictions.values())))
        high=float(np.ceil(max(v.max() for v in predictions.values())))
        color_limits[period]=[low,high]
        for method,ax in zip(methods,axes):
            grid=np.full(lon.shape,np.nan); grid[mask]=predictions[method]
            im=ax.pcolormesh(lon,lat,grid,cmap='coolwarm',vmin=low,vmax=high,shading='auto',rasterized=True)
            ax.set_title(method.replace(' ','\n',1),fontsize=13); ax.set_xlabel('Longitude',fontsize=12); ax.set_ylabel('Latitude',fontsize=12); ax.tick_params(labelsize=12); ax.set_aspect(1/np.cos(np.deg2rad(37.5)))
        colorbar=fig.colorbar(im,ax=axes,label='Temperature (Celsius)')
        colorbar.ax.tick_params(labelsize=12); colorbar.set_label('Temperature (Celsius)',size=12)
        fig.suptitle(f'{period.title()} climatology, 2015-2024',fontsize=15)
        fig.savefig(ROOT/cfg['figures_dir']/'temperature_maps'/f'{period}.png',dpi=140); plt.close(fig)
        prefix='map_' if dem is not None else 'map_dem_independent_'
        np.savez_compressed(ROOT/cfg['results_dir']/f'{prefix}{period}.npz',longitude=ll[:,0],latitude=ll[:,1],**predictions)
        print('Maps complete:',period,flush=True)
    manifest=dict(methods=methods,external_dem_used=dem is not None,color_limits_celsius=color_limits,
                  dem_provenance=dem.attrs['provenance'] if dem is not None else None,
                  figures={f'{period}.png':hashlib.sha256((ROOT/cfg['figures_dir']/'temperature_maps'/f'{period}.png').read_bytes()).hexdigest() for period in periods})
    (ROOT/cfg['results_dir']/'map_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')

def diagnostic_figures(pred,rank,cfg):
    import matplotlib.pyplot as plt
    for period,g in pred.groupby('period',sort=False):
        fig,ax=plt.subplots(figsize=(9,4),layout='constrained')
        data=[abs(g[g.method==m].error.to_numpy()) for m in cfg['methods']]
        ax.boxplot(data,tick_labels=cfg['methods'],showmeans=True,whis=(10,90)); ax.tick_params(axis='x',labelrotation=15); ax.set_ylabel('Absolute prediction error (Celsius)'); ax.set_title(period)
        fig.savefig(ROOT/cfg['figures_dir']/'error_analysis'/f'{period}.png',dpi=140); plt.close(fig)
        fig,panel_axes=plt.subplots(2,3,figsize=(11,8),layout='constrained')
        axes=panel_axes.ravel()[:5]; panel_axes.ravel()[5].set_visible(False)
        for method,ax in zip(cfg['methods'],axes):
            s=g[g.method==method]; ax.scatter(s.reference,s.prediction,s=18); lo=min(s.reference.min(),s.prediction.min()); hi=max(s.reference.max(),s.prediction.max()); ax.plot([lo,hi],[lo,hi],'k--'); ax.set_title(method.replace(' ','\n',1)+f'\nrho={spearmanr(s.reference,s.prediction).statistic:.3f}',fontsize=12); ax.set_xlabel('Reference (Celsius)',fontsize=12); ax.set_ylabel('Prediction (Celsius)',fontsize=12); ax.tick_params(labelsize=12)
        fig.savefig(ROOT/cfg['figures_dir']/'correlation'/f'{period}.png',dpi=140); plt.close(fig)
        fig,ax=plt.subplots(figsize=(9,4),layout='constrained')
        for name in cfg['topsis_criteria']:
            s=rank[(rank.period==period)&(rank.configuration==name)].set_index('method').loc[cfg['methods']]
            ax.plot(cfg['methods'],s['rank'],'o-',label=name)
        ax.invert_yaxis(); ax.set_ylabel('TOPSIS rank'); ax.tick_params(axis='x',labelrotation=15); ax.legend(); ax.set_title(period)
        fig.savefig(ROOT/cfg['figures_dir']/'rankings'/f'{period}.png',dpi=140); plt.close(fig)

def dependency_versions():
    """Report distributions without requiring every optional xarray backend."""
    versions={}
    for name in ('numpy','pandas','scipy','xarray','netCDF4','h5netcdf','shapely','pyproj','matplotlib','joblib'):
        try:
            versions[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name]=None
    return versions

def run(maps=True, dem_independent=False):
    cfg=configuration(); output_dirs(cfg); periods,boundary,dem,audit=load_data(cfg,require_dem=maps and not dem_independent)
    (ROOT/cfg['results_dir']/'data_audit.json').write_text(json.dumps(audit,indent=2,default=str),encoding='utf-8')
    pred,aggregate=cross_validate(periods,cfg)
    pred.to_csv(ROOT/cfg['results_dir']/'cross_validation_predictions.csv',index=False)
    aggregate.to_csv(ROOT/cfg['results_dir']/'cross_validation_metrics.csv',index=False)
    rank,weights=rankings(aggregate,cfg)
    rank.to_csv(ROOT/cfg['results_dir']/'topsis_rankings.csv',index=False); weights.to_csv(ROOT/cfg['results_dir']/'topsis_weights.csv',index=False)
    diagnostic_figures(pred,rank,cfg)
    if maps: make_maps(periods,boundary,dem,cfg)
    versions=dependency_versions()
    manifest=dict(python=platform.python_version(),dependencies=versions,configuration=cfg,input_sha256=audit['input_sha256'],maps_executed=maps,
                  map_methods=([m for m in cfg['methods'] if m!='IDW'] if dem_independent else cfg['methods']) if maps else [],
                  dem_provenance=audit['dem_provenance'],replacement_dem_status='used' if dem is not None else 'not required by this run')
    (ROOT/cfg['results_dir']/'run_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    return pred,aggregate,rank

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(); group=parser.add_mutually_exclusive_group()
    group.add_argument('--no-maps',action='store_true'); group.add_argument('--maps-without-dem',action='store_true',help='Reconstruct the four methods that do not use external elevations')
    args=parser.parse_args(); run(not args.no_maps,args.maps_without_dem)
