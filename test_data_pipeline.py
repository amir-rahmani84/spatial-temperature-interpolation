"""Synthetic fixtures test raster mathematics and DEM-free evaluation, not SRTM observations."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
import xarray as xr
import rasterio
from rasterio.transform import from_origin
import analysis as a
import prepare_dem as d

class DataPipelineChecks(unittest.TestCase):
    def test_optional_backend_version_reporting(self):
        actual=a.importlib.metadata.version
        def version(name):
            if name in ('netCDF4','h5netcdf'):
                raise a.importlib.metadata.PackageNotFoundError(name)
            return actual(name)
        with patch.object(a.importlib.metadata,'version',side_effect=version):
            versions=a.dependency_versions()
        self.assertIsNone(versions['netCDF4']); self.assertIsNone(versions['h5netcdf'])
        self.assertIsInstance(versions['numpy'],str)
        json.dumps(versions)  # Missing metadata remains serializable in run manifests.
    def test_bilinear_plane_and_nodata(self):
        tr=from_origin(49,38,d.PIXEL,d.PIXEL)
        raster=np.arange(16,dtype=float).reshape(4,4)
        points=np.array([tr*(1,1),tr*(2,2)])
        np.testing.assert_allclose(d.bilinear(raster,tr,points,-32768),[2.5,7.5])
        raster[0,0]=-32768
        with self.assertRaisesRegex(ValueError,'No silent filling'): d.bilinear(raster,tr,points,-32768)
        # A zero-weight invalid neighbor cannot invalidate an exact posting.
        np.testing.assert_allclose(d.bilinear(raster,tr,np.array([tr*(1.5,1.5)]),-32768),[5])
    def test_merge_seam_preserves_lattice_and_order(self):
        with tempfile.TemporaryDirectory() as folder:
            paths=[]
            for i in range(2):
                p=Path(folder)/f'tile{i}.tif'; paths.append(p)
                tr=from_origin(49+i*3*d.PIXEL,38,d.PIXEL,d.PIXEL)
                data=np.tile(np.arange(i*3,i*3+3,dtype=np.float32),(3,1))
                with rasterio.open(p,'w',driver='GTiff',width=3,height=3,count=1,dtype='float32',crs='EPSG:4326',transform=tr,nodata=-32768) as out: out.write(data,1)
            base=from_origin(49,38,d.PIXEL,d.PIXEL)
            points=np.array([base*(4.5,1.5),base*(3,1.5),base*(1.5,1.5)])
            values,_=d.sample_rasters(paths,points)
            np.testing.assert_allclose(values,[4,2.5,1],atol=1e-7)
            with rasterio.open(paths[0],'r+') as src: src.crs='EPSG:3857'
            with self.assertRaisesRegex(ValueError,'EPSG:4326'): d.sample_rasters(paths,points)
    def test_comparison_summary_and_idw_delta(self):
        with tempfile.TemporaryDirectory() as folder:
            old=pd.DataFrame(dict(lon=[49,49.1,49.2],lat=[37,37.1,37.2],elevation=[100,200,300]))
            path=Path(folder)/'old.csv'; old.to_csv(path,index=False)
            new=old.copy(); new.elevation+=10
            summary=d.compare(path,new,Path(folder)/'private',-.0065)
            self.assertAlmostEqual(summary['elevation']['mean_signed_difference'],10)
            self.assertAlmostEqual(summary['idw_temperature_change']['rmse'],.065)
    def test_documented_dem_required_for_full_maps(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'gilan_dem_masked.csv'; p.write_text('lon,lat,elevation\n49,37,100\n')
            with self.assertRaises(FileNotFoundError): a.load_dem(p)
    def test_no_maps_cli_without_dem(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); data=root/'data'; data.mkdir()
            shutil.copyfile(a.ROOT/'analysis.py',root/'analysis.py')
            cfg=a.configuration(); (root/'config.json').write_text(json.dumps(cfg),encoding='utf-8')
            times=pd.date_range(cfg['start'],cfg['end'],freq='MS')
            vals=np.broadcast_to(np.array([[280.,281.,282.],[282.,284.,286.]]),(len(times),2,3)).copy()
            ds=xr.Dataset(dict(t2m=(('valid_time','latitude','longitude'),vals),z=(('valid_time','latitude','longitude'),np.ones_like(vals)*980.665)),coords=dict(valid_time=times,latitude=[37.,37.25],longitude=[49.,49.25,49.5]))
            ds.t2m.attrs['units']='K'; ds.z.attrs['units']='m**2 s**-2'; ds.to_netcdf(data/'data.nc')
            geo=dict(type='FeatureCollection',features=[dict(type='Feature',properties=dict(NAME_1='Gilan'),geometry=dict(type='Polygon',coordinates=[[[48,36],[51,36],[51,39],[48,39],[48,36]]]))])
            (data/cfg['boundary']).write_text(json.dumps(geo))
            with patch.dict(os.environ,{'SPATIAL_DATA_DIR':str(data)}),patch.object(a,'load_dem',side_effect=AssertionError('DEM touched')):
                _,_,dem,audit=a.load_data(cfg,require_dem=False)
                self.assertIsNone(dem); self.assertIsNone(audit['dem_provenance'])
            env=os.environ.copy(); env.update(SPATIAL_DATA_DIR=str(data),MPLBACKEND='Agg')
            result=subprocess.run([sys.executable,str(root/'analysis.py'),'--no-maps'],env=env,capture_output=True,text=True,timeout=180)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertFalse((data/cfg['dem']).exists())
            manifest=json.loads((root/'results/run_manifest.json').read_text())
            self.assertFalse(manifest['maps_executed']); self.assertIsNone(manifest['dem_provenance'])
            self.assertEqual(len(pd.read_csv(root/'results/cross_validation_predictions.csv')),150)

if __name__=='__main__': unittest.main()
