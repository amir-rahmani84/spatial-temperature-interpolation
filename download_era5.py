"""Explicit replacement CDS request, not a recovered acquisition receipt.

Install cdsapi>=0.7.7 and configure your account outside this repository.
Run with --request-only to inspect without authentication or downloading.
"""
from pathlib import Path
import argparse
import json

DATASET='reanalysis-era5-single-levels-monthly-means'
REQUEST={
    'product_type':['monthly_averaged_reanalysis'],
    'variable':['2m_temperature','geopotential'],
    'year':[str(y) for y in range(2015,2025)],
    'month':[f'{m:02d}' for m in range(1,13)],
    'time':['00:00'],
    'area':[38.5,48.5,36.5,50.5],
    'data_format':'netcdf',
    'download_format':'unarchived'
}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--request-only',action='store_true')
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'data'/'data.nc')
    args=parser.parse_args()
    if args.request_only:
        print(json.dumps(dict(dataset=DATASET,request=REQUEST,target=args.output.name),indent=2)); return
    import cdsapi
    if args.output.exists(): raise FileExistsError('Choose another output to preserve the verified existing input')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    cdsapi.Client().retrieve(DATASET,REQUEST,str(args.output))

if __name__=='__main__': main()
