"""Download the six SRTM GL1 tiles from OpenTopography's public mirror.

No account, API key or user profile is required. This is the EGM96 GL1 data,
not the GL1 ellipsoidal variant. Source acknowledgement is required.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import urllib.request

TILES=['N36E049','N36E050','N37E048','N37E049','N37E050','N38E048']
BASE='https://opentopography.s3.sdsc.edu/raster/SRTM_GL1/SRTM_GL1_srtm/'
PLATFORM='https://portal.opentopography.org/raster?opentopoID=OTSRTM.082015.4326.1'

def download(folder):
    folder=Path(folder); folder.mkdir(parents=True,exist_ok=True)
    records=[]
    for tile in TILES:
        path=folder/(tile+'.tif'); url=BASE+path.name
        # Refuse to invent acquisition dates for files from another download.
        if path.exists(): raise FileExistsError(f'{path}: preserve its existing source record or download into a new empty directory')
        temporary=path.with_suffix('.tif.part'); digest=hashlib.sha256()
        with urllib.request.urlopen(url,timeout=180) as response,temporary.open('wb') as output:
            expected=response.headers.get('Content-Length'); size=0
            while block:=response.read(1024*1024):
                output.write(block); digest.update(block); size+=len(block)
        if expected is not None and size!=int(expected): raise IOError(f'Incomplete download: {tile}')
        temporary.replace(path)
        records.append(dict(tile_id=tile,filename=path.name,url=url,sha256=digest.hexdigest(),verified_at_utc=datetime.now(timezone.utc).isoformat()))
        print(f'{tile}: {size} bytes',flush=True)
    record=dict(product='NASA SRTM GL1 v3, distributed by OpenTopography',product_doi='10.5069/G9445JDF',platform=PLATFORM,acquired_on=datetime.now(timezone.utc).date().isoformat(),tiles=records)
    (folder/'source.json').write_text(json.dumps(record,indent=2),encoding='utf-8')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tile-dir',type=Path,default=Path(__file__).resolve().parent/'data'/'srtm')
    args=parser.parse_args(); download(args.tile_dir)
