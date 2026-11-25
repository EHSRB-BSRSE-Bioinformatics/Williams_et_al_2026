"""Freeze the existing locked R 4.2 environment plus three cached dependencies.

No package is updated in the parent environment. All conda artifacts are recorded
as explicit URLs with checksums; fgsea is compiled from the preserved source.
"""
import json
import os
import pathlib
import shutil
import subprocess

root = pathlib.Path.cwd()
out = root / 'docs/revision/gsea_unique_ids_v1'
lib = root / 'data/large/gsea_unique_ids_v1/library'
base = root / '.pixi/envs/deg-study'
cache = pathlib.Path(os.environ['CONDA_PACKAGE_CACHE'])
out.mkdir(parents=True, exist_ok=True)
lib.mkdir(parents=True, exist_ok=True)
extras = ['r-data.table-1.14.6-r42h06615bd_0',
          'r-fastmatch-1.1_4-r42h57805ef_0', 'r-cowplot-1.1.3-r42hc72bb7e_0']
records = [json.loads(p.read_text()) for p in sorted((base/'conda-meta').glob('*.json'))]
for package in extras:
    folder = cache/package
    records.append(json.loads((folder/'info/repodata_record.json').read_text()))
    for src in (folder/'lib/R/library').iterdir():
        dest = lib/src.name
        if not dest.exists():
            shutil.copytree(src, dest, symlinks=True)
lines = ['@EXPLICIT']
for r in records:
    url = r.get('url') or r['channel'].rstrip('/')+'/'+r['subdir']+'/'+r['fn']
    lines.append(url+'#'+r['md5'])
(out/'conda-explicit.txt').write_text('\n'.join(lines)+'\n')
shutil.copy2('pixi.lock', out/'base-pixi.lock')
env = dict(os.environ, R_LIBS_USER=str(lib), MAKEFLAGS='-j4',
           PATH=str(base/'bin')+os.pathsep+os.environ['PATH'])
subprocess.run([str(base/'bin/R'), 'CMD', 'INSTALL', '-l', str(lib),
                'data/large/reactome_ora_v1/annotation/fgsea_1.26.0.tar.gz'],
               env=env, check=True)
