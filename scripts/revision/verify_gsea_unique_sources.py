"""Verify frozen DEG source identity, without extracting/refitting any results."""
import csv
import hashlib
from pathlib import Path

archive = Path('../2023_Meier_FDR/bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10')
rows = []
for record in csv.DictReader(Path('data/reference/deg_source_archives.tsv').open(), delimiter='\t'):
    p = archive/record['archive_file']
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with p.open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            sha.update(block)
            md5.update(block)
    assert sha.hexdigest() == record['sha256'], f'Changed frozen source: {p}'
    rows.append(dict(archive=p.name,sha256=sha.hexdigest(),md5=md5.hexdigest(),bytes=p.stat().st_size))
    print('PASS: frozen archive SHA-256', p.name, flush=True)
out = Path('docs/revision/gsea_unique_ids_v1/source_archives.tsv')
with out.open('w') as f:
    writer = csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
