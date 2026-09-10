#!/usr/bin/env python3
"""Package large exact ledgers/transcripts without cluttering the source tree."""
import hashlib
import json
from pathlib import Path
import zipfile

HERE=Path(__file__).resolve().parent


def main():
    files=sorted((HERE/'models').glob('*.json'))
    files += [HERE/n for n in ('ledger.json','six_gate_checkpoint.json','run.log',
              'additional_gates.log','point_hunt_results.json','point_hunt.log','seed_validation.log')]
    manifest={}
    archive=HERE/'arithmetic_evidence.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for p in files:
            data=p.read_bytes();name=str(p.relative_to(HERE))
            manifest[name]=dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
            info=zipfile.ZipInfo(name,date_time=(2026,9,10,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o100644<<16
            z.writestr(info,data)
    out=dict(archive=archive.name,archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
             files=manifest)
    (HERE/'evidence_manifest.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(files=len(files),archive_bytes=archive.stat().st_size,
                          uncompressed_bytes=sum(r['bytes'] for r in manifest.values()))))


if __name__=='__main__':main()
