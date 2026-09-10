#!/usr/bin/env python3
"""Restore and SHA256-check retained data; never overwrite different files."""
import hashlib
import json
from pathlib import Path
import zipfile

HERE=Path(__file__).resolve().parent


def main():
    m=json.loads((HERE/'evidence_manifest.json').read_text())
    p=HERE/m['archive'];raw=p.read_bytes()
    assert hashlib.sha256(raw).hexdigest()==m['archive_sha256']
    checked=0
    with zipfile.ZipFile(p) as z:
        assert set(z.namelist())==set(m['files'])
        for name,r in m['files'].items():
            dest=(HERE/name).resolve()
            assert dest.is_relative_to(HERE.resolve())
            data=z.read(name)
            assert len(data)==r['bytes'] and hashlib.sha256(data).hexdigest()==r['sha256']
            if dest.exists():
                if dest.read_bytes()!=data:
                    raise SystemExit(f'Refusing to overwrite changed local data: {name}')
            else:
                dest.parent.mkdir(parents=True,exist_ok=True)
                dest.write_bytes(data)
            checked+=1
    print(json.dumps(dict(result='PASS',files_checked_or_restored=checked)))


if __name__=='__main__':main()
