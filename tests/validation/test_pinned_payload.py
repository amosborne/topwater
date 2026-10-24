from io import BytesIO
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile
import json
import numpy as np
import pytest
from .common_output_sources import pinned_bytes


@pytest.mark.parametrize('nested',[False,True])
def test_pinned_reader_accepts_both_actual_archive_layouts_and_rejects_changed_arrays(tmp_path,nested):
    path=tmp_path/'sample.json';npz=path.with_suffix('.npz');np.savez(npz,basis=np.arange(12).reshape(2,6))
    report={'arrays_sha256':sha256(npz.read_bytes()).hexdigest()};payload=BytesIO()
    if nested:
        with ZipFile(payload,'w') as z:z.writestr('arrays.npz',npz.read_bytes())
        report['payload_package']={'sha256':sha256(payload.getvalue()).hexdigest()}
    raw=json.dumps(report).encode();path.write_bytes(raw);bundle=tmp_path/'input.zip'
    with ZipFile(bundle,'x') as z:
        z.writestr('validation/sample.json',raw)
        z.writestr('validation/sample.payload.zip' if nested else 'validation/sample.npz',payload.getvalue() if nested else npz.read_bytes())
    pinned_bytes(bundle,'sample',path)
    npz.write_bytes(npz.read_bytes()+b'changed')
    with pytest.raises(ValueError):pinned_bytes(bundle,'sample',path)
