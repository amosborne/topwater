from hashlib import sha256
import numpy as np
from .ainf_sources import canonical_matrix_hash
from .fixed_source_export import dense_plan


def test_row_stream_operator_hash_preserves_canonical_bytes_for_layouts():
    values=np.arange(120).reshape(10,12)+1j*np.arange(120).reshape(10,12)/7
    for matrix in (values,np.asfortranarray(values),values[:,::2]):
        assert canonical_matrix_hash(matrix)==sha256(np.ascontiguousarray(matrix).tobytes()).hexdigest()


def test_largest_export_plan_counts_cache_pristine_lu_and_transient():
    plan=dense_plan(7700)
    assert plan['one_complex_matrix_bytes']==948640000 and plan['dense_matrices']==6
    assert plan['estimated_peak_bytes']==8113672704 and plan['estimated_peak_bytes']<8*1024**3
