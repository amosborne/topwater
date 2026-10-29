from pathlib import Path
from types import SimpleNamespace
import capytaine as cpt
import numpy as np
from capytaine.green_functions.delhommeau import Delhommeau
from topwater.solve.source_cache import SourceCache


def hydro():
    mesh=cpt.Mesh(vertices=np.array([[-1.,-1.,-.4],[1.,-1.,-.4],[1.,1.,-.4],[-1.,1.,-.4]]),
                  faces=np.array([[0,1,2,3]]),auto_clean=False,auto_check=False)
    body=cpt.FloatingBody(mesh=mesh,dofs=cpt.rigid_body_dofs(rotation_center=(0.,0.,0.)))
    solver=cpt.BEMSolver(green_function=Delhommeau(tabulation_cache_dir=str(Path.home()/'.cache/topwater/capytaine-3.0.0')))
    return SimpleNamespace(body=body,solver=solver,gravity=9.81,density=1025.,depth=np.inf)


def test_identity_uses_numeric_dofs_and_invalidates_physics(tmp_path):
    cache=SourceCache(tmp_path)
    a,b=hydro(),hydro()
    w=np.array([1.,2.])
    assert cache.path(a,w)==cache.path(b,w)
    b.gravity=9.8
    assert cache.path(a,w)!=cache.path(b,w)
    assert cache.path(a,w)!=cache.path(a,w+.1)


def test_atomic_numeric_cache_and_corruption_rebuild(tmp_path):
    cache=SourceCache(tmp_path)
    path=tmp_path/'case.npz'
    values=np.arange(12).reshape(2,1,6).astype(complex)
    assert cache.load(path,values.shape) is None
    assert cache.save(path,values)
    actual=cache.load(path,values.shape)
    np.testing.assert_array_equal(actual,values)
    assert not actual.flags.writeable
    path.write_bytes(b'broken')
    assert cache.load(path,values.shape) is None
    assert not list(tmp_path.glob('.sources-*'))


def test_cache_write_failure_is_nonfatal_and_removes_staging(tmp_path, monkeypatch):
    def denied(*args):
        raise PermissionError('cache destination is read-only')
    monkeypatch.setattr('topwater.solve.source_cache.os.replace', denied)
    cache = SourceCache(tmp_path)
    assert not cache.save(tmp_path/'case.npz', np.zeros((1, 1, 6), dtype=complex))
    assert not list(tmp_path.iterdir())


def test_cache_staging_cleanup_failure_is_nonfatal(tmp_path, monkeypatch, capsys):
    def denied(*args, **kwargs):
        raise PermissionError('read-only cache')
    with monkeypatch.context() as patch:
        patch.setattr('topwater.solve.source_cache.os.replace', denied)
        patch.setattr(Path, 'unlink', denied)
        assert not SourceCache(tmp_path).save(tmp_path/'case.npz', np.zeros((1, 1, 6), complex))
    assert 'staging cleanup skipped' in capsys.readouterr().err
    for path in tmp_path.iterdir():
        path.unlink()
