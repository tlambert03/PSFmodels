import platform

import pytest


@pytest.mark.skipif(platform.system() == "Linux", reason="annoying on ubuntu CI")
def test_napari_widget(monkeypatch):
    from magicgui import magicgui
    from psfmodels._napari import make_psf

    with monkeypatch.context() as m:
        # avoid no module named 'napari' error
        m.setitem(make_psf.__annotations__, "return", None)
        wdg = magicgui(make_psf)
    assert wdg().shape == (101, 101, 101)
