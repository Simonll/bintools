import numpy as np
import pandas as pd
import pytest

pytest.importorskip("sklearn")

from bintools.cabc.utils import transform  # noqa: E402


@pytest.mark.xfail(reason="fix pending")
def test_log10_inverse_is_power_of_ten():
    t = transform("log10", pd.DataFrame({"a": [1.0, 10.0, 100.0]}))
    out = t.inverse_transform(np.array([[2.0], [10.0]]))
    assert out[0, 0] == pytest.approx(100.0)
    assert out[1, 0] == pytest.approx(1e10)


@pytest.mark.xfail(reason="fix pending")
def test_log10_roundtrip():
    df = pd.DataFrame({"a": [1.0, 10.0, 100.0]})
    t = transform("log10", df)
    back = t.inverse_transform(t.transform(df))
    assert np.allclose(np.asarray(back).ravel(), [1.0, 10.0, 100.0])
