import pytest

from bintools.align.utils import compute_XpY
from bintools.align.utils import compute_XpY_oe

FIX = pytest.mark.xfail(reason="fix pending")


def test_compute_XpY_value():
    assert compute_XpY("ACGCG", "C", "G") == {"CpG": pytest.approx(0.5)}


def test_compute_XpY_absent_base():
    assert compute_XpY("AAAA", "C", "G") == {"CpG": 0}


@FIX
def test_compute_XpY_single_base_no_zero_division():
    compute_XpY("A", "C", "G")


@FIX
def test_compute_XpY_empty_no_zero_division():
    compute_XpY("", "C", "G")


def test_compute_XpY_oe_value():
    # acgt=6, CpG=2, C=3, G=2 -> (2/5) / (3*2/36)
    r = compute_XpY_oe("ACGCGC", "C", "G")
    assert r["CpGo/e"] == pytest.approx((2 / 5) / (3 * 2 / 36))


@FIX
def test_compute_XpY_oe_x_absent_no_zero_division():
    compute_XpY_oe("ATATAT", "C", "G")


@FIX
def test_compute_XpY_oe_y_absent_no_zero_division():
    compute_XpY_oe("ACACAC", "C", "G")
