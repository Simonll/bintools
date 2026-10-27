import datetime
from types import SimpleNamespace

import pytest

from bintools.utils import utils as bt_utils

BT_T0 = datetime.datetime(2026, 1, 1, 12, 0, 0)


def bt_patch(monkeypatch, delta_us):
    t1 = BT_T0 + datetime.timedelta(microseconds=delta_us)

    class BtDT(datetime.datetime):
        @classmethod
        def now(cls, tz=None):
            return t1

    class BtDate(datetime.date):
        @classmethod
        def today(cls):
            return BT_T0.date()

    monkeypatch.setattr(
        bt_utils, "datetime", SimpleNamespace(datetime=BtDT, date=BtDate)
    )


def test_delta_time_microseconds(monkeypatch):
    bt_patch(monkeypatch, 50)
    assert bt_utils.compute_delta_time(BT_T0) == "50.00us"


def test_delta_time_milliseconds(monkeypatch):
    bt_patch(monkeypatch, 5000)
    assert bt_utils.compute_delta_time(BT_T0) == "5.00ms"


def test_delta_time_seconds(monkeypatch):
    bt_patch(monkeypatch, 10_000_000)
    assert bt_utils.compute_delta_time(BT_T0) == "10.00s"


def test_delta_time_minutes(monkeypatch):
    bt_patch(monkeypatch, 240_000_000)
    assert bt_utils.compute_delta_time(BT_T0) == "4.00min"
