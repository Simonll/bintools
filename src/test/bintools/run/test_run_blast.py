import pytest

from bintools.run.run import generate_blast_cmd


def test_tblastx_calls_tblastx_binary():
    kwargs = dict([("-query", "q.fa")])
    cmd = generate_blast_cmd("tblastx", "-v /a:/b", **kwargs)
    assert "tblastx" in cmd.split()
    assert "blastx" not in cmd.split()
