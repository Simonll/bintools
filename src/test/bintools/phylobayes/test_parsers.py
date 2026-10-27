import pytest

from bintools.phylobayes.gtrgaa_parser import posterior_GTRGAA
from bintools.phylobayes.mcmc_parser import chunck_chain
from bintools.phylobayes.mcmc_parser import posterior_M0_GTR

FIX = pytest.mark.xfail(reason="fix pending")


@FIX
def test_chunck_chain_exact_multiple():
    lines = list(range(30))
    chunks = list(chunck_chain(lines, 10))
    assert len(chunks) == 3
    assert chunks[0] == list(range(10))
    assert chunks[2] == list(range(20, 30))


def test_chunck_chain_leading_chunks():
    chunks = list(chunck_chain(list(range(35)), 10))
    assert len(chunks) >= 3
    assert chunks[0] == list(range(10))
    assert chunks[1] == list(range(10, 20))
    assert chunks[2] == list(range(20, 30))


BT_M0 = (
    "(a,b);\n0.1\t0.2\t0.3\t0.4\n1\t2\t3\t4\t5\t6\n0.5\n"
    "(a,c);\n0.2\t0.2\t0.2\t0.4\n6\t5\t4\t3\t2\t1\n0.7\n"
)


def test_m0_gtr_read_valid(tmp_path):
    f = tmp_path / "chain.trace"
    f.write_text(BT_M0)
    p = posterior_M0_GTR.read(f)
    assert p is not None
    assert p.list_of_trees == ["(a,b);", "(a,c);"]


@FIX
def test_m0_gtr_read_trailing_blank_line(tmp_path):
    f = tmp_path / "chain.trace"
    f.write_text(BT_M0 + "\n")
    p = posterior_M0_GTR.read(f)
    assert p is not None
    assert p.list_of_trees[:2] == ["(a,b);", "(a,c);"]
    assert p.list_of_phi[0] == [0.1, 0.2, 0.3, 0.4]


@FIX
def test_gtrgaa_sample_uses_rho():
    probe = posterior_GTRGAA([], [], [], [], [])
    phi = [float(i) for i in range(len(probe.list_of_phiID))]
    rho = [100.0 + i for i in range(len(probe.list_of_rhoID))]
    p = posterior_GTRGAA(
        list_of_sampleID=[0],
        list_of_trees=["(a,b);"],
        list_of_phi=[phi],
        list_of_rho=[rho],
        list_of_alpha=[0.5],
    )
    d = p.sample()
    for i, r in enumerate(probe.list_of_rhoID):
        assert d["rho_" + r] == rho[i]
    assert d["phi_" + probe.list_of_phiID[3]] == phi[3]
