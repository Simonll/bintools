import io

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from bintools.align.align import ali
from bintools.align.align import find_missing_col
from bintools.align.align import read_fasta
from bintools.align.align import read_phylip
from bintools.align.align import write_fasta_from_align

FIX = pytest.mark.xfail(reason="fix pending")


def bt_make_ali(seqs):
    d = {k: {i: b for i, b in enumerate(s)} for k, s in seqs.items()}
    n_site = len(next(iter(seqs.values())))
    return ali(n_taxa=len(seqs), n_site=n_site, dict_of_seq=d, genetic_code="Standard")


def bt_seq_str(seq):
    return "".join(seq.values()) if isinstance(seq, dict) else "".join(seq)


BT_SEQS = {"a": "ACGT", "b": "TGCA", "c": "AAAA"}


@FIX
def test_get_site_i():
    a = bt_make_ali(BT_SEQS)
    assert a.get_site_i(1) == {"a": "C", "b": "G", "c": "A"}
    assert a.get_site_i(0) == {"a": "A", "b": "T", "c": "A"}


def test_get_site_i_negative():
    with pytest.raises(AssertionError):
        bt_make_ali(BT_SEQS).get_site_i(-1)


@FIX
def test_get_site_i_too_large():
    with pytest.raises(AssertionError):
        bt_make_ali(BT_SEQS).get_site_i(4)


@FIX
def test_iter_yields_taxon_seq_pairs():
    a = bt_make_ali(BT_SEQS)
    got = {sp: bt_seq_str(seq) for sp, seq in a}
    assert got == BT_SEQS


@FIX
def test_find_missing_col_with_gaps():
    a = bt_make_ali({"s1": "ATGAAA---CCC", "s2": "ATGAAAGGGCCC"})
    kept = find_missing_col(a, ["---"], ratio=1)
    assert kept == [0, 1, 3]
    assert all(0 <= k < 4 for k in kept)


@FIX
def test_find_missing_col_no_gaps_keeps_all_codons():
    a = bt_make_ali({"s1": "ATGAAACCC", "s2": "ATGAAAGGG"})
    assert find_missing_col(a, ["---"], ratio=1) == [0, 1, 2]


def test_read_fasta_simple():
    a = read_fasta(io.StringIO(">x\nACGT\n>y\nTGCA\n"))
    assert a.get_n_taxa() == 2
    assert a.get_n_site() == 4
    assert bt_seq_str(a.get_dict_of_seq()["y"]) == "TGCA"


@FIX
def test_read_fasta_wrapped():
    a = read_fasta(io.StringIO(">x\nACGT\nAAAA\n>y\nTGCA\nCCCC\n"))
    d = a.get_dict_of_seq()
    assert bt_seq_str(d["x"]) == "ACGTAAAA"
    assert bt_seq_str(d["y"]) == "TGCACCCC"
    assert a.get_n_site() == 8


def bt_msa():
    return MultipleSeqAlignment(
        [SeqRecord(Seq("ACGT"), id="a"), SeqRecord(Seq("TGCA"), id="b")]
    )


@FIX
def test_write_fasta_twice_no_duplicates(tmp_path):
    f = str(tmp_path / "out.fa")
    assert write_fasta_from_align(f, bt_msa())
    assert write_fasta_from_align(f, bt_msa())
    with open(f) as fh:
        assert fh.read() == ">a\nACGT\n>b\nTGCA\n"


def test_write_fasta_once(tmp_path):
    f = str(tmp_path / "out.fa")
    assert write_fasta_from_align(f, bt_msa())
    with open(f) as fh:
        assert fh.read() == ">a\nACGT\n>b\nTGCA\n"


def test_read_phylip_simple():
    a = read_phylip(io.StringIO("2\t4\na  ACGT\nb  TGCA\n"))
    assert a.get_n_taxa() == 2
    assert a.get_n_site() == 4
    assert bt_seq_str(a.get_dict_of_seq()["a"]) == "ACGT"


@FIX
def test_read_phylip_empty_raises_clean_error():
    with pytest.raises((ValueError, RuntimeError, EOFError)):
        read_phylip(io.StringIO(""))
