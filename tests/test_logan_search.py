"""Unit tests for Logan Search (job ID, parquet / TSV export) support (no network required)."""
import pytest

from logan_blaster import extract_accessions_from_parquet, extract_accessions_from_tsv, is_kmviz_session


def test_is_kmviz_session():
    assert is_kmviz_session("kmviz-b2bce461-ca13-4a45-b0b4-6c894eacf103")
    assert not is_kmviz_session("2f780abfeae3f54e")


def test_extract_accessions_from_tsv_with_header(tmp_path):
    tsv = tmp_path / "job.tsv"
    tsv.write_text(
        "score\taccession\tshard\n"
        "1.0\tSRR12825347\tTRANSCRIPTOMIC_INV_26\n"
        "\n"
        "0.9\tERR333791\tTRANSCRIPTOMIC_INV_28\n"
    )
    out = tmp_path / "acc.txt"
    assert extract_accessions_from_tsv(str(tsv), str(out)) == 2
    assert out.read_text() == "SRR12825347\nERR333791\n"


def test_extract_accessions_from_tsv_without_accession_header(tmp_path):
    tsv = tmp_path / "job.tsv"
    tsv.write_text("run\tscore\nSRR1\t1.0\n")
    out = tmp_path / "acc.txt"
    assert extract_accessions_from_tsv(str(tsv), str(out)) == 1
    assert out.read_text() == "SRR1\n"


def test_extract_accessions_from_parquet(tmp_path):
    pa = pytest.importorskip("pyarrow")
    pq = pytest.importorskip("pyarrow.parquet")
    parquet_file = tmp_path / "job.parquet"
    table = pa.table({"accession": ["SRR12825347", None, "ERR333791"], "score": [1.0, 0.9, 0.8]})
    pq.write_table(table, str(parquet_file))
    out = tmp_path / "acc.txt"
    assert extract_accessions_from_parquet(str(parquet_file), str(out)) == 2
    assert out.read_text() == "SRR12825347\nERR333791\n"
