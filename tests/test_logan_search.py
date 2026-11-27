"""Unit tests for brc-analytics (parquet) and logan-search.org (session.json) support (no network required)."""
import json
import zipfile

import pytest

from logan_blaster import (
    DEFAULT_LOGAN_SEARCH_API,
    extract_accessions_from_parquet,
    extract_kmviz_session_json,
    fasta_sequences,
    is_kmviz_session,
    parse_session_argument,
    read_kmviz_session,
    read_query_from_parquet,
)


def test_is_kmviz_session():
    assert is_kmviz_session("kmviz-b2bce461-ca13-4a45-b0b4-6c894eacf103")
    assert not is_kmviz_session("2f780abfeae3f54e")


def test_extract_accessions_from_parquet(tmp_path):
    pa = pytest.importorskip("pyarrow")
    pq = pytest.importorskip("pyarrow.parquet")
    parquet_file = tmp_path / "job.parquet"
    table = pa.table({"accession": ["SRR12825347", None, "ERR333791"], "score": [1.0, 0.9, 0.8]})
    pq.write_table(table, str(parquet_file))
    out = tmp_path / "acc.txt"
    assert extract_accessions_from_parquet(str(parquet_file), str(out)) == 2
    assert out.read_text() == "SRR12825347\nERR333791\n"


def test_parse_session_argument():
    assert parse_session_argument("58831987628c3b84") == ("58831987628c3b84", None)
    assert parse_session_argument("kmviz-b2bce461") == ("kmviz-b2bce461", None)
    assert parse_session_argument("https://brc-analytics.org/logan-search?job=58831987628c3b84") == (
        "58831987628c3b84", "https://brc-analytics.org/api/v1")
    assert parse_session_argument("https://dev.brc-analytics.org/logan-search?job=58831987628c3b84") == (
        "58831987628c3b84", "https://dev.brc-analytics.org/api/v1")


def test_default_api_is_prod():
    assert DEFAULT_LOGAN_SEARCH_API == "https://brc-analytics.org/api/v1"


def test_parse_session_argument_url_without_job():
    with pytest.raises(SystemExit):
        parse_session_argument("https://dev.brc-analytics.org/logan-search")


def test_fasta_sequences():
    assert fasta_sequences(">a desc\nacgt\nAC\n>b\nTT\n") == ["ACGTAC", "TT"]


def test_read_query_from_parquet(tmp_path):
    pa = pytest.importorskip("pyarrow")
    pq = pytest.importorskip("pyarrow.parquet")
    table = pa.table({"accession": ["SRR1"]})
    without_query = tmp_path / "old.parquet"
    pq.write_table(table, str(without_query))
    assert read_query_from_parquet(str(without_query)) is None

    with_query = tmp_path / "new.parquet"
    table = table.replace_schema_metadata({"logan.query_fasta": ">JK976550.1 some description\nACGT\n",
                                           "logan.kmer_size": "31"})
    pq.write_table(table, str(with_query))
    assert read_query_from_parquet(str(with_query)) == ">JK976550.1 some description\nACGT\n"


KMVIZ_SESSION = {
    "kmviz-test": {
        "JK976550.1": {
            "Logan": {
                "_metadata": {
                    "ID": {"0": "ERR1", "1": "SRR2", "10": "SRR3"},
                    "kmer_coverage": {"0": 0.5, "1": 0.9, "10": 0.7},
                },
                "_query": {"_name": "JK976550.1", "_seq": "ACGTACGT"},
            }
        }
    }
}


def test_read_kmviz_session(tmp_path):
    session_json = tmp_path / "session.json"
    session_json.write_text(json.dumps(KMVIZ_SESSION))
    query_fasta, accessions = read_kmviz_session(str(session_json))
    assert query_fasta == ">JK976550.1\nACGTACGT\n"
    assert accessions == ["SRR2", "SRR3", "ERR1"]  # decreasing k-mer coverage


def test_read_kmviz_session_several_queries(tmp_path):
    session = json.loads(json.dumps(KMVIZ_SESSION))
    session["kmviz-test"]["other"] = {"Logan": {"_metadata": {"ID": {"0": "SRR9"}},
                                                "_query": {"_name": "other", "_seq": "TTTT"}}}
    session_json = tmp_path / "session.json"
    session_json.write_text(json.dumps(session))
    with pytest.raises(SystemExit):
        read_kmviz_session(str(session_json))


def test_extract_kmviz_session_json(tmp_path):
    zip_file = tmp_path / "kmviz-test.zip"
    with zipfile.ZipFile(zip_file, "w") as archive:
        archive.writestr("session.json", json.dumps(KMVIZ_SESSION))
        archive.writestr("README.md", "readme")
    out = tmp_path / "out.json"
    extract_kmviz_session_json(str(zip_file), str(out))
    assert json.loads(out.read_text()) == KMVIZ_SESSION
