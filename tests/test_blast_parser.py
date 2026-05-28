"""Unit tests for blast parser functions (no external tools or network required)."""
import io
import contextlib
import pytest

from logan_blaster import (
    get_query_ACGT,
    get_query_name,
    get_query_length,
    parse_blastn,
    run_blast_parser,
)

QUERY_LENGTH = 963
QUERY_NAME = "my_query"
QUERY_PREFIX = "ATGATATTTTCAACTTTAGA"
QUERY_SUFFIX = "TCCTAATTGA"


class TestGetQueryACGT:
    def test_length(self, query_fa):
        assert len(get_query_ACGT(query_fa)) == QUERY_LENGTH

    def test_prefix(self, query_fa):
        assert get_query_ACGT(query_fa).startswith(QUERY_PREFIX)

    def test_suffix(self, query_fa):
        assert get_query_ACGT(query_fa).endswith(QUERY_SUFFIX)

    def test_only_nucleotides(self, query_fa):
        seq = get_query_ACGT(query_fa)
        assert set(seq.upper()) <= set("ACGT"), "Sequence should contain only ACGT"

    def test_returns_first_sequence_only(self, tmp_path):
        fa = tmp_path / "multi.fa"
        fa.write_text(">seq1\nAAAA\n>seq2\nCCCC\n")
        assert get_query_ACGT(str(fa)) == "AAAA"

    def test_handles_multiline_fasta(self, tmp_path):
        fa = tmp_path / "multiline.fa"
        fa.write_text(">seq1\nATGC\nATGC\n")
        assert get_query_ACGT(str(fa)) == "ATGCATGC"


class TestGetQueryName:
    def test_reads_name_from_blast_output(self, self_blast_txt):
        assert get_query_name(self_blast_txt) == QUERY_NAME

    def test_returns_none_when_missing(self, tmp_path):
        f = tmp_path / "empty.txt"
        f.write_text("no query here\n")
        assert get_query_name(str(f)) is None


class TestGetQueryLength:
    def test_reads_length_from_blast_output(self, self_blast_txt):
        assert get_query_length(self_blast_txt) == QUERY_LENGTH

    def test_returns_none_when_missing(self, tmp_path):
        f = tmp_path / "empty.txt"
        f.write_text("no length here\n")
        assert get_query_length(str(f)) is None


class TestParseBLASTN:
    def test_returns_correct_name(self, self_blast_txt):
        name, _, _ = parse_blastn(self_blast_txt)
        assert name == QUERY_NAME

    def test_returns_correct_length(self, self_blast_txt):
        _, length, _ = parse_blastn(self_blast_txt)
        assert length == QUERY_LENGTH

    def test_positions_vector_size(self, self_blast_txt):
        _, length, positions = parse_blastn(self_blast_txt)
        assert len(positions) == length

    def test_all_positions_covered_by_perfect_match(self, self_blast_txt):
        _, _, positions = parse_blastn(self_blast_txt)
        assert all(p >= 1 for p in positions), "Every position should be covered at least once"

    def test_secondary_alignments_increase_count(self, self_blast_txt):
        # The Plus/Minus alignment covers query 54-75 → those positions get count >= 2
        _, _, positions = parse_blastn(self_blast_txt)
        for i in range(54, 76):
            assert positions[i - 1] >= 2, f"Position {i} should have count >= 2"

    def test_spot_check_triple_overlap(self, self_blast_txt):
        # Positions 54-61: covered by perfect match + secondary Plus/Minus + secondary Plus/Plus
        _, _, positions = parse_blastn(self_blast_txt)
        for i in range(54, 62):
            assert positions[i - 1] == 3, f"Position {i} should have count == 3"

    def test_spot_check_single_coverage(self, self_blast_txt):
        # Position 1: only the perfect match, count must be exactly 1
        _, _, positions = parse_blastn(self_blast_txt)
        assert positions[0] == 1

    def test_max_count(self, self_blast_txt):
        _, _, positions = parse_blastn(self_blast_txt)
        assert max(positions) == 3

    def test_graceful_on_empty_blast(self, tmp_path):
        f = tmp_path / "empty.txt"
        f.write_text("BLASTN 2.x\n\nQuery= test\n\nLength=5\n\n")
        name, length, positions = parse_blastn(str(f))
        assert name == "test"
        assert length == 5
        assert positions == [0, 0, 0, 0, 0]

    def test_minimal_blast_positions(self, tmp_path):
        blast = tmp_path / "minimal.txt"
        # Length=15 so positions 3-12 are within bounds
        blast.write_text(
            "Query= minimal\n\nLength=15\n\nQuery  3  ATGATATTTT  12\n"
        )
        _, length, positions = parse_blastn(str(blast))
        assert length == 15
        assert positions[0:2] == [0, 0]           # positions 1-2: not covered
        assert positions[2:12] == [1] * 10        # positions 3-12: covered once
        assert positions[12:15] == [0, 0, 0]      # positions 13-15: not covered

    def test_returns_empty_positions_on_missing_length(self, tmp_path):
        blast = tmp_path / "nolength.txt"
        blast.write_text("Query= test\n\nQuery  1  ATGAT  5\n")
        name, length, positions = parse_blastn(str(blast))
        assert length is None
        assert positions == []


class TestRunBlastParser:
    def test_output_matches_expected_reference(self, query_fa, self_blast_txt, expected_self_synth):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            run_blast_parser(query_fa, self_blast_txt, abundance=True)
        assert buf.getvalue() == expected_self_synth

    def test_raises_on_sequence_length_mismatch(self, tmp_path, self_blast_txt):
        wrong_fa = tmp_path / "short.fa"
        wrong_fa.write_text(">wrong\nATGAT\n")
        with pytest.raises(AssertionError, match="length"):
            run_blast_parser(str(wrong_fa), self_blast_txt, abundance=True)

    def test_output_contains_query_name(self, query_fa, self_blast_txt):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            run_blast_parser(query_fa, self_blast_txt, abundance=True)
        assert QUERY_NAME in buf.getvalue()

    def test_output_contains_sequence(self, query_fa, self_blast_txt):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            run_blast_parser(query_fa, self_blast_txt, abundance=True)
        assert QUERY_PREFIX in buf.getvalue()
