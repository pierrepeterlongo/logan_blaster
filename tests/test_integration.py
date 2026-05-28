"""Integration tests for the LoganBlaster pipeline.

- TestRunBlast: tests _run_blast() with local files only (blastn required)
- TestFullPipelineLocal: tests _process_accessions() with a pre-placed .fa.zst
  (blastn + back_to_sequences required, no network)
- TestFullPipelineNetwork: full end-to-end pipeline with download from Logan
  (all tools + network required — only runs with pytest --network)
"""
import io
import os
import shutil
import subprocess
import contextlib
import pytest

from logan_blaster import LoganBlaster, run_blast_parser

REQUIRED_TOOLS = ["blastn", "back_to_sequences", "zstd"]

requires_tools = pytest.mark.skipif(
    not all(shutil.which(t) for t in REQUIRED_TOOLS),
    reason=f"requires {REQUIRED_TOOLS}",
)


def _make_blaster(tmp_path, query_fa, accession_file=None, unitigs=False):
    b = LoganBlaster(
        session_id=None,
        accession_file=str(accession_file) if accession_file else None,
        query_file=str(query_fa),
        delete=False,
        unitigs=unitigs,
        kmer_size=17,
        limit=0,
        output_dir=str(tmp_path),
    )
    b.failed_accession_list = str(tmp_path / "failed_accessions.txt")
    return b


@requires_tools
class TestRunBlast:
    """Tests for _run_blast using query vs itself — no download needed."""

    def test_creates_blast_and_synth_files(self, tmp_path, query_fa):
        aln_dir = tmp_path / LoganBlaster.ALIGNEMENT_DIR_NAME
        aln_dir.mkdir()
        blaster = _make_blaster(tmp_path, query_fa)

        orig = os.getcwd()
        os.chdir(str(tmp_path))
        try:
            blaster._run_blast(str(query_fa), str(query_fa))
        finally:
            os.chdir(orig)

        assert list(aln_dir.glob("*.txt")), "blastn output file should be created"
        assert list(aln_dir.glob("synth_*.txt")), "synth file should be created"

    def test_synth_content_matches_reference(self, tmp_path, query_fa, expected_self_synth):
        aln_dir = tmp_path / LoganBlaster.ALIGNEMENT_DIR_NAME
        aln_dir.mkdir()
        blaster = _make_blaster(tmp_path, query_fa)

        orig = os.getcwd()
        os.chdir(str(tmp_path))
        try:
            blaster._run_blast(str(query_fa), str(query_fa))
        finally:
            os.chdir(orig)

        synth_files = list(aln_dir.glob("synth_*.txt"))
        assert synth_files
        assert synth_files[0].read_text() == expected_self_synth

    def test_synth_contains_query_name(self, tmp_path, query_fa):
        aln_dir = tmp_path / LoganBlaster.ALIGNEMENT_DIR_NAME
        aln_dir.mkdir()
        blaster = _make_blaster(tmp_path, query_fa)

        orig = os.getcwd()
        os.chdir(str(tmp_path))
        try:
            blaster._run_blast(str(query_fa), str(query_fa))
        finally:
            os.chdir(orig)

        synth_files = list(aln_dir.glob("synth_*.txt"))
        assert "my_query" in synth_files[0].read_text()


@requires_tools
class TestFullPipelineLocal:
    """Full pipeline test using a pre-placed .fa.zst — no download needed."""

    ACCESSION = "LOCAL_TEST_ACCESSION"

    def _setup_pipeline_dir(self, tmp_path, query_fa):
        """Create directory structure and pre-place the .zst file."""
        input_dir = tmp_path / LoganBlaster.INPUT_DATA_DIR_NAME
        logan_dir = tmp_path / LoganBlaster.LOGAN_DIR_NAME
        aln_dir = tmp_path / LoganBlaster.ALIGNEMENT_DIR_NAME
        for d in (input_dir, logan_dir, aln_dir):
            d.mkdir()

        # Copy query into input_data
        query_dest = input_dir / os.path.basename(query_fa)
        shutil.copy(query_fa, str(query_dest))

        # Minimal accession file with a single fake accession
        acc_file = input_dir / "accessions.txt"
        acc_file.write_text(f"{self.ACCESSION}\n")

        # Pre-place the "downloaded" contigs as a zstd-compressed copy of the query
        zst_path = logan_dir / f"{self.ACCESSION}.contigs.fa.zst"
        subprocess.run(["zstd", "-q", query_fa, "-o", str(zst_path)], check=True)

        return query_dest, acc_file, aln_dir

    def test_pipeline_creates_synth_file(self, tmp_path, query_fa):
        query_dest, acc_file, aln_dir = self._setup_pipeline_dir(tmp_path, query_fa)
        blaster = _make_blaster(tmp_path, str(query_dest), accession_file=acc_file)
        blaster.accession_file = str(acc_file)

        orig = os.getcwd()
        os.chdir(str(tmp_path))
        try:
            blaster._process_accessions()
        finally:
            os.chdir(orig)

        assert list(aln_dir.glob("synth_*.txt")), "synth file should be produced"

    def test_pipeline_synth_contains_query(self, tmp_path, query_fa):
        query_dest, acc_file, aln_dir = self._setup_pipeline_dir(tmp_path, query_fa)
        blaster = _make_blaster(tmp_path, str(query_dest), accession_file=acc_file)
        blaster.accession_file = str(acc_file)

        orig = os.getcwd()
        os.chdir(str(tmp_path))
        try:
            blaster._process_accessions()
        finally:
            os.chdir(orig)

        synth_files = list(aln_dir.glob("synth_*.txt"))
        assert synth_files
        content = synth_files[0].read_text()
        assert "my_query" in content
        assert len(content) > 100

    def test_pipeline_synth_matches_reference(self, tmp_path, query_fa, expected_self_synth):
        """Since query == contigs, the recruited file is identical to the query,
        so blastn produces the same alignment as query vs itself."""
        query_dest, acc_file, aln_dir = self._setup_pipeline_dir(tmp_path, query_fa)
        blaster = _make_blaster(tmp_path, str(query_dest), accession_file=acc_file)
        blaster.accession_file = str(acc_file)

        orig = os.getcwd()
        os.chdir(str(tmp_path))
        try:
            blaster._process_accessions()
        finally:
            os.chdir(orig)

        synth_files = list(aln_dir.glob("synth_*.txt"))
        assert synth_files
        assert synth_files[0].read_text() == expected_self_synth

    def test_no_failed_accessions(self, tmp_path, query_fa):
        query_dest, acc_file, aln_dir = self._setup_pipeline_dir(tmp_path, query_fa)
        blaster = _make_blaster(tmp_path, str(query_dest), accession_file=acc_file)
        blaster.accession_file = str(acc_file)

        orig = os.getcwd()
        os.chdir(str(tmp_path))
        try:
            blaster._process_accessions()
        finally:
            os.chdir(orig)

        failed_file = tmp_path / "failed_accessions.txt"
        assert not failed_file.exists() or failed_file.read_text().strip() == ""


@pytest.mark.network
class TestFullPipelineNetwork:
    """End-to-end pipeline tests downloading from Logan.
    Only run with: pytest --network
    """

    def test_first_example_accession_produces_output(self, tmp_path, query_fa, accessions_txt):
        """Run the pipeline on the first accession of example/accessions.txt."""
        with open(accessions_txt) as f:
            first_accession = f.readline().strip().split()[0]

        input_dir = tmp_path / LoganBlaster.INPUT_DATA_DIR_NAME
        logan_dir = tmp_path / LoganBlaster.LOGAN_DIR_NAME
        aln_dir = tmp_path / LoganBlaster.ALIGNEMENT_DIR_NAME
        for d in (input_dir, logan_dir, aln_dir):
            d.mkdir()

        query_dest = input_dir / os.path.basename(query_fa)
        shutil.copy(query_fa, str(query_dest))

        acc_file = input_dir / "accessions.txt"
        acc_file.write_text(f"{first_accession}\n")

        blaster = _make_blaster(tmp_path, str(query_dest), accession_file=acc_file)
        blaster.accession_file = str(acc_file)

        orig = os.getcwd()
        os.chdir(str(tmp_path))
        try:
            blaster._process_accessions()
        finally:
            os.chdir(orig)

        # Structural checks (content depends on Logan data)
        blast_files = list(aln_dir.glob("*.txt"))
        assert blast_files, f"No output for accession {first_accession}"

    def test_first_example_accession_synth_format(self, tmp_path, query_fa, accessions_txt):
        """Check synth output format on the first example accession."""
        with open(accessions_txt) as f:
            first_accession = f.readline().strip().split()[0]

        input_dir = tmp_path / LoganBlaster.INPUT_DATA_DIR_NAME
        logan_dir = tmp_path / LoganBlaster.LOGAN_DIR_NAME
        aln_dir = tmp_path / LoganBlaster.ALIGNEMENT_DIR_NAME
        for d in (input_dir, logan_dir, aln_dir):
            d.mkdir()

        query_dest = input_dir / os.path.basename(query_fa)
        shutil.copy(query_fa, str(query_dest))
        acc_file = input_dir / "accessions.txt"
        acc_file.write_text(f"{first_accession}\n")

        blaster = _make_blaster(tmp_path, str(query_dest), accession_file=acc_file)
        blaster.accession_file = str(acc_file)

        orig = os.getcwd()
        os.chdir(str(tmp_path))
        try:
            blaster._process_accessions()
        finally:
            os.chdir(orig)

        synth_files = list(aln_dir.glob("synth_*.txt"))
        if not synth_files:
            pytest.skip(f"No matches found for accession {first_accession} — cannot check synth format")

        content = synth_files[0].read_text()
        assert content.startswith("Query: my_query"), "synth should start with query name"
        assert "query" in content, "synth should contain 'query' position lines"
