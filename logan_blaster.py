#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import shutil
from urllib.request import urlopen
from urllib.parse import urlparse, parse_qs
from pathlib import Path
import ssl
import json
import zipfile
from concurrent.futures import ThreadPoolExecutor

__author__ = 'Pierre Peterlongo'

def _resolve_version():
    # Primary: read live from git tags — no extra dependency
    try:
        result = subprocess.run(
            ["git", "describe", "--tags", "--abbrev=0"],
            capture_output=True, text=True, check=True,
            cwd=os.path.dirname(os.path.abspath(__file__)) or ".",
        )
        return result.stdout.strip().lstrip("v")
    except Exception:
        pass
    # Fallback: installed package metadata (pip install without a git repo)
    try:
        from importlib.metadata import version as _pkg_version
        return _pkg_version("logan_blaster")
    except Exception:
        pass
    return "0+unknown"

__version__ = _resolve_version()

# Message colors
GREEN = "\033[0;32m"
RED = "\033[0;31m"
YELLOW = "\033[0;33m"
BLUE = "\033[1;34m"
CYAN = "\033[1;36m"
NOCOLOR = "\033[0m"


# --- Blast parser utilities ---
def get_query_ACGT(file_path):
    """Returns the first sequence from a fasta file. Possibly multiline"""
    seq = ""
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(">"):
                if seq != "":
                    return seq
                continue
            seq += line.strip()
    return seq


def get_query_name(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('Query='):
                return line.split('=')[1].strip()
    return None


def get_query_length(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('Length='):
                return int(line.split('=')[1].strip())
    return None


def parse_blastn(file_path):
    query_name = get_query_name(file_path)
    query_length = get_query_length(file_path)
    if query_length is None:
        return None, None, []
    query_positions = [0] * query_length
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('Query '):
                parts = line.split()
                try:
                    start = int(parts[1])
                    end = int(parts[3])
                    for pos in range(start, end + 1):
                        query_positions[pos - 1] += 1
                except (IndexError, ValueError):
                    continue
    return query_name, query_length, query_positions


def print_value(v, matched_positions):
    if matched_positions:
        if v == 0:
            print('-', end='')
        else:
            if v <= 26:
                print(chr(ord('a') + v - 1), end='')
            if v > 26:
                print('Z', end='')
    else:
        if v == 0:
            print('-', end='')
        else:
            print('|', end='')


def print_spaces(n):
    for _ in range(n):
        print(' ', end='')


def visualize_matches(query_ACGT, query_name, query_length, matched_positions, print_abundance=False):
    nb_chars_before_line = 9
    nb_chars_for_len = len(str(query_length))
    nb_chars_before_line += nb_chars_for_len
    len_line = 80
    print(f"Query: {query_name}")
    pos = 0
    while True:
        diff = len(str(query_length)) - len(str(pos + 1))
        print("query  ", end='')
        print_spaces(diff)
        print(f"{pos + 1}  ", end='')
        if pos + len_line < query_length:
            print(query_ACGT[pos: pos + len_line])
            print_spaces(nb_chars_before_line)
            for v in matched_positions[pos: pos + len_line]:
                print_value(v, print_abundance)
            print("\n")
            pos += len_line
        else:
            print(query_ACGT[pos:])
            print_spaces(nb_chars_before_line)
            for v in matched_positions[pos:]:
                print_value(v, print_abundance)
            print("\n")
            break


def run_blast_parser(fasta_file, blastn_file, abundance=False):
    query_ACGT = get_query_ACGT(fasta_file)
    query_name, query_length, matched_positions = parse_blastn(blastn_file)
    assert len(query_ACGT) == query_length, (
        f"Error, query given in {fasta_file} is of length {len(query_ACGT)}, "
        f"while the blastn result indicates a sequence of length {query_length}"
    )
    visualize_matches(query_ACGT, query_name, query_length, matched_positions, print_abundance=abundance)


def download_file(url, destination, hint=None):
    try:
        context = ssl._create_unverified_context()
        with urlopen(url, context=context) as response, open(destination, "wb") as f:
            while chunk := response.read(8192):
                f.write(chunk)
    except Exception as e:
        print(f"{RED}Error: Failed to download {url}.{NOCOLOR}")
        print(e)
        if hint:
            print(f"{YELLOW}[HINT] {hint}{NOCOLOR}")
        sys.exit(1)


# Logan Search (kmindex jobs served through Galaxy by BRC Analytics)
DEFAULT_LOGAN_SEARCH_API = "https://brc-analytics.org/api/v1"
KMVIZ_SESSION_PREFIX = "kmviz-"
# Key of the parquet schema metadata storing the query (fasta) used for the job
PARQUET_QUERY_FASTA_KEY = b"logan.query_fasta"


def parse_session_argument(session):
    """Split a --session value into (session_id, api_url).

    session may be a bare job / kmviz session ID, or a Logan Search URL such as
    https://brc-analytics.org/logan-search?job=<job_id>. For a URL, the
    API base URL is derived from its host (<scheme>://<host>/api/v1); otherwise api_url is None.
    """
    session = session.strip()
    if not session.lower().startswith(("http://", "https://")):
        return session, None
    parsed = urlparse(session)
    job_ids = parse_qs(parsed.query).get("job")
    if not job_ids or not job_ids[0].strip():
        print(f"{RED}Error: no 'job' parameter found in URL {session}.{NOCOLOR}")
        sys.exit(1)
    return job_ids[0].strip(), f"{parsed.scheme}://{parsed.netloc}/api/v1"


def is_kmviz_session(session_id):
    """True for logan-search.org sessions (kmviz-...), False for brc-analytics Logan Search job IDs."""
    return session_id.startswith(KMVIZ_SESSION_PREFIX)


def read_kmviz_session(session_json):
    """Parse a logan-search.org (kmviz) session.json file.

    Returns (query_fasta, accessions): the query as a fasta string, and the accessions
    sorted by decreasing k-mer coverage (as in the logan-search.org result table).
    A session containing several distinct queries is not supported.
    """
    with open(session_json, "r") as f:
        session = json.load(f)
    queries = {}
    best_coverage = {}
    for per_query in session.values():
        for per_index in per_query.values():
            for block in per_index.values():
                if not isinstance(block, dict):
                    continue
                query = block.get("_query") or {}
                if query.get("_name") and query.get("_seq"):
                    queries[query["_name"]] = query["_seq"].strip()
                metadata = block.get("_metadata") or {}
                ids = metadata.get("ID") or {}
                coverages = metadata.get("kmer_coverage") or {}
                if isinstance(ids, list):
                    ids = {str(i): accession for i, accession in enumerate(ids)}
                    if isinstance(coverages, list):
                        coverages = {str(i): coverage for i, coverage in enumerate(coverages)}
                for row, accession in ids.items():
                    if not accession or not str(accession).strip():
                        continue
                    accession = str(accession).strip()
                    coverage = coverages.get(row)
                    coverage = float(coverage) if coverage is not None else float("-inf")
                    if accession not in best_coverage or coverage > best_coverage[accession]:
                        best_coverage[accession] = coverage
    if len(queries) > 1:
        print(f"{RED}Error: {session_json} contains several queries ({', '.join(queries)}). "
              f"This is not supported: use --accessions (-a) and --query (-q) instead.{NOCOLOR}")
        sys.exit(1)
    query_fasta = None
    if queries:
        name, sequence = next(iter(queries.items()))
        query_fasta = f">{name}\n{sequence}\n"
    accessions = sorted(best_coverage, key=lambda accession: -best_coverage[accession])
    return query_fasta, accessions


def extract_kmviz_session_json(zip_file, output_file):
    """Extract the session.json file of a logan-search.org session zip archive into output_file."""
    with zipfile.ZipFile(zip_file) as archive:
        if "session.json" not in archive.namelist():
            print(f"{RED}Error: {zip_file} does not contain a session.json file.{NOCOLOR}")
            sys.exit(1)
        with archive.open("session.json") as fin, open(output_file, "wb") as fout:
            shutil.copyfileobj(fin, fout)


def import_pyarrow_parquet():
    """Return the pyarrow.parquet module, or exit with an explicit message if it cannot be imported."""
    try:
        import pyarrow.parquet as pq
    except ImportError as error:
        print(f"{RED}Error: reading Logan Search parquet exports requires pyarrow, which could not be imported "
              f"by {sys.executable}:{NOCOLOR}")
        print(f"{RED}  {error}{NOCOLOR}")
        print(f"{RED}Install it for this interpreter ({sys.executable} -m pip install pyarrow, "
              f"or conda install -c conda-forge pyarrow).{NOCOLOR}")
        sys.exit(1)
    return pq


def extract_accessions_from_parquet(parquet_file, output_file):
    """Write the accessions of a Logan Search parquet export (one per line) into output_file.

    Only the 'accession' column is read. Requires pyarrow.
    Returns the number of written accessions.
    """
    pq = import_pyarrow_parquet()
    nb_accessions = 0
    parquet = pq.ParquetFile(parquet_file)
    with open(output_file, "w") as fout:
        for batch in parquet.iter_batches(columns=["accession"]):
            for accession in batch.column(0).to_pylist():
                if accession and accession.strip():
                    fout.write(f"{accession.strip()}\n")
                    nb_accessions += 1
    return nb_accessions


def read_query_from_parquet(parquet_file):
    """Return the query fasta stored in the metadata of a Logan Search parquet export, or None.

    Recent Logan Search exports store the query in the 'logan.query_fasta' schema metadata.
    """
    pq = import_pyarrow_parquet()
    metadata = pq.ParquetFile(parquet_file).schema_arrow.metadata or {}
    query_fasta = metadata.get(PARQUET_QUERY_FASTA_KEY)
    if not query_fasta:
        return None
    query_fasta = query_fasta.decode("utf-8").strip()
    if not query_fasta.startswith(">"):
        return None
    return query_fasta + "\n"


def fasta_sequences(fasta_text):
    """Return the list of (upper-case) sequences of a fasta text, headers being ignored."""
    sequences = []
    for record in fasta_text.split(">")[1:]:
        lines = record.splitlines()
        sequences.append("".join(line.strip() for line in lines[1:]).upper())
    return sequences


class LoganBlaster:
    LOGAN_DIR_NAME = "logan_data"
    ALIGNEMENT_DIR_NAME = "alignments"
    INPUT_DATA_DIR_NAME = "input_data"
    PREFETCH = 2  # number of accessions downloaded in advance

    def __init__(self, session_id, accession_file, query_file, delete, unitigs, kmer_size, limit, output_dir,
                 api_url=DEFAULT_LOGAN_SEARCH_API):
        self.session_id = session_id
        self.api_url = api_url.rstrip("/")
        self.accession_file = accession_file
        self.query_file = query_file
        self.delete = delete
        self.unitigs = unitigs
        self.kmer_size = kmer_size
        self.limit = limit
        self.main_dir_name = output_dir
        self.type = "unitig" if unitigs else "contig"
        self.failed_accession_list = ""
        self.cli_installed = shutil.which("aws") is not None
        self.stopping = False

    def _setup_directories(self):
        if not self.main_dir_name:
            if self.session_id:
                self.main_dir_name = f"session_{self.session_id}"
            else:
                naming_file = self.query_file if self.query_file else self.accession_file
                query_basename = os.path.basename(naming_file).split(".")[0]
                found_free_name = False
                for i in range(1, 1001):
                    candidate = f"{query_basename}_{i}"
                    if not os.path.exists(candidate):
                        self.main_dir_name = candidate
                        found_free_name = True
                        break
                if not found_free_name:
                    print(f"{RED}Error: Could not find a free directory name based on {query_basename} "
                          f"after 1000 attempts. Please specify an output directory with --output.{NOCOLOR}")
                    sys.exit(1)

        os.makedirs(self.main_dir_name, exist_ok=True)
        os.chdir(self.main_dir_name)
        os.makedirs(self.LOGAN_DIR_NAME, exist_ok=True)
        os.makedirs(self.ALIGNEMENT_DIR_NAME, exist_ok=True)

        if not self.unitigs:
            self.failed_accession_list = "failed_accessions.txt"
            Path(self.failed_accession_list).touch()

    @staticmethod
    def _write_query(source_fasta, abs_query_file, output_file, source_label):
        """Write the query to use into output_file (in the current directory) and return its name.

        source_fasta is the query provided by the input (logan-search.org session or
        brc-analytics parquet export), or None. A user query (abs_query_file) has priority,
        with a warning if its sequence differs from source_fasta.
        """
        if abs_query_file:
            if source_fasta is not None:
                with open(abs_query_file, "r") as f:
                    user_sequences = fasta_sequences(f.read())
                if user_sequences != fasta_sequences(source_fasta):
                    print(f"{YELLOW}[WARNING] The provided query file ({abs_query_file}) differs from the query "
                          f"of {source_label}. The provided query is used.{NOCOLOR}")
            output_file = os.path.basename(abs_query_file)
            if os.path.abspath(output_file) != abs_query_file:
                shutil.copy(abs_query_file, output_file)
            return output_file
        if source_fasta is None:
            print(f"{RED}Error: {source_label} does not contain the query sequence. "
                  f"Please provide it with --query (-q).{NOCOLOR}")
            sys.exit(1)
        with open(output_file, "w") as f:
            f.write(source_fasta)
        print(f"{YELLOW}[INFO] Query of {source_label}: {source_fasta.lstrip('>').split(chr(10), 1)[0]}{NOCOLOR}")
        return output_file

    @staticmethod
    def _write_accessions(accessions, output_file):
        with open(output_file, "w") as f:
            for accession in accessions:
                f.write(f"{accession}\n")

    def _setup_local_files(self, abs_query_file, abs_accession_file):
        os.makedirs(self.INPUT_DATA_DIR_NAME, exist_ok=True)
        os.chdir(self.INPUT_DATA_DIR_NAME)
        input_file = os.path.basename(abs_accession_file)
        if os.path.abspath(input_file) != abs_accession_file:
            shutil.copy(abs_accession_file, input_file)
        basename, extension = os.path.splitext(input_file)
        extension = extension.lower()
        accession_file = input_file
        query_file = None

        if extension == ".parquet":
            print(f"{YELLOW}[INFO] Extracting accession IDs from {input_file}...{NOCOLOR}")
            accession_file = f"{basename}_acc.txt"
            extract_accessions_from_parquet(input_file, accession_file)
            query_file = self._write_query(read_query_from_parquet(input_file), abs_query_file,
                                           f"{basename}_query.fa", input_file)
        elif extension in (".json", ".zip"):
            session_json = input_file
            if extension == ".zip":
                session_json = f"{basename}.json"
                extract_kmviz_session_json(input_file, session_json)
            print(f"{YELLOW}[INFO] Extracting accession IDs and query from {session_json}...{NOCOLOR}")
            query_fasta, accessions = read_kmviz_session(session_json)
            accession_file = f"{basename}_acc.txt"
            self._write_accessions(accessions, accession_file)
            query_file = self._write_query(query_fasta, abs_query_file, f"{basename}_query.fa", input_file)
        else:
            if extension == ".csv":
                print(f"{YELLOW}[INFO] Extracting accession IDs from {input_file}...{NOCOLOR}")
                accession_file = f"{basename}_acc.txt"
                with open(input_file, "r") as fin:
                    fin.readline()
                    accessions = [line.split("\t")[0].strip() for line in fin if line.split("\t")[0].strip()]
                self._write_accessions(accessions, accession_file)
            query_file = self._write_query(None, abs_query_file, None, input_file)

        os.chdir("..")
        self.accession_file = os.path.join(self.INPUT_DATA_DIR_NAME, accession_file)
        self.query_file = os.path.join(self.INPUT_DATA_DIR_NAME, query_file)

    def _setup_session(self, abs_query_file=None):
        if is_kmviz_session(self.session_id):
            self._setup_kmviz_session(abs_query_file)
        else:
            self._setup_search_job(abs_query_file)

    def _setup_search_job(self, abs_query_file):
        """Retrieve the accessions and the query of a brc-analytics Logan Search job.

        Accessions are read from the parquet export of the job, and the query from its
        metadata ('logan.query_fasta'). A user query (-q) has priority.
        """
        job_url = f"{self.api_url}/galaxy/kmindex/jobs/{self.session_id}"
        os.makedirs(self.INPUT_DATA_DIR_NAME, exist_ok=True)
        os.chdir(self.INPUT_DATA_DIR_NAME)

        parquet_file = f"{self.session_id}.parquet"
        if not os.path.exists(parquet_file):
            print(f"{YELLOW}[INFO] Downloading Logan Search results for job {self.session_id} "
                  f"from {self.api_url}...{NOCOLOR}")
            download_file(f"{job_url}/export?format=parquet", f"{parquet_file}.part",
                          hint=f"Job {self.session_id} was looked up on {self.api_url}. A job run on another "
                               f"brc-analytics instance (e.g. dev) is only visible there: pass the full job URL "
                               f"(e.g. https://dev.brc-analytics.org/logan-search?job={self.session_id}) "
                               f"or --api-url.")
            os.replace(f"{parquet_file}.part", parquet_file)
        else:
            print(f"{YELLOW}[INFO] Using existing local version of {parquet_file}...{NOCOLOR}")

        accession_file = f"{self.session_id}_acc.txt"
        nb_accessions = extract_accessions_from_parquet(parquet_file, accession_file)
        print(f"{YELLOW}[INFO] {nb_accessions} accessions extracted from {parquet_file}{NOCOLOR}")
        query_file = self._write_query(read_query_from_parquet(parquet_file), abs_query_file,
                                       f"{self.session_id}_query.fa", f"Logan Search job {self.session_id}")

        os.chdir("..")
        self.accession_file = os.path.join(self.INPUT_DATA_DIR_NAME, accession_file)
        self.query_file = os.path.join(self.INPUT_DATA_DIR_NAME, query_file)

    def _setup_kmviz_session(self, abs_query_file):
        """Retrieve the accessions and the query of a logan-search.org (kmviz) session.

        The session zip archive is downloaded, and its session.json file is parsed.
        """
        os.makedirs(self.INPUT_DATA_DIR_NAME, exist_ok=True)
        os.chdir(self.INPUT_DATA_DIR_NAME)
        session_json = f"{self.session_id}.json"
        if not os.path.exists(session_json):
            zip_file = f"{self.session_id}.zip"
            if not os.path.exists(zip_file):
                print(f"{YELLOW}[INFO] Downloading session data for session ID {self.session_id}...{NOCOLOR}")
                download_file(f"https://logan-search.org/api/download/{self.session_id}", f"{zip_file}.part")
                os.replace(f"{zip_file}.part", zip_file)
            else:
                print(f"{YELLOW}[INFO] Using existing local version of {zip_file}...{NOCOLOR}")
            extract_kmviz_session_json(zip_file, session_json)
        else:
            print(f"{YELLOW}[INFO] Using existing local version of {session_json}...{NOCOLOR}")

        print(f"{YELLOW}[INFO] Extracting accession IDs and query from {session_json}...{NOCOLOR}")
        query_fasta, accessions = read_kmviz_session(session_json)
        accession_file = f"{self.session_id}_acc.txt"
        self._write_accessions(accessions, accession_file)
        print(f"{YELLOW}[INFO] {len(accessions)} accessions extracted from {session_json}{NOCOLOR}")
        query_file = self._write_query(query_fasta, abs_query_file, f"{self.session_id}_query.fa",
                                       f"logan-search.org session {self.session_id}")

        os.chdir("..")
        self.accession_file = os.path.join(self.INPUT_DATA_DIR_NAME, accession_file)
        self.query_file = os.path.join(self.INPUT_DATA_DIR_NAME, query_file)

    def _run_coverage_stats(self, fasta_file, label):
        if shutil.which("count_logan_tig_coverage") is None:
            return
        print(f"{YELLOW}[INFO] Coverage statistics for {label}:{NOCOLOR}")
        cmd = ["count_logan_tig_coverage", "--in", fasta_file]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"{YELLOW}[WARNING] count_logan_tig_coverage failed: {result.stderr.strip()}{NOCOLOR}")
        else:
            for line in result.stdout.splitlines():
                print(f"  {line}")

    def _download_command(self, accession, destination):
        prefix = "c" if not self.unitigs else "u"
        file_name = f"{accession}.{self.type}s.fa.zst"
        if self.cli_installed:
            return ["aws", "s3", "cp", f"s3://logan-pub/{prefix}/{accession}/{file_name}", destination,
                    "--no-sign-request", "--region", "us-east-1", "--only-show-errors"]
        return ["wget", "--quiet", "-O", destination, f"https://s3.amazonaws.com/logan-pub/{prefix}/{accession}/{file_name}"]

    def _fetch_accession(self, accession):
        """Download the Logan file of an accession into LOGAN_DIR_NAME (runs in a worker thread).

        The file is first written to a '.part' file then atomically renamed, so the final path
        exists only when the download is complete. Never prints: sys.stdout may be redirected
        by the main thread (see _run_blast). Returns (status, messages) with status in
        {"existing", "downloaded", "failed"}.
        """
        local_file = os.path.join(self.LOGAN_DIR_NAME, f"{accession}.{self.type}s.fa.zst")
        if os.path.exists(local_file):
            return "existing", []

        part_file = f"{local_file}.part"
        cmd_dl = self._download_command(accession, part_file)
        messages = [f"{GREEN}Running command: {' '.join(cmd_dl)}{NOCOLOR}"]
        for attempt in range(3):
            if self.stopping:
                break
            result = subprocess.run(cmd_dl, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
            if result.returncode == 0 and os.path.exists(part_file):
                os.replace(part_file, local_file)
                return "downloaded", messages
            messages.append(f"{YELLOW}[WARNING] Attempt {attempt + 1} download failed for {accession}.{self.type}s.fa.zst. "
                            f"{result.stderr.strip()}{NOCOLOR}")

        try:
            os.remove(part_file)
        except FileNotFoundError:
            pass
        return "failed", messages

    def _run_blast(self, query_fasta, target_fasta):
        query_basename = os.path.basename(query_fasta).split(".")[0]
        target_basename = os.path.basename(target_fasta).split(".")[0]
        with open(query_fasta, "r") as f:
            query_id = f.readline().strip().lstrip(">").split()[0]

        output_name = f"{query_id}_vs_{target_basename}.txt"
        print(f"{YELLOW}[INFO] Aligning {target_basename} vs {query_basename}...{NOCOLOR}")

        cmd = [
            "blastn",
            "-query", query_fasta,
            "-subject", target_fasta,
            "-out", os.path.join(self.ALIGNEMENT_DIR_NAME, output_name),
            "-outfmt", "0",
            "-sorthits", "0",
            "-word_size", "11",
            "-gapextend", "2",
            "-gapopen", "5",
            "-reward", "2",
            "-penalty", "-3"
        ]

        print(f"{GREEN}Running command: {' '.join(cmd)}{NOCOLOR}")
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=open("error.log", "w"))
        except subprocess.CalledProcessError:
            print(f"{RED}Error: blastn failed{NOCOLOR}")
            with open("error.log", "r") as f:
                print(f.read())
            return

        print(f"{YELLOW}[INFO] Synthesize blast results{NOCOLOR}")
        blastn_file = os.path.join(self.ALIGNEMENT_DIR_NAME, output_name)
        synth_file = os.path.join(self.ALIGNEMENT_DIR_NAME, f"synth_{output_name}")
        with open(synth_file, "w") as f:
            sys.stdout = f
            run_blast_parser(self.query_file, blastn_file, abundance=True)
            sys.stdout = sys.__stdout__

    def _read_accessions(self):
        accessions = []
        seen = set()
        with open(self.accession_file, "r") as f:
            for line in f:
                if not line.strip():
                    continue
                accession = line.strip().split()[0]
                if accession in seen:
                    print(f"{YELLOW}[WARNING] Duplicated accession {accession} ignored.{NOCOLOR}")
                    continue
                seen.add(accession)
                accessions.append(accession)
        if self.limit != 0:
            accessions = accessions[:self.limit]
        return accessions

    def _process_accessions(self):
        accessions = self._read_accessions()
        futures = {}
        executor = ThreadPoolExecutor(max_workers=self.PREFETCH)
        try:
            for idx, accession in enumerate(accessions):
                # Keep the current accession and the next PREFETCH ones scheduled for download
                for next_idx in range(idx, min(idx + 1 + self.PREFETCH, len(accessions))):
                    if next_idx not in futures:
                        futures[next_idx] = executor.submit(self._fetch_accession, accessions[next_idx])
                self._process_accession(accession, futures.pop(idx))
        finally:
            self.stopping = True
            for future in futures.values():
                future.cancel()
            executor.shutdown(wait=True)

    def _process_accession(self, accession, download_future):
        print(f"\n{BLUE}=========================================={NOCOLOR}")
        print(f"{CYAN}>>> Processing accession: {accession} <<<{NOCOLOR}")
        print(f"{BLUE}=========================================={NOCOLOR}")

        local_file = os.path.join(self.LOGAN_DIR_NAME, f"{accession}.{self.type}s.fa.zst")
        if not download_future.done():
            print(f"{YELLOW}[INFO] Waiting for download of {accession}.{self.type}s.fa.zst...{NOCOLOR}")
        status, messages = download_future.result()
        for message in messages:
            print(message)

        if status == "existing":
            print(f"{YELLOW}[INFO] Using existing local version of {local_file}...{NOCOLOR}")
        elif status == "downloaded":
            print(f"{YELLOW}[INFO] Downloaded {local_file}{NOCOLOR}")
        else:
            print(f"{RED}Error: Failed to download {accession}.{self.type}s.fa.zst after 3 attempts.{NOCOLOR}")
            if not self.unitigs:
                print(f"{YELLOW}[INFO] Contigs do not exist for accession {accession}. Adding {accession} to failed accession list.{NOCOLOR}")
                with open(self.failed_accession_list, "a") as f:
                    f.write(f"{accession}\n")
            return

        self._run_coverage_stats(local_file, f"downloaded {self.type}s ({accession})")

        recruited_file = os.path.join(self.LOGAN_DIR_NAME, f"{accession}.recruited_{self.type}s.fa")
        print(f"{YELLOW}[INFO] Recruiting sequences from {accession}.{self.type}s.fa.zst with a match with {self.query_file}...{NOCOLOR}")
        cmd_recruit = [
            "back_to_sequences",
            "--kmer-size", str(self.kmer_size),
            "--in-kmers", self.query_file,
            "--in-sequences", local_file,
            "--out-sequences", recruited_file
        ]
        print(f"{GREEN}Running command: {' '.join(cmd_recruit)}{NOCOLOR}")
        try:
            subprocess.run(cmd_recruit, check=True, stdout=subprocess.DEVNULL, stderr=open("error.log", "w"))
        except subprocess.CalledProcessError:
            print(f"{RED}Error: back_to_sequences failed for accession {accession}.{NOCOLOR}")
            try:
                os.remove(recruited_file)
                os.remove(local_file)
            except FileNotFoundError:
                pass
            with open("error.log", "r") as f:
                print(f.read())
            with open(self.failed_accession_list, "a") as f:
                f.write(f"{accession}\n")
            return

        if os.path.getsize(recruited_file) == 0:
            print(f"{YELLOW}[INFO]\tNo sequences were recruited from {accession}.{self.type}s.fa.zst. Skipping BLAST step.{NOCOLOR}")
            if self.delete:
                print(f"{YELLOW}[INFO] Deleting {recruited_file} and {local_file}...{NOCOLOR}")
                try:
                    os.remove(recruited_file)
                    os.remove(local_file)
                except FileNotFoundError:
                    pass
            if self.type == "contig":
                with open(self.failed_accession_list, "a") as f:
                    f.write(f"{accession}\n")
            return

        self._run_coverage_stats(recruited_file, f"recruited {self.type}s ({accession})")

        print(f"{YELLOW}[INFO] Aligning recruited sequences from {accession}.{self.type}s.fa.zst with {self.query_file}...{NOCOLOR}")
        self._run_blast(self.query_file, recruited_file)

        if self.delete:
            print(f"{YELLOW}[INFO] Deleting {recruited_file} and {local_file}...{NOCOLOR}")
            try:
                os.remove(recruited_file)
                os.remove(local_file)
            except FileNotFoundError:
                pass

    def run(self, abs_query_file=None, abs_accession_file=None):
        self._setup_directories()

        if not self.session_id:
            self._setup_local_files(abs_query_file, abs_accession_file)
        else:
            self._setup_session(abs_query_file)

        if not os.path.exists(self.query_file):
            print(f"{RED}Error: Query file '{self.query_file}' does not exist.{NOCOLOR}")
            shutil.rmtree(self.main_dir_name, ignore_errors=True)
            sys.exit(1)

        if not os.path.exists(self.accession_file):
            print(f"{RED}Error: Accessions file '{self.accession_file}' does not exist.{NOCOLOR}")
            sys.exit(1)

        self._process_accessions()


def main():
    parser = argparse.ArgumentParser(description="Process Logan session or accession/query files.")
    parser.add_argument("-s", "--session", type=str,
                        help="brc-analytics Logan Search job ID or job URL "
                             "(e.g. https://brc-analytics.org/logan-search?job=<job_id>), "
                             "or logan-search.org session ID (e.g. kmviz-4563ab1a-82d7-4463-bedd-bd1604a5c981)")
    parser.add_argument("-a", "--accessions", type=str,
                        help="Path to accessions.txt file, or .csv file (containing a first header line, ignored, "
                             "and storing accessions as the first column), or brc-analytics Logan Search .parquet "
                             "export, or logan-search.org session (.zip archive or its session.json file)")
    parser.add_argument("-q", "--query", type=str,
                        help="Path to query fasta file. Optional with a session / job ID, a .parquet, .json or .zip "
                             "file, that provide the query (if given, it is used instead of the provided query)")
    parser.add_argument("-o", "--output", type=str, default=None, help="Output directory name (default: based on query name if using --accessions and --query or session ID if using --session)")
    parser.add_argument("-u", "--unitigs", action="store_true", help="Use unitigs instead of contigs")
    parser.add_argument("-k", "--kmer-size", type=int, default=17, help="K-mer size for sequence recruitment")
    parser.add_argument("-l", "--limit", type=int, default=0, help="Limit number of accessions to process")
    parser.add_argument("--api-url", type=str, default=None,
                        help="Logan Search API base URL, used with a Logan Search job ID "
                             f"(default: derived from the job URL if given, else {DEFAULT_LOGAN_SEARCH_API})")
    parser.add_argument("-d", "--delete", action="store_true", help="Delete intermediate files after processing")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args()

    # Argument checks
    api_url = args.api_url
    if args.session:
        args.session, url_api = parse_session_argument(args.session)
        if api_url is None:
            api_url = url_api
    if api_url is None:
        api_url = DEFAULT_LOGAN_SEARCH_API

    query_in_accessions = bool(args.accessions) and args.accessions.lower().endswith((".parquet", ".json", ".zip"))
    if not args.session and (not args.accessions or (not args.query and not query_in_accessions)):
        print(f"{RED}Error: You must provide either --session (-s), or both --accessions (-a) and --query (-q) "
              f"(--query is optional with a .parquet, .json or .zip file).{NOCOLOR}")
        sys.exit(1)

    if args.session and args.accessions:
        print(f"{RED}Error: --session (-s) cannot be combined with --accessions (-a).{NOCOLOR}")
        sys.exit(1)

    if args.accessions and not os.path.isfile(args.accessions):
        print(f"{RED}Error: Accessions file '{args.accessions}' does not exist.{NOCOLOR}")
        sys.exit(1)

    if args.query and not os.path.isfile(args.query):
        print(f"{RED}Error: Query file '{args.query}' does not exist.{NOCOLOR}")
        sys.exit(1)

    if args.kmer_size <= 0:
        print(f"{RED}Error: K-mer size must be a positive integer.{NOCOLOR}")
        sys.exit(1)

    if args.limit < 0:
        print(f"{RED}Error: Limit must be a non-negative integer.{NOCOLOR}")
        sys.exit(1)

    if shutil.which("count_logan_tig_coverage") is None:
        print(f"{YELLOW}[WARNING] 'count_logan_tig_coverage' not found — coverage statistics will be skipped.{NOCOLOR}")

    for cmd in ["back_to_sequences", "blastn"]:
        if shutil.which(cmd) is None:
            print(f"{RED}Error: '{cmd}' could not be found. Please install it and ensure it's in your PATH.{NOCOLOR}")
            sys.exit(1)

    # Resolve absolute paths before any os.chdir
    abs_query_file = os.path.abspath(args.query) if args.query else None
    abs_accession_file = os.path.abspath(args.accessions) if args.accessions else None

    blaster = LoganBlaster(
        session_id=args.session,
        accession_file=args.accessions,
        query_file=args.query,
        delete=args.delete,
        unitigs=args.unitigs,
        kmer_size=args.kmer_size,
        limit=args.limit,
        output_dir=args.output,
        api_url=api_url,
    )
    blaster.run(abs_query_file=abs_query_file, abs_accession_file=abs_accession_file)

    print(f"\n{BLUE}================")
    print(f"{CYAN}>>> All done <<<")
    print(f"{BLUE}================\n")

    if not args.delete:
        print(f"{YELLOW}[INFO] You did not use --delete option. So you can manually remove all intermediate files "
              f"(recruited {blaster.type}s and {blaster.type}s files) by running:{NOCOLOR}")
        print(f"rm -rf {blaster.main_dir_name}/{LoganBlaster.LOGAN_DIR_NAME}")

    os.chdir("..")
    print(f"{YELLOW}[INFO] Results can be found in directory {CYAN}{blaster.main_dir_name}{NOCOLOR}")

    if not args.unitigs and os.path.getsize(f"{blaster.main_dir_name}/{blaster.failed_accession_list}") > 0:
        nb_failed = len(open(f"{blaster.main_dir_name}/{blaster.failed_accession_list}").readlines())
        print(f"{YELLOW}[INFO] {nb_failed} accession{'s' if nb_failed > 1 else ''} failed to download contigs or had no recruited sequences.{NOCOLOR}")
        print(f"{YELLOW}[INFO] List of failed accessions: {CYAN}{blaster.main_dir_name}/{blaster.failed_accession_list}{NOCOLOR}")
        print(f"{YELLOW}[INFO] You can try to re-run the script with --unitigs option and this accession list.{NOCOLOR}")
        print(f"{YELLOW}[INFO] Command example:{NOCOLOR}")
        delete_flag = "-d" if args.delete else ""
        print(f"logan_blaster -a {blaster.main_dir_name}/{blaster.failed_accession_list} -q {blaster.main_dir_name}/{blaster.query_file} --unitigs -k {args.kmer_size} {delete_flag} {NOCOLOR}")


if __name__ == "__main__":
    main()
