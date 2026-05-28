#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import shutil
from urllib.request import urlopen
from pathlib import Path
import ssl

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


def download_file(url, destination):
    try:
        context = ssl._create_unverified_context()
        with urlopen(url, context=context) as response, open(destination, "wb") as f:
            while chunk := response.read(8192):
                f.write(chunk)
    except Exception as e:
        print(f"{RED}Error: Failed to download {url}.{NOCOLOR}")
        print(e)
        sys.exit(1)


class LoganBlaster:
    LOGAN_DIR_NAME = "logan_data"
    ALIGNEMENT_DIR_NAME = "alignments"
    INPUT_DATA_DIR_NAME = "input_data"

    def __init__(self, session_id, accession_file, query_file, delete, unitigs, kmer_size, limit, output_dir):
        self.session_id = session_id
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

    def _setup_directories(self):
        if not self.main_dir_name:
            if self.session_id:
                self.main_dir_name = f"session_{self.session_id}"
            else:
                query_basename = os.path.basename(self.query_file).split(".")[0]
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

    def _setup_local_files(self, abs_query_file, abs_accession_file):
        os.makedirs(self.INPUT_DATA_DIR_NAME, exist_ok=True)
        shutil.copy(abs_query_file, self.INPUT_DATA_DIR_NAME)
        shutil.copy(abs_accession_file, self.INPUT_DATA_DIR_NAME)
        self.query_file = os.path.join(self.INPUT_DATA_DIR_NAME, os.path.basename(abs_query_file))
        self.accession_file = os.path.join(self.INPUT_DATA_DIR_NAME, os.path.basename(abs_accession_file))

        if self.accession_file.endswith(".csv"):
            accession_basename = self.accession_file.rstrip(".csv")
            print(f"{YELLOW}[INFO] Extracting accession IDs from {accession_basename}.csv...{NOCOLOR}")
            cmd_jq_acc = f"cat {accession_basename}.csv | tail -n +2 | cut -f 1 > {accession_basename}_acc.txt"
            subprocess.run(cmd_jq_acc, shell=True, check=True)
            self.accession_file = os.path.join(
                self.INPUT_DATA_DIR_NAME, os.path.basename(f"{accession_basename}_acc.txt")
            )

    def _setup_session(self):
        os.makedirs(self.INPUT_DATA_DIR_NAME, exist_ok=True)
        os.chdir(self.INPUT_DATA_DIR_NAME)
        accession_file = f"accessions_{self.session_id}.txt"
        if not os.path.exists(accession_file):
            zip_file = f"{self.session_id}.zip"
            if not os.path.exists(zip_file):
                print(f"{YELLOW}[INFO] Downloading session data for session ID {self.session_id}...{NOCOLOR}")
                download_file(f"https://logan-search.org/api/download/{self.session_id}", zip_file)
            else:
                print(f"{YELLOW}[INFO] Using existing local version of {zip_file}...{NOCOLOR}")

            subprocess.run(["unzip", "-o", zip_file, "session.json"], check=True)
            os.rename("session.json", f"{self.session_id}.json")

            print(f"{YELLOW}[INFO] Extracting accession IDs from {self.session_id}.json...{NOCOLOR}")
            cmd_jq_acc = f"jq -r '.. | ._metadata?.ID? // empty | .[]' {self.session_id}.json > {self.session_id}_acc.txt"
            subprocess.run(cmd_jq_acc, shell=True, check=True)

            print(f"{YELLOW}[INFO] Extracting query name and sequence from {self.session_id}.json...{NOCOLOR}")
            with open(f"{self.session_id}_query.fa", "w") as f:
                cmd_jq_name = f"jq -r '.. | ._query?._name? // empty' {self.session_id}.json"
                query_name = subprocess.check_output(cmd_jq_name, shell=True, text=True).strip()
                f.write(f">{query_name}\n")

                cmd_jq_seq = f"jq -r '.. | ._query?._seq? // empty' {self.session_id}.json"
                query_seq = subprocess.check_output(cmd_jq_seq, shell=True, text=True).strip()
                f.write(f"{query_seq}\n")

        os.chdir("..")
        self.accession_file = os.path.join(self.INPUT_DATA_DIR_NAME, f"{self.session_id}_acc.txt")
        self.query_file = os.path.join(self.INPUT_DATA_DIR_NAME, f"{self.session_id}_query.fa")

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

    def _process_accessions(self):
        counter = 0

        with open(self.accession_file, "r") as f:
            for accession in f:
                accession = accession.strip().split()[0]

                if self.limit != 0 and counter >= self.limit:
                    print(f"\n{YELLOW}[INFO] Reached limit of {self.limit} accessions. Stopping further processing.{NOCOLOR}")
                    break
                counter += 1

                print(f"\n{BLUE}=========================================={NOCOLOR}")
                print(f"{CYAN}>>> Processing accession: {accession} <<<{NOCOLOR}")
                print(f"{BLUE}=========================================={NOCOLOR}")

                local_file = os.path.join(self.LOGAN_DIR_NAME, f"{accession}.{self.type}s.fa.zst")
                print(f"{YELLOW}[INFO] Checking for local file {local_file}...{NOCOLOR}")
                if not os.path.exists(local_file):
                    print(f"{YELLOW}[INFO] Downloading {accession}.{self.type}s.fa.zst...{NOCOLOR}")
                    if self.cli_installed:
                        if not self.unitigs:
                            cmd_dl = f"aws s3 cp s3://logan-pub/c/{accession}/{accession}.contigs.fa.zst . --no-sign-request"
                        else:
                            cmd_dl = f"aws s3 cp s3://logan-pub/u/{accession}/{accession}.unitigs.fa.zst . --no-sign-request"
                    else:
                        if not self.unitigs:
                            cmd_dl = f"wget https://s3.amazonaws.com/logan-pub/c/{accession}/{accession}.contigs.fa.zst"
                        else:
                            cmd_dl = f"wget https://s3.amazonaws.com/logan-pub/u/{accession}/{accession}.unitigs.fa.zst"

                    print(f"{GREEN}Running command: {cmd_dl}{NOCOLOR}")
                    for attempt in range(3):
                        try:
                            subprocess.run(cmd_dl, shell=True, check=True)
                            break
                        except subprocess.CalledProcessError:
                            print(f"{YELLOW}[WARNING] Attempt {attempt + 1} download failed for {accession}.{self.type}s.fa.zst. {NOCOLOR}")

                    if not os.path.exists(f"{accession}.{self.type}s.fa.zst"):
                        print(f"{RED}Error: Failed to download {accession}.{self.type}s.fa.zst after 3 attempts.{NOCOLOR}")
                        if not self.unitigs:
                            print(f"{YELLOW}[INFO] Contigs do not exist for accession {accession}. Adding {accession} to failed accession list.{NOCOLOR}")
                            with open(self.failed_accession_list, "a") as f:
                                f.write(f"{accession}\n")
                        continue
                    shutil.move(f"{accession}.{self.type}s.fa.zst", self.LOGAN_DIR_NAME)
                else:
                    print(f"{YELLOW}[INFO] Using existing local version of {local_file}...{NOCOLOR}")

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
                    continue

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
                    continue

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
            self._setup_session()

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
    parser.add_argument("-s", "--session", type=str, help="Logan session ID")
    parser.add_argument("-a", "--accessions", type=str, help="Path to accessions.txt file or .csv file (containing a first header line, ignored, and storing accessions as the first column)")
    parser.add_argument("-q", "--query", type=str, help="Path to query fasta file")
    parser.add_argument("-o", "--output", type=str, default=None, help="Output directory name (default: based on query name if using --accessions and --query or session ID if using --session)")
    parser.add_argument("-u", "--unitigs", action="store_true", help="Use unitigs instead of contigs")
    parser.add_argument("-k", "--kmer-size", type=int, default=17, help="K-mer size for sequence recruitment")
    parser.add_argument("-l", "--limit", type=int, default=0, help="Limit number of accessions to process")
    parser.add_argument("-d", "--delete", action="store_true", help="Delete intermediate files after processing")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args()

    # Argument checks
    if not args.session and (not args.accessions or not args.query):
        print(f"{RED}Error: You must provide either --session (-s) or both --accessions (-a) and --query (-q).{NOCOLOR}")
        sys.exit(1)

    if args.session and (args.accessions or args.query):
        print(f"{RED}Error: --session (-s) cannot be combined with --accessions (-a) or --query (-q).{NOCOLOR}")
        sys.exit(1)

    if args.kmer_size <= 0:
        print(f"{RED}Error: K-mer size must be a positive integer.{NOCOLOR}")
        sys.exit(1)

    if args.limit < 0:
        print(f"{RED}Error: Limit must be a non-negative integer.{NOCOLOR}")
        sys.exit(1)

    for cmd in ["back_to_sequences", "blastn", "jq"]:
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
