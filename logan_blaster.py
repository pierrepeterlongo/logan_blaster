#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import shutil
from urllib.request import urlopen
from pathlib import Path
import ssl

__version__ = '0.1.1'
__author__ = 'Pierre Peterlongo'

# --- Using blast_parser ---
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
    """Run the blast parser on the given fasta and blastn result files."""
    query_ACGT = get_query_ACGT(fasta_file)
    query_name, query_length, matched_positions = parse_blastn(blastn_file)
    assert len(query_ACGT) == query_length, f"Error, query given in {fasta_file} is of length {len(query_ACGT)}, while the blastn result indicates a sequence of length {query_length}"
    visualize_matches(query_ACGT, query_name, query_length, matched_positions, print_abundance=abundance)


# Message colors
GREEN = "\033[0;32m"
RED = "\033[0;31m"
YELLOW = "\033[0;33m"
BLUE = "\033[1;34m"
CYAN = "\033[1;36m"
NOCOLOR = "\033[0m"

# (Horrible) global variables
SESSION_ID = ""
ACCESSION_FILE = ""
QUERY_FILE = ""
DELETE = False
UNITIGS = False
LIMIT = 0
KMER_SIZE = 17
MAIN_DIR_NAME = ""
LOGAN_DIR_NAME = "logan_data"
ALIGNEMENT_DIR_NAME = "alignments"
INPUT_DATA_DIR_NAME = "input_data"
TYPE = ""
cli_installed = shutil.which("aws") is not None
FAILED_ACCESSION_LIST = ""


def run_blast(query_fasta, target_fasta):
    """Run blastn to align query_fasta against target_fasta."""
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
        "-out", os.path.join(ALIGNEMENT_DIR_NAME, output_name),
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

    # Run the integrated blast_parser on blast results
    print(f"{YELLOW}[INFO] Synthesize blast results{NOCOLOR}")
    blastn_file = os.path.join(ALIGNEMENT_DIR_NAME, output_name)
    synth_file = os.path.join(ALIGNEMENT_DIR_NAME, f"synth_{output_name}")
    with open(synth_file, "w") as f:
        sys.stdout = f
        run_blast_parser(QUERY_FILE, blastn_file, abundance=True)
        sys.stdout = sys.__stdout__


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

def process_accessions():
    counter = 0

    with open(ACCESSION_FILE, "r") as f:
        for accession in f:
            accession = accession.strip()
            if LIMIT != 0 and counter >= LIMIT:
                print(f"\n{YELLOW}[INFO] Reached limit of {LIMIT} accessions. Stopping further processing.{NOCOLOR}")
                break
            counter += 1

            print(f"\n{BLUE}=========================================={NOCOLOR}")
            print(f"{CYAN}>>> Processing accession: {accession} <<<{NOCOLOR}")
            print(f"{BLUE}=========================================={NOCOLOR}")

            local_file = os.path.join(LOGAN_DIR_NAME, f"{accession}.{TYPE}s.fa.zst")
            print(f"{YELLOW}[INFO] Checking for local file {local_file}...{NOCOLOR}")
            if not os.path.exists(local_file):
                print(f"{YELLOW}[INFO] Downloading {accession}.{TYPE}s.fa.zst...{NOCOLOR}")
                if cli_installed:
                    if not UNITIGS:
                        cmd_dl = f"aws s3 cp s3://logan-pub/c/{accession}/{accession}.contigs.fa.zst . --no-sign-request"
                    else:
                        cmd_dl = f"aws s3 cp s3://logan-pub/u/{accession}/{accession}.unitigs.fa.zst . --no-sign-request"
                else:
                    if not UNITIGS:
                        cmd_dl = f"wget https://s3.amazonaws.com/logan-pub/c/{accession}/{accession}.contigs.fa.zst"
                    else:
                        cmd_dl = f"wget https://s3.amazonaws.com/logan-pub/u/{accession}/{accession}.unitigs.fa.zst"

                print(f"{GREEN}Running command: {cmd_dl}{NOCOLOR}")
                # try to download the file, 3 attempts
                for attempt in range(3):
                    try:
                        subprocess.run(cmd_dl, shell=True, check=True)
                        break
                    except subprocess.CalledProcessError:
                        print(f"{YELLOW}[WARNING] Attempt {attempt + 1} download failed for {accession}.{TYPE}s.fa.zst. {NOCOLOR}")
                    
                if not os.path.exists(f"{accession}.{TYPE}s.fa.zst"):
                    print(f"{RED}Error: Failed to download {accession}.{TYPE}s.fa.zst after 3 attempts.{NOCOLOR}")
                    if not UNITIGS:
                        print(f"{YELLOW}[INFO] Contigs do not exist for accession {accession}. Adding {accession} to failed accession list.{NOCOLOR}")
                        with open(FAILED_ACCESSION_LIST, "a") as f:
                            f.write(f"{accession}\n")
                    continue      
                shutil.move(f"{accession}.{TYPE}s.fa.zst", LOGAN_DIR_NAME)
            else:
                print(f"{YELLOW}[INFO] Using existing local version of {local_file}...{NOCOLOR}")

            recruited_file = os.path.join(LOGAN_DIR_NAME, f"{accession}.recruited_{TYPE}s.fa")
            print(f"{YELLOW}[INFO] Recruiting sequences from {accession}.{TYPE}s.fa.zst with a match with {QUERY_FILE}...{NOCOLOR}")
            cmd_recruit = [
                "back_to_sequences",
                "--kmer-size", str(KMER_SIZE),
                "--in-kmers", QUERY_FILE,
                "--in-sequences", local_file,
                "--out-sequences", recruited_file
            ]
            print(f"{GREEN}Running command: {' '.join(cmd_recruit)}{NOCOLOR}")
            # Run back_to_sequences
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
                with open(FAILED_ACCESSION_LIST, "a") as f:
                    f.write(f"{accession}\n")
                continue

            if os.path.getsize(recruited_file) == 0:
                print(f"{YELLOW}[INFO]\tNo sequences were recruited from {accession}.{TYPE}s.fa.zst. Skipping BLAST step.{NOCOLOR}")
                if DELETE:
                    print(f"{YELLOW}[INFO] Deleting {recruited_file} and {local_file}...{NOCOLOR}")
                    try: 
                        os.remove(recruited_file)
                        os.remove(local_file)
                    except FileNotFoundError:
                        pass
                if TYPE == "contig":
                    with open(FAILED_ACCESSION_LIST, "a") as f:
                        f.write(f"{accession}\n")
                continue

            print(f"{YELLOW}[INFO] Aligning recruited sequences from {accession}.{TYPE}s.fa.zst with {QUERY_FILE}...{NOCOLOR}")
            run_blast(QUERY_FILE, recruited_file)

            if DELETE:
                print(f"{YELLOW}[INFO] Deleting {recruited_file} and {local_file}...{NOCOLOR}")
                try:
                    os.remove(recruited_file)
                    os.remove(local_file)
                except FileNotFoundError:
                    pass

def main():
    global SESSION_ID, ACCESSION_FILE, QUERY_FILE, DELETE, UNITIGS, LIMIT, KMER_SIZE, MAIN_DIR_NAME, TYPE, FAILED_ACCESSION_LIST

    parser = argparse.ArgumentParser(description="Process Logan session or accession/query files.")
    parser.add_argument("-s", "--session", type=str, help="Logan session ID")
    parser.add_argument("-a", "--accessions", type=str, help="Path to accessions.txt file")
    parser.add_argument("-q", "--query", type=str, help="Path to query fasta file")
    parser.add_argument("-o", "--output", type=str, default=None, help="Output directory name (default: based on query name if using --accessions and --query or session ID if using --session)")
    parser.add_argument("-u", "--unitigs", action="store_true", help="Use unitigs instead of contigs")
    parser.add_argument("-k", "--kmer-size", type=int, default=17, help="K-mer size for sequence recruitment")
    parser.add_argument("-l", "--limit", type=int, default=0, help="Limit number of accessions to process")
    parser.add_argument("-d", "--delete", action="store_true", help="Delete intermediate files after processing")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args()

    SESSION_ID = args.session
    ACCESSION_FILE = args.accessions
    QUERY_FILE = args.query
    DELETE = args.delete
    UNITIGS = args.unitigs
    KMER_SIZE = args.kmer_size
    LIMIT = args.limit
    # get the full absolute path of QUERY_FILE and ACCESSION_FILE

    # Argument checks
    if not SESSION_ID and (not ACCESSION_FILE or not QUERY_FILE):
        print(f"{RED}Error: You must provide either --session (-s) or both --accessions (-a) and --query (-q).{NOCOLOR}")
        # print_help()
        sys.exit(1)

    if SESSION_ID and (ACCESSION_FILE or QUERY_FILE):
        print(f"{RED}Error: --session (-s) cannot be combined with --accessions (-a) or --query (-q).{NOCOLOR}")
        # print_help()
        sys.exit(1)

    # Define the working directory
    if args.output:
        MAIN_DIR_NAME = args.output
        ABS_QUERY_FILE = os.path.abspath(QUERY_FILE)
        ABS_ACCESSION_FILE = os.path.abspath(ACCESSION_FILE)
            
    else:
        if SESSION_ID:
            MAIN_DIR_NAME = f"session_{SESSION_ID}"
        else:
            ABS_QUERY_FILE = os.path.abspath(QUERY_FILE)
            ABS_ACCESSION_FILE = os.path.abspath(ACCESSION_FILE)
            query_basename = os.path.basename(QUERY_FILE).split(".")[0]
            found_free_name = False
            for i in range(1, 1001):
                MAIN_DIR_NAME = f"{query_basename}_{i}"
                if not os.path.exists(MAIN_DIR_NAME):
                    found_free_name = True
                    break
            if not found_free_name:
                print(f"{RED}Error: Could not find a free directory name based on {query_basename} after 1000 attempts. Please specify an output directory with --output.{NOCOLOR}")
                sys.exit(1)
    
    os.makedirs(MAIN_DIR_NAME, exist_ok=True)
    os.chdir(MAIN_DIR_NAME)
    os.makedirs(LOGAN_DIR_NAME, exist_ok=True)
    os.makedirs(ALIGNEMENT_DIR_NAME, exist_ok=True)
    if not UNITIGS:
        FAILED_ACCESSION_LIST = "failed_accessions.txt"
        Path(FAILED_ACCESSION_LIST).touch()
    
    # If no session ID, copy input files to input_data directory
    if not SESSION_ID:
        os.makedirs(INPUT_DATA_DIR_NAME, exist_ok=True)
        shutil.copy(f"{ABS_QUERY_FILE}", INPUT_DATA_DIR_NAME)
        shutil.copy(f"{ABS_ACCESSION_FILE}", INPUT_DATA_DIR_NAME)
        QUERY_FILE = os.path.join(INPUT_DATA_DIR_NAME, os.path.basename(QUERY_FILE))
        ACCESSION_FILE = os.path.join(INPUT_DATA_DIR_NAME, os.path.basename(ACCESSION_FILE))

    # If session ID is provided, download and extract data
    if SESSION_ID:
        os.makedirs(INPUT_DATA_DIR_NAME, exist_ok=True)
        os.chdir(INPUT_DATA_DIR_NAME)
        accession_file = f"accessions_{SESSION_ID}.txt"
        if not os.path.exists(accession_file):
            zip_file = f"{SESSION_ID}.zip"
            if not os.path.exists(zip_file):
                print(f"{YELLOW}[INFO] Downloading session data for session ID {SESSION_ID}...{NOCOLOR}")
                download_file(f"https://logan-search.org/api/download/{SESSION_ID}", zip_file)
            else:
                print(f"{YELLOW}[INFO] Using existing local version of {zip_file}...{NOCOLOR}")

            subprocess.run(["unzip", "-o", zip_file, "session.json"], check=True)
            os.rename("session.json", f"{SESSION_ID}.json")

            # Using jq to extract accession IDs
            print(f"{YELLOW}[INFO] Extracting accession IDs from {SESSION_ID}.json...{NOCOLOR}")
            cmd_jq_acc = f"jq -r '.. | ._metadata?.ID? // empty | .[]' {SESSION_ID}.json > {SESSION_ID}_acc.txt"
            subprocess.run(cmd_jq_acc, shell=True, check=True)

            # Using jq to extract the query name and sequence
            print(f"{YELLOW}[INFO] Extracting query name and sequence from {SESSION_ID}.json...{NOCOLOR}")
            with open(f"{SESSION_ID}_query.fa", "w") as f:
                cmd_jq_name = f"jq -r '.. | ._query?._name? // empty' {SESSION_ID}.json"
                query_name = subprocess.check_output(cmd_jq_name, shell=True, text=True).strip()
                f.write(f">{query_name}\n")

                cmd_jq_seq = f"jq -r '.. | ._query?._seq? // empty' {SESSION_ID}.json"
                query_seq = subprocess.check_output(cmd_jq_seq, shell=True, text=True).strip()
                f.write(f"{query_seq}\n")

            ACCESSION_FILE = f"{SESSION_ID}_acc.txt"
            QUERY_FILE = f"{SESSION_ID}_query.fa"
        os.chdir("..")
        ACCESSION_FILE = os.path.join(INPUT_DATA_DIR_NAME, f"{SESSION_ID}_acc.txt")
        QUERY_FILE = os.path.join(INPUT_DATA_DIR_NAME, f"{SESSION_ID}_query.fa")

    # Additional checks
    if not os.path.exists(QUERY_FILE):
        print(f"{RED}Error: Query file '{QUERY_FILE}' does not exist.{NOCOLOR}")
        shutil.rmtree(MAIN_DIR_NAME, ignore_errors=True)
        sys.exit(1)

    if not os.path.exists(ACCESSION_FILE):
        print(f"{RED}Error: Accessions file '{ACCESSION_FILE}' does not exist.{NOCOLOR}")
        sys.exit(1)

    if KMER_SIZE <= 0:
        print(f"{RED}Error: K-mer size must be a positive integer.{NOCOLOR}")
        sys.exit(1)

    if LIMIT < 0:
        print(f"{RED}Error: Limit must be a non-negative integer.{NOCOLOR}")
        sys.exit(1)

    # Check for required commands
    for cmd in ["back_to_sequences", "blastn", "jq"]:
        if shutil.which(cmd) is None:
            print(f"{RED}Error: '{cmd}' could not be found. Please install it and ensure it's in your PATH.{NOCOLOR}")
            sys.exit(1)

    # Main loop
    
    TYPE = "unitig" if UNITIGS else "contig"
    process_accessions()

    print(f"\n{BLUE}================")
    print(f"{CYAN}>>> All done <<<")
    print(f"{BLUE}================\n")

    if not DELETE:
        print(f"{YELLOW}[INFO] You did not use --delete option. So you can manually remove all intermediate files (recruited {TYPE}s and {TYPE}s files) by running:{NOCOLOR}")
        print(f"rm -rf {MAIN_DIR_NAME}/{LOGAN_DIR_NAME}")

    os.chdir("..")
    print(f"{YELLOW}[INFO] Results can be found in directory {CYAN}{MAIN_DIR_NAME}{NOCOLOR}")

    if not UNITIGS and os.path.getsize(f"{MAIN_DIR_NAME}/{FAILED_ACCESSION_LIST}") > 0:
        nb_failed = len(open(f"{MAIN_DIR_NAME}/{FAILED_ACCESSION_LIST}").readlines())
        print(f"{YELLOW}[INFO] {nb_failed} accession{'s' if nb_failed > 1 else ''} failed to download contigs or had no recruited sequences.{NOCOLOR}")
        print(f"{YELLOW}[INFO] List of failed accessions: {CYAN}{MAIN_DIR_NAME}/{FAILED_ACCESSION_LIST}{NOCOLOR}")
        print(f"{YELLOW}[INFO] You can try to re-run the script with --unitigs option and this accession list.{NOCOLOR}")
        print(f"{YELLOW}[INFO] Command example:{NOCOLOR}")
        ## indicate how to re-run the script, with all the used arguments replacing the accession list by the FAILED_ACCESSION_LIST and using --unitigs
        ## if SESSION_ID is used, keep it, else use the FAILED_ACCESSION_LIST
        ## indicate all other used arguments (-l, -k, -d)
        if args.delete:
            delete_flag = "-d"
        else:
            delete_flag = ""
        
        print(f"logan_blaster -a {MAIN_DIR_NAME}/{FAILED_ACCESSION_LIST} -q {MAIN_DIR_NAME}/{QUERY_FILE} --unitigs -k {args.kmer_size} {delete_flag} {NOCOLOR}")  
        

if __name__ == "__main__":
    main()
