#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import shutil
import urllib.request
import urllib.error
from pathlib import Path

__version__ = '0.1.0'
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
    query_ACGT = get_query_ACGT(fasta_file)
    query_name, query_length, matched_positions = parse_blastn(blastn_file)
    assert len(query_ACGT) == query_length, f"Error, query given in {fasta_file} is of length {len(query_ACGT)}, while the blastn result indicates a sequence of length {query_length}"
    visualize_matches(query_ACGT, query_name, query_length, matched_positions, print_abundance=abundance)

# --- Using blast_parser ---

# Couleurs pour les messages
GREEN = "\033[0;32m"
RED = "\033[0;31m"
YELLOW = "\033[0;33m"
BLUE = "\033[1;34m"
CYAN = "\033[1;36m"
NOCOLOR = "\033[0m"

# Global variables
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
cli_installed = shutil.which("aws") is not None
failed_accession_list = ""


def run_blast(query_fasta, target_fasta):
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
        os.remove("error.log")
        sys.exit(1)
    os.remove("error.log")

    # Run the integrated blast_parser on blast results
    print(f"{YELLOW}[INFO] Synthesize blast results{NOCOLOR}")
    blastn_file = os.path.join(ALIGNEMENT_DIR_NAME, output_name)
    synth_file = os.path.join(ALIGNEMENT_DIR_NAME, f"synth_{output_name}")
    with open(synth_file, "w") as f:
        sys.stdout = f
        run_blast_parser(QUERY_FILE, blastn_file, abundance=True)
        sys.stdout = sys.__stdout__

def download_file(url, destination):
    """
    Download a file from a URL and save it to a local destination.
    
    This function downloads a file from the specified URL using urllib and saves it
    to the given destination path. The download is performed in chunks of 8192 bytes
    to handle large files efficiently.
    
        url (str): The URL of the file to download.
        destination (str): The local file path where the downloaded file will be saved.
    
    Raises:
        SystemExit: Exits with code 1 if the download fails due to HTTP errors,
                   URL errors, or other exceptions.
    
    Returns:
        None
    
    Examples:
        >>> download_file("https://example.com/file.zip", "./downloads/file.zip")
    """
    try:
        with urllib.request.urlopen(url) as response:
            with open(destination, "wb") as f:
                # Read and write in chunks
                while True:
                    chunk = response.read(8192)
                    if not chunk:
                        break
                    f.write(chunk)
    except urllib.error.HTTPError as e:
        print(f"{RED}Error: Failed to download {url} (HTTP {e.code}).{NOCOLOR}")
        sys.exit(1)
    except urllib.error.URLError as e:
        print(f"{RED}Error: Failed to download {url}.{NOCOLOR}")
        print(e.reason)
        sys.exit(1)
    except Exception as e:
        print(f"{RED}Error: Failed to download {url}.{NOCOLOR}")
        print(e)
        sys.exit(1)
        
# def download_file(url, destination):
#     try:
#         response = requests.get(url, stream=True)
#         response.raise_for_status()
#         with open(destination, "wb") as f:
#             for chunk in response.iter_content(chunk_size=8192):
#                 f.write(chunk)
#     except requests.RequestException as e:
#         print(f"{RED}Error: Failed to download {url}.{NOCOLOR}")
#         print(e)
#         sys.exit(1)

def main():
    global SESSION_ID, ACCESSION_FILE, QUERY_FILE, DELETE, UNITIGS, LIMIT, KMER_SIZE, MAIN_DIR_NAME, failed_accession_list

    parser = argparse.ArgumentParser(description="Process Logan session or accession/query files.")
    parser.add_argument("-s", "--session", help="Logan session ID")
    parser.add_argument("-a", "--accessions", help="Path to accessions.txt file")
    parser.add_argument("-q", "--query", help="Path to query fasta file")
    parser.add_argument("-u", "--unitigs", action="store_true", help="Use unitigs instead of contigs")
    parser.add_argument("-k", "--kmer-size", type=int, default=17, help="K-mer size for sequence recruitment")
    parser.add_argument("-l", "--limit", type=int, default=0, help="Limit number of accessions to process")
    parser.add_argument("-d", "--delete", action="store_true", help="Delete intermediate files after processing")
    args = parser.parse_args()

    SESSION_ID = args.session
    ACCESSION_FILE = args.accessions
    QUERY_FILE = args.query
    DELETE = args.delete
    UNITIGS = args.unitigs
    KMER_SIZE = args.kmer_size
    LIMIT = args.limit

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
    if SESSION_ID:
        MAIN_DIR_NAME = f"session_{SESSION_ID}"
        os.makedirs(MAIN_DIR_NAME, exist_ok=True)
        os.chdir(MAIN_DIR_NAME)
        if not UNITIGS:
            failed_accession_list = "failed_accessions.txt"
            Path(failed_accession_list).touch()
    else:
        query_basename = os.path.basename(QUERY_FILE).split(".")[0]
        # get the full absolute path of QUERY_FILE and ACCESSION_FILE
        ABS_QUERY_FILE = os.path.abspath(QUERY_FILE)
        ABS_ACCESSION_FILE = os.path.abspath(ACCESSION_FILE)
        for i in range(1, 1001):
            MAIN_DIR_NAME = f"{query_basename}_{i}"
            if not os.path.exists(MAIN_DIR_NAME):
                os.makedirs(MAIN_DIR_NAME)
                os.chdir(MAIN_DIR_NAME)
                os.makedirs(INPUT_DATA_DIR_NAME)
                shutil.copy(f"{ABS_QUERY_FILE}", INPUT_DATA_DIR_NAME)
                shutil.copy(f"{ABS_ACCESSION_FILE}", INPUT_DATA_DIR_NAME)
                QUERY_FILE = os.path.join(INPUT_DATA_DIR_NAME, os.path.basename(QUERY_FILE))
                ACCESSION_FILE = os.path.join(INPUT_DATA_DIR_NAME, os.path.basename(ACCESSION_FILE))
                if not UNITIGS:
                    failed_accession_list = "failed_accessions.txt"
                    Path(failed_accession_list).touch()
                break

    os.makedirs(LOGAN_DIR_NAME, exist_ok=True)
    os.makedirs(ALIGNEMENT_DIR_NAME, exist_ok=True)

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
    type_ = "unitig" if UNITIGS else "contig"
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

            local_file = os.path.join(LOGAN_DIR_NAME, f"{accession}.{type_}s.fa.zst")
            if not os.path.exists(local_file):
                print(f"{YELLOW}[INFO] Downloading {accession}.{type_}s.fa.zst...{NOCOLOR}")
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
                        print(f"{YELLOW}[WARNING] Attempt {attempt + 1} download failed for {accession}.{type_}s.fa.zst. {NOCOLOR}")
                    
                if not os.path.exists(f"{accession}.{type_}s.fa.zst"):
                    print(f"{RED}Error: Failed to download {accession}.{type_}s.fa.zst after 3 attempts.{NOCOLOR}")
                    if not UNITIGS:
                        with open(failed_accession_list, "a") as f:
                            f.write(f"{accession}\n")
                    continue      
                shutil.move(f"{accession}.{type_}s.fa.zst", LOGAN_DIR_NAME)
            else:
                print(f"{YELLOW}[INFO] Using existing local version of {local_file}...{NOCOLOR}")

            recruited_file = os.path.join(LOGAN_DIR_NAME, f"{accession}.recruited_{type_}s.fa")
            print(f"{YELLOW}[INFO] Recruiting sequences from {accession}.{type_}s.fa.zst with a match with {QUERY_FILE}...{NOCOLOR}")
            cmd_recruit = [
                "back_to_sequences",
                "--kmer-size", str(KMER_SIZE),
                "--in-kmers", QUERY_FILE,
                "--in-sequences", local_file,
                "--out-sequences", recruited_file
            ]
            print(f"{GREEN}Running command: {' '.join(cmd_recruit)}{NOCOLOR}")
            try:
                subprocess.run(cmd_recruit, check=True, stdout=subprocess.DEVNULL, stderr=open("error.log", "w"))
            except subprocess.CalledProcessError:
                print(f"{RED}Error: back_to_sequences failed for accession {accession}.{NOCOLOR}")
                with open("error.log", "r") as f:
                    print(f.read())
                os.remove("error.log")
                sys.exit(1)
            os.remove("error.log")

            if os.path.getsize(recruited_file) == 0:
                print(f"{YELLOW}[INFO]\tNo sequences were recruited from {accession}.{type_}s.fa.zst. Skipping BLAST step.{NOCOLOR}")
                if DELETE:
                    print(f"{YELLOW}[INFO] Deleting {recruited_file} and {local_file}...{NOCOLOR}")
                    os.remove(recruited_file)
                    os.remove(local_file)
                with open(failed_accession_list, "a") as f:
                    f.write(f"{accession}\n")
                continue

            print(f"{YELLOW}[INFO] Aligning recruited sequences from {accession}.{type_}s.fa.zst with {QUERY_FILE}...{NOCOLOR}")
            run_blast(QUERY_FILE, recruited_file)

            if DELETE:
                print(f"{YELLOW}[INFO] Deleting {recruited_file} and {local_file}...{NOCOLOR}")
                os.remove(recruited_file)
                os.remove(local_file)

    print(f"\n{BLUE}================")
    print(f"{CYAN}>>> All done <<<")
    print(f"{BLUE}================\n")

    if not DELETE:
        print(f"{YELLOW}[INFO] You did not use --delete option. So you can manually remove all intermediate files (recruited {type_}s and {type_}s files) by running:{NOCOLOR}")
        print(f"rm -rf {MAIN_DIR_NAME}/{LOGAN_DIR_NAME}")

    os.chdir("..")
    print(f"{YELLOW}[INFO] Results can be found in directory {CYAN}{MAIN_DIR_NAME}{NOCOLOR}")

    if not UNITIGS and os.path.getsize(f"{MAIN_DIR_NAME}/{failed_accession_list}") > 0:
        print(f"{YELLOW}[INFO] Some accessions failed to download contigs or had no recruited sequences.{NOCOLOR}")
        print(f"{YELLOW}[INFO] List of failed accessions: {CYAN}{MAIN_DIR_NAME}/{failed_accession_list}{NOCOLOR}")
        print(f"{YELLOW}[INFO] You can try to re-run the script with --unitigs option and this accession list.{NOCOLOR}")
        print(f"{YELLOW}[INFO] Command example:{NOCOLOR}")
        ## indicate how to re-run the script, with all the used arguments replacing the accession list by the failed_accession_list and using --unitigs
        ## if SESSION_ID is used, keep it, else use the failed_accession_list
        ## indicate all other used arguments (-l, -k, -d)
        if args.delete:
            delete_flag = "-d"
        else:
            delete_flag = ""
        
        print(f"logan_blaster -a {MAIN_DIR_NAME}/{failed_accession_list} -q {MAIN_DIR_NAME}/{QUERY_FILE} --unitigs -k {args.kmer_size} {delete_flag} {NOCOLOR}")  
        

if __name__ == "__main__":
    main()
