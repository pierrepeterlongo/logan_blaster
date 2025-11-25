import argparse

def get_query_ACGT(file_path):
    """Returns the first sequence from a fasta file. Possibly multiline

    Args:
        file_path (String): fasta file

    Returns:
        String: ACGT sequence
    """
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
        if v == 0: print('-', end='')
        else:
            # 1 - 26 -> a-z
            if v <= 26:
                print(chr(ord('a') + v - 1), end='')
            if v > 26: 
                print('Z', end='')
                
    else: 
        if v == 0: print('-', end='')
        else: print('|', end='')

def print_spaces(n):
    for _ in range(n):
        print(' ', end='')

def visualize_matches(query_ACGT, query_name, query_length, matched_positions, print_abundance=False):
    nb_chars_before_line = 9 #Query  123  (here 9 is for "query" + 4)
    # count the number of chars needed to print 123
    nb_chars_for_len = len(str(query_length))
    nb_chars_before_line += nb_chars_for_len
    len_line = 80
    print(f"Query: {query_name}")
    pos = 0
    while True: 
        diff = len(str(query_length))-len(str(pos+1))
        print("query  ", end='')
        print_spaces(diff)
        print("{pos}  ".format(pos=pos+1), end='')
        if pos + len_line < query_length:
            print(query_ACGT[pos: pos+len_line])
            print_spaces(nb_chars_before_line)
            for v in matched_positions[pos: pos+len_line]:
                print_value(v, print_abundance)
            print("\n") # two new lines
            pos += len_line
        else:
            print(query_ACGT[pos:])
            print_spaces(nb_chars_before_line)
            for v in matched_positions[pos:]:
                print_value(v, print_abundance)
            print("\n") # two new lines
            break


def main():
    """Analyses a blastn output
    """
    parser = argparse.ArgumentParser(description="Pretty print blastn results. \
        Given a query fasta file and a blastn output file, prints the query sequence \
        and indicates which positions were matched in the blastn results.\n\
        If --abundance is given, indicates how many times each position was matched.\
        Otherwise, just indicates whether a position was matched or not.\n\
        Without --abundance, matched positions are indicated with '|', unmatched with '-'.\n\
        With --abundance, matched positions are indicated with letters a-z (1-26 matches), \
        and 'Z' for more than 26 matches. Unmatched positions are indicated with '-'.")
    parser.add_argument("--fasta_file", type=str, required=True, help="Input query file")
    parser.add_argument("--blastn_file", type=str, required=True, help="Input blastn file")
    parser.add_argument("--abundance", action="store_true", help="Prints how much time each position was matched", default=False)
    args = parser.parse_args()
    
    query_fasta_file = args.fasta_file
    blastn_file = args.blastn_file
    
    query_ACGT = get_query_ACGT(query_fasta_file)
    query_name, query_length, matched_positions = parse_blastn(blastn_file)
    assert len(query_ACGT) == query_length, f"Error, query given in {query_fasta_file} \
        is of length {len(query_ACGT)}, while the blastn result indicates a sequence \
            of length {query_length}"
    visualize_matches(query_ACGT, query_name, query_length, matched_positions, print_abundance=args.abundance)
    
if __name__ == "__main__":
    main()