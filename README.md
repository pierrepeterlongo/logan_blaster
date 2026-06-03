# logan_blaster: Blast a query against Logan data

[![GitHub release](https://img.shields.io/github/v/tag/pierrepeterlongo/logan_blaster?label=version)](https://github.com/pierrepeterlongo/logan_blaster/releases)

Align genomic sequences with [Logan](https://github.com/IndexThePlanet/Logan/) contigs or unitigs.

## Two main modes
1. From a query and a list of Logan accessions.  
2. From a [Logan-Search](https://logan-search.org/) session id. The query and Logan accessions are then automaticaly retreived.

In any case, for each accession, `logan_blaster` 
1. Downloads the Logan contigs,
2. Computes coverage statistics (mean, median) of the downloaded contigs,
3. Recruits contigs that contain at least one shared k-mer (k=17 by default) with the query (uses [`back_to_sequences`](https://github.com/pierrepeterlongo/back_to_sequences)),
4. Computes coverage statistics (mean, median) of the recruited contigs (uses [`count_logan_tig_coverage`](https://github.com/pierrepeterlongo/count_logan_tig_coverage)),
5. Runs a local blast between the query and this subset of contigs.
6. Analyses the blast results: prints the portion(s) of the query matched by at least one contig of the accession

## Install

### From conda / mamba

```bash
# Clone logan_blaster repository
git clone https://github.com/pierrepeterlongo/logan_blaster
cd logan_blaster

# Create the environment with dependencies
mamba env create -f environment.yml

# Activate the environment
mamba activate logan_blaster

# Install back_to_sequences (Rust tool) if not already installed
cargo install back_to_sequences

# Install count_tig_coverage (Rust tool)
cargo install --git https://github.com/pierrepeterlongo/count_logan_tig_coverage
```

### From sources

#### Requires 

- blast: *on mac* `brew install blast` or look at [blast installation web page](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- back_to_sequences: see back_to_sequences installation [web page](https://b2s-doc.readthedocs.io/en/latest/usage.html#installation])
- jq: see [jq installation web page](https://jqlang.org/)
- count_tig_coverage: `cargo install --git https://github.com/pierrepeterlongo/count_logan_tig_coverage`

#### Clone the repository

```bash
git clone https://github.com/pierrepeterlongo/logan_blaster
cd logan_blaster

# Use the py script: 
python logan_blaster.py --help
```

## Usage

```bash
logan_blaster -h
usage: logan_blaster [-h] [-s SESSION] [-a ACCESSIONS] [-q QUERY] [-u] [-k KMER_SIZE] [-l LIMIT] [-d]

Process Logan session or accession/query files.

options:
  -h, --help            show this help message and exit
  -s, --session SESSION
                        Logan session ID
  -a, --accessions ACCESSIONS
                        Path to accessions.txt file or a .csv file 
                        (containing a first header line, ignored, 
                        and storing accessions as the first column)
  -q, --query QUERY     Path to query fasta file
  -u, --unitigs         Use unitigs instead of contigs
  -k, --kmer-size KMER_SIZE
                        K-mer size for sequence recruitment
  -l, --limit LIMIT     Limit number of accessions to process
  -d, --delete          Delete intermediate files after processing
```

## Examples

### Example running from session.

```bash
logan_blaster -s kmviz-b2bce461-ca13-4a45-b0b4-6c894eacf103
```

### Example running from accessions and query files.

This usage enables to select specific accessions to process, also ordering them, and to provide any custom query file.

```bash
logan_blaster  -a example/accessions.txt -q example/query.fa
```

Note, the `-a` option also accepts .csv files directly downloaded from the logan-search interface

```bash
logan_blaster  -a my_data.csv -q example/query.fa
```

## Output

### Created files and directories

The program creates a directory named either with the session name (if run from a session id) or with the query file name (if run from accessions and query files).
Here this the structure of this directory (here with two successful accessions on which the query was aligned):

```
.
|-- alignments
|   |-- my_query_vs_SRR1608527.txt
|   |-- my_query_vs_SRR1608810.txt
|   |-- synth_my_query_vs_SRR1608527.txt
|   |-- synth_my_query_vs_SRR1608810.txt
|-- failed_accessions.txt
|-- input_data
|   |-- num_accession.txt
|   `-- seq.fa
`-- logan_data
    |-- SRR1608527.contigs.fa.zst
    |-- SRR1608527.recruited_contigs.fa
    |-- SRR1608810.contigs.fa.zst
    |-- SRR1608810.recruited_contigs.fa
```

### Interpretation of the results

- The `failed_accessions.txt` file contains the list of accessions on which the query was not aligned (or not existing on logan data).
- In the `input_data` directory, 
  - `seq.fa` is the query fasta file,
  - `num_accession.txt` is the accessions file.
- In the `logan_data` directory, 
  - files named `<ACCESSION>.contigs.fa.zst` are the downloaded Logan contigs,
  - files named `<ACCESSION>.recruited_contigs.fa` are the contigs that were recruited because they share at least one k-mer with the query (found thanks to back_to_sequences).
- In the `alignments` directory, 
  - files named `my_query_vs_<ACCESSION>.txt` contain the blast alignments between the query and the recruited contigs from accession `<ACCESSION>`.
  - files named `synth_my_query_vs_<ACCESSION>.txt` contain a synthesis of these alignments, indicating for each position of the query how many contigs aligned to it (see below).

#### Synthesis of the alignments

In the `alignments` directory, files named `synth*` contain this piece of information: 

```bash
query    1  CTTCCCTCTAGAACGGGACGAGGTGATGCCCCCACCCATCTCGCACCCCCCATGTGAAAGGCTCACATTCCCGAAGGGGC
            ---aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbb

query   81  TTCCCTGGGCTCCGAAGGTCAGGGAGAAGGATATTGAGATGTTCCTTGAAAACAGTCGCAGCAAATTCATTGGCTACACG
            aaaaaaaaaaaabbbbbbbbbcccccccccbbbbbcdddccccccdddddddccbbaaaaaaaaaaaaaaaaaaaaaaaa
```

In this case, except for the first three nucleotides, the full query was aligned to at least one contig. Some portions were aligned to two contigs (positions labeled 'b'), or three contigs (positions labeled 'c'), and so on.
All positions aligned to more than 26 contigs are labeled 'Z'.


## Tests

The test suite uses [pytest](https://pytest.org) and is organised in two files:

| File | Content | External tools required |
|---|---|---|
| `tests/test_blast_parser.py` | Unit tests for all blast-parser functions | none |
| `tests/test_integration.py` | Pipeline integration tests with local and remote data | `blastn`, `back_to_sequences`, `zstd` |

### Running the tests

```bash
# Unit tests + local integration tests (no network, ~3 s)
pytest

# Include network tests (download one Logan accession)
pytest --network
```

### Test categories

**Unit tests** (`test_blast_parser.py`) — always fast, no external tools:
- `TestGetQueryACGT` — FASTA reading (single sequence, multiline, multi-record)
- `TestGetQueryName` / `TestGetQueryLength` — blast output header parsing
- `TestParseBLASTN` — position-coverage vector: spot-checks on known overlaps (single, double, triple coverage computed from `tests/data/self_blast.txt`)
- `TestRunBlastParser` — byte-exact comparison of the full visualisation output against `tests/data/expected_self_synth.txt`

**Local integration tests** (`test_integration.py`, no network):
- `TestRunBlast` — calls `_run_blast()` with the query aligned against itself; verifies that the blastn and synth files are created and match the reference
- `TestFullPipelineLocal` — runs the complete `_process_accessions()` loop with a pre-placed `.fa.zst` file (the query compressed with `zstd`); checks file creation, synth content, and absence of failed accessions

**Network integration tests** (`test_integration.py`, `--network` flag required):
- `TestFullPipelineNetwork` — downloads the first accession from `example/accessions.txt` and verifies output structure and synth format

### Reference test data

```
tests/
├── data/
│   ├── self_blast.txt           blastn output: example/query.fa vs itself
│   └── expected_self_synth.txt  reference synth visualisation (byte-exact)
├── conftest.py                  shared fixtures and --network option
├── test_blast_parser.py
└── test_integration.py
```

`self_blast.txt` and `expected_self_synth.txt` are committed and serve as the non-regression baseline. Regenerate them if the query file or blastn parameters change:

```bash
# Regenerate self_blast.txt
blastn -query example/query.fa -subject example/query.fa \
       -out tests/data/self_blast.txt \
       -outfmt 0 -sorthits 0 -word_size 11 -gapextend 2 -gapopen 5 -reward 2 -penalty -3

# Regenerate expected_self_synth.txt
python3 -c "
import sys, io, contextlib
sys.path.insert(0, '.')
from logan_blaster import run_blast_parser
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    run_blast_parser('example/query.fa', 'tests/data/self_blast.txt', abundance=True)
open('tests/data/expected_self_synth.txt', 'w').write(buf.getvalue())
"
```

## Versioning

Versions follow [Semantic Versioning](https://semver.org/) and are driven by **git tags**.
The version displayed by `logan_blaster --version` is derived automatically from the latest tag.


## Authors

- [Pierre Peterlongo](https://people.rennes.inria.fr/Pierre.Peterlongo/)
- [Téo Lemane](https://tlemane.github.io/)
