# logan_blaster: Blast a query against Logan data

[![GitHub release](https://img.shields.io/github/v/release/pierrepeterlongo/blast_logan_search_results?label=version)](https://github.com/pierrepeterlongo/blast_logan_search_results/releases)

Align genomic sequences with [Logan](https://github.com/IndexThePlanet/Logan/) contigs or unitigs.

## Two main modes
1. From a query and a list of Logan accessions.  
2. From a [Logan-Search](https://logan-search.org/) session id. The query and Logan accessions are then automaticaly retreived.

In any case, for each accession, `logan_blaster` 
1. Downloads the Logan contigs,
2. Recruits contigs that contain at least one shared k-mer (k=17 by default) with the query (uses `[back_to_sequences](https://github.com/pierrepeterlongo/back_to_sequences)`), 
3. Runs a local blast between the query and this subset of contigs.
4. Analyses the blast results: prints the portion(s) of the query matched by at least contig of the accession

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
```

### From sources

#### Requires 

- blast: *on mac* `brew install blast` or look at [blast installation web page](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- back_to_sequences: see back_to_sequences installation [web page](https://b2s-doc.readthedocs.io/en/latest/usage.html#installation])
- jq: see [jq installation web page](https://jqlang.org/)

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


## Versioning

Versions follow [Semantic Versioning](https://semver.org/) and are driven by **git tags**.
The version displayed by `logan_blaster --version` is derived automatically from the latest tag.

### Creating a new release

```bash
# 1. Commit all changes
git commit -am "Release vX.Y.Z"

# 2. Create an annotated tag
git tag -a vX.Y.Z -m "Release vX.Y.Z"

# 3. Push the tag
git push origin vX.Y.Z
```

The tag is the single source of truth — no need to edit the version manually anywhere in the code.

## Authors

- [Pierre Peterlongo](https://people.rennes.inria.fr/Pierre.Peterlongo/)
- [Téo Lemane](https://tlemane.github.io/)
