# vclust-dev

vclust is an alignment-based tool for fast and accurate calculation of Average Nucleotide Identity (ANI) between complete or metagenomically-assembled viral genomes. The tool also performs ANI-based clustering of genomes according to standards recommended by international virus consortia, including *International Committee on Taxonomy of Viruses* (ICTV), *Minimum Information about an Uncultivated Virus Genome* (MIUViG), and 
*Integrated Microbial Genomes - Virus Resource* (IMG/VR). 

## Features

#### 1. Accurate ANI calculations

vclust uses a Lempel-Ziv-based pairwise sequence [aligner](https://github.com/refresh-bio/LZ-ANI) for ANI calculation. It is magnitudes faster than BLAST-based tools and equally accurate as the most sensitive BLASTn searches.

#### 2. Multiple similarity measures

vclust offers multiple similarity measures between two genome sequences, whereas other tools typically provide only one or two.
- **ANI**: number of identical bases across local alignments divided by the total length of the alignments.
- **Global ANI (gANI)**: number of identical bases across local alignments divided by the length of the query/target genome.
- **Total ANI (tANI)**: number of identical bases between query-target and target-query genomes divided by the sum length of both genomes. tANI is equivalent to VIRIDIC's intergenomic similarity.
- **Coverage (alignment fraction)**: proportion of query sequence that is aligned with target sequence.
- Number of alignments
- Length ratio

#### 3. Multiple clustering algorithms

vclust offers support for six clustering algorithms tailored to various scenarios, including taxonomic classification and dereplication of viral genomes.
- Single-linkage
- Complete-linkage
- UCLUST
- CD-HIT (Greedy incremental)
- Greedy set cover (adopted from MMseqs2)
- Leiden algorithm

#### 4. Speed and efficiency

vclust uses three efficient C++ tools - [Kmer-db](https://github.com/refresh-bio/kmer-db), [LZ-ANI](https://github.com/refresh-bio/LZ-ANI), [rapid-cluster](https://github.com/refresh-bio/rapid-cluster) - for prefiltering, aligning, calculating ANI, and clustering viral genomes. This combination enables the processing of millions of virus genomes within a few hours on a mid-range workstation.

#### 5. Web service

For datasets of fewer than 1000 viral genomes, vclust is available [on-line](www.google.pl).


## Installation

To install vclust you can compile dependencies from source or download statically compiled binaries.

### Option 1: Compile from source

Clone this repository recursively with submodules and compile it.

```bash
git clone --recurse-submodules https://github.com/refresh-bio/vclust-dev
cd vclust-dev
make -j
```

## Quick usage

1. **Prefilter** similar genome sequence pairs before conducting pairwise alignments.

```bash
./vclust.py prefilter -i example/multifasta.fna -o fltr.txt
```

2. **Align** similar genome sequence pairs and calculate pairwise ANI measures.

```bash
./vclust.py align -i example/multifasta.fna -o ani.tsv --filter fltr.txt
```

3. **Cluster** genome sequences based on given ANI measure and minimum threshold.

```bash
./vclust.py cluster -i ani.tsv -o clusters.tsv --ids ani.data.tsv --metric ani --ani 0.95
```

## Getting started

vclust comprises three commands: `prefilter`, `align`, and `cluster`. Calls to these commands follow the structure:

```
./vclust.py command -i <input> -o <output> [options]
```

### Prefilter

> For small datasets, using prefilter is optional. However, note that without it, vclust will perform all-versus-all pairwise sequence alignments.

The `prefilter` command creates a pre-alignment filter, which eliminates dissimilar genome pairs before calculating pairwise alignments. This process reduces the number of potential genome pairs to only those with sufficient *k*-mer-based sequence similarity. The *k*-mer-based sequence similarity between two genomes is controlled by two options: minimum number of common *k*-mers (`--min-kmers`) and minimum sequence identity of the shorter sequence (`--min-ident`). Both *k*-mer-based similarities are computed using [Kmer-db](https://github.com/refresh-bio/kmer-db) which evaluates the entire set of *k*-mers, overcoming the sampling constraints typical of sketching methods like fastANI and Mash. In addittion, the use of sparse distance matrices enables memory-efficient processing of millions of genomes.

```bash
# Prefilter creates a pre-alignment filter with genome sequence pairs that have
# at least 30 common k-mers and 50% of identity over the shorter sequence.
./vclust.py prefilter -i genomes.fna -o fltr.txt --min-kmers 30 --min-ident 0.5 
```

Large sets of genome sequences can be processed in smaller, equally-sized, batches of sequences to reduce memory consumption. 

```bash
# Prefilter works with batches of 5 million sequences each.
./vclust.py prefilter -i genomes.fna -o fltr.txt --batch-size 5000000
```

### Align

The `align` command conducts parwise sequence alignments among viral genomes and rovides similarity measures like ANI and coverage. For this purpose, vclust uses the sensitive and efficent [LZ-ANI](https://github.com/refresh-bio/LZ-ANI) aligner, which is based on the Lempel-Ziv parsing.

```bash
# Align genome pairs filtered by the prefiltering command will be aligned.
./vclust.py align -i genomes.fna -o ani.tsv --filter fltr.txt
```

The `align` command allows for further filtering based on filter's minimum threshold.

```bash
# Align only the genome pairs with sequence identity greater than or equal to 80%.
./vclust.py align -i genomes.fna -o ani.tsv --filter fltr.txt --filter-threshold 0.8
```

Without `--filter` vclust aligns all possible genome pairs.

```bash
./vclust.py align -i genomes.fna -o ani.tsv
```

#### Align output

The `align` command creates two TSV files: one contains ANI values for genome pairs, and the other lists genome identifiers sorted by decreasing sequence length. Both TSV files are used as an input for vclust's clustering. For example, the following command will create two TSV files: [ani.tsv](./example/output/ani.tsv) and [ani.ids.tsv](./example/output/ani.ids.tsv):

```bash
./vclust.py align -i genomes.fna -o ani.tsv --filter fltr.txt
```

For the TSV file with ANI values, vclust has three output formats that differ in the number of columns (`standard`, `lite`, `complete`).

```bash
./vclust.py align -i genomes.fna -o ani.tsv --filter fltr.txt --outfmt lite
```

| Column | Standard | Lite | Complete | Description |
| --- | :---: |:---: | :---: | --- |
| idx1 | + | + | +  | index of sequence 1 |
| idx2 | + | + | +  | index of sequence 2 |
| id1 | + | - | +  | identifier (name) of sequence 1 |
| id2 | + | - | +  | identifier (name) of sequence 2 |
| tani | + | + | +  | total ANI [0-1] |
| gani | + | + | +  | global ANI [0-1] |
| ani | + | + | +  | ANI [0-1] |
| cov | + | + | +  | Coverage (alignment fraction) [0-1] |
| num_alns | + | + | +  | Number of alignments |
| len_ratio | + | + | +  | Length ratio between sequence 1 and sequence 2 |
| len1 | - | - | +  | Length of sequence 1 |
| len2 | - | - | +  | Length of sequence 2|
| nt_match | - | - | +  | Number of matching nucleotides across alignments |
| nt_mismatch | - | - | +  | Number of mismatching nucleotides across alignments |

#### Align output filtering

The `align` command enables filtering output by setting minimum thresholds for similarity measures. This allows reporting only genome pairs meeting the specified criteria, significantly reducing the output TSV file size:

```bash
# Output only genome pairs with ANI ≥ 0.75 and coverage ≥ 0.75.
./vclust.py align -i genomes.fna -o ani.tsv --filter fltr.txt --out-ani 0.75 --out-cov 0.75
```

### Cluster

The `cluster` command takes as input two TSV files (outputs of vclust align) and clusters viral genomes using a user-selected clustering algorithm and similarity measures. Internally, vclust uses [rapid-cluster](https://github.com/refresh-bio/rapid-cluster), enabling large-scale clustering of millions of genome sequences by leveraging sparse distance matrices.

vclustr supports six clustering algorithms (`--algorithm`):

1. `single`: Single-linkage clustering
2. `complete`: Complete-linkage clustering
3. `uclust`: UCLUST-like clustering
4. `cd-hit`: Greedy incremental clustering
5. `set-cover`: Set-Cover (greedy)
6. `leiden`: the Leiden algorithm

Clustering can use one of three similarity measures (`--metric`): `tani`, `gani`, or `ani`. The user-provided value of the selected similarity measure serves as the minimum similarity threshold for connecting genomes.

```bash
# Cluster genomes based on ANI and connect genomes if their ANI ≥ 95%.
./vclust cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm single \
--metric ani --ani 0.95
```

Additionally, vclust may use extra threshold values for various combinations of similarity measures (`--ani`, `--gani`, `--tani`, `--cov`, `--num_alns`). If a genome pair fails to meet the specified thresholds, it is excluded from clustering.

```bash
# Cluster genomes based on ANI similarity measure, but connect them only
# if ANI ≥ 95% and coverage ≥ 85%.
./vclust cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm single \
--metric ani --ani 0.95 --cov 0.85
```

#### Cluster output

The `cluster` command creates a TSV file listing genome identifiers followed by 0-based numberical cluster identifiers. The file contains identifiers of all genomes, both clustered and singleton ones. Genome identifiers are listed in the same order as in the ids file (sorted by decreasing sequence length).

```
object	cluster
NC_025457.alt2	3
NC_005091.alt2	0
NC_005091.alt1	0
NC_005091.ref	0
NC_002486.alt	1
NC_002486.ref	1
NC_025457.ref	4
NC_025457.alt1	5
NC_010807.alt2	2
NC_010807.alt3	2
NC_010807.alt1	2
NC_010807.ref	2
```

Alternatively, instead of numerical cluster identifiers, vclust can output a representative genome for each cluster using the `--out-repr` option. The representative genome within each cluster is determined as the one with the longest sequence.

```
object	cluster
NC_025457.alt2	NC_025457.alt2
NC_005091.alt2	NC_005091.alt2
NC_005091.alt1	NC_005091.alt2
NC_005091.ref	NC_005091.alt2
NC_002486.alt	NC_002486.alt
NC_002486.ref	NC_002486.alt
NC_025457.ref	NC_025457.ref
NC_025457.alt1	NC_025457.alt1
NC_010807.alt2	NC_010807.alt2
NC_010807.alt3	NC_010807.alt2
NC_010807.alt1	NC_010807.alt2
NC_010807.ref	NC_010807.alt2
```

## Examples

### Classify complete viral genomes into putative species and genera following ICTV recommendations

```bash
# Create a pre-alignment filter for genome pairs with a minimum of 10 common k-mers.
./vclust prefilter -i genomes.fna -o fltr.txt --min_kmers 10
```

```bash
# Calculate ANI measures for genome pairs specified in the filter.
./vclust align -i genomes -o ani.tsv --filter fltr.txt
```

```bash
# Assign viruses into putative species (tANI ≥ 95%).
./vclust cluster -i ani.tsv -o species.tsv --ids ani.ids.tsv --algorithm single \
--metric tani --tani 0.95

# Assign viruses into putative genus (tANI ≥ 70%).
./vclust cluster -i ani.tsv -o genus.tsv --ids ani.ids.tsv --algorithm single \
--metric tani --tani 0.70
```

### Assign metagenomic viral contigs into vOTUs following MIUViG standards

```bash
# Create a pre-alignment filter by processing batches of 5 million genomes.
./vclust prefilter -i genomes.fna -o fltr.txt --min_kmers 30 --min_ident 0.75 \
--batch-size 5000000
```

```bash
# Calculate ANI measures for genome pairs specified in the filter. Keep the output TSV
# file relatively small-sized: use lite output format and report only genome pairs 
# with ANI ≥ 90% and coverage ≥ 80%.
./vclust align -i genomes -o ani.tsv --filter fltr.txt --outfmt lite --ani 0.9 --cov 0.8
```

```bash
# Cluster contigs into vOTUs using the MIUVIG recommended-thresholds.
./vclust cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm single \
--metric ani --ani 0.95 --cov 0.85
```

### Dereplicate all IMG/VR contigs

```bash
# Create a pre-alignment filter by processing batches of 5 million genomes.
./vclust prefilter -i genomes.fna -o fltr.txt --min_kmers 30 --min_ident 0.75 \
--batch-size 5000000
```

```bash
# Calculate ANI measures for genome pairs specified in the filter. Keep the output TSV
# file relatively small-sized: use lite output format and report only genome pairs 
# with ANI ≥ 90% and coverage ≥ 80%.
./vclust align -i genomes.fna -o ani.tsv --filter fltr.txt --outfmt lite \
--ani 0.9 --cov 0.8
```

```bash
# Cluster contigs into vOTUs using the Leiden algorithm.
./vclust cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm leiden \
--leiden-resolution 0.7 --metric ani --ani 0.95 --cov 0.85
```

### Calculate pairwise similarities between all-versus-all genomes

```
./vclust.py align -i genomes.fna -o ani.tsv --outfmt complete
```



## Test

To ensure that vclust works as expected, you can run tests using [pytest](https://docs.pytest.org/).

```bash
pytest test.py
```

## Cite


