# <img src="./images/logo.svg" alt="Vclust logo" /> vclust-dev

Vclust is an alignment-based tool for fast and accurate calculation of Average Nucleotide Identity (ANI) between complete or metagenomically-assembled viral genomes. The tool also performs ANI-based clustering of genomes according to standards recommended by international virus consortia, including *International Committee on Taxonomy of Viruses* (ICTV) and *Minimum Information about an Uncultivated Virus Genome* (MIUViG). 


## Table of contents

1. [Features](#1-features)
2. [Installation](#2-installation)
3. [Quick Start](#3-quick-start)
4. [Input data](#4-input-data)
5. [Usage](#5-usage)
   1. [Prefilter](#51-prefilter)
   2. [Align](#52-align)
      * [Align output](#align-output)
      * [Align output filtering](#align-output-filtering)
   3. [Cluster](#53-cluster)
      * [Cluster output](#cluster-output)
6. [Examples](#examples)
   1. [Classify viruses into species and genera using the ICTV standards](#61-classify-viruses-into-species-and-genera-using-the-ictv-standards)
   2. [Assign viral contigs into vOTUs using the MIUViG standards](#62-assign-viral-contigs-into-votus-using-the-miuvig-standards)
   3. [Dereplicate genomes](#63-dereplicate-genomes)
   4. [Calculate pairwise similarities between all-versus-all genomes](#64-calculate-pairwise-similarities-between-all-versus-all-genomes)
   5. [Process large datasets](#65-process-large-datasets)
7. [Tests](#7-test)
8. [Cite](#8-cite)
9. [License](#9-license)


## 1. Features

#### :gem: Accurate ANI calculations

Vclust uses a Lempel-Ziv-based pairwise sequence aligner ([LZ-ANI](https://github.com/refresh-bio/LZ-ANI)) for ANI calculation. LZ-ANI achieves high sensitivity in detecting matched and mismatched nucleotides, ensuring accurate ANI determination. Its efficiency comes from a simplified indel handling model, making LZ-ANI  magnitudes faster than alignment-based tools (e.g., BLASTn, MegaBLAST) while maintaining comparable accuracy to the most sensitive BLASTn searches.

#### :triangular_ruler: Multiple similarity measures

Vclust offers multiple similarity measures between two genome sequences, whereas other tools typically provide only one or two.
- **ANI**: The number of identical nucleotides across local alignments divided by the total length of the alignments.
- **Global ANI (gANI)**: The number of identical nucleotides across local alignments divided by the length of the query/target genome.
- **Total ANI (tANI)**: The number of identical nucleotides between query-target and target-query genomes divided by the sum length of both genomes. tANI is equivalent to the VIRIDIC's intergenomic similarity.
- **Coverage (alignment fraction)**: The proportion of the query sequence aligned with the target sequence.
- **Number of local alignments**: The count of individual alignments found between the sequences.
- **Ratio between query and target genome lengths**: A measure comparing the lengths of the two genomes.

#### :star2: Multiple clustering algorithms 

Vclust provides six clustering algorithms tailored to various scenarios, including taxonomic classification and dereplication of viral genomes.
- Single-linkage
- Complete-linkage
- UCLUST
- CD-HIT (Greedy incremental)
- Greedy set cover (adopted from MMseqs2)
- Leiden algorithm

#### :fire: Speed and efficiency 

Vclust uses three efficient C++ tools - [Kmer-db](https://github.com/refresh-bio/kmer-db), [LZ-ANI](https://github.com/refresh-bio/LZ-ANI), [Clusty](https://github.com/refresh-bio/clusty) - for prefiltering, aligning, calculating ANI, and clustering viral genomes. This combination enables the processing of millions of virus genomes within a few hours on a mid-range workstation.

#### :earth_americas: Web service

For datasets containing up to 1000 viral genomes, Vclust is available at [http://www.vclust.org](http://www.vclust.org).


## 2. Installation

To install Vclust you can compile dependencies from source or download statically compiled binaries.

### Option 1: Compile from source

Clone this repository recursively with submodules and compile it.

```bash
git clone --recurse-submodules https://github.com/refresh-bio/vclust-dev
cd vclust-dev
make -j
```

### Option 2: Download precompiled binaries

TODO

## 3. Quick start

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

## 4. Input data

Vclust accepts a single FASTA file containing viral genomic sequences ([example](./example/multifasta.fna)) or a directory of FASTA files (one genome per file) ([example](./example/fna/)). The input file(s) can be gzipped.

## 5. Usage

Vclust provides three commands: `prefilter`, `align`, and `cluster`. Calls to these commands follow the structure:

```
./vclust.py command -i <input> -o <output> [options]
```


### 5.1. Prefilter

The `prefilter` command creates a pre-alignment filter, which eliminates dissimilar genome pairs before calculating pairwise alignments. This process reduces the number of potential genome pairs to only those with sufficient *k*-mer-based sequence similarity. The *k*-mer-based sequence similarity between two genomes is controlled by two options: minimum number of common *k*-mers (`--min-kmers`) and minimum sequence identity of the shorter sequence (`--min-ident`). Both *k*-mer-based similarities are computed using [Kmer-db 2](https://github.com/refresh-bio/kmer-db) which evaluates the entire set of *k*-mers, overcoming the sampling constraints typical of sketching methods like FastANI and Mash. In addittion, the use of sparse distance matrices enables memory-efficient processing of millions of genomes.

```bash
# Create a pre-alignment filter with genome sequence pairs that have
# at least 30 common 25-mers and 70% of identity over the shorter sequence.
./vclust.py prefilter -i genomes.fna -o fltr.txt --k 25 --min-kmers 30 --min-ident 0.7 
```

To reduce memory consumption, large sets of genome sequences can be processed in smaller, equally-sized, batches of sequences.

```bash
# Process genomes in batches of 5 million sequences each.
./vclust.py prefilter -i genomes.fna -o fltr.txt --batch-size 5000000
```

### 5.2. Align

The `align` command conducts pairwise sequence alignments among viral genomes and provides similarity measures like ANI and coverage (alignment fraction). For this purpose, Vclust uses the fast and efficient [LZ-ANI](https://github.com/refresh-bio/LZ-ANI) aligner.

```bash
# Align genome pairs filtered by the prefiltering command.
./vclust.py align -i genomes.fna -o ani.tsv --filter fltr.txt
```

The `align` command allows for further filtering based on filter's minimum threshold.

```bash
# Align only the genome pairs with sequence identity greater than or equal to 80%.
./vclust.py align -i genomes.fna -o ani.tsv --filter fltr.txt --filter-threshold 0.8
```

Without `--filter` Vclust aligns all possible genome pairs.

```bash
./vclust.py align -i genomes.fna -o ani.tsv
```

#### Align output

The `align` command creates two TSV files: one with ANI values for genome pairs and another with genome identifiers sorted by decreasing sequence length. Both TSV files can be used as input for Vclust's clustering. 

The following command will create two TSV files: [ani.tsv](./example/output/ani.tsv) and [ani.ids.tsv](./example/output/ani.ids.tsv):

```bash
./vclust.py align -i genomes.fna -o ani.tsv --filter fltr.txt
```

Sample output:

```
idx1 idx2 id1             id2             tani      gani      ani       cov num_alns len_rat
8    9    NC_010807.alt2  NC_010807.alt3  0.969793  0.980939  0.986608  0.994255  51  1.0166
9    8    NC_010807.alt3  NC_010807.alt2  0.969793  0.958462  0.985184  0.972876  48  0.9836
8   10    NC_010807.alt2  NC_010807.alt1  0.986229  0.986266  0.986363  0.999901  29  1.0448
10   8    NC_010807.alt1  NC_010807.alt2  0.986229  0.986191  0.986216  0.999974  28  0.9570
8   11    NC_010807.alt2  NC_010807.ref   0.989467  0.989693  0.989693  1.000000   8  1.0448
11   8    NC_010807.ref   NC_010807.alt2  0.989467  0.989231  0.989231  1.000000   8  0.9570
9   10    NC_010807.alt3  NC_010807.alt1  0.971502  0.958863  0.992810  0.965807  58  1.0277
10   9    NC_010807.alt1  NC_010807.alt3  0.971502  0.984491  0.993113  0.991318  56  0.9730
9   11    NC_010807.alt3  NC_010807.ref   0.982555  0.968865  0.995800  0.972951  44  1.0277
11   9    NC_010807.ref   NC_010807.alt3  0.982555  0.996625  0.996676  0.999948  44  0.9730
10  11    NC_010807.alt1  NC_010807.ref   0.996664  0.996548  0.996548  1.000000  22  1.0000
11  10    NC_010807.ref   NC_010807.alt1  0.996664  0.996780  0.996780  1.000000  22  1.0000
```

To manage the output TSV file size, Vclust offers three formats: `standard`, `lite`, `complete`.

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

The `align` command enables filtering output by setting minimum thresholds for similarity measures. This ensures only genome pairs meeting the thresholds are reported, significantly reducing the output TSV file size.

```bash
# Output only genome pairs with ANI ≥ 0.75 and coverage ≥ 0.75.
./vclust.py align -i genomes.fna -o ani.tsv --filter fltr.txt --out-ani 0.75 --out-cov 0.75
```

### 5.3. Cluster

The `cluster` command takes as input two TSV files (outputs of `vclust align`) and clusters viral genomes using a user-selected clustering algorithm and similarity measures. Internally, Vclust uses [Clusty](https://github.com/refresh-bio/clusty), enabling large-scale clustering of millions of genome sequences by leveraging sparse distance matrices.

Vclust supports six clustering algorithms (`--algorithm`):

1. `single`: Single-linkage clustering
2. `complete`: Complete-linkage clustering
3. `uclust`: UCLUST-like clustering
4. `cd-hit`: Greedy incremental clustering
5. `set-cover`: Set-Cover (greedy)
6. `leiden`: the Leiden algorithm

Clustering uses one of three similarity measures (`--metric`): `tani`, `gani`, or `ani`. The user-provided threshold for the selected measure sets the minimum similarity required for connecting genomes.

```bash
# Cluster genomes based on ANI and connect genomes if their ANI ≥ 95%.
./vclust cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm single \
--metric ani --ani 0.95
```

Additionally, Vclust may use extra threshold values for various combinations of similarity measures (`--ani`, `--gani`, `--tani`, `--cov`, `--num_alns`). If a genome pair fails to meet the specified thresholds, it is excluded from clustering.

```bash
# Cluster genomes based on ANI similarity measure, but connect them only
# if ANI ≥ 95% and coverage ≥ 85%.
./vclust cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm single \
--metric ani --ani 0.95 --cov 0.85
```

#### Cluster output

The cluster command generates a TSV file with genome identifiers followed by 0-based cluster identifiers. The file includes all input genomes, both clustered and singletons, listed in the same order as the IDs file (sorted by decreasing sequence length). The first genome in each cluster is the representative. For example, cluster `0` has `NC_005091.alt2` as its representative and includes `NC_005091.alt1` and `NC_005091.ref`.

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

Alternatively, instead of numerical cluster identifiers, Vclust can output a representative genome for each cluster using the `--out-repr` option. The representative genome within each cluster is determined as the one with the longest sequence.

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

## 6. Examples

### 6.1. Classify viruses into species and genera using the ICTV standards

The following commands perform VIRIDIC-like analysis by calculating the total ANI (tANI) between complete virus genomes and classifying these viruses into species and genera based on 95% and 70% tANI, respectively.

```bash
# Create a pre-alignment filter for genome pairs with a minimum of 10 common k-mers
# and a minimum sequence identity of 70% (relative to the shortest sequence).
./vclust prefilter -i genomes.fna -o fltr.txt --min_kmers 10 --min-ident 0.7
```

```bash
# Calculate ANI measures for genome pairs specified in the filter.
./vclust align -i genomes -o ani.tsv --filter fltr.txt
```

```bash
# Assign viruses into putative species (tANI ≥ 95%).
./vclust cluster -i ani.tsv -o species.tsv --ids ani.ids.tsv --algorithm complete \
--metric tani --tani 0.95

# Assign viruses into putative genera (tANI ≥ 70%).
./vclust cluster -i ani.tsv -o genus.tsv --ids ani.ids.tsv --algorithm complete \
--metric tani --tani 0.70
```

### 6.2. Assign viral contigs into vOTUs using the MIUViG standards

The following commands assign contigs into viral operational taxonomic units (vOTUs) based on the MIUViG thresholds (ANI ≥ 95% and aligned fraction ≥ 85%).

```bash
# Create a pre-alignment filter.
./vclust prefilter -i genomes.fna -o fltr.txt --min_kmers 30 --min_ident 0.9
```

```bash
# Calculate ANI measures for genome pairs specified in the filter.
./vclust align -i genomes -o ani.tsv --filter fltr.txt
```

```bash
# Cluster contigs into vOTUs using the MIUVIG thresholds and the Leiden algorithm.
./vclust cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm leiden \
--metric ani --ani 0.95 --cov 0.85
```

### 6.3. Dereplicate genomes

The following commands reduce the sequence dataset to representative genomes.

```bash
# Create a pre-alignment filter.
./vclust prefilter -i genomes.fna -o fltr.txt --min_kmers 30 --min_ident 0.90
```

```bash
# Calculate ANI measures for genome pairs specified in the filter.
./vclust align -i genomes.fna -o ani.tsv --filter fltr.txt --out-ani 0.9
```

```bash
# Cluster contigs using the UCLUST algorithm and show representative genome.
./vclust cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm uclust \
--metric ani --ani 0.95 --cov 0.85 --out-repr
```

### 6.4. Calculate pairwise similarities between all-versus-all genomes

The following command calculates ANI measures between all genome pairs in the dataset. For small datasets, using `prefilter` is optional. However, note that without it, Vclust will perform all-versus-all pairwise sequence alignments.

```bash
./vclust.py align -i genomes.fna -o ani.tsv --outfmt complete
```

### 6.5. Process large datasets

The following commands reduce RAM usage and hard disk storage for processing metagenomic contigs from the [IMG/VR](https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html) database.

```bash
# Create a pre-alignment filter by processing batches of 5 million genomes.
./vclust prefilter -i genomes.fna -o fltr.txt --min_kmers 30 --min_ident 0.90 \
--batch-size 5000000
```

```bash
# Calculate ANI measures for genome pairs specified in the filter. Keep the output TSV
# file relatively small-sized: use the lite output format and report only genome pairs 
# with ANI ≥ 90% and coverage ≥ 80%.
./vclust align -i genomes -o ani.tsv --filter fltr.txt --outfmt lite \ 
--out-ani 0.9 --out-cov 0.8
```

```bash
# Cluster contigs into vOTUs using the MIUVIG thresholds and the Leiden algorithm.
./vclust cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm leiden \
--metric ani --ani 0.95 --cov 0.85
```


## 7. Test

To ensure that Vclust works as expected, you can run tests using [pytest](https://docs.pytest.org/).

```bash
pytest test.py
```

## 8. Cite

Zielezinski A, Gudyś A, Barylski J, Siminski K, Rozwalak P, Dutilh BE, Deorowicz S. *Ultrafast and accurate sequence alignment and clustering of viral genomes*. bioRxiv [[doi](https://google.pl)][[pubmed](https://google.pl)].

## 9. License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
