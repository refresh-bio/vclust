# <img src="./images/logo.svg" alt="Vclust logo" /> Vclust

![version](https://img.shields.io/badge/version-1.2.3-blue.svg)
[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/vclust/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/vclust/releases)
[![GitHub Actions CI](../../workflows/GitHub%20Actions%20CI/badge.svg)](../../actions/workflows/main.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

![x86-64](https://img.shields.io/static/v1?label=%E2%80%8B&message=x86-64&color=yellow&logo=PCGamingWiki&logoColor=white)
![ARM](https://img.shields.io/static/v1?label=%E2%80%8B&message=ARM&color=yellow&logo=Raspberry%20Pi&logoColor=white)
![Apple M](https://img.shields.io/static/v1?label=%E2%80%8B&message=Apple%20M&color=yellow&logo=Apple&logoColor=white)
![Linux](https://img.shields.io/static/v1?label=%E2%80%8B&message=Linux&color=00A98F&logo=linux&logoColor=white)
![macOS](https://img.shields.io/badge/%E2%80%8B-macOS-00A98F?logo=apple)

Vclust is an alignment-based tool for fast and accurate calculation of Average Nucleotide Identity (ANI) between complete or metagenomically-assembled viral genomes. The tool also performs ANI-based clustering of genomes according to standards recommended by international virus consortia, including *International Committee on Taxonomy of Viruses* (ICTV) and *Minimum Information about an Uncultivated Virus Genome* (MIUViG). 


## Table of contents

1. [Features](#1-features)
2. [Requirements](#2-requirements)
3. [Installation](#3-installation)
4. [Quick Start](#4-quick-start)
5. [Input data](#5-input-data)
6. [Usage](#6-usage)
   1. [Prefilter](#61-prefilter)
   2. [Align](#62-align)
   3. [Cluster](#63-cluster)
7. [Use cases](#7-use-cases)
   1. [Classify viruses into species and genera following ICTV standards](#71-classify-viruses-into-species-and-genera-following-ictv-standards)
   2. [Assign viral contigs into vOTUs following MIUViG standards](#72-assign-viral-contigs-into-votus-following-miuvig-standards)
   3. [Dereplicate viral contigs into representative genomes](#73-dereplicate-viral-contigs-into-representative-genomes)
   4. [Calculate pairwise similarities between all-versus-all genomes](#74-calculate-pairwise-similarities-between-all-versus-all-genomes)
   5. [Process large dataset of diverse virus genomes (IMG/VR)](#75-process-large-dataset-of-diverse-virus-genomes-imgvr)
   6. [Process large dataset of highly redundant virus genomes](#76-process-large-dataset-of-highly-redundant-virus-genomes)
   7. [Cluster plasmid genomes into pOTUs](#77-cluster-plasmid-genomes-into-potus)
8. [FAQ](#8-faq)
9. [Tests](#9-tests)
10. [Citation](#10-citation)
11. [License](#11-license)


## 1. Features

#### :gem: Accurate ANI calculations

Vclust uses a Lempel-Ziv-based pairwise sequence aligner ([LZ-ANI](https://github.com/refresh-bio/LZ-ANI)) for ANI calculation. LZ-ANI achieves high sensitivity in detecting matched and mismatched nucleotides, ensuring accurate ANI determination. Its efficiency comes from a simplified indel handling model, making LZ-ANI  magnitudes faster than alignment-based tools (e.g., BLASTn, MegaBLAST) while maintaining comparable accuracy to the most sensitive BLASTn searches.

#### :triangular_ruler: Multiple similarity measures

Vclust offers multiple similarity measures between two genome sequences:
- **ANI**: The number of identical nucleotides across local alignments divided by the total length of the alignments.
- **Global ANI (gANI)**: The number of identical nucleotides across local alignments divided by the length of the query/reference genome.
- **Total ANI (tANI)**: The number of identical nucleotides between query-reference and reference-query genomes divided by the sum length of both genomes. tANI is equivalent to the VIRIDIC's intergenomic similarity.
- **Coverage (alignment fraction)**: The proportion of the query/reference sequence aligned with the reference/query sequence.
- **Number of local alignments**: The count of individual alignments found between the sequences.
- **Ratio between genome lengths**: Ratio by dividing the length of the shorter sequence by the length of the longer sequence.

#### :star2: Multiple clustering algorithms 

Vclust provides six clustering algorithms tailored to various scenarios, including taxonomic classification and dereplication of viral genomes.
- Single-linkage
- Complete-linkage
- UCLUST
- CD-HIT (Greedy incremental)
- Greedy set cover (adopted from MMseqs2)
- Leiden algorithm [optional]

#### :fire: Speed and efficiency 

Vclust uses three efficient C++ tools - [Kmer-db](https://github.com/refresh-bio/kmer-db), [LZ-ANI](https://github.com/refresh-bio/LZ-ANI), [Clusty](https://github.com/refresh-bio/clusty) - for prefiltering, aligning, calculating ANI, and clustering viral genomes. This combination enables the processing of millions of virus genomes within a few hours on a mid-range workstation.

#### :earth_americas: Web service

For datasets containing up to 1000 viral genomes, Vclust is available at [http://www.vclust.org](http://www.vclust.org).

## 2. Requirements

Vclust requires Python 3.7 or higher.

## 3. Installation

To install Vclust, you can either download the pre-compiled binaries, install from Bioconda, or compile the dependencies from source. The compilation process typically takes a few minutes.

### 3.1. Download precompiled binaries

The quickest way to get started is by downloading prebuilt binaries from the [Releases tab](https://github.com/refresh-bio/vclust/releases). These binaries include the Leiden algorithm by default. Select your platform and download the tool.

### 3.2. Installation via Bioconda

Vclust is also available on Bioconda.

```bash
TODO
```

### 3.3. Compile from source

#### 3.3.1. Requirements

1. `make`
2. `g++` version 11 or higher
3. `cmake` version 3.12 or higher

#### 3.3.2. Default Installation

The default installation of Vclust includes all functionalities except for the Leiden clustering algorithm.

```bash
git clone --recurse-submodules https://github.com/refresh-bio/vclust
cd vclust
make -j
```

#### 3.3.3. Installation with Leiden algorithm support

Vclust provides igraph's implementation of the Leiden algorithm. However, because igraph requires several external dependencies (CMake 3.18, Flex, Bison), it is not integrated with Vclust by default. To install these dependencies under Debian/Ubuntu Linux, use the following command:

```bash
sudo apt-get install cmake flex bison
```

Then, build Vclust with Leiden algorithm support:

```bash
git clone --recurse-submodules https://github.com/refresh-bio/vclust
cd vclust
make -j LEIDEN=true
```

## 4. Quick start

Follow these steps to quickly run Vclust on the provided example genomes. The process takes just a few seconds.

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
./vclust.py cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --metric ani --ani 0.95
```

## 5. Input data

Vclust accepts a single FASTA file containing viral genomic sequences ([example](./example/multifasta.fna)) or a directory of FASTA files (one genome per file) ([example](./example/fna/)). The input file(s) can be gzipped.

## 6. Usage

Vclust provides three commands: `prefilter`, `align`, and `cluster`. Calls to these commands follow the structure:

```
./vclust.py command -i <input> -o <output> [options]
```


### 6.1. Prefilter

The `prefilter` command creates a pre-alignment filter that reduces the number of genome pairs to be aligned by filtering out dissimilar sequences before the alignment step. This process siginificantly lowers the computational demands by keeping only genome pairs with sufficient *k*-mer-based sequence similarity.

Two main options control the prefiltering process:

* `--min-kmers`: Minimum number of common k-mers required between two sequences (default: 30).
* `--min-ident`: Minimum sequence identity, calculated over the shorter sequence [0-1] (default: 0.7).

Vclust computes these similarities using the entire set of *k*-mers with [Kmer-db 2](https://github.com/refresh-bio/kmer-db), overcoming the sampling limitations of sketching methods like FastANI and Mash. It also uses sparse distance matrices to efficiently process large datasets.

By default, the `prefilter` command selects genome pairs with at least 30 common 25-mers and 70% identity over the shorter sequence:

```bash
# Create pre-alignment filter with 30 common 25-mers
# and 70% identity over the shorter sequence.
./vclust.py prefilter -i genomes.fna -o fltr.txt 
```

The sequence identity computed during prefiltering is **higher** than the ANI value from the alignment step. Thus, you can safely set the `--min-ident` value close to the final ANI threshold without losing relevant genome pairs. For example, if you're targeting an ANI threshold of 95% or higher, set `--min-ident` to around 0.95:

```bash
# Filter genome pairs with at least 30 common 21-mers and 90% identity
./vclust.py prefilter -i genomes.fna -o fltr.txt --k 21 --min-kmers 30 --min-ident 0.90 
```

#### 6.1.1. Optimizing RAM usage and speed

To reduce memory consumption or improve processing speed, you can use the following three methods, either individually or in combination:

##### 1. Process genomes in smaller batches

Vclust can lower memory usage by processing the genome dataset in smaller, equally-sized batches. While this approach may slightly increase runtime, it significantly reduces memory consumption:

```bash
# Process genomes in batches of 5 million sequences.
./vclust.py prefilter -i genomes.fna -o fltr.txt --batch-size 5000000
```

##### 2. Analyze only a fraction of *k*-mers

To speed up filtering, Vclust can analyze only a fraction of the *k*-mers for each genome. The `--kmers-fraction` option controls the proportion [0-1] of *k*-mers used in comparisons:

```bash
# Process genomes in batches and analyze 10% of k-mers in each genome sequence.
./vclust.py prefilter -i genomes.fna -o fltr.txt --batch-size 5000000 --kmers-fraction 0.1
```

##### 3. Limit the number of target sequences per query

For highly redundant datasets (e.g., hundreds of thousands of nearly identical genomes), the `prefilter` step may still pass a large number of genome pairs, increasing both memory usage and runtime. The `--max-seqs` option limits the number of target sequences reported for each query genome, reducing the overall number of genome pairs passed to alignment. For each query, `--max-seqs` returns up to *n* sequences that have passed the `--min-kmers` and `--min-ident` filters, and have the highest sequence identity to query sequence. For example, in a dataset containing 1 million nearly identical genomes, the total number of possible genome pairs is almost 500 billion. Setting `--max-seqs 1000` reduces this to 1 billion genome pairs, significantly decreasing memory usage and runtime.

```bash
# Limit the number of target sequences to 1000 per query genome.
./vclust.py prefilter -i genomes.fna -o fltr.txt --batch-size 100000 --max-seqs 1000
```


### 6.2. Align

The `align` command performs pairwise sequence alignments of viral genomes and provides similarity measures like ANI and coverage (alignment fraction). It uses the [LZ-ANI](https://github.com/refresh-bio/LZ-ANI) aligner, which, like BLAST-based methods, finds multiple local alignments (similar to HSPs in BLAST) between two genomic sequences to estimate ANI and other sequence similarity measures for that genome pair.

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

#### 6.2.1. ANI output

The `align` command creates two TSV files: one with ANI values for genome pairs and another with genome identifiers sorted by decreasing sequence length. Both TSV files can be used as input for Vclust's clustering. 

The following command will create two TSV files: [ani.tsv](./example/output/ani.tsv) and [ani.ids.tsv](./example/output/ani.ids.tsv):

```bash
./vclust.py align -i genomes.fna -o ani.tsv --filter fltr.txt
```

Sample output:

```
qidx  ridx  query reference   tani  gani  ani   qcov  rcov  num_alns len_ratio
9  8  NC_010807.alt3 NC_010807.alt2 0.972839 0.960192 0.986657 0.973177 0.997608 60 0.9836
8  9  NC_010807.alt2 NC_010807.alt3 0.972839 0.985279 0.987642 0.997608 0.973177 67 0.9836
10 8  NC_010807.alt1 NC_010807.alt2 0.987250 0.987041 0.987117 0.999923 0.999901 34 0.9571
8  10 NC_010807.alt2 NC_010807.alt1 0.987250 0.987449 0.987547 0.999901 0.999923 36 0.9571
11 8  NC_010807.ref  NC_010807.alt2 0.989807 0.989540 0.989617 0.999923 1.000000 14 0.9571
8  11 NC_010807.alt2 NC_010807.ref  0.989807 0.990063 0.990063 1.000000 0.999923 14 0.9571
10 9  NC_010807.alt1 NC_010807.alt3 0.979963 0.993250 0.994557 0.998686 0.972575 71 0.9730
9  10 NC_010807.alt3 NC_010807.alt1 0.979963 0.967035 0.994304 0.972575 0.998686 70 0.9730
11 9  NC_010807.ref  NC_010807.alt3 0.983839 0.997166 0.997217 0.999948 0.974230 52 0.9730
9  11 NC_010807.alt3 NC_010807.ref  0.983839 0.970871 0.996552 0.974230 0.999948 52 0.9730
11 10 NC_010807.ref  NC_010807.alt1 0.997462 0.997475 0.997475 1.000000 1.000000 23 1.0000
10 11 NC_010807.alt1 NC_010807.ref  0.997462 0.997449 0.997449 1.000000 1.000000 23 1.0000
```

To manage the output TSV file size, Vclust offers three formats: `standard`, `lite`, `complete`.

```bash
./vclust.py align -i genomes.fna -o ani.tsv --filter fltr.txt --outfmt lite
```

| Column | Standard | Lite | Complete | Description |
| --- | :---: |:---: | :---: | --- |
| qidx | + | + | +  | Index of query sequence |
| ridx | + | + | +  | Index of reference sequence |
| query | + | - | +  | Identifier (name) of query sequence |
| reference | + | - | +  | Identifier (name) of reference sequence |
| tani | + | + | +  | total ANI [0-1] |
| gani | + | + | +  | global ANI [0-1] |
| ani | + | + | +  | ANI [0-1] |
| qcov | + | + | +  | Query coverage (aligned fraction) [0-1] |
| rcov | + | + | +  | Reference coverage (aligned fraction) [0-1] |
| num_alns | + | + | +  | Number of local alignments |
| len_ratio | + | + | +  | Length ratio between shorter and longer sequence [0-1] |
| qlen | - | - | +  | Query sequence length |
| rlen | - | - | +  | Reference sequence length |
| nt_match | - | - | +  | Number of matching nucleotides across alignments |
| nt_mismatch | - | - | +  | Number of mismatching nucleotides across alignments |

#### 6.2.2. ANI output filtering

The `align` command enables filtering output by setting minimum thresholds for similarity measures. This ensures only genome pairs meeting the thresholds are reported, significantly reducing the output TSV file size.

```bash
# Output only genome pairs with ANI ≥ 0.95 and query coverage ≥ 0.85.
./vclust.py align -i genomes.fna -o ani.tsv --filter fltr.txt --out-ani 0.95 --out-qcov 0.85
```

#### 6.2.3. Alignment output

Vclust can output alignment details in a separate TSV file. This output format is similar to the BLASTn tabular output and includes information on each local alignment between two genomes, such as the coordinates in both the query and reference sequences, strand orientation, the number of matched and mismatched nucleotides, and the percentage of sequence identity.

```bash
./vclust.py align -i genomes.fna -o ani.tsv --aln ani.aln.tsv
```

Sample output:

```
query reference   pident   alnlen   qstart   qend  rstart   rend  nt_match nt_mismatch
NC_025457.alt2 NC_025457.ref  89.2893  999   22119 23117 14207 15163 892   107
NC_025457.alt2 NC_025457.ref  89.8305  826   3373  4198  2202  3020  742   84
NC_025457.alt2 NC_025457.ref  91.0804  796   41697 42492 27680 28475 725   71
NC_025457.alt2 NC_025457.ref  87.2483  745   38039 38783 24969 25688 650   95
NC_025457.alt2 NC_025457.ref  89.8860  702   7269  7970  5077  5778  631   71
NC_025457.alt2 NC_025457.ref  93.2081  692   62572 63263 41329 42020 645   47
NC_025457.alt2 NC_025457.ref  90.9565  575   31121 31695 20438 21003 523   52
NC_025457.alt2 NC_025457.ref  90.6195  565   11476 12040 7999  8563  512   53
NC_025457.alt2 NC_025457.ref  91.6211  549   10905 11453 7455  8003  503   46
NC_025457.alt2 NC_025457.ref  86.7041  534   29624 30157 19067 19586 463   71
NC_025457.alt2 NC_025457.ref  93.5673  513   10149 10661 6915  7427  480   33
NC_025457.alt2 NC_025457.ref  89.3701  508   34017 34524 22188 22695 454   54
NC_025457.alt2 NC_025457.ref  88.0240  501   18330 18830 11549 12049 441   60
```

| Column | Description |
| --- | --- |
| query | Identifier (name) of query sequence |
| reference | Identifier (name) of reference sequence |
| pident | Percent identity of local alignment |
| alnlen | Alignment length |
| qstart | Start of alignment in query |
| qend | End of alignment in query |
| rstart | Start of alignment in reference |
| rend | End of alignment in reference |
| nt_match | Number of matched (identical) nucleotides  |
| nt_mismatch | Number of mismatching nucleotides |


### 6.3. Cluster

The `cluster` command takes as input two TSV files (outputs of `vclust align`) and clusters viral genomes using a user-selected clustering algorithm and similarity measures. Internally, Vclust uses [Clusty](https://github.com/refresh-bio/clusty), enabling large-scale clustering of millions of genome sequences by leveraging sparse distance matrices.

Vclust supports six clustering algorithms (`--algorithm`):

1. `single`: Single-linkage clustering
2. `complete`: Complete-linkage clustering
3. `uclust`: UCLUST-like clustering
4. `cd-hit`: Greedy incremental clustering
5. `set-cover`: Set-Cover (greedy)
6. `leiden`: the Leiden algorithm


#### 6.3.1. Clustering with similarity measures and thresholds

Vclust performs threshold-based clustering by assigning a genome sequence to a cluster if its similarity (e.g., ANI) to the cluster meets or exceeds a user-defined threshold. The specific similarity criterion depends on the clustering algorithm and may refer to the nearest genome, the furthest genome, or the centroid within the cluster.

Users can choose from three similarity measures for clustering using the `--metric` option: `tani`, `gani`, or `ani`. The threshold provided for the selected measure determines the minimum similarity required to connect two genomes. For asymmetric measures like ANI and gANI (where the similarity between genome A and B differs from that between B and A), Vclust uses the maximum value as the weight for clustering.

```bash
# Cluster genomes based on tANI and connect genomes if their ANI ≥ 95%.
./vclust.py cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm single \
--metric ani --ani 0.95
```

Additionally, Vclust may use extra threshold values for various combinations of similarity measures (`--ani`, `--gani`, `--tani`, `--qcov`, `--rcov`, `--num_alns`, `--len_ratio`). If a genome pair fails to meet the specified thresholds, it is excluded from clustering.

```bash
# Cluster genomes based on ANI similarity measure, but connect them only
# if ANI ≥ 95% and query coverage ≥ 85% and reference coverage ≥ 85%.
./vclust.py cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm single \
--metric ani --ani 0.95 --qcov 0.85 --rcov 0.85
```

#### 6.3.2. Clustering with the Leiden algorithm

Vclust leverages the [Leiden algorithm](https://doi.org/10.1038/s41598-019-41695-z) via the [igraph](https://igraph.org) implementation integrated into [Clusty](https://github.com/refresh-bio/clusty). As with other Vclust algorithms, you can select from three similarity measures: tANI, gANI, or ANI. These measures define the edge weights used by the Leiden algorithm.

```bash
# Cluster genomes based on ANI similarity measure, but connect them only
# if ANI ≥ 95% and query coverage ≥ 85%.
./vclust.py cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm leiden \
--metric ani --ani 0.95 --qcov 0.85
```

If you wish to include coverage information in the edge weight, following the [IMG/PR's approach](https://github.com/apcamargo/pyleiden) (where edge weight = ANI x max(qcov, rcov)), use the `gani` measure for clustering.

```bash
./vclust.py cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm leiden \
--metric gani --gani 0.5 --ani 0.95 --qcov 0.85
```

The granularity of clustering can be adjusted with the `--leiden-resolution` parameter [range: 0-1, default: `0.7`]. Higher resolutions produce more, smaller clusters, while lower resolutions result in fewer, larger clusters.

```bash
./vclust.py cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm leiden \
--metric ani --ani 0.95 --qcov 0.85 --leiden-resolution 1
```

The clustering can be further refined using two additional parameters. The `--leiden-beta` parameter (default: `0.01`) controls the randomness in the Leiden algorithm, allowing for fine-tuning of the clustering process. The `--leiden-iterations` parameter (default: `2`) specifies the number of iterations for the Leiden algorithm, where increasing the number may lead to improved clustering quality.

```bash
./vclust.py cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm leiden \
--metric ani --ani 0.95 --qcov 0.85 --leiden-resolution 0.8 --leiden-beta 0.02 \
--leiden-iterations 3
```


#### 6.3.3. Cluster output

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

## 7. Use cases

### 7.1. Classify viruses into species and genera following ICTV standards

The following commands perform VIRIDIC-like analysis by calculating the total ANI (tANI) between complete virus genomes and classifying these viruses into species and genera based on 95% and 70% tANI, respectively.

```bash
# Create a pre-alignment filter for genome pairs with a minimum of 10 common k-mers
# and a minimum sequence identity of 70% (relative to the shortest sequence).
./vclust.py prefilter -i genomes.fna -o fltr.txt --min-kmers 10 --min-ident 0.7
```

```bash
# Calculate ANI measures for genome pairs specified in the filter.
./vclust.py align -i genomes -o ani.tsv --filter fltr.txt
```

```bash
# Assign viruses into putative species (tANI ≥ 95%).
./vclust.py cluster -i ani.tsv -o species.tsv --ids ani.ids.tsv --algorithm complete \
--metric tani --tani 0.95

# Assign viruses into putative genera (tANI ≥ 70%).
./vclust.py cluster -i ani.tsv -o genus.tsv --ids ani.ids.tsv --algorithm complete \
--metric tani --tani 0.70
```

### 7.2. Assign viral contigs into vOTUs following MIUViG standards

The following commands assign contigs into viral operational taxonomic units (vOTUs) based on the MIUViG thresholds (ANI ≥ 95% and aligned fraction ≥ 85%).

```bash
# Create a pre-alignment filter.
./vclust.py prefilter -i genomes.fna -o fltr.txt --min-kmers 30 --min-ident 0.9
```

```bash
# Calculate ANI measures for genome pairs specified in the filter.
./vclust.py align -i genomes -o ani.tsv --filter fltr.txt
```

```bash
# Cluster contigs into vOTUs using the MIUVIG thresholds and the Leiden algorithm.
./vclust.py cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm leiden \
--metric ani --ani 0.95 --qcov 0.85
```

### 7.3. Dereplicate viral contigs into representative genomes

The following commands reduce the sequence dataset to representative genomes.

```bash
# Create a pre-alignment filter.
./vclust.py prefilter -i genomes.fna -o fltr.txt --min-kmers 30 --min-ident 0.90
```

```bash
# Calculate ANI measures for genome pairs specified in the filter.
./vclust.py align -i genomes.fna -o ani.tsv --filter fltr.txt --out-ani 0.9
```

```bash
# Cluster contigs using the CD-HIT algorithm and show representative genome.
./vclust.py cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm cd-hit \
--metric ani --ani 0.95 --qcov 0.85 --out-repr
```

### 7.4. Calculate pairwise similarities between all-versus-all genomes

The following command calculates ANI measures between all genome pairs in the dataset. For small datasets, using `prefilter` is optional. However, note that without it, Vclust will perform all-versus-all pairwise sequence alignments.

```bash
./vclust.py align -i genomes.fna -o ani.tsv
```

### 7.5. Process large dataset of diverse virus genomes (IMG/VR)

Vclust is optimized for efficient comparison of large viral genome and contig datasets, especially those containing diverse viruses across a broad range of sequence identities. Below is an example of processing over 15 million metagenomic contigs from the [IMG/VR](https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html) database.

```bash
# Create a pre-alignment filter by processing batches of 5 million genomes.
./vclust.py prefilter -i genomes.fna -o fltr.txt --min-kmers 30 --min-ident 0.90 \
--batch-size 5000000
```

```bash
# Calculate ANI values for genome pairs specified in the filter. To keep the output
# file compact, use the `lite` format and only report genome pairs with ANI ≥ 95%
# and query coverage ≥ 85%.
./vclust.py align -i genomes -o ani.tsv --filter fltr.txt --outfmt lite \ 
--out-ani 0.95 --out-qcov 0.85
```

```bash
# Cluster contigs into vOTUs using the MIUVIG standards and the Leiden algorithm.
./vclust.py cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm leiden \
--metric ani --ani 0.95 --qcov 0.85
```

### 7.6. Process large dataset of highly redundant virus genomes

When working with large datasets containing highly redundant sequences (e.g., hundreds of thousands of nearly identical genomes), prefiltering may still pass a large number of genome pairs for alignment, even when using high thresholds for `--min-kmers` and `--min-ident`. Since most sequences in these datasets are almost identical, this can lead to increased memory usage and longer runtimes for all three Vclust commands (`prefilter`, `align`, `cluster`). To address this, Vclust offers three additional options in the `prefilter` step to reduce RAM usage and improve performance (as detailed in: [6.1.1. Optimizing RAM usage and speed](#611-optimizing-ram-usage-and-speed)). The example below shows the use of all three options simultaneously:

```bash
# Create a pre-alignment filter by processing batches of 100,000 genomes,
# analyzing only 10% of k-mers, and limiting each query genome to 1,000
# target sequences with the highest sequence identity.
./vclust.py prefilter -i genomes.fna -o fltr.txt --min-kmers 30 --min-ident 0.95 \
--batch-size 100000 --kmers-fraction 0.1 --max-seqs 1000
```

```bash
# Calculate ANI values for genome pairs specified in the filter. To keep the output
# file compact, use the `lite` format and only report genome pairs with ANI ≥ 95%
# and query coverage ≥ 95%.
./vclust.py align -i genomes -o ani.tsv --filter fltr.txt --outfmt lite \ 
--out-ani 0.95 --out-qcov 0.95
```

```bash
# Cluster highly similar contigs using the Leiden algorithm.
./vclust.py cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm leiden \
--metric ani --ani 0.97 --qcov 0.95
```

### 7.7. Cluster plasmid genomes into pOTUs

The following commands cluster plasmid genomes into plasmid taxonomic units (PTUs).

```bash
# Create a pre-alignment filter passing genome pairs with at least common 30 k-mers
# and have minimum sequence identity of 30%.
./vclust.py prefilter -i genomes.fna -o fltr.txt --min-kmers 30 --min-ident 0.30
```

```bash
# Calculate ANI values for genome pairs specified in the filter. Report genome pairs
# that meet the IMG/PR thresholds (ANI ≥ 70% and query coverage ≥ 50%)
./vclust.py align -i genomes -o ani.tsv --filter fltr.txt --out-ani 0.70 --out-qcov 0.50
```

```bash
# Cluster the genomes using the Leiden algorithm. The cluster edges are weighted by
# the gani metric (ANI x coverage), with a Leiden resolution parameter set to 0.9.
./vclust.py cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm leiden \
--metric gani --gani 0.35 --leiden-resolution 0.9
```

## 8. FAQ

**Q: Does Vclust handle circularly permuted bacteriophage genomes?**

**A:** Yes, Vclust handles circularly permuted bacteriophage genomes by being robust to sequence rearrangements (e.g., translocations and circular permutations). It calculates ANI and alignment fraction (coverage) across all local alignments between two genomes, even when homologous segments are reordered. In tests with circularly permuted genomes, Vclust showed minimal inaccuracies in ANI and coverage, with a mean absolute error of 0.04% compared to non-permuted genomes. These small discrepancies are due to short alignment breaks at the breakpoint positions in circular genomes.


**Q: How does Vclust's sensitivity compare to BLASTn and MegaBLAST?**

**A:** Vclust is designed to match the sensitivity of BLASTn, which is considered highly reliable for estimating ANI. Like BLASTn, Vclust uses an anchor length of 11 nucleotides to align sequences with high precision. MegaBLAST, in comparison, uses a larger word size of 28 nucleotides, making it less sensitive.


**Q: Can I increase the default minimum sequence identity (0.7) in prefilter if I'm aiming for a higher ANI threshold (0.95)?**

**A:** Yes, you can safely increase the default minimum sequence identity (`--min-ident`) in the prefilter step to target a higher ANI threshold. We designed the sequence identity calculation in the `prefilter` command to be higher than the ANI derived from the subsequent `align` step. Specifically, while the sequence identity is calculated similarly to ANI in Mash, Vclust's calculation is based on the shorter sequence. As a result, the default `--min-ident` of `0.7` can be raised to values closer to the final alignment-based ANI threshold.

In our tests for vOTU clustering (ANI ≥ 95% and AF ≥ 85%), even increasing `--min-ident` to `0.95` during prefiltering did not exclude any genome pairs with an alignment-based ANI of ≥ 95%. Additionally, raising the default `--min-ident` from `0.7` to `0.95` significantly reduces the number of genome pairs requiring alignment, thereby speeding up the alignment step.


## 9. Tests

To ensure that Vclust works as expected, you can run tests using [pytest](https://docs.pytest.org/).

```bash
pytest test.py
```

## 10. Citation

Zielezinski A, Gudyś A, Barylski J, Siminski K, Rozwalak P, Dutilh BE, Deorowicz S. *Ultrafast and accurate sequence alignment and clustering of viral genomes*. bioRxiv [[doi:10.1101/2024.06.27.601020](https://www.biorxiv.org/content/10.1101/2024.06.27.601020)].

## 11. License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
