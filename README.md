# <img src="./images/logo.svg" alt="Vclust logo" /> Vclust

![version](https://img.shields.io/badge/version-1.2.7-blue.svg)
[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/vclust/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/vclust/releases)
[![Build and tests](../../workflows/Build%20and%20tests/badge.svg)](../../actions/workflows/main.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

![x86-64](https://img.shields.io/static/v1?label=%E2%80%8B&message=x86-64&color=yellow&logo=PCGamingWiki&logoColor=white)
![ARM](https://img.shields.io/static/v1?label=%E2%80%8B&message=ARM&color=yellow&logo=Raspberry%20Pi&logoColor=white)
![Apple M](https://img.shields.io/static/v1?label=%E2%80%8B&message=Apple%20M&color=yellow&logo=Apple&logoColor=white)
![Linux](https://img.shields.io/static/v1?label=%E2%80%8B&message=Linux&color=00A98F&logo=linux&logoColor=white)
![macOS](https://img.shields.io/badge/%E2%80%8B-macOS-00A98F?logo=apple)

Vclust is an alignment-based tool for fast and accurate calculation of Average Nucleotide Identity (ANI) between complete or metagenomically-assembled viral genomes. The tool also performs ANI-based clustering of genomes according to standards recommended by international virus consortia, including *International Committee on Taxonomy of Viruses* (ICTV) and *Minimum Information about an Uncultivated Virus Genome* (MIUViG). 

## Features

#### :gem: Accurate ANI calculations

Vclust uses a Lempel-Ziv-based pairwise sequence aligner ([LZ-ANI](https://github.com/refresh-bio/LZ-ANI)) for ANI calculation. LZ-ANI achieves high sensitivity in detecting matched and mismatched nucleotides, ensuring accurate ANI determination. Its efficiency comes from a simplified indel handling model, making LZ-ANI  magnitudes faster than alignment-based tools (e.g., BLASTn, MegaBLAST) while maintaining comparable accuracy to the most sensitive BLASTn searches.

#### :triangular_ruler: Multiple similarity measures

Vclust offers multiple similarity measures between two genome sequences:
- **ANI**: The number of identical nucleotides across local alignments divided by the total length of the alignments.
- **Global ANI (gANI)**: The number of identical nucleotides across local alignments divided by the length of the query/reference genome.
- **Total ANI (tANI)**: The number of identical nucleotides between query-reference and reference-query genomes divided by the sum length of both genomes. tANI is equivalent to the [VIRIDIC's intergenomic similarity](https://doi.org/10.3390/v12111268).
- **Coverage (alignment fraction)**: The proportion of the query/reference sequence aligned with the reference/query sequence.
- **Number of local alignments**: The number of local alignments between the two genome sequences.
- **Ratio between genome lengths**: The length of the shorter genome divided by the longer one.

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

## Quick start

```bash
# Clone repository and build Vclust
git clone --recurse-submodules https://github.com/refresh-bio/vclust
cd vclust && make -j

# Prefilter similar genome sequence pairs before conducting pairwise alignments.
./vclust.py prefilter -i example/multifasta.fna -o fltr.txt

# Align similar genome sequence pairs and calculate pairwise ANI measures.
./vclust.py align -i example/multifasta.fna -o ani.tsv --filter fltr.txt

# Cluster genome sequences based on given ANI measure and minimum threshold.
./vclust.py cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --metric ani --ani 0.95
```
## Documentation

The Vclust documentation is available on the [GitHub Wiki](https://github.com/refresh-bio/vclust/wiki) and includes the following sections:

1. [Features](https://github.com/refresh-bio/vclust/wiki/1-Features)
2. [Installation](https://github.com/refresh-bio/vclust/wiki/2-Installation)
3. [Quick Start](https://github.com/refresh-bio/vclust/wiki/3-Quick-start)
4. [Usage](https://github.com/refresh-bio/vclust/wiki/4-Usage)
   1. Input data
   2. Prefilter
   3. Align
   4. Cluster
5. [Optimizing sensitivity and resource usage](https://github.com/refresh-bio/vclust/wiki/5-Optimizing-sensitivity-and-resource-usage)
6. [Use cases](https://github.com/refresh-bio/vclust/wiki/6-Use-cases)
   1. Classify viruses into species and genera following ICTV standards
   2. Assign viral contigs into vOTUs following MIUViG standards
   3. Dereplicate viral contigs into representative genomes
   4. Calculate pairwise similarities between all-versus-all genomes
   5. Process large dataset of diverse virus genomes (IMG/VR)
   6. Process large dataset of highly redundant virus genomes
   7. Cluster plasmid genomes into pOTUs
7. [FAQ: Frequently Asked Questions](https://github.com/refresh-bio/vclust/wiki/7-FAQ:-Frequently-Asked-Questions)


## Citation

Zielezinski A, Gudyś A, Barylski J, Siminski K, Rozwalak P, Dutilh BE, Deorowicz S. *Ultrafast and accurate sequence alignment and clustering of viral genomes*. bioRxiv [[doi:10.1101/2024.06.27.601020](https://www.biorxiv.org/content/10.1101/2024.06.27.601020)].

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
