# vclust-dev

## Klonowanie i kompilacja
Klonowanie wraz z zależnościami:

`git clone --recurse-submodules https://github.com/refresh-bio/vclust-dev`

Kompilacja

```bash
cd vclust-dev
make -j
```

## Wymagania (tymczasaowe)
Zainstalowany moduł `sklearn` (scikit-learn)

## Przetwarzanie danych z zestawu viridica
```bash
bin/kmer-db build -k 18 -t 32 data/viridic.txt data/viridic.kmer-db
bin/kmer-db all2all -sparse -above 3 data/viridic.kmer-db data/viridic.a2a
bin/lz-ani all2all -t 32 --verbose 2 -bs -cd 128 -reg 80 -mml 5 -mdl 11 -mlrim 32 -aw 16 -am 6 -ar 2 --in-file-names data/viridic.txt -out data/viridic.lz-ani -filter data/viridic.a2a 5
python3 src/cluster.py data/viridic.task data/
```


## Quick start

1. Align virus genome sequences:

```bash
./vlust.py align test/fna/ test/alignment.txt
```

2. calculate ANI from alignment:

```bash
./vclust.py test/alignment.txt test/ani.csv
```

3. cluster virus genomic sequences (TODO):

```bash
./vclust.py cluster -h
```

## Usage

```bash
./vclust.py -h
./vclust.py align -h
./vclust.py calcani -h
./vclust.py cluster -h      # TODO
```

## Test

You can run tests to ensure that the scripts work as expected.

```bash
./test.py
```

