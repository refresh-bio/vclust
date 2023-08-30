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
```