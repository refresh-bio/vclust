all: prep

prep:
	cd src/3rd-party/kmer-db; make -j;
	cd src/3rd-party/lz-ani; make -j;
	cp src/3rd-party/kmer-db/kmer-db ./bin/kmer-db
	cp src/3rd-party/lz-ani/lz-ani-0.1 ./bin/lz-ani