all: prep

prep:
	cd src/3rd-party/kmer-db; make -j;
	cd src/3rd-party/lz-ani; make -j;
