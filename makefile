all: prep

prep:
	cd 3rd_party/kmer-db; make -j;
	cd 3rd_party/lz-ani; make -j;
	cd 3rd_party/rapid-cluster; make -j;
	cd 3rd_party/ref-utils; make -j;
	mkdir -p bin
	cp 3rd_party/kmer-db/kmer-db ./bin/
	cp 3rd_party/lz-ani/lz-ani ./bin/
	cp 3rd_party/rapid-cluster/rapid-cluster ./bin/
	cp 3rd_party/ref-utils/multi-fasta-split/multi-fasta-split ./bin/
	
clean:
	cd 3rd_party/kmer-db; make clean;
	cd 3rd_party/lz-ani; make clean;
	cd 3rd_party/rapid-cluster; make clean;
	cd 3rd_party/ref-utils; make clean;
	rm bin/*
