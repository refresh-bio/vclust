all: prep

prep:
	cd 3rd_party/kmer-db; make -j;
	cd 3rd_party/lz-ani; make -j;
	cd 3rd_party/rapid-cluster; make -j;
	mkdir -p bin
	cp 3rd_party/kmer-db/kmer-db ./bin/kmer-db
	cp 3rd_party/lz-ani/lz-ani ./bin/lz-ani
	cp 3rd_party/rapid-cluster/rapid-cluster ./bin/rapid-cluster
	
clean:
	cd 3rd_party/kmer-db; make clean;
	cd 3rd_party/lz-ani; make clean;
	cd 3rd_party/rapid-cluster; make clean;
	rm bin/*
