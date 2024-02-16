all: prep

prep:
	cd 3rd-party/kmer-db; make -j;
	cd 3rd-party/lz-ani; make -j;
	cd 3rd-party/rapid-cluster; make -j;
	mkdir -p bin
	cp 3rd-party/kmer-db/kmer-db ./bin/kmer-db
	cp 3rd-party/lz-ani/lz-ani ./bin/lz-ani
	cp 3rd-party/rapid-cluseter/lz-ani ./bin/rapid-cluseter
	
clean:
	cd src/3rd-party/kmer-db; make clean;
	cd src/3rd-party/lz-ani; make clean;
	cd src/3rd-party/rapid-cluster; make clean;
	rm bin/*
