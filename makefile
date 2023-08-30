all: prep

prep:
	cd src/3rd-party/kmer-db; make -j;
	cd src/3rd-party/lz-ani; make -j;
	mkdir -p bin
	cp src/3rd-party/kmer-db/kmer-db ./bin/kmer-db
	cp src/3rd-party/lz-ani/lz-ani-0.1 ./bin/lz-ani
	
clean:
	cd src/3rd-party/kmer-db; make clean;
	cd src/3rd-party/lz-ani; make clean;
	rm bin/*
	rm data/*.a2a
	rm data/*.csv
	rm data/*.kmer-db
	rm data/*.lz-ani
