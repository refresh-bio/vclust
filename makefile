all: prep

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
	uname_M := "x86_64"
else                          # If uname not available => 'not'
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
	uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')
endif

prep:
	cd 3rd_party/kmer-db && $(MAKE) -j 
	cd 3rd_party/lz-ani && $(MAKE) -j 
	cd 3rd_party/clusty && $(MAKE) -j 
	cd 3rd_party/ref-utils && $(MAKE) -j 
	mkdir -p bin
	cp 3rd_party/kmer-db/kmer-db ./bin/
	cp 3rd_party/lz-ani/lz-ani ./bin/
	cp 3rd_party/clusty/clusty ./bin/
	cp 3rd_party/ref-utils/multi-fasta-split/multi-fasta-split ./bin/
	
clean:
	cd 3rd_party/kmer-db && $(MAKE) clean
	cd 3rd_party/lz-ani && $(MAKE) clean
	cd 3rd_party/clusty && $(MAKE) clean
	cd 3rd_party/ref-utils && $(MAKE) clean
	rm bin/*
