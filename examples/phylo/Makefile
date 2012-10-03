PHYLOC = phylomain.cc phylofunc.cc 
PHYLOH = phylofunc.hh

all: phylo

phylo: phylomain.cc phylofunc.cc 
	g++ -std=c++0x -O0 -g -I ../../include -L. -I. $(PHYLOC) -lsmctc -lgsl -lgslcblas -o phylo

clean:
	rm -f *.o phylo

