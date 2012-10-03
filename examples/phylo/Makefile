PHYLOC = phylomain.cc phylofunc.cc
PHYLOH = phylofunc.hh

CCFLAGS := $(CCFLAGS) -I../../include -L../../lib
LFLAGS := -I../../include -L../../lib $(LFLAGS)

all: phylo

phylo: phylomain.cc phylofunc.cc
	g++ -std=c++0x -O0 -g $(CCFLAGS) $(PHYLOC) -lsmctc -lgsl -lgslcblas -o phylo

clean:
	rm -f *.o phylo

