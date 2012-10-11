PHYLOC = phylomain.cc phylofunc.cc
PHYLOH = phylofunc.hh

CC ?= gcc
CXX ?= g++

CCFLAGS := $(CCFLAGS) -I../../include -L../../lib -g
LFLAGS := -I../../include -L../../lib $(LFLAGS)

all: phylo

phylo: phylomain.cc phylofunc.cc phylofunc.hh hmsbeagle.h
	$(CXX) -std=c++0x -O0 -g $(CCFLAGS) $(PHYLOC) -lsmctc -lgsl -lgslcblas -o phylo `pkg-config --libs --cflags hmsbeagle-1`

clean:
	rm -f *.o phylo

