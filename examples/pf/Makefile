include ../../Makefile.in

C = main.cc rng.cc history.cc simfunctions.cc
H = smctc.hh sampler.hh particle.hh moveset.hh history.hh rng.hh

PFC = pfexample.cc pffuncs.cc 
PFO = pfexample.o pffuncs.o
PFH = pffuncs.hh

CCFLAGS := $(CCFLAGS) -I../../include -L../../lib 
LFLAGS := -I../../include -L../../lib $(LFLAGS)

all: pf

pf: $(PFC) $(PFH)
	$(CC) $(CCFLAGS) -c $(PFC)
	$(CC) $(PFO) -lsmctc $(LFLAGS) -opf
	cp pf ../../bin
	cp data.csv ../../bin

clean:
	-rm *.o
	-rm *~
	-rm pf