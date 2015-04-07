include ../../Makefile.in

C = main.cc rng.cc history.cc simfunctions.cc
H = smctc.hh sampler.hh particle.hh moveset.hh history.hh rng.hh

CLUSTERC = cluster.cc
CLUSTERO = cluster.o
PFH = 

CCFLAGS := $(CCFLAGS) -I../../include -L../../lib 
LFLAGS := -I../../include -L../../lib $(LFLAGS)

all: cluster

cluster: $(CLUSTERC) $(CLUSTERH)
	$(CC) $(CCFLAGS) -c $(CLUSTERC)
	$(CC) $(CLUSTERO) -lsmctc $(LFLAGS) -ocluster
	cp cluster ../../bin

clean:
	-rm *.o
	-rm *~
	-rm cluster
