include ../Makefile.in

SMCC = rng.cc history.cc smc-exception.cc
SMCO = rng.o history.o smc-exception.o

all: libsmctc.a

.PHONY: clean

clean:
	-rm *.o
	-rm *.a
	-rm *~

libsmctc.a: $(SMCC) Makefile
	$(CC) $(CCFLAGS) -c -I../include $(SMCC)
	ar rcs libsmctc.a $(SMCO)
	cp libsmctc.a ../lib
