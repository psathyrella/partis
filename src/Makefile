CFLAGS = -Wall -Wextra -Isamtools -std=c99 -g -O2 -pthread
CXXFLAGS = -Wall -std=c++0x
LDLIBS = -lz

ALL = sw-align sort-fq-by-seq
DEPS = kstring.c ksw.c
OBJ = $(DEPS:.c=.o)

PREFIX=..

all: $(ALL)

clean:
	rm -f $(ALL) $(ALL:=.o) $(OBJ)

sw-align: LDLIBS+=-lpthread

sw-align: $(OBJ)
sw-align.o: kseq.h ksort.h kstring.h ksw.h kvec.h

install: $(ALL)
	mkdir -p $(PREFIX)/bin
	install -m 775 $< $(PREFIX)/bin

.PHONY: all clean


		
