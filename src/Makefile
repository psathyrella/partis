CFLAGS = -Wall -Wextra -Isamtools -std=c99 -g -O2 -pthread
LDLIBS = -lpthread -lz

ALL = sw-align
DEPS = kstring.c ksw.c
OBJ = $(DEPS:.c=.o)

PREFIX=..

all: $(ALL)

clean:
	rm -f $(ALL) $(ALL:=.o) $(OBJ)

sw-align: $(OBJ)

sw-align.o: kseq.h ksort.h kstring.h ksw.h kvec.h

install: sw-align
	mkdir -p $(PREFIX)/bin
	install -m 775 $< $(PREFIX)/bin

.PHONY: all clean


		
