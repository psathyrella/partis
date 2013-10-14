CFLAGS = -Wall -Wextra -Isamtools -std=c99 -g -O2 -pthread
LDLIBS = -lpthread -lz

ALL = sw-align
DEPS = kstring.c ksw.c
OBJ = $(DEPS:.c=.o)
DFILES = $(DEPS:.c=.d)

PREFIX=..

all: $(ALL) sw-align.d

-include $(DFILES)

clean:
	rm -f $(ALL) $(ALL:=.o) $(OBJ)

sw-align: $(OBJ)

install: sw-align
	mkdir -p $(PREFIX)/bin
	install -m 775 $< $(PREFIX)/bin

%.d: %.c
	$(CC) -MM $(CFLAGS) $*.c > $*.d

%.o: %.c %.d
	$(CC) -c $(CFLAGS) $*.c -o $*.o
	$(CC) -MM $(CFLAGS) $*.c > $*.d

.PHONY: all clean


		
