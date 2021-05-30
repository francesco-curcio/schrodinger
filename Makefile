libdir=/usr/local/lib/
incdir=/usr/local/include/
CC=g++
CFLAGS=-Wall -O3 -DHAVE_INLINE
LFLAGS=-L$(libdir)
IFLAGS=-I$(incdir)
libs=-lm -lgsl -lgslcblas

main: main.o
	$(CC) -o $@ $(CFLAGS) $(LFLAGS) $(libs) $<

main.o: main.cpp
	$(CC) -c $< $(CFLAGS) $(IFLAGS)
	
.PHONY: clean

clean: 
	rm -f *.o 
