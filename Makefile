CC = gcc
IMAGER_DIR = $(shell pwd)/src
CFLAGS=-I$(IMAGER_DIR)
LDFLAGS=-L$(IMAGER_DIR) -L/usr/local/lib
LIBS=-lhdf5 -lbitshuffle -llz4 -ljpeg

all: sporta

sporta: main.o
	$(CC) $(LDFLAGS) main.o -o sporta $(LIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c -o main.o main.c

clean:
	rm -f sporta main.o
