CC     = gcc
CFLAGS = -fopenmp -O3

all: scascade

scascade: source/scascade.c source/queue.c source/prelim.c
	$(CC) $(CFLAGS) -o bin/scascade source/scascade.c

clean:
	rm -f bin/scascade
