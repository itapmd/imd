#
# Makefile 2dvis
#

CC = gcc
CFLAGS = -O6 -w -L/usr/X11/lib -DTWOD
LIBS = -lvogle -lX11 -lm

all:	2dvis

2dvis:	2dvis.o client.o sockutil.o
	$(CC) $(CFLAGS) -o 2dvis 2dvis.o client.o sockutil.o $(LIBS)

2dvis.o:	2dvis.c
	$(CC) $(CFLAGS) -c 2dvis.c

client.o:	client.c client.h sockutil.h tokens.h
	$(CC) $(CFLAGS) -c client.c

sockutil.o:	sockutil.c sockutil.h
	$(CC) $(CFLAGS) -c sockutil.c

clean:
	rm -f *.o 














