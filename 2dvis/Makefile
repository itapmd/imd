#
# Makefile 2dvis/3dvis
#

# Target
MAKETARGET = $(@)

# Compiler
CC_iris4d = cc
CC_i386-linux = gcc
CC_alpha = gcc
CC = $(CC_$(HOSTTYPE))

# Compiler flags
CFLAGS = -g

ifneq (,$(strip $(findstring 3d,${MAKETARGET})))
CFLAGS += -DTWOD
endif

ifneq (,$(strip $(findstring debug,${MAKETARGET})))
CFLAGS += -g
endif

# Linker
LIBP = -L/usr/X11/lib
LIBS = -lvogle -lX11 -lm

all:	2dvis

3dvis:	2dvis.o client.o sockutil.o param.o
	$(CC) $(LIBP) $(CFLAGS) -o 2dvis 2dvis.o client.o sockutil.o param.o $(LIBS)
	mv 2dvis $(HOME)/bin/$(HOSTTYPE)

2dvis:	2dvis.o client.o sockutil.o param.o
	$(CC) $(LIBP) $(CFLAGS) -o 2dvis 2dvis.o client.o sockutil.o param.o $(LIBS)
	mv 2dvis $(HOME)/bin/$(HOSTTYPE)

2dvis.o:	2dvis.c globals.h prototypes.h
	$(CC) $(CFLAGS) -c 2dvis.c

client.o:	client.c client.h sockutil.h tokens.h globals.h prototypes.h 
	$(CC) $(CFLAGS) -c client.c

sockutil.o:	sockutil.c sockutil.h globals.h prototypes.h
	$(CC) $(CFLAGS) -c sockutil.c

param.o:	param.c  globals.h prototypes.h
	$(CC) $(CFLAGS) -c param.c

clean:
	rm -f *.o *~ \#*  














