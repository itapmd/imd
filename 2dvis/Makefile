# Makefile beavis
#

# Compiler
CC_iris4d = cc
CC_i386-linux = gcc
CC_alpha = gcc
CC = 

# Linker
LIBP = -L/usr/X11/lib
LIBS = -lvogle -lm -lXaw

# Compiler
CC_iris4d = cc
CC_i386-linux = gcc
CC_alpha = gcc
CC = $(CC_$(HOSTTYPE))

# Headers
HEADERS 	= sockutil.h tokens.h globals.h prototypes.h
SOURCES		= main.c client.c sockutil.c param.c menu.c

# Compiler flags
#CFLAGS = -g

ifneq (,$(strip $(findstring 2d,${MAKETARGET})))
CFLAGS += -DTWOD
endif

ifneq (,$(strip $(findstring debug,${MAKETARGET})))
CFLAGS += -g
endif

# Linker
LIBP = -L/usr/X11/lib
LIBS = -lvogle -lXaw -lm

# Objects
# Substitute .o for .c to get the names of the object files
OBJECTS := $(subst .c,.o,${SOURCES})

# How to compile *.c files
%.o: %.c ${HEADERS}
	${CC} ${CFLAGS} -c $(subst .o,.c,$@)

# How to link
$(MAKETARGET): ${OBJECTS} ${HEADERS}
	$(CC) $(LIBP) $(CFLAGS) -o $@ ${OBJECTS} ${LIBS}
	mv $@ $(HOME)/bin/$(HOSTTYPE)

.DEFAULT:
	${MAKE} MAKETARGET='$@' STAGE2

# Second recursion sets all other Variables (depending on MAKETARGET) and compiles
# An empty MAKETARGET Variable would create an infinite recursion, so we check
STAGE2:
ifneq (,${MAKETARGET})
	${MAKE} OBJECTS='${OBJECTS}'\
                CC='${CC}'\
                CFLAGS='${CFLAGS}'\
                LIBS='${LIBS}'\
                MAKETARGET='${MAKETARGET}' ${MAKETARGET}
else
	@echo 'No TARGET specified.'
endif

clean:
	rm -f *.o *~ \#*  














