##############################################################################
#
# Makefile for imd package
#
# This only works with gnumake.
# 
# You need to set the HOSTTYPE envirorment variable according to your machine.
#
# Usage: make <target>
#
# <target> has the form imd_<ensemble>_<architecture>[_<suffix>...]
#
# <ensemble> is one of nve, nvt, and, mik, frac, pull, npt_iso, npt_axial
#
# <architecture> is either risc or mpi or gcc
#
# <suffix> is one or several of 2d, debug, mono, app, prof
#
# This makefile analyzes the target name. It sets the 
# preprocessor flags and chooses the source files accordingly.
#
# Examples:
#
# make imd_nve_risc            ; Make single cpu 3d nve simulation
# make imd_frac_mpi_debug_2d   ; 2d fractures simlation with debugging support
# make imd_and_risc_prof       ; 3d simluation with Andersen thermostat and 
#                                profiling
#
# initially set up by Joerg Stadler (1995-1997)
# revised and extended by Martin Hohl (1997-1998)
#
# $RCSfile$
# $Revision$
# $Date$
##############################################################################

# On some machines EXTRAFLAGS might work only if the shell used by make
# is the same as the shell from which make is invoked (T3E, e.g.).
# In those cases it is possible to set the shell used by make to your
# working shell, by uncommenting and adapting the following line. 
# SHELL         = tcsh

GET           = co

IMDHEADERS    = config.h globals.h imd.h makros.h prototypes.h types.h 

SOURCES       = imd.c imd_maxwell.c imd_integrate.c imd_forces.c \
	        imd_misc.c imd_param.c imd_alloc.c
SOURCES2D     = ${SOURCES} imd_main_2d.c imd_io_2d.c imd_geom_2d.c \
		imd_pictures_2d.c
SOURCES3D     = ${SOURCES} imd_main_3d.c imd_io_3d.c imd_geom_3d.c \
		imd_pictures_3d.c

RISCSOURCES2D = imd_main_risc_2d.c
RISCSOURCES3D = imd_main_risc_3d.c

MPISOURCES    = imd_mpi_util.c
MPISOURCES2D  = ${MPISOURCES} imd_main_mpi_2d.c imd_geom_mpi_2d.c
MPISOURCES3D  = ${MPISOURCES} imd_main_mpi_3d.c imd_geom_mpi_3d.c \
		imd_savemem_3d.c

PACXSOURCES2D   = ${MPISOURCES2D} imd_pacx.c
PACXSOURCES3D   = ${MPISOURCES3D} imd_pacx.c

SOCKHEADERS   = socket_tokens.h sockets.h sockutil.h
SOCKSOURCES   = socket_io.c sockutil.c

FRACSOURCES   = imd_load.c

SHEARSOURCES  = imd_mikshear.c

QUASISOURCES  = imd_qc.c

EAMSOURCES    = imd_forces_eam.c
TTBPSOURCES   = imd_forces_ttbp.c

MCSOURCES     = imd_mc.c

CORRSOURCES   = imd_correl.c

# Collection of my favourite compiler flags

# RISC Version Flags

CFLAGS_iris4d			= -Dsgi -DRISC -O2 -n32 -mips3 -xansi
CFLAGS_iris4d_DEBUG		= -Dsgi -DRISC -g -fullwarn -xansi  
CFLAGS_decstation		= -DRISC -O -float -D_POSIX_SOURCE
CFLAGS_decstation_DEBUG		= -DRISC -g -float -D_POSIX_SOURCE
CFLAGS_alpha			= -DRISC -O3 -float -fp_reorder -tune host
CFLAGS_alpha_DEBUG		= -DRISC -g  -float
CFLAGS_speedy			= -DRISC -O3 -float -fp_reorder -arch host -tune host
CFLAGS_speedy_DEBUG		= -DRISC -g -float
CFLAGS_rs6000			= -DRISC -O3 -qarch=pwr2
CFLAGS_c90			= -DRISC -O3 -htask0,report=fv -DPVPCRAY
CFLAGS_c90_DEBUG		= -DRISC -g
# tcsh-Version
CFLAGS_i386-linux		= -DRISC -O3 -m486
CFLAGS_i386-linux_DEBUG		= -DRISC -g -m486
CFLAGS_i386-FreeBSD		= -DRISC -O3 -m486
CFLAGS_i386-FreeBSD_DEBUG	= -DRISC -g -m486
# Linux/FreeBSD mit bash
CFLAGS_i386			= -DRISC -O3 -m486
CFLAGS_i386_DEBUG		= -DRISC -g -m486
CFLAGS_sx4			= -DRISC -Oinfomsg,div,move,nooverlap -T1 -pvctl,fullmsg -DSX4
CFLAGS_sx4_DEBUG		= -DRISC -g -T1 -DSX4
CFLAGS_t3e			= -Dt3e -DRISC -O3 -haggress,report=isf -X1
CFLAGS_cray-t3e			= -Dt3e -DRISC -O3 -haggress,report=isf -X1

# OS/2
CFLAGS_os2			= -O3 -DRISC

# generic GCC

CFLAGS_GCC_alpha		= -DRISC -O3
CFLAGS_GCC_alpha_DEBUG		= -DRISC -g
CFLAGS_GCC_iris4d		= -DRISC -mcpu=r4000 -mips2 -O2
CFLAGS_GCC_iris4d_DEBUG		= -DRISC -g
CFLAGS_GCC_i386-linux		= -DRISC -O3 -m486
CFLAGS_GCC_i386-linux_DEBUG	= -DRISC -m486 -g
# Pentium GCC (pgcc) mit bash
CFLAGS_GCC_i386			= -DRISC -O3 -mpentium
CFLAGS_GCC_i386_DEBUG		= -DRISC -mpentium -g

# Generic MPI

CFLAGS_MPI_paragon		= -DMPI -O4 -DAR
CFLAGS_MPI_paragon_DEBUG	= -DMPI -g  -DAR
CFLAGS_MPI_decstation		= -DMPI -O -float -D_POSIX_SOURCE
CFLAGS_MPI_decstation_DEBUG	= -DMPI -O -float -D_POSIX_SOURCE
CFLAGS_MPI_rs6000		= -DMPI -O3 -qarch=pwr2
CFLAGS_MPI_sx4			= -DMPI -hpthread_specific -Onodiv,bomove -pvctl -DSX4
CFLAGS_MPI_sx4_DEBUG		= -DMPI -h2,pthread_specific,novector -DSX4
CFLAGS_MPI_t3e			= -Dt3e -DMPI -DAR -O3 -htolerant,aggress,report=isf
CFLAGS_MPI_t3e_DEBUG		= -Dt3e -DMPI -DAR -g
CFLAGS_MPI_cray-t3e		= ${CFLAGS_MPI_t3e}
CFLAGS_MPI_cray-t3e_DEBUG	= ${CFLAGS_MPI_t3e_DEBUG}
CFLAGS_MPI_HITACHI		= -DMPI +O4 +Ol +Oprefetch -I/usr/local/mpi/include  -DMPI_hmpp2
CFLAGS_MPI_HITACHI_DEBUG	= -DMPI -g -I/usr/local/mpi/include  -DMPI_hmpp2
CFLAGS_MPI_i386                 = -DMPI -O3 -m486
CFLAGS_MPI_i386_DEBUG           = -DMPI -O -g -Wall -m486

# PACX MPI

PACXDIR				= ${HOME}/WORK/PACX

CFLAGS_PACX_t3e			= -Dt3e -DPACX -DMPI -I${PACXDIR}/include -DAR -O3 -htolerant,aggress,report=isf
CFLAGS_PACX_t3e_DEBUG		= -Dt3e -DPACX -DMPI -I${PACXDIR}/include -DAR -g
CFLAGS_PACX_cray-t3e		= ${CFLAGS_PACX_t3e}
CFLAGS_PACX_cray-t3e_DEBUG	= ${CFLAGS_PACX_t3e_DEBUG}
CFLAGS_PACX_paragon		= -DPACX -DMPI -O4 -I${PACXDIR} -DAR
CFLAGS_PACX_paragon_DEBUG	= -DPACX -DMPI -g  -I${PACXDIR} -DAR

# Files save to vectorize on sx4
SAVETOVECTORIZE			= imd_forces.o imd_integrate.o

# Compilers
CC_RISC_decstation		= cc
CC_RISC_rs6000			= cc
CC_RISC_sx4			= cc
CC_RISC_os2			= gcc
CC_RISC_iris4d			= cc
CC_RISC_alpha			= cc
CC_RISC_t3e			= cc
CC_RISC_cray-t3e		= cc
# Linux mit tcsh
CC_RISC_i386-linux		= gcc
CC_RISC_i386-linux_DEBUG	= gcc
# Linux mit bash
CC_RISC_i386 			= gcc
CC_RISC_i386_DEBUG		= gcc
# FreeBSD mit tcsh
CC_RISC_i386-FreeBSD		= gcc
CC_RISC_i386-FreeBSD_DEBUG	= gcc
#
CC_GCC_alpha			= gcc
CC_GCC_alpha_debug		= gcc
CC_GCC_iris4d			= gcc
CC_GCC_i386-linux		= /app/bin/gcc
CC_GCC_i386			= /usr/local/bin/gcc

#CC_MPI_paragon			= mpicc
CC_MPI_paragon			= cc
CC_MPI_decstation		= mpicc
CC_MPI_rs6000			= mpcc
CC_MPI_sx4			= cc
CC_MPI_t3e			= cc
CC_MPI_cray-t3e			= cc
CC_MPI_c90			= cc -Tcray_t3d
CC_MPI_HITACHI			= cc
CC_MPI_i386			= hcc

CC_PACX_paragon			= cc
CC_PACX_t3e			= cc
CC_PACX_cray-t3e		= cc

#Profiling
CFLAGS_PROF_t3e			:= -happrentice
CFLAGS_PROF_cray-t3e		:= -happrentice
CFLAGS_PROF_sx4			:= -hacct

# Libraries
MPILIBS_t3e     		= -lmpi 
MPILIBS_paragon			= -lmpi -nx
MPILIBS_decstation		= 
MPILIBS_sx4			= -lmpi -lpthread
MPILIBS_HITACHI			= -L/usr/local/mpi/lib/hmpp2/cml -lmpi
MPILIBS_i386                    = -lmpi

PACXLIBS_t3e			= -L ${PACXDIR}/lib -lpacx -llzo -lmpi
PACXLIBS_cray-t3e		= -L ${PACXDIR}/lib -lpacx -llzo -lmpi
PACXLIBS_paragon		= -lpacx -llzo -lmpi -nx

LIBS_PROF_t3e 			:= -lapp
LIBS_PROF_cray-t3e 		:= -lapp

# set flags accoring to MAKETARGET
# this uses the HOSTTYPE envirorment variable
# currently, there are CLFAGS for the following HOSTTYPES:
# 
# iris4d		SGI Irix 5.3		RISC
# decstation            DEC Ultirx 4.4          RISC MPI
# alpha                 DEC UNIX 3.2            RISC
# rs6000                AIX (Power2 arch)       RISC MPI
# paragon               Intel Paragon                MPI
# c90                   Cray Unicos             RISC MPI(T3D hosted by c90)
# sx4                   NEC Super UX            RISC
# t3e			Cray Unicos/MP		     MPI
# HITACHI		HI-UX/MPP		     MPI

CC_MPI		  = ${CC_MPI_${HOSTTYPE}}
CC_PACX		  = ${CC_PACX_${HOSTTYPE}}
CC_GCC		  = ${CC_GCC_${HOSTTYPE}}
CC_RISC           = ${CC_RISC_${HOSTTYPE}}

MPILIBS		  = ${MPILIBS_${HOSTTYPE}}	
PACXLIBS	  = ${PACXLIBS_${HOSTTYPE}}	

# to link or not to link

LIBS     	:= -lm
RISCLIBS 	:= 
T3DLIBS  	:= -lmpi -lsci 

# Executeables are put into $HOME/$BINDIR$/HOSTTYPE
BINDIR   = bin

# Local releases are done by coping the source to /app/itapsoft
RELEASEPATH = /app/itapsoft/simulation/imd


###########################################################################
#
#	 IMD Configuration rules
#
###########################################################################


# Sockets are not included by default
ifneq	(,$(strip $(findstring sock,${MAKETARGET})))
CFLAGS += -DUSE_SOCKETS
endif

# Erst machen wir RISC
# Architecture & Dimension
ifneq (,$(strip $(findstring risc,${MAKETARGET})))
CC      := ${CC_RISC}

# twod?
ifneq (,$(strip $(findstring 2d,${MAKETARGET})))
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES2D} ${RISCSOURCES2D}
CFLAGS  += -DTWOD
else
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES3D} ${RISCSOURCES3D}
endif

#general compiler flags
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += ${CFLAGS_${HOSTTYPE}_DEBUG}
else
CFLAGS += ${CFLAGS_${HOSTTYPE}}
endif
endif

#extra preprocessor flags
ifneq (,$(strip $(EXTRAFLAGS)))
CFLAGS += ${EXTRAFLAGS}
endif

# Dann machen wir MPI
# Architecture & Dimension
ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
CC      := ${CC_MPI}
LIBS    += ${MPILIBS_${HOSTTYPE}}

# twod?
ifneq (,$(strip $(findstring 2d,${MAKETARGET})))
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES2D} ${MPISOURCES2D}
CFLAGS  += -DTWOD
else
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES3D} ${MPISOURCES3D}
endif

#general compiler flags
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += ${CFLAGS_MPI_${HOSTTYPE}_DEBUG}
else
CFLAGS += ${CFLAGS_MPI_${HOSTTYPE}}
endif
endif

# Und zuletzt machen wir PACX
# Architecture & Dimension
ifneq (,$(strip $(findstring pacx,${MAKETARGET})))
CC      := ${CC_PACX}
LIBS    += ${PACXLIBS_${HOSTTYPE}}

# twod?
ifneq (,$(strip $(findstring 2d,${MAKETARGET})))
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES2D} ${PACXSOURCES2D}
CFLAGS  += -DTWOD
else
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES3D} ${PACXSOURCES3D}
endif

#general compiler flags
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += ${CFLAGS_PACX_${HOSTTYPE}_DEBUG}
else
CFLAGS += ${CFLAGS_PACX_${HOSTTYPE}}
endif
endif

# WAS SOLL DAS HIER ???
ifneq (,$(strip $(findstring gcc,${MAKETARGET})))
CC      := ${CC_GCC}

# twod?
ifneq (,$(strip $(findstring 2d,${MAKETARGET})))
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES2D} ${RISCSOURCES2D}
CFLAGS  += -DTWOD
else
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES3D} ${RISCSOURCES3D}
endif

#general compiler flags
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += ${CFLAGS_GCC_${HOSTTYPE}_DEBUG}
else
CFLAGS += ${CFLAGS_GCC_${HOSTTYPE}}
endif
endif

#socket interface
ifneq (,$(strip $(findstring USE_SOCKETS,${CFLAGS})))
HEADERS += ${SOCKHEADERS}
SOURCES += ${SOCKSOURCES}
endif

#profiling
ifneq (,$(findstring prof,${MAKETARGET}))
CFLAGS += ${CFLAGS_PROF_${HOSTTYPE}}
LIBS   += ${LIBS_PROF_${HOSTTYPE}}
endif

# Ensembles
ifneq (,$(findstring nve,${MAKETARGET}))
CFLAGS += -DNVE
endif

ifneq (,$(findstring nvt,${MAKETARGET}))
CFLAGS += -DNVT
endif

ifneq (,$(findstring and,${MAKETARGET}))
CFLAGS += -DAND
endif

ifneq (,$(findstring npt_axial,${MAKETARGET}))
CFLAGS += -DNPT -DNPT_axial
endif

ifneq (,$(findstring npt_iso,${MAKETARGET}))
CFLAGS += -DNPT -DNPT_iso
endif

ifneq (,$(findstring frac,${MAKETARGET}))
CFLAGS += -DFRAC
SOURCES += ${FRACSOURCES}
endif

ifneq (,$(findstring pull,${MAKETARGET}))
CFLAGS += -DPULL
endif

ifneq (,$(findstring shock,${MAKETARGET}))
CFLAGS += -DSHOCK
endif

ifneq (,$(findstring mik,${MAKETARGET}))
CFLAGS += -DMIK
endif

ifneq (,$(findstring mikshear,${MAKETARGET}))
CFLAGS += -DMIKSHEAR
SOURCES += ${SHEARSOURCES}
endif

ifneq (,$(findstring quasi,${MAKETARGET}))
CFLAGS += -DQUASI
SOURCES += ${QUASISOURCES}
endif

ifneq (,$(findstring disloc,${MAKETARGET}))
CFLAGS += -DDISLOC
endif

# EAM
ifneq (,$(strip $(findstring eam,${MAKETARGET})))
CFLAGS  += -DEAM
SOURCES += ${EAMSOURCES}
endif

# TTBP
ifneq (,$(strip $(findstring ttbp,${MAKETARGET})))
CFLAGS  += -DTTBP
SOURCES += ${TTBPSOURCES}
endif

# Monte Carlo
ifneq (,$(strip $(findstring mc,${MAKETARGET})))
CFLAGS  += -DMC
SOURCES += ${MCSOURCES}
endif

# Correlation
ifneq (,$(strip $(findstring corr,${MAKETARGET})))
CFLAGS  += -DCORRELATE -DCORR_PBC
SOURCES += ${CORRSOURCES}
endif

# mean square displacement
ifneq (,$(strip $(findstring msqd,${MAKETARGET})))
CFLAGS  += -DMSQD -DCORR_PBC
SOURCES += ${CORRSOURCES}
endif

# MONOLJ Case
ifneq (,$(findstring mono,${MAKETARGET}))
CFLAGS += -DMONOLJ
endif

# Single precision
ifneq (,$(findstring single,${MAKETARGET}))
CFLAGS += -DSINGLE
endif

# Substitute .o for .c to get the names of the object files
OBJECTS := $(subst .c,.o,${SOURCES})

###########################################################################
#
#	 Rules
#
###########################################################################

# How to compile *.c files
%.o: %.c ${HEADERS}
	${CC} ${CFLAGS} -c $(subst .o,.c,$@)

# How to link
${MAKETARGET}: ${OBJECTS} ${HEADERS}
	${CC} ${CFLAGS} -o $@ ${OBJECTS} ${LIBS}
#	mv $@ ${HOME}/${BINDIR}/${HOSTTYPE}

# First recursion only set the MAKETARGET Variable
.DEFAULT:
	${MAKE} MAKETARGET='$@' STAGE2

# Second recursion sets all other Variables (depening on MAKETARGET) and compiles
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

###########################################################################
#
#	 Misc. TARGETs
#
###########################################################################

clean:
	rm -f *.o *.u *~ *.V *.T *.O

release:	Makefile ${IMDHEADERS} ${SOURCES3D} ${RISCSOURCES3D} ${MPISOURCES3D} ${PACXSOURCES3D}  ${SOURCES2D} ${RISCSOURCES2D} ${MPISOURCES2D} ${PACXSOURCES2D} ${MISCSOURCES} ${FRACSOURCES} ${SOCKHEADERS} ${SOCKSOURCES}
	cp Makefile         ${RELEASEPATH}
	cp ${IMDHEADERS}    ${RELEASEPATH}
	cp ${SOURCES3D}     ${RELEASEPATH}
	cp ${RISCSOURCES3D} ${RELEASEPATH}
	cp ${MPISOURCES3D}  ${RELEASEPATH}
	cp ${PACXSOURCES3D} ${RELEASEPATH}
	cp ${SOURCES2D}     ${RELEASEPATH}
	cp ${RISCSOURCES2D} ${RELEASEPATH}
	cp ${MPISOURCES2D}  ${RELEASEPATH}
	cp ${PACXSOURCES2D} ${RELEASEPATH}
	cp ${MISCSOURCES}   ${RELEASEPATH}
	cp ${SOCKHEADERS}   ${RELEASEPATH}
	cp ${SOCKSOURCES}   ${RELEASEPATH}
	cp ${FRACSOURCES}   ${RELEASEPATH}
	cp ${QUASISOURCES}   ${RELEASEPATH}

archive:
	tar -zcvf imd-arch.tar.gz Makefile Parameters README *.[ch] *.perl

count:
	wc -l Makefile *.[ch]


