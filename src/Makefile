##############################################################################
#
# IMD Makefile - works only with GNU make
#
# $Revision$
# $Date$
# 
# Usage: make <target>
#
# <target> has the form imd_<arch>_<ensemble>[_<option>[_<option>...]]
#
# <arch> must be one of:
#
# risc	compile for serial execution, with cc where available
# gcc	compile for serial execution, with gcc
# mpi	compile for parallel execution, using MPI
# omp   compile for parallel execution, using OpenMP
# ompi  compile for parallel execution, using OpenMP and MPI
# pacx  compile for parallel execution, using PACX
#
# Machine and compiler specific flags, libraries, etc., are chosen as
# a function of <arch> and the shell variable HOSTTYPE, which must
# be set to a suitable value. Currently supported are:
#
#   HOSTTYPE			OS			<arch>
#
#   i?86-linux, i686		Linux    		risc gcc mpi
#   alpha			True64 Unix             risc gcc mpi
#   iris4d			SGI Irix		risc gcc omp mpi
#   t3e				Cray Unicos/mk		risc mpi pacx
#   SR8000                      HI-UX/MPP               mpi omp ompi
#   SR2201, HITACHI		HI-UX/MPP		mpi
#
# <ensemble> is one of nve, mik, nvt, npt_iso, npt_axial, frac, nvx, stm
#
# <option> is one or several of 2d, debug, prof, mono, binary, ....
#
# A comprehensive list of options can be found below in this makefile.
# Further compiler flags can be handed over with the EXTRAFLAGS variable.
# The makefile analyzes the target name. It sets the preprocessor flags 
# and chooses the source files accordingly.
#
# Examples:
#
# make imd_risc_nve                     ; single CPU 3d nve simulation
# make imd_mpi_frac_debug_2d            ; 2d fracture simlation with debugging
# make imd_risc_nvt_prof                ; 3d NVT simulation with profiling
# make EXTRAFLAGS=-DMYSYMB imd_gcc_mik  ; define preprocessor flag MYSYMB
#
# On some machines EXTRAFLAGS might work only if the shell used by make
# is the same as the shell from which make is invoked (T3E, e.g.).
# In those cases it is possible to set the shell used by make to your
# working shell, by uncommenting and adapting the following line. 
# SHELL         = tcsh


###########################################################################
#
#  The IMD sources
#
###########################################################################

IMDHEADERS      = config.h globals.h imd.h makros.h potaccess.h \
                  prototypes.h types.h 

SOURCES         = imd.c imd_maxwell.c imd_integrate.c imd_misc.c \
	          imd_param.c imd_alloc.c imd_io.c imd_potential.c \
                  imd_time.c imd_generate.c imd_histogram.c
SOURCES2D       = ${SOURCES} imd_main_2d.c imd_io_2d.c imd_geom_2d.c \
		  imd_pictures_2d.c
SOURCES3D       = ${SOURCES} imd_main_3d.c imd_io_3d.c imd_geom_3d.c \
		  imd_pictures_3d.c

RISCSOURCES2D   = imd_main_risc_2d.c
RISCSOURCES3D   = imd_main_risc_3d.c

MPISOURCES      = imd_mpi_util.c
MPISOURCES2D    = ${MPISOURCES} imd_main_mpi_2d.c imd_geom_mpi_2d.c \
		  imd_comm_force_2d.c
MPISOURCES3D    = ${MPISOURCES} imd_main_mpi_3d.c imd_geom_mpi_3d.c \
		  imd_comm_force_3d.c

PAIRSOURCES     = imd_forces.c
EAMSOURCES      = imd_forces_eam.c
EAM2SOURCES     = imd_forces_eam2.c
COVALENTSOURCES = imd_forces.c imd_forces_covalent.c
UNIAXSOURCES    = imd_forces_uniax.c imd_gay_berne.c

DEFORMSOURCES   = imd_deform.c

SOCKHEADERS     = socket_tokens.h sockets.h sockutil.h
SOCKSOURCES     = socket_io.c sockutil.c

FRACSOURCES     = imd_load.c

QUASISOURCES    = imd_qc.c

CORRSOURCES     = imd_correl.c

TRANSSOURCES    = imd_transport.c

SAVEMEMSOURCES  = imd_savemem_3d.c


###########################################################################
#
#  Machine specific flags, libraries, etc.
#
###########################################################################

# PACX MPI
PACXDIR				= ${HOME}/WORK/PACX

# flags for Linux
CC_RISC_i586		= gcc
CFLAGS_RISC_i586	= -DRISC -mcpu=i586 -O3
CFLAGS_RISC_i586_DEBUG	= -DRISC -mcpu=i586 -g
CFLAGS_RISC_i586_PROF	= -g3 -pg

# flags for Linux with gcc (same as RISC)
CC_GCC_i586		= gcc
CFLAGS_GCC_i586		= -DRISC -mcpu=i586 -O3
CFLAGS_GCC_i586_DEBUG	= -DRISC -mcpu=i586 -g
CFLAGS_GCC_i586_PROF	= -g3 -pg

# flags for Linux with LAM/MPI
CC_MPI_i586		= hcc
CFLAGS_MPI_i586		= -DMPI -mcpu=i586 -O3
CFLAGS_MPI_i586_DEBUG	= -DMPI -mcpu=i586 -g
CFLAGS_MPI_i586_PROF	= -g3 -pg
MPILIBS_i586		= -lmpi

# flags for alpha with cc
CC_RISC_alpha		= cc
CFLAGS_RISC_alpha	= -Dalpha -DRISC -O3 -float -fp_reorder -arch ev56 -tune host
CFLAGS_RISC_alpha_DEBUG	= -Dalpha -DRISC -g
CFLAGS_RISC_alpha_PROF	= -g3 -pg

# flags for alpha with gcc
CC_GCC_alpha		= gcc
CFLAGS_GCC_alpha	= -Dalpha -DRISC -O3
CFLAGS_GCC_alpha_DEBUG	= -Dalpha -DRISC -g
CFLAGS_GCC_alpha_PROF	= -g3 -pg

# flags for alpha with OpenMP
#CC_OMP_alpha		= cc.alt
#CFLAGS_OMP_alpha	= -Dalpha -DOMP -mp -O3 -float -fp_reorder -arch ev56 -tune host
#CFLAGS_OMP_alpha_DEBUG	= -Dalpha -DOMP -mp -g
#CFLAGS_OMP_alpha_PROF	= -g3 -pg

# flags for alpha with LAM/MPI
CC_MPI_alpha		= hcc
CFLAGS_MPI_alpha	= -Dalpha -DMPI -O3
CFLAGS_MPI_alpha_DEBUG	= -Dalpha -DMPI -g
MPILIBS_alpha		= -lmpi -lots

# flags for IRIX with cc
CC_RISC_iris4d			= cc
CFLAGS_RISC_iris4d		= -Dsgi -DRISC -O3 -n32 -mips3 -xansi -woff 1174
CFLAGS_RISC_iris4d_DEBUG	= -Dsgi -DRISC -g -xansi -woff 1174
CFLAGS_RISC_iris4d_PROF		= -g3

# flags for IRIX with gcc
CC_GCC_iris4d			= gcc
CFLAGS_GCC_iris4d		= -Dsgi -DRISC -O3
CFLAGS_GCC_iris4d_DEBUG		= -Dsgi -DRISC -g
CFLAGS_GCC_iris4d_PROF		= -g3 -pg

# flags for IRIX with OpenMP
CC_OMP_iris4d			= cc
CFLAGS_OMP_iris4d		= -Dsgi -DOMP -mp -O3 -n32 -mips3 -xansi -woff 1174
CFLAGS_OMP_iris4d_DEBUG		= -Dsgi -DOMP -mp -g -xansi -woff 1174
CFLAGS_OMP_iris4d_PROF		= -g3

# flags for IRIX with LAM/MPI
CC_MPI_iris4d			= hcc
CFLAGS_MPI_iris4d		= -Dsgi -DMPI -O2
MPILIBS_iris4d			= -lmpi

# flags for T3E RISC
CC_RISC_t3e		= cc
CFLAGS_RISC_t3e		= -Dt3e -DRISC -O3 -haggress,report=isf -X1
CFLAGS_RISC_t3e_DEBUG	= -Dt3e -DRISC -g -X1
CFLAGS_RISC_t3e_PROF	= -Gf -happrentice

# flags for T3E MPI
CC_MPI_t3e		= cc
CFLAGS_MPI_t3e		= -Dt3e -DMPI -DAR -O3 -htolerant,aggress,report=isf
CFLAGS_MPI_t3e_DEBUG	= -Dt3e -DMPI -DAR -g
CFLAGS_MPI_t3e_PROF	= -Gf -happrentice
MPILIBS_t3e    	        = -lmpi 

# flags for T3E PACX
CC_PACX_t3e		= cc
CFLAGS_PACX_t3e		= -Dt3e -DPACX -DMPI -I${PACXDIR}/include -DAR -O3 -htolerant,aggress,report=isf
CFLAGS_PACX_t3e_DEBUG	= -Dt3e -DPACX -DMPI -I${PACXDIR}/include -DAR -g
CFLAGS_PACX_t3e_PROF	= -Gf -happrentice
PACXLIBS_t3e		= -L ${PACXDIR}/lib -lpacx -llzo -lmpi

LIBS_PROF_t3e 		:= -lapp

# flags for Hitachi SR2201
CC_MPI_SR2201		 = cc
CFLAGS_MPI_SR2201	 = -DMPI +O4 +Ol +Oprefetch -I/usr/local/mpi/include  -DMPI_hmpp2
CFLAGS_MPI_SR2201_DEBUG	 = -DMPI -g -I/usr/local/mpi/include  -DMPI_hmpp2
CFLAGS_MPI_SR2201_PROF   = -G
MPILIBS_SR2201		 = -L/usr/local/mpi/lib/hmpp2/cml -lmpi

# flags for Hitachi SR8000 MPI
CC_MPI_SR8000		 = cc
CFLAGS_MPI_SR8000	 = -DMPI +O4 +Ol -I/usr/mpi/include
CFLAGS_MPI_SR8000_DEBUG	 = -DMPI -g  -I/usr/mpi/include
CFLAGS_MPI_SR8000_PROF   = -G
MPILIBS_SR8000		 = -L/usr/mpi/lib -lmpi

# flags for Hitachi SR8000 OpenMP
CC_OMP_SR8000		 = cc
CFLAGS_OMP_SR8000	 = -DOMP -omp -par +O4 +Ol
CFLAGS_OMP_SR8000_DEBUG	 = -DOMP -omp -par=1 -O2 -g
CFLAGS_OMP_SR8000_PROF   = -G

# flags for Hitachi SR8000 MPI + OpenMP
CC_OMPI_SR8000		 = cc
CFLAGS_OMPI_SR8000	 = -DOMPI -omp -par +O4 +Ol -I/usr/mpi/include
CFLAGS_OMPI_SR8000_DEBUG = -DOMPI -omp -par +O2 -g -I/usr/mpi/include
CFLAGS_OMPI_SR8000_PROF  = -G


###########################################################################
#
#  IMD host type
#
###########################################################################

# the default for IMDHOST is the shell variable HOSTTYPE
IMDHOST=${HOSTTYPE}

# deal with some special cases:
ifneq (,$(strip $(findstring linux,${HOSTTYPE})))
IMDHOST = i586
endif
ifneq (,$(strip $(findstring i686,${HOSTTYPE})))
IMDHOST = i586
endif
ifneq (,$(strip $(findstring cray,${HOSTTYPE})))
IMDHOST = t3e
endif
ifneq (,$(strip $(findstring HITACHI,${HOSTTYPE})))
IMDHOST = SR2201
endif

CC_RISC           = ${CC_RISC_${IMDHOST}}
CC_GCC		  = ${CC_GCC_${IMDHOST}}
CC_MPI		  = ${CC_MPI_${IMDHOST}}
CC_OMP		  = ${CC_OMP_${IMDHOST}}
CC_OMPI		  = ${CC_OMPI_${IMDHOST}}
CC_PACX		  = ${CC_PACX_${IMDHOST}}

# to link or not to link
LIBS    	  := -lm
MPILIBS		  = ${MPILIBS_${IMDHOST}}	
PACXLIBS	  = ${PACXLIBS_${IMDHOST}}	

# Executeables are put into $HOME/$BINDIR/$HOSTTYPE
BINDIR   = bin


#########################################################
#
# IMD Configuration rules
#
#########################################################

# Architecture & Dimension

###  RISC  ##############################################

ifneq (,$(strip $(findstring risc,${MAKETARGET})))
CC      := ${CC_RISC}

# twod or not twod
ifneq (,$(strip $(findstring 2d,${MAKETARGET})))
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES2D} ${RISCSOURCES2D}
CFLAGS  += -DTWOD
else
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES3D} ${RISCSOURCES3D}
endif

# general compiler flags
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += ${CFLAGS_RISC_${IMDHOST}_DEBUG}
else
CFLAGS += ${CFLAGS_RISC_${IMDHOST}}
endif
ifneq (,$(findstring prof,${MAKETARGET}))
CFLAGS += ${CFLAGS_RISC_${IMDHOST}_PROF}
endif

endif


###  GCC  ###############################################

ifneq (,$(strip $(findstring gcc,${MAKETARGET})))
CC      := ${CC_GCC}

# twod or not twod
ifneq (,$(strip $(findstring 2d,${MAKETARGET})))
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES2D} ${RISCSOURCES2D}
CFLAGS  += -DTWOD
else
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES3D} ${RISCSOURCES3D}
endif

# general compiler flags
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += ${CFLAGS_GCC_${IMDHOST}_DEBUG}
else
CFLAGS += ${CFLAGS_GCC_${IMDHOST}}
endif
ifneq (,$(findstring prof,${MAKETARGET}))
CFLAGS += ${CFLAGS_GCC_${IMDHOST}_PROF}
endif

endif


###  MPI and OpenMP #####################################

ifneq (,$(strip $(findstring ompi,${MAKETARGET})))
CC      := ${CC_OMPI}
LIBS    += ${MPILIBS_${IMDHOST}}

# twod or not twod
ifneq (,$(strip $(findstring 2d,${MAKETARGET})))
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES2D} ${MPISOURCES2D}
CFLAGS  += -DTWOD
else
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES3D} ${MPISOURCES3D}
endif

# general compiler flags
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += ${CFLAGS_OMPI_${IMDHOST}_DEBUG}
else
CFLAGS += ${CFLAGS_OMPI_${IMDHOST}}
endif
ifneq (,$(findstring prof,${MAKETARGET}))
CFLAGS += ${CFLAGS_OMPI_${IMDHOST}_PROF}
LIBS   += ${LIBS_PROF_${IMDHOST}}
endif

else


###  MPI  ###############################################

ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
CC      := ${CC_MPI}
LIBS    += ${MPILIBS_${IMDHOST}}

# twod or not twod
ifneq (,$(strip $(findstring 2d,${MAKETARGET})))
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES2D} ${MPISOURCES2D}
CFLAGS  += -DTWOD
else
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES3D} ${MPISOURCES3D}
endif

# general compiler flags
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += ${CFLAGS_MPI_${IMDHOST}_DEBUG}
else
CFLAGS += ${CFLAGS_MPI_${IMDHOST}}
endif
ifneq (,$(findstring prof,${MAKETARGET}))
CFLAGS += ${CFLAGS_MPI_${IMDHOST}_PROF}
LIBS   += ${LIBS_PROF_${IMDHOST}}
endif

endif


###  OpenMP  ############################################

ifneq (,$(strip $(findstring omp,${MAKETARGET})))
CC      := ${CC_OMP}

# twod or not twod
ifneq (,$(strip $(findstring 2d,${MAKETARGET})))
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES2D} ${RISCSOURCES2D}
CFLAGS  += -DTWOD
else
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES3D} ${RISCSOURCES3D}
endif

# general compiler flags
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += ${CFLAGS_OMP_${IMDHOST}_DEBUG}
else
CFLAGS += ${CFLAGS_OMP_${IMDHOST}}
endif
ifneq (,$(findstring prof,${MAKETARGET}))
CFLAGS += ${CFLAGS_OMP_${IMDHOST}_PROF}
endif

endif

endif


###  PACX  ##############################################

ifneq (,$(strip $(findstring pacx,${MAKETARGET})))
CC      := ${CC_PACX}
LIBS    += ${PACXLIBS_${IMDHOST}}

# twod or not twod
ifneq (,$(strip $(findstring 2d,${MAKETARGET})))
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES2D} ${MPISOURCES2D}
CFLAGS  += -DTWOD
else
HEADERS := ${IMDHEADERS}
SOURCES	:= ${SOURCES3D} ${MPISOURCES3D}
endif

# general compiler flags
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += ${CFLAGS_PACX_${IMDHOST}_DEBUG}
else
CFLAGS += ${CFLAGS_PACX_${IMDHOST}}
endif
ifneq (,$(findstring prof,${MAKETARGET}))
CFLAGS += ${CFLAGS_PACX_${IMDHOST}_PROF}
endif

endif


###  INTERACTION  #######################################

# pair interaction is the default
FORCESOURCES = ${PAIRSOURCES}

# EAM2 or EAM
ifneq (,$(strip $(findstring eam2,${MAKETARGET})))
CFLAGS  += -DEAM2
FORCESOURCES = ${EAM2SOURCES}
else
ifneq (,$(strip $(findstring eam,${MAKETARGET})))
CFLAGS  += -DEAM
CFLAGS  += -UAR
FORCESOURCES = ${EAMSOURCES}
endif
endif

# TTBP
ifneq (,$(strip $(findstring ttbp,${MAKETARGET})))
CFLAGS  += -DTTBP
FORCESOURCES = ${COVALENTSOURCES}
endif

# TERSOFF
ifneq (,$(strip $(findstring tersoff,${MAKETARGET})))
CFLAGS  += -DTERSOFF
FORCESOURCES = ${COVALENTSOURCES}
endif

# UNIAX
ifneq (,$(strip $(findstring uniax,${MAKETARGET})))
CFLAGS  += -DUNIAX
FORCESOURCES = ${UNIAXSOURCES}
endif

SOURCES += ${FORCESOURCES}


###  ENSEMBLES  #########################################

ifneq (,$(findstring nve,${MAKETARGET}))
CFLAGS += -DNVE
endif

ifneq (,$(findstring mik,${MAKETARGET}))
CFLAGS += -DMIK
endif

ifneq (,$(findstring nvt,${MAKETARGET}))
CFLAGS += -DNVT
endif

ifneq (,$(findstring npt_iso,${MAKETARGET}))
CFLAGS += -DNPT -DNPT_iso
endif

ifneq (,$(findstring npt_axial,${MAKETARGET}))
CFLAGS += -DNPT -DNPT_axial
endif

ifneq (,$(findstring frac,${MAKETARGET}))
CFLAGS += -DFRAC
SOURCES += ${FRACSOURCES}
endif

ifneq (,$(findstring nvx,${MAKETARGET}))
CFLAGS += -DNVX
SOURCES += ${TRANSSOURCES}
endif

ifneq (,$(findstring stm,${MAKETARGET}))
CFLAGS += -DSTM
endif

###  OPTIONS  ############################################

# extra preprocessor flags
ifneq (,$(strip $(EXTRAFLAGS)))
CFLAGS += ${EXTRAFLAGS}
endif

# socket interface
ifneq (,$(strip $(findstring sock,${MAKETARGET})))
HEADERS += ${SOCKHEADERS}
SOURCES += ${SOCKSOURCES}
CFLAGS  += -DUSE_SOCKETS
endif

# timing
ifneq (,$(findstring timing,${MAKETARGET}))
CFLAGS += -DTIMING
endif

ifneq (,$(findstring and,${MAKETARGET}))
CFLAGS += -DAND
endif

ifneq (,$(findstring fbc,${MAKETARGET}))
CFLAGS += -DFBC
endif

ifneq (,$(findstring sr,${MAKETARGET}))
CFLAGS += -DSR
endif

ifneq (,$(findstring fnorm,${MAKETARGET}))
CFLAGS += -DFNORM
endif

ifneq (,$(findstring glok,${MAKETARGET}))
CFLAGS += -DGLOK -DFNORM
endif

ifneq (,$(findstring efilter,${MAKETARGET}))
CFLAGS += -DEFILTER
endif

ifneq (,$(findstring deform,${MAKETARGET}))
CFLAGS += -DDEFORM
SOURCES += ${DEFORMSOURCES}
endif

ifneq (,$(findstring homdef,${MAKETARGET}))
CFLAGS += -DHOMDEF
SOURCES += ${DEFORMSOURCES}
endif

ifneq (,$(findstring shock,${MAKETARGET}))
CFLAGS += -DSHOCK
endif

ifneq (,$(findstring stress,${MAKETARGET}))
CFLAGS += -DSTRESS_TENS
endif

ifneq (,$(findstring quasi,${MAKETARGET}))
CFLAGS += -DQUASI
SOURCES += ${QUASISOURCES}
endif

ifneq (,$(findstring disloc,${MAKETARGET}))
CFLAGS += -DDISLOC
endif

ifneq (,$(findstring sllod,${MAKETARGET}))
CFLAGS += -DSLLOD -DNVT
endif

ifneq (,$(findstring avpos,${MAKETARGET}))
CFLAGS += -DAVPOS
endif

ifneq (,$(findstring ordpar,${MAKETARGET}))
CFLAGS += -DORDPAR
endif

# Correlation
ifneq (,$(strip $(findstring corr,${MAKETARGET})))
CFLAGS  += -DCORRELATE
SOURCES += ${CORRSOURCES}
endif

# RNEMD heat transport
ifneq (,$(strip $(findstring rnemd,${MAKETARGET})))
CFLAGS  += -DRNEMD
SOURCES += ${TRANSSOURCES}
endif

# mean square displacement
ifneq (,$(strip $(findstring msqd,${MAKETARGET})))
CFLAGS  += -DMSQD
SOURCES += ${CORRSOURCES}
endif

# MONOLJ Case
ifneq (,$(findstring monolj,${MAKETARGET}))
CFLAGS += -DMONOLJ
ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
SOURCES += ${SAVEMEMSOURCES}
endif
endif

# Single precision
ifneq (,$(findstring single,${MAKETARGET}))
CFLAGS += -DSINGLE
endif

# monoatomic system (performance tweak)
ifneq (,$(findstring mono,${MAKETARGET}))
CFLAGS += -DMONO
endif

# binary system (performance tweak)
ifneq (,$(findstring binary,${MAKETARGET}))
CFLAGS += -DBINARY
endif

# high precision output (checkpoints)
ifneq (,$(findstring hpo,${MAKETARGET}))
CFLAGS += -DHPO
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
	${CC} ${CFLAGS} -o $@ ${OBJECTS} ${LIBS} ${ONLINELIBS}
	mv $@ ${HOME}/${BINDIR}/${HOSTTYPE}


# First recursion only set the MAKETARGET Variable
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

###########################################################################
#
#	 Misc. TARGETs
#
###########################################################################

clean:
	rm -f *.o *.u *~ \#* *.V *.T *.O

archive:
	tar -zcvf imd-arch.tar.gz Makefile Parameters README *.[ch] *.perl

count:
	wc -l Makefile *.[ch]

help:
	echo "Usage: gmake imd_<ensemble>_<architeture>[_<suffix>...]"
