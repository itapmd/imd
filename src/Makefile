#############################################################################
#
# IMD -- The ITAP Molecular Dynamics Program
#
# Copyright 1996-2001 Institute for Theoretical and Applied Physics,
# University of Stuttgart, D-70550 Stuttgart
#
#
# IMD Makefile - works only with GNU make
#
# $Revision$
# $Date$
# 
# Usage: make <target>
#
# <target> has the form 
#
#    imd[_<comp>][_<parallel>]_<ensemble>[_<option>[_<option>...]]
#
# The parallelization method <parallel> can be one of:
#
#    mpi   compile for parallel execution, using MPI
#    omp   compile for parallel execution, using OpenMP
#    ompi  compile for parallel execution, using OpenMP and MPI
#    pacx  compile for parallel execution, using PACX
#
# If no parallelization method is specified, the default is
# "serial". For each machine-parallelization method combination
# (including serial), a default compiler is defined. If a different
# compiler is desired, its name can be specified as the <comp> 
# componend in the make target. For each compiler-parallelization
# method-machine combination, a suitable set of compiler flags can
# be defined. This is already done for many such combinations.
# The machine type is derived from the shell variable HOSTTYPE.
# Currently supported machine-parallelization combinations are:
#
#   HOSTTYPE			OS		Parallelization
#
#   i?86-linux,i686,i386	Linux    	serial mpi omp
#   alpha		 	True64 Unix     serial mpi omp
#   iris4d		 	SGI Irix	serial mpi omp
#   t3e			 	Cray Unicos/mk	serial mpi pacx
#   SR8000               	HI-UX/MPP       mpi omp ompi
#   SR2201, HITACHI	 	HI-UX/MPP	mpi
#   rs6000			AIX             serial mpi
#
# <ensemble> is one or several of nve, mik, nvt, npt_iso, npt_axial, 
# frac, nvx, ftg, stm
#
# <option> is one or several of 2d, debug, prof, mono, binary, ....
#
# The executable, which has the name <target>, is moved to the directory
# $HOME/bin/$HOSTTYPE, which must exist, and which should be in the
# search $PATH.
#
# A comprehensive list of options can be found below in this makefile.
# Further compiler flags can be handed over with the EXTRAFLAGS variable.
# The makefile analyzes the target name. It sets the preprocessor flags 
# and chooses the source files accordingly.
#
# Examples:
#
# make imd_nve                          ; single CPU 3d nve simulation
# make imd_mpi_frac_debug_2d            ; 2d fracture simlation with debugging
# make imd_nvt_prof                     ; 3d NVT simulation with profiling
# make EXTRAFLAGS=-DMYSYMB imd_gcc_mik  ; define preprocessor flag MYSYMB,
#                                       ; compile with gcc
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
EAM2SOURCES     = imd_forces_eam2.c
CGSOURCES	= imd_main_cg_3d.c
COVALENTSOURCES = imd_forces_covalent.c
UNIAXSOURCES    = imd_forces_uniax.c imd_gay_berne.c
EWALDSOURCES    = imd_forces_ewald.c

EPITAXSOURCES   = imd_epitax.c

DEFORMSOURCES   = imd_deform.c

SOCKHEADERS     = socket_tokens.h sockets.h sockutil.h
SOCKSOURCES     = socket_io.c sockutil.c

QUASISOURCES    = imd_qc.c

CORRSOURCES     = imd_correl.c

TRANSSOURCES    = imd_transport.c


###########################################################################
#
# flags for Linux
#
###########################################################################

# default compilers and libraries
# if available, we take the Intel compiler
ifneq (,$(strip $(findstring ia32,${PATH})))
CC_DEFAULT_SERIAL_i586		= icc
CC_DEFAULT_OMP_i586		= icc
CC_DEFAULT_OMPI_i586		= mpicc
else
CC_DEFAULT_SERIAL_i586		= gcc
CC_DEFAULT_OMP_i586		= guidec
endif
CC_DEFAULT_MPI_i586		= mpicc

# gcc, serial
#CFLAGS_gcc_SERIAL_i586		= -O6 -funroll-all-loops -fomit-frame-pointer -march=pentiumpro -static -mfpmath=sse -msse -mmmx
CFLAGS_gcc_SERIAL_i586		= -O -funroll-all-loops -fomit-frame-pointer -march=pentiumpro -static -mfpmath=sse -msse #-mmmx
#CFLAGS_gcc_SERIAL_i586		= -mcpu=i686 -O3
CFLAGS_gcc_SERIAL_i586_DEBUG	= -mcpu=i686 -g
CFLAGS_gcc_SERIAL_i586_PROF	= -g3 -pg

# pgcc, serial
CFLAGS_pgcc_SERIAL_i586        	= -fast -tp p6 -Mvect=sse
CFLAGS_pgcc_SERIAL_i586_DEBUG  	= -g -Xa
CFLAGS_pgcc_SERIAL_i586_PROF   	= -Mprof=func

# icc, serial
CFLAGS_icc_SERIAL_i586        	= -O -ip -xK
CFLAGS_icc_SERIAL_i586_PROF	= -prof_gen

# icc, OpenMP
CFLAGS_icc_OMP_i586        	= -O -ip -openmp -static 
CFLAGS_icc_OMP_i586_PROF	= -prof_gen

# kcc, serial
CFLAGS_KCC_SERIAL_i586       	= +K3 -O3 --diag_suppress 177,550

# MPICH
# there are many mpicc; if available, we take the Intel compiler
ifneq (,$(strip $(findstring ia32,${PATH})))
#CFLAGS_mpicc_MPI_i586		= -O -ip -tpp7 -static -L. -lgm
CFLAGS_mpicc_MPI_i586		= -O -ip -static 
CFLAGS_mpicc_MPI_i586_DEBUG	= -g
CFLAGS_mpicc_MPI_i586_PROF	= -prof_gen
else
CFLAGS_mpicc_MPI_i586		= -O3
CFLAGS_mpicc_MPI_i586_DEBUG	= -g
CFLAGS_mpicc_MPI_i586_PROF	= -g3 -pg
endif

# MPICH + OpenMP
# there are many mpicc; if available, we take the Intel compiler
ifneq (,$(strip $(findstring ia32,${PATH})))
CFLAGS_mpicc_OMPI_i586		= -O -ip -openmp -static
CFLAGS_mpicc_OMPI_i586_DEBUG	= -g
CFLAGS_mpicc_OMPI_i586_PROF	= -prof_gen
#else
#CFLAGS_mpicc_OMPI_i586		= -O3
#CFLAGS_mpicc_OMPI_i586_DEBUG	= -g
#CFLAGS_mpicc_OMPI_i586_PROF	= -g3 -pg
endif

# guidec, OpenMP
CFLAGS_guidec_OMP_i586         	= +K3 -O3 --diag_suppress 177,550
CFLAGS_guidec_OMP_i586_DEBUG   	= -g +K0 -WGstrict --diag_suppress 177,550
CFLAGS_guidec_OMP_i586_PROF    	= -p

# omni, OpenMP
CFLAGS_omcc_OMP_i586         	= -O
CFLAGS_omcc_OMP_i586_DEBUG   	= -g

# pgcc, OpenMP
CFLAGS_pgcc_OMP_i586         	= -fast
CFLAGS_pgcc_OMP_i586_DEBUG   	= -g -Xa
CFLAGS_pgcc_OMP_i586_PROF    	= -fast -Mprof=func

###########################################################################
#
# flags for Linux IA64
#
###########################################################################

# default compilers and libraries
CC_DEFAULT_SERIAL_ia64		= ecc
CC_DEFAULT_OMP_ia64		= ecc
#CC_DEFAULT_OMPI_ia64		= mpicc

# ecc, serial
CFLAGS_ecc_SERIAL_ia64        	= -O -ipo -static
CFLAGS_ecc_SERIAL_ia64_PROF	= -prof_gen

# ecc, OpenMP
CFLAGS_ecc_OMP_ia64        	= -O -ipo -static -openmp
CFLAGS_ecc_OMP_ia64_PROF	= -prof_gen

# gcc, serial
CFLAGS_gcc_SERIAL_ia64        	= -O3
CFLAGS_gcc_SERIAL_ia64_PROF	= -pg

###########################################################################
#
# flags for alpha
#
###########################################################################

# default compilers and libraries
CC_DEFAULT_SERIAL_alpha		= cc
CC_DEFAULT_OMP_alpha		= cc.alt
CC_DEFAULT_MPI_alpha		= hcc
MPILIBS_alpha			= -lmpi -lots

# cc, serial
CFLAGS_cc_SERIAL_alpha		= -DALPHA -O3 -float -fp_reorder -arch ev56 -tune host
CFLAGS_cc_SERIAL_alpha_DEBUG	= -DALPHA -g
CFLAGS_cc_SERIAL_alpha_PROF	= -g3 -pg

# gcc, serial
CFLAGS_gcc_SERIAL_alpha		= -DALPHA -O3
CFLAGS_gcc_SERIAL_alpha_DEBUG	= -DALPHA -g
CFLAGS_gcc_SERIAL_alpha_PROF	= -g3 -pg

# OpenMP
CFLAGS_cc.alt_OMP_alpha		= -DALPHA -mp -O3 -float -fp_reorder -arch ev56 -tune host
CFLAGS_cc.alt_OMP_alpha_DEBUG	= -DALPHA -mp -g
CFLAGS_cc.alt_OMP_alpha_PROF	= -g3 -pg

# LAM/MPI
CFLAGS_hcc_MPI_alpha		= -DALPHA -O3
CFLAGS_hcc_MPI_alpha_DEBUG	= -DALPHA -g


###########################################################################
#
# flags for IRIX
#
###########################################################################

# default compilers and libraries
CC_DEFAULT_SERIAL_iris4d	= cc
CC_DEFAULT_OMP_iris4d		= cc

# cc
CFLAGS_cc_SERIAL_iris4d		= -Dsgi -O3 -n32 -mips3 -xansi -woff 1174
CFLAGS_cc_SERIAL_iris4d_DEBUG	= -Dsgi -g -xansi -woff 1174
CFLAGS_cc_SERIAL_iris4d_PROF	= -g3

# gcc
CFLAGS_gcc_SERIAL_iris4d	= -Dsgi -O3
CFLAGS_gcc_SERIAL_iris4d_DEBUG	= -Dsgi -g
CFLAGS_gcc_SERIAL_iris4d_PROF	= -g3 -pg

# OpenMP
CFLAGS_cc_OMP_iris4d		= -Dsgi -mp -O3 -n32 -mips3 -xansi -woff 1174
CFLAGS_cc_OMP_iris4d_DEBUG	= -Dsgi -mp -g -xansi -woff 1174
CFLAGS_cc_OMP_iris4d_PROF	= -g3


###########################################################################
#
# flags for sparc
#
###########################################################################

# default compilers and libraries
CC_DEFAULT_SERIAL_sparc		= cc
CC_DEFAULT_OMP_sparc		= cc
CC_DEFAULT_MPI_sparc		= mpcc
MPILIBS_sparc			= -lmpi

# cc
CFLAGS_cc_SERIAL_sparc		= -fast -xtarget=ultra3 -xarch=v9b
CFLAGS_cc_SERIAL_sparc_DEBUG	= -g -xO3
CFLAGS_cc_SERIAL_sparc_PROF	= -p

# gcc - slow
CFLAGS_gcc_SERIAL_sparc		= -O3
CFLAGS_gcc_SERIAL_sparc_DEBUG	= -g
CFLAGS_gcc_SERIAL_sparc_PROF	= -g3 -pg

# cc
CFLAGS_cc_OMP_sparc		= -fast -xopenmp -xtarget=ultra3 -xarch=v9b
CFLAGS_cc_OMP_sparc_DEBUG	= -g -xopenmp -xO3
CFLAGS_cc_OMP_sparc_PROF	= -p

# MPI
CFLAGS_mpcc_MPI_sparc		= -fast -xtarget=ultra3 -xarch=v9b
CFLAGS_mpcc_MPI_sparc_DEBUG	= -g -xO3
CFLAGS_mpcc_MPI_sparc_PROF	= -p


###########################################################################
#
# flags for T3E
#
###########################################################################

# default compilers and libraries
CC_DEFAULT_SERIAL_t3e		= cc
CC_DEFAULT_MPI_t3e		= cc
CC_DEFAULT_PACX_t3e		= cc
MPILIBS_t3e	    	        = -lmpi 
PACXDIR				= ${HOME}/WORK/PACX
PACXLIBS_t3e			= -L ${PACXDIR}/lib -lpacx -llzo -lmpi
LIBS_PROF_t3e	 		= -lapp

# serial
CFLAGS_cc_SERIAL_t3e		= -Dt3e -O3 -haggress,report=isf -X1
CFLAGS_cc_SERIAL_t3e_DEBUG	= -Dt3e -g -X1
CFLAGS_cc_SERIAL_t3e_PROF	= -Gf -happrentice

# MPI
CFLAGS_cc_MPI_t3e		= -Dt3e -O3 -htolerant,aggress,report=isf
CFLAGS_cc_MPI_t3e_DEBUG		= -Dt3e -g
CFLAGS_cc_MPI_t3e_PROF		= -Gf -happrentice

# PACX
CFLAGS_cc_PACX_t3e		= -Dt3e -I${PACXDIR}/include -O3 -htolerant,aggress,report=isf
CFLAGS_cc_PACX_t3e_DEBUG	= -Dt3e -I${PACXDIR}/include -g
CFLAGS_cc_PACX_t3e_PROF		= -Gf -happrentice

###########################################################################
#
# flags for Hitachi SR2201
#
###########################################################################

# default compilers and libraries
CC_DEFAULT_MPI_SR2201	 	= cc
MPILIBS_SR2201		 	= -L/usr/local/mpi/lib/hmpp2/cml -lmpi

# MPI
CFLAGS_cc_MPI_SR2201	 	= +O4 +Ol +Oprefetch -I/usr/local/mpi/include  -DMPI_hmpp2
CFLAGS_cc_MPI_SR2201_DEBUG   	= -g -I/usr/local/mpi/include -DMPI_hmpp2
CFLAGS_cc_MPI_SR2201_PROF  	= -G

###########################################################################
#
# flags for Hitachi SR8000
#
###########################################################################

# default compilers and libraries
CC_DEFAULT_MPI_SR8000		= mpicc
CC_DEFAULT_OMP_SR8000		= cc
CC_DEFAULT_OMPI_SR8000		= mpicc
LIBS_PROF_SR8000		= -lpl

# MPI
CFLAGS_mpicc_MPI_SR8000        	= -O4 +Op -noparallel -msg=e
CFLAGS_mpicc_MPI_SR8000_DEBUG  	= -g
CFLAGS_mpicc_MPI_SR8000_PROF   	= -Xfuncmonitor

# OpenMP
CFLAGS_cc_OMP_SR8000	 	= -omp -par -O4 +Op -msg=e
CFLAGS_cc_OMP_SR8000_DEBUG	= -omp -par=1 -O2 -g
CFLAGS_cc_OMP_SR8000_PROF   	= -Xfuncmonitor

# MPI + OpenMP
CFLAGS_mpicc_OMPI_SR8000	= -omp -par -O4 +Op -msg=e
CFLAGS_mpicc_OMPI_SR8000_DEBUG 	= -omp -par -O2 -g
CFLAGS_mpicc_OMPI_SR8000_PROF  	= -Xfuncmonitor

# MPI (crosscompiler linux -> SR8000)
CFLAGS_xmpicc_MPI_i586        	= -O4 +Op -noparallel -msg=e
CFLAGS_xmpicc_MPI_i586_DEBUG  	= -g
CFLAGS_xmpicc_MPI_i586_PROF   	= -Xfuncmonitor

# OpenMP (crosscompiler linux -> SR8000)
CFLAGS_xcc_OMP_i586	 	= -omp -par -O4 +Op -msg=e
CFLAGS_xcc_OMP_i586_DEBUG	= -omp -par=1 -O2 -g
CFLAGS_xcc_OMP_i586_PROF   	= -Xfuncmonitor

# MPI + OpenMP (crosscompiler linux -> SR8000)
CFLAGS_xmpicc_OMPI_i586		= -omp -par -O4 +Op -msg=e
CFLAGS_xmpicc_OMPI_i586_DEBUG 	= -omp -par -O2 -g
CFLAGS_xmpicc_OMPI_i586_PROF  	= -Xfuncmonitor

###########################################################################
#
# flags for IBM SP-SMP
#
###########################################################################

# default compilers and libraries
CC_DEFAULT_SERIAL_rs6000	= cc
CC_DEFAULT_MPI_rs6000		= mpcc
CC_DEFAULT_OMP_rs6000		= xlc_r
CC_DEFAULT_OMPI_rs6000		= mpcc_r #?
LIBS_PROF_rs6000		= -lpl   #?

# SERIAL
CFLAGS_cc_SERIAL_rs6000        	= -O4 -qtune=auto  # -O4 automatically invokes -qarch=auto

# MPI
CFLAGS_mpcc_MPI_rs6000        	= -O4 -qtune=auto  # -O4 automatically invokes -qarch=auto
CFLAGS_mpcc_MPI_rs6000_DEBUG  	= -g 
CFLAGS_mpcc_MPI_rs6000_PROF   	= -p 

# OpenMP
#CFLAGS_xlc_OMP_rs6000	 	= -O4 -qtune=auto -qsmp=omp 
#CFLAGS_cc_OMP_SR8000_DEBUG	= -omp -par=1 -O2 -g
#CFLAGS_cc_OMP_SR8000_PROF   	= -Xfuncmonitor

# MPI + OpenMP
#CFLAGS_mpicc_OMPI_SR8000	= -omp -par -O4 +Op -msg=e
#CFLAGS_mpicc_OMPI_SR8000_DEBUG 	= -omp -par -O2 -g
#CFLAGS_mpicc_OMPI_SR8000_PROF  	= -Xfuncmonitor


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
ifneq (,$(strip $(findstring i386,${HOSTTYPE})))
IMDHOST = i586
endif
ifneq (,$(strip $(findstring cray,${HOSTTYPE})))
IMDHOST = t3e
endif
ifneq (,$(strip $(findstring HITACHI,${HOSTTYPE})))
IMDHOST = SR2201
endif


###########################################################################
#
#  Parallelization method
#
###########################################################################

# default is serial
PARALLEL = SERIAL
# MPI
ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
PARALLEL = MPI
CFLAGS  += -DMPI
LIBS    += ${MPILIBS_${IMDHOST}}
endif
# OpenMP
ifneq (,$(strip $(findstring omp,${MAKETARGET})))
PARALLEL = OMP
CFLAGS  += -DOMP
endif
# MPI + OpenMP
ifneq (,$(strip $(findstring ompi,${MAKETARGET})))
PARALLEL = OMPI
endif
# PACX
ifneq (,$(strip $(findstring pacx,${MAKETARGET})))
PARALLEL = PACX
CFLAGS  += -DMPI -DPACX
LIBS    += PACXLIBS_${IMDHOST}
endif


###########################################################################
#
#  Compiler
#
###########################################################################

CC = ${CC_DEFAULT_${PARALLEL}_${IMDHOST}}
# gcc
ifneq (,$(strip $(findstring gcc,${MAKETARGET})))
CC = gcc
endif
# pgcc
ifneq (,$(strip $(findstring pgcc,${MAKETARGET})))
CC = pgcc
endif
# icc
ifneq (,$(strip $(findstring icc,${MAKETARGET})))
CC = icc
endif
# ecc
ifneq (,$(strip $(findstring ecc,${MAKETARGET})))
CC = ecc
endif
# kcc
ifneq (,$(strip $(findstring kcc,${MAKETARGET})))
CC = KCC
endif
# guidec
ifneq (,$(strip $(findstring guidec,${MAKETARGET})))
CC = guidec
endif
# omni
ifneq (,$(strip $(findstring omcc,${MAKETARGET})))
CC = omcc
endif
# SR8000 crosscompiler
ifneq (,$(strip $(findstring sr8k,${MAKETARGET})))
ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
CC = xmpicc
else
CC = xcc
endif
endif


###########################################################################
#
#  Libraries
#
###########################################################################

LIBS    	 += -lm
MPILIBS		  = ${MPILIBS_${IMDHOST}}	
PACXLIBS	  = ${PACXLIBS_${IMDHOST}}	

# where to put executable
ifneq (,$(strip $(findstring sr8k,${MAKETARGET})))
BINDIR = hwwfs1:sr8k/bin/SR8000
MV = scp
else
BINDIR   = ${HOME}/bin/${HOSTTYPE}
MV =mv
endif


#########################################################
#
# IMD Configuration rules
#
#########################################################

HEADERS := ${IMDHEADERS}
# twod or not twod
ifneq (,$(strip $(findstring 2d,${MAKETARGET})))
# serial or mpi
ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
SOURCES	:= ${SOURCES2D} ${MPISOURCES2D}
else
SOURCES	:= ${SOURCES2D} ${RISCSOURCES2D}
endif
CFLAGS  += -DTWOD
else
# serial or mpi
ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
SOURCES	:= ${SOURCES3D} ${MPISOURCES3D}
else
SOURCES	:= ${SOURCES3D} ${RISCSOURCES3D}
endif
endif

# general compiler flags, debug and profile
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += ${CFLAGS_${CC}_${PARALLEL}_${IMDHOST}_DEBUG}
else
CFLAGS += ${CFLAGS_${CC}_${PARALLEL}_${IMDHOST}}
endif
ifneq (,$(findstring prof,${MAKETARGET}))
CFLAGS += ${CFLAGS_${CC}_${PARALLEL}_${IMDHOST}_PROF}
LIBS   += ${LIBS_PROF_${IMDHOST}}
endif


###  INTERACTION  #######################################

# pair interaction is the default
FORCESOURCES = ${PAIRSOURCES}

# EAM2 or EAM  -  this is now the same
ifneq (,$(strip $(findstring eam,${MAKETARGET})))
CFLAGS  += -DEAM2
ifneq (,$(strip $(findstring asympot,${MAKETARGET})))
CFLAGS += -DASYMPOT
FORCESOURCES = ${EAM2SOURCES}
else
FORCESOURCES += ${EAM2SOURCES}
endif
endif


# TTBP
ifneq (,$(strip $(findstring ttbp,${MAKETARGET})))
CFLAGS  += -DTTBP
FORCESOURCES += ${COVALENTSOURCES}
endif

# STIWEB
ifneq (,$(strip $(findstring stiweb,${MAKETARGET})))
CFLAGS  += -DSTIWEB
FORCESOURCES += ${COVALENTSOURCES}
endif

# TERSOFF
ifneq (,$(strip $(findstring tersoff,${MAKETARGET})))
CFLAGS  += -DTERSOFF
FORCESOURCES += ${COVALENTSOURCES}
endif

# KEATING
ifneq (,$(strip $(findstring keating,${MAKETARGET})))
CFLAGS  += -DKEATING
FORCESOURCES += ${COVALENTSOURCES}
endif

# UNIAX
ifneq (,$(strip $(findstring uniax,${MAKETARGET})))
CFLAGS  += -DUNIAX
FORCESOURCES = ${UNIAXSOURCES}
endif

# EWALD
ifneq (,$(strip $(findstring ewald,${MAKETARGET})))
CFLAGS  += -DEWALD
FORCESOURCES  += ${EWALDSOURCES}
endif

# SPRING
ifneq (,$(strip $(findstring spring,${MAKETARGET})))
CFLAGS  += -DSPRING
endif

SOURCES += ${FORCESOURCES}


###  ENSEMBLES  #########################################

ifneq (,$(findstring nve,${MAKETARGET}))
CFLAGS += -DNVE
endif

ifneq (,$(findstring mik,${MAKETARGET}))
CFLAGS += -DMIK
endif

# CG
ifneq (,$(strip $(findstring cg,${MAKETARGET})))
SOURCES += ${CGSOURCES}
IMDHEADERS += cg_util.h  
CFLAGS  += -DCG
ifneq (,$(strip $(findstring cgef,${MAKETARGET})))
CFLAGS  += -DCGEF
else
ifneq (,$(strip $(findstring cge,${MAKETARGET})))
CFLAGS  += -DCGE
endif
endif
ifneq (,$(strip $(findstring cgf,${MAKETARGET})))
CFLAGS  += -DCGF
endif
ifneq (,$(strip $(findstring fr,${MAKETARGET})))
CFLAGS  += -DFR
endif
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
# SOURCES += ${FRACSOURCES}
endif

ifneq (,$(findstring nvx,${MAKETARGET}))
CFLAGS += -DNVX
SOURCES += ${TRANSSOURCES}
endif

ifneq (,$(findstring ftg,${MAKETARGET}))
CFLAGS += -DFTG
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

ifneq (,$(findstring sendrec,${MAKETARGET}))
CFLAGS += -DSR
endif

ifneq (,$(findstring einstein,${MAKETARGET}))
CFLAGS += -DEINSTEIN
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

ifneq (,$(findstring writef,${MAKETARGET}))
CFLAGS += -DWRITEF
endif

ifneq (,$(findstring deform,${MAKETARGET}))
CFLAGS += -DDEFORM
SOURCES += ${DEFORMSOURCES}
endif

ifneq (,$(findstring glokdeform,${MAKETARGET}))
CFLAGS += -DGLOKDEFORM
endif

ifneq (,$(findstring snap,${MAKETARGET}))
CFLAGS += -DSNAPSHOT
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

ifneq (,$(findstring atdist,${MAKETARGET}))
CFLAGS += -DATDIST
endif

ifneq (,$(findstring ordpar,${MAKETARGET}))
CFLAGS += -DORDPAR
endif

# EPITAX
ifneq (,$(strip $(findstring epitax,${MAKETARGET})))
CFLAGS  += -DEPITAX
SOURCES += ${EPITAXSOURCES}
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

# use 4-point-interpolation
ifneq (,$(findstring 4point,${MAKETARGET}))
CFLAGS += -DFOURPOINT
endif

# debug flag
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += -DDEBUG
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
	${MV} $@ ${BINDIR}; rm -f $@


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
