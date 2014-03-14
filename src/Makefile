############################################################################
#
# Imd -- The ITAP Molecular Dynamics Program
#
# Copyright 1996-2012 Institute for Theoretical and Applied Physics,
# University of Stuttgart, D-70550 Stuttgart
#
# $Revision$
# $Date$
#
############################################################################
#
# Beware: This Makefile works only with GNU make (gmake)!
#
# Usage:  gmake <target>
#
# <target> has the form
#
#    imd[_<parallel>][_<option>[_<option>...]]
#
# The parallelization method <parallel> can be one of:
#
#    mpi   compile for parallel execution, using MPI
#    omp   compile for parallel execution, using OpenMP
#    ompi  compile for parallel execution, using OpenMP and MPI
#    pacx  compile for parallel execution, using PACX
#
# The compilation options, whose names are parsed from the make target,
# must include at least one ensemble (integrator).
#
###########################################################################
#
# Customizing this Makefile
#
# As IMD supports a large number of compile options, you will have to
# compile IMD freqently. Before doing so, however, you must check whether
# the settings in this Makefile fit your needs. You possibly have to
# customize these setttings. Before you can do that, we have to explain
# a bit how the compilation process works.
#
# The compilation process requires the environment variable IMDSYS to
# be set to a recognized value. It specifies what system you have, and
# what compiler you are using. The flags for the compiler and the linker
# are then selected as a function of this variable. It is also possible
# to pass the value of IMDSYS on the command line, e.g.:
#
#   gmake IMDSYS=P4-icc imd_mpi_nve
#
# Another important ingredient is the parallelization method, which is
# determined from the make target. The parallelization method is stored
# in the variable PARALLEL, which takes as value one of SERIAL, MPI,
# OMP, OMPI, or PACX.
#
# Depending on the value of ${IMDSYS}, a number of variables must be
# set, from which everything else is constructed.
#
# CC_${PARALLEL} defines the compiler to be used for parallelization
# method ${PARALLEL}. If not defined, the parallelization method
# ${PARALLEL} is not available.
#
# BIN_DIR defines the directory where the IMD binary is put. Note that
# this directory must exist.
#
# MV defines the program used to move the IMD binary to ${BIN_DIR}.
# The default is mv, which is usually ok.
#
# The compilation options are stored in the variable CFLAGS.
# The initial value of CFLAGS is set to the variable FLAGS,
# which can be given on the command line as explained above for
# IMDSYS, although this is usually not necessary.
#
# If the option debug was specified, ${DEBUG_FLAGS} is then appended
# to ${CFLAGS}, otherwise ${OPT_FLAGS}. If the option prof was specified
# (for profiling), ${PROF_FLAGS} is also appended to ${CFLAGS}. However,
# before appending ${OPT_FLAGS} or ${DEBUG_FLAGS} to ${CFLAGS}, some
# parallelization specific flags are appended to them:
#
#   OPT_FLAGS   += ${${PARALLEL}_FLAGS} ${OPT_${PARALLEL}_FLAGS}
#   DEBUG_FLAGS += ${${PARALLEL}_FLAGS} ${DEBUG_${PARALLEL}_FLAGS}
#
# If any of these variables is not defined, it is assumed to be empty.
# This setup should provide sufficient flexibility to set one's favorite
# flags, depending on parallelization, profiling, and optimization/debugging.
#
# Similarly, the link libraries are stored in the variable LIBS,
# to which ${${PARALLEL}_LIBS} and possibly ${PROF_LIBS} (for profiling)
# is appended.
#
# preprocessor flags controlling compile options of IMD are stored in the
# variable ${PP_FLAGS}.
#
# You may have to change the setting for an existing value of IMDSYS.
# or you have to add support for a new value of IMDSYS. The latter is
# best done by using the folloing template for IMDSYS=sys-cc:
#
# ifeq (sys-cc,${IMDSYS})
#   CC_SERIAL		= serial-compiler
#   CC_OMP		= OpenMP-compiler
#   CC_MPI		= MPI-compiler
#   CC_OMPI		= OpenMP+MPI-compiler
#   BIN_DIR		= ${HOME}/bin
#   OPT_FLAGS		+= generic flags for optimization
#   OPT_MPI_FLAGS	+= MPI-specific flags for optimization
#                          similar variables for other parallelizations
#   MPI_FLAGS		+= MPI-specific flags
#                          similar variables for other parallelizations
#   DEBUG_FLAGS		+= generic flags for debugging
#   DEBUG_MPI_FLAGS	+= MPI-specific flags for debugging
#                          similar variables for other parallelizations
#   PROF_FLAGS		+= flags for profiling
#   LIBS		+= generically needed libraries
#   MPI_LIBS		+= MPI-specific libraries
#                          similar variables for other parallelizations
#   PROF_LIBS		+= libraries for profiling
# endif
#
# Variables remaining empty need not be mentioned.

###########################################################################
#
#  Defaults for some variables
#
###########################################################################

MV		= mv      # program to move imd binary to ${BIN_DIR}
LIBS		+= -lm
MPI_FLAGS	+= #-DMPI
OMP_FLAGS	+= #-DOMP
OMPI_FLAGS	+= #-DMPI -DOMP
PACX_FLAGS	+= #-DMPI -DPACX
DEBUG_FLAGS	+= #-DDEBUG # -Wall # very noisy

# directory of the openkim API
KIM_DIR 	= /data/daniel/openkim/openkim-api

###########################################################################
#
#  flags for IA32-Linux
#
###########################################################################

# lima
ifeq (lima,${IMDSYS})
  CC_SERIAL     = icc
  CC_MPI        = mpicc
  MPICH_CC      = icc
  MPICH_CLINKER = icc
  OMPI_MPICC    = icc
  BIN_DIR       = ${HOME}/bin/${HOSTTYPE}
#  OPT_FLAGS     += -Ofast -xHost -ipo                # second best performance, but worse than current "deprecated" binaries
  OPT_FLAGS     += -Ofast -xHost -ipo -fno-alias     # best performance - only applicable if no pointer aliases are used!!!
  DEBUG_FLAGS   += -g
  PROF_FLAGS    += -g3 -pg
endif

# woody
ifeq (woody,${IMDSYS})
  CC_SERIAL	= icc
  CC_MPI	= mpicc
  MPICH_CC      = icc
  MPICH_CLINKER = icc
  OMPI_MPICC    = icc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  OPT_FLAGS	+= -O -funroll-all-loops -fomit-frame-pointer -ipo -xT
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
#  LFLAGS        +=  -static
#  export        OMPI_MPICC MPICH_CC # MPICH_CLINKER
endif


# intel core 2 duo gcc
ifeq (x86_64-c2_gcc,${IMDSYS})
  CC_SERIAL	= gcc
  CC_MPI	= mpicc
  MPICH_CC      = gcc
#  MPICH_CLINKER = icc
  OMPI_MPICC    = gcc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
#  OPT_FLAGS	+= -O -m64 -Wno-unused -march=core2   -mmmx -msse -ffast-math  -funroll-all-loops -fomit-frame-pointer
  OPT_FLAGS	+= -O3  -march=core2     -funroll-all-loops -fomit-frame-pointer
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
  LFLAGS        +=  -static
  export        OMPI_MPICC MPICH_CC # MPICH_CLINKER
endif

ifeq (x86_64-c2_gcc4,${IMDSYS})
  CC_SERIAL	= gcc
  CC_MPI	= mpicc
  MPICH_CC      = gcc
#  MPICH_CLINKER = icc
  OMPI_MPICC    = gcc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
#  OPT_FLAGS	+= -O -m64 -Wno-unused -march=core2   -mmmx -msse -ffast-math  -funroll-all-loops -fomit-frame-pointer
#  OPT_FLAGS	+= -O3 -m64 -march=core2  -mtune=core2  -funroll-all-loops -fomit-frame-pointer -ffast-math
  OPT_FLAGS	+= -O -m64 -march=core2  -mtune=core2  -funroll-all-loops -fomit-frame-pointer -ffast-math
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
# static seems to have problems with the mpi implemenation?
#  LFLAGS        +=  -static
  export        OMPI_MPICC MPICH_CC # MPICH_CLINKER
endif

# intel core 2 duo, icc
ifeq (x86_64-c2_icc,${IMDSYS})
  CC_SERIAL	= icc
  CC_MPI	= mpicc
  MPICH_CC      = icc
  MPICH_CLINKER = icc
  OMPI_MPICC    = icc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  OPT_FLAGS	+= -O3 -ip  -fno-builtin -axT -march=core2  -mtune=core2  # remove -axP for Opteron
  MPI_FLAGS     +=
  OMP_FLAGS	+= -openmp
  OMPI_FLAGS	+= -openmp
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -prof_gen
  RCD_FLAGS	+= # -DRCD -rcd
  MPI_LIBS      +=
  LFLAGS        += # -i-static
  export        OMPI_MPICC MPICH_CC MPICH_CLINKER
endif



# AMD Opteron or Intel EM46T, gcc3
ifeq (x86_64-gcc3,${IMDSYS})
  CC_SERIAL	= gcc
  CC_OMP	= gcc
  CC_MPI	= mpicc
  MPICH_CC      = gcc
#  MPICH_CLINKER = icc
  OMPI_MPICC    = gcc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  OPT_FLAGS	+= -O -m64 -Wno-unused
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
  LFLAGS        += # -static
  export        OMPI_MPICC MPICH_CC # MPICH_CLINKER
endif

# AMD Opteron or Intel EM64T, icc
ifeq (x86_64-icc,${IMDSYS})
  CC_SERIAL	= icc
  CC_OMP	= icc
  CC_MPI	= mpicc
  MPICH_CC      = icc
  MPICH_CLINKER = icc
  OMPI_MPICC    = icc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  OPT_FLAGS	+= -O3 -ip # -fno-builtin # -axP # remove -axP for Opteron
  MPI_FLAGS     +=
  OMP_FLAGS	+= -openmp
  OMPI_FLAGS	+= -openmp
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -prof_gen
  RCD_FLAGS	+= # -DRCD -rcd
  MPI_LIBS      +=
  LFLAGS        += # -i-static
  export        OMPI_MPICC MPICH_CC MPICH_CLINKER
endif

# generic x86, gcc3 or gcc4
ifeq (x86-gcc3,${IMDSYS})
  CC_SERIAL	= gcc
  CC_MPI	= mpicc
  MPICH_CC      = gcc
  OMPI_MPICC    = gcc
#  MPICH_CLINKER = gcc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  OPT_FLAGS	+= -O
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
  PAPI_INC      =  -I /opt/include
  PAPI_LIBS     =  -L /opt/lib -lpapi -lperfctr
  LFLAGS        += # -static
  export        OMPI_MPICC MPICH_CC # MPICH_CLINKER
endif

# generic x86, icc
ifeq (x86-icc,${IMDSYS})
  CC_SERIAL	= icc
  CC_OMP	= icc
  CC_MPI	= mpicc
  MPICH_CC      = icc
  OMPI_MPICC    = icc
  MPICH_CLINKER = icc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  OPT_FLAGS	+= -O -ip
  OMP_FLAGS	+= -openmp
  OMPI_FLAGS	+= -openmp
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -prof_gen
  PAPI_INC      =  -I /opt/include
  PAPI_LIBS     =  -L /opt/lib -lpapi -lperfctr
  RCD_FLAGS	+= -DRCD -rcd
  LFLAGS        += # -static
  export        OMPI_MPICC MPICH_CC MPICH_CLINKER
endif

# Athlon MP or XP, icc
ifeq (AthlonMP-icc,${IMDSYS})
  CC_SERIAL	= icc
  CC_OMP	= icc
  CC_MPI	= mpicc
  CC_OMPI	= mpicc
  MPICH_CC      = icc
  MPICH_CLINKER = icc
  OMPI_MPICC    = icc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  OPT_FLAGS	+= -O -ip -tpp6 -axK
  OMP_FLAGS	+= -openmp
  OMPI_FLAGS	+= -openmp
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -prof_gen
  RCD_FLAGS	+= -DRCD -rcd
  LFLAGS        += -i-static
  export OMPI_MPICC MPICH_CC  MPICH_CLINKER
endif

# Athlon MP or XP, gcc3
ifeq (AthlonMP-gcc3,${IMDSYS})
  CC_SERIAL	= gcc
  CC_MPI	= mpicc
  MPICH_CC      = gcc
#  MPICH_CLINKER = gcc
  OMPI_MPICC    = gcc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  OPT_FLAGS	+= -O -march=athlon-mp
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
  LFLAGS        += #-static
  export        OMPI_MPICC MPICH_CC # MPICH_CLINKER
endif

# Pentium 4 or Xeon, icc
ifeq (P4-icc,${IMDSYS})
  CC_SERIAL	= icc
  CC_OMP	= icc
  CC_MPI	= mpicc
  CC_OMPI	= mpicc
  MPICH_CC      = icc
  MPICH_CLINKER = icc
  OMPI_MPICC    = icc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  OPT_FLAGS	+= -O -ip
  OMP_FLAGS	+= -openmp
  OMPI_FLAGS	+= -openmp
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -prof_gen
  RCD_FLAGS	+= -DRCD -rcd
  LFLAGS        += -i-static
  export        OMPI_MPICC MPICH_CC MPICH_CLINKER
endif

# woody
ifeq (woody-icc,${IMDSYS})
  CC_SERIAL	= icc
  CC_OMP	= icc
  CC_MPI	= mpicc
  CC_OMPI	= mpicc
  MPICH_CC      = icc
  MPICH_CLINKER = icc
  OMPI_MPICC    = icc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
#  OPT_FLAGS	+= -O -ip -m64 -mkl  -msse2
#  OPT_FLAGS	+= -O -ip   -O3 -xP
  OPT_FLAGS	+= -O -ip  -msse2
  OMP_FLAGS	+= -openmp
  OMPI_FLAGS	+= -openmp
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -prof_gen
  RCD_FLAGS	+= -DRCD -rcd
  LFLAGS        += -i-static
  export        OMPI_MPICC MPICH_CC MPICH_CLINKER
endif

# Pentium 4 or Xeon, gcc3
ifeq (P4-gcc3,${IMDSYS})
  CC_SERIAL	= gcc
  CC_MPI	= mpicc
  MPICH_CC      = gcc
  OMPI_MPICC    = gcc
#  MPICH_CLINKER = gcc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
#  OPT_FLAGS	+= -O -march=pentium4
  OPT_FLAGS += -mtune=native
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
  PAPI_INC      =  -I /opt/include
  PAPI_LIBS     =  -L /opt/lib -lpapi -lperfctr
  LFLAGS        += # -static
  export        OMPI_MPICC MPICH_CC # MPICH_CLINKER
endif

# Pentium 4 or Xeon, gcc3
ifeq (woody_gcc,${IMDSYS})
  CC_SERIAL	= gcc
  CC_MPI	= mpicc
  MPICH_CC      = gcc
  OMPI_MPICC    = gcc
#  MPICH_CLINKER = gcc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  OPT_FLAGS	+= -O -march=nocona
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
  PAPI_INC      =  -I /opt/include
  PAPI_LIBS     =  -L /opt/lib -lpapi -lperfctr
  LFLAGS        += # -static
  export        OMPI_MPICC MPICH_CC # MPICH_CLINKER
endif


# Pentium III, icc
ifeq (P3-icc,${IMDSYS})
  CC_SERIAL	= icc
  CC_OMP	= icc
  CC_MPI	= mpicc
  CC_OMPI	= mpicc
  MPICH_CC      = icc
  MPICH_CLINKER = icc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  FFTW_DIR 	= /common/linux/paket/fftw-3.0.1
  OPT_FLAGS	+= -O -ip -tpp6 -axK # -static
  OMP_FLAGS	+= -openmp
  OMPI_FLAGS	+= -openmp
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -prof_gen
  RCD_FLAGS	+= -DRCD -rcd
  LFLAGS        += -i-static
  export        MPICH_CC  MPICH_CLINKER
endif

# Pentium III, gcc3
ifeq (P3-gcc3,${IMDSYS})
  CC_SERIAL	= gcc
  CC_MPI	= mpicc
  MPICH_CC      = gcc
#  MPICH_CLINKER = gcc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  FFTW_DIR 	= /common/linux/paket/fftw-3.0.1
  OPT_FLAGS	+= -O -mcpu=pentium3 # -static
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
  export        MPICH_CC # MPICH_CLINKER
endif

# generic ia32 CPU, gcc 2.95 - slow!
ifeq (ia32-gcc2,${IMDSYS})
  CC_SERIAL	= gcc
  CC_MPI	= mpicc
  MPICH_CC      = gcc
#  MPICH_CLINKER = gcc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  FFTW_DIR 	= /common/linux/paket/fftw-3.0.1
  OPT_FLAGS	+= -O3
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
  export        MPICH_CC # MPICH_CLINKER
endif

###########################################################################
#
#  flags for IA64-Linux
#
###########################################################################

# Itanium 2, icc 8.0 or higher
ifeq (ia64-icc,${IMDSYS})
  CC_SERIAL	= icc
  CC_OMP	= icc
  CC_MPI	= mpicc
  CC_OMPI	= mpicc
  BIN_DIR	= ${HOME}/bin/ia64
  OPT_FLAGS	+= -O -ip -ivdep_parallel -wd1572 # remove -wd1572 for icc 8.0
  NOALIAS       = -fno-fnalias
  MPI_FLAGS     +=
  OMP_FLAGS	+= -openmp
  OMPI_FLAGS	+= -openmp -quiet
  PROF_FLAGS	+= -prof_gen
  DEBUG_FLAGS	+= -g
  LFLAGS        += -i-static # -p
  PAPI_INC      =  -I /opt/papi/3.0/include
  PAPI_LIBS     =  -L /opt/papi/3.0/lib -lpapi
endif

# Itanium 2, ecc
ifeq (ia64-ecc,${IMDSYS})
  CC_SERIAL	= ecc
  CC_OMP	= ecc
  CC_MPI	= mpicc
  CC_OMPI	= mpicc
  BIN_DIR	= ${HOME}/bin/ia64
  OPT_FLAGS	+= -O -ip -ivdep_parallel
  OMP_FLAGS	+= -openmp
  OMPI_FLAGS	+= -openmp
  PROF_FLAGS	+= -p
  DEBUG_FLAGS	+= -g
  LFLAGS        += -i-static -p -static
endif

###########################################################################
#
#  flags for alpha
#
###########################################################################

# alpha ev56 or higher, cc
ifeq (alpha-cc,${IMDSYS})
  CC_SERIAL	= cc
  CC_OMP	= cc
  BIN_DIR	= ${HOME}/bin/alpha
  OPT_FLAGS	+= -DALPHA -O3 -float -fp_reorder -arch ev56 -tune host
  OMP_FLAGS	+= -mp
  OMPI_FLAGS	+= -mp
  PROF_FLAGS	+= -prof_gen
  DEBUG_FLAGS	+= -g3 -pg
endif

# alpha, gcc2 - slow!
ifeq (alpha-gcc2,${IMDSYS})
  CC_SERIAL	= gcc
  BIN_DIR	= ${HOME}/bin/$alpha
  OPT_FLAGS	+= -DALPHA -O3
  PROF_FLAGS	+= -g3 -pg
  DEBUG_FLAGS	+= -g
endif

###########################################################################
#
#  flags for IRIX
#
###########################################################################

# irix-mips3, cc
ifeq (irix-cc,${IMDSYS})
  CC_SERIAL	= cc
  CC_OMP	= cc
  BIN_DIR	= ${HOME}/bin/iris4d
  OPT_FLAGS	+= -Dsgi -O3 -n32 -mips3 -xansi -woff 1174
  OMP_FLAGS	+= -mp
  OMPI_FLAGS	+= -mp
  PROF_FLAGS	+= -g3
  DEBUG_FLAGS	+= -Dsgi -g  -n32 -mips3 -xansi -woff 1174
endif

# irix, gcc2 - slow!
ifeq (irix-gcc2,${IMDSYS})
  CC_SERIAL	= gcc
  BIN_DIR	= ${HOME}/bin/iris4d
  OPT_FLAGS	+= -Dsgi -O3
  PROF_FLAGS	+= -g3 -pg
  DEBUG_FLAGS	+= -Dsgi -g
endif

###########################################################################
#
#  flags for sparc
#
###########################################################################

# UltraSparc VI, cc
ifeq (USparc4-cc,${IMDSYS})
  CC_SERIAL	= cc
  CC_OMP	= cc
  CC_MPI	= mpcc
  BIN_DIR	= ${HOME}/bin/sparc
  MPI_LIBS	+= -lmpi
  OPT_FLAGS	+= -fast -xarch=sparcvis2 # -m64
  OMP_FLAGS	+= -xopenmp
  OMPI_FLAGS	+= -xopenmp
  PROF_FLAGS	+= -p
  DEBUG_FLAGS	+= -g -xO3
  LFLAGS        += -fast -xarch=sparcvis2 # -m64
endif

# UltraSparc III, cc
ifeq (USparc3-cc,${IMDSYS})
  CC_SERIAL	= cc
  CC_OMP	= cc
  CC_MPI	= mpcc
  BIN_DIR	= ${HOME}/bin/sparc
  MPI_LIBS	+= -lmpi
  OPT_FLAGS	+= -fast -xarch=sparcvis2 # -m64
  OMP_FLAGS	+= -xopenmp
  OMPI_FLAGS	+= -xopenmp
  PROF_FLAGS	+= -p
  DEBUG_FLAGS	+= -g -xO3
  LFLAGS        += -fast -xarch=sparcvis2 # -m64
endif

# sparc, gcc3 - slow!
ifeq (sparc-gcc3,${IMDSYS})
  CC_SERIAL	= gcc
  BIN_DIR	= ${HOME}/bin/sparc
  OPT_FLAGS	+= -O
  PROF_FLAGS	+= -g3 -pg
  DEBUG_FLAGS	+= -g
endif

# sparc, gcc2 - slow!
ifeq (sparc-gcc2,${IMDSYS})
  CC_SERIAL	= gcc
  BIN_DIR	= ${HOME}/bin/sparc
  OPT_FLAGS	+= -O3
  PROF_FLAGS	+= -g3 -pg
  DEBUG_FLAGS	+= -g
endif

###########################################################################
#
#  flags for T3E
#
###########################################################################

# Cray T3E, cc
ifeq (T3E-cc,${IMDSYS})
  CC_MPI	= cc
  CC_PACX	= cc
  BIN_DIR	= ${HOME}/bin/t3e
  MPI_LIBS	+= -lmpi
  PROF_LIBS	+= -lapp
  OPT_FLAGS	+= -Dt3e -O3 -htolerant,aggress,report=isf
  PROF_FLAGS	+= -Gf -happrentice
  DEBUG_FLAGS	+= -Dt3e -g
  PACX_DIR	= ${HOME}/WORK/PACX
  PACX_LIBS	+= -L ${PACX_DIR}/lib -lpacx -llzo -lmpi
  PACX_FLAGS	+= -I${PACX_DIR}/include
endif

###########################################################################
#
#  flags for Hitachi SR8000
#
###########################################################################

# Hitachi sr8k, cc
ifeq (sr8k-cc,${IMDSYS})
  CC_MPI		= mpicc
  CC_OMP		= cc
  CC_OMPI		= mpicc
  BIN_DIR		= ${HOME}/bin/SR8000
  OPT_FLAGS		+= -O4 +Op -msg=e
  OPT_MPI_FLAGS		+= -noparallel
  OPT_OMP_FLAGS		+= -omp -par
  OPT_OMPI_FLAGS	+= -omp -par
  DEBUG_FLAGS		+= -g
  DEBUG_OMP_FLAGS	+= -omp -par=1 -O2
  DEBUG_OMPI_FLAGS	+= -omp -par   -O2
  PROF_FLAGS		+= -Xfuncmonitor
  PROF_LIBS		+= -lpl
endif

# Hitachi sr8k, xcc cross compiler
ifeq (sr8k-xcc,${IMDSYS})
  CC_MPI		= xmpicc
  CC_OMP		= xcc
  CC_OMPI		= xmpicc
  BIN_DIR		= hwwfs1:sr8k/bin/SR8000
  MV			= scp   # we move binary to different machine
  OPT_FLAGS		+= -O4 +Op -msg=e
  OPT_MPI_FLAGS		+= -noparallel
  OPT_OMP_FLAGS		+= -omp -par
  OPT_OMPI_FLAGS	+= -omp -par
  DEBUG_FLAGS		+= -g
  DEBUG_OMP_FLAGS	+= -omp -par=1 -O2
  DEBUG_OMPI_FLAGS	+= -omp -par   -O2
  PROF_FLAGS		+= -Xfuncmonitor
  PROF_LIBS		+= -lpl
endif

###########################################################################
#
#  flags for NEC SX crosscompiler
#
###########################################################################

ifeq (SX8,${IMDSYS})
  CC_SERIAL             = sxcc
  CC_MPI		= sxmpicc
#  CC_OMP		=
#  CC_OMPI		=
  BIN_DIR		= ${HOME}/bin/sx8
  OPT_FLAGS		+= -Xa -O -sx8 -C hopt -pi auto -pi file=imd_geom_3d.c:imd_geom_mpi_3d.c -pvctl fullmsg -w off=177:550 -O extendreorder # -gv
#  OPT_MPI_FLAGS	+=
#  OPT_OMP_FLAGS	+=
#  OPT_OMPI_FLAGS	+=
  DEBUG_FLAGS		+= -Xa -sx8 -gv
#  DEBUG_OMP_FLAGS	+=
#  DEBUG_OMPI_FLAGS	+=
  PROF_FLAGS		+= -ftrace -DFTRACE
  PROF_LIBS		+= -lftrace
endif

ifeq (SX6,${IMDSYS})
  CC_SERIAL             = sxcc
#  CC_MPI		=
#  CC_OMP		=
#  CC_OMPI		=
  BIN_DIR		= ${HOME}/bin/sx6
  OPT_FLAGS		+= -Xa -O -sx6 -C hopt -pi auto -pi file=imd_geom_3d.c:imd_geom_mpi_3d.c -pvctl fullmsg -w off=177:550 -O extendreorder # -gv
#  OPT_MPI_FLAGS	+=
#  OPT_OMP_FLAGS	+=
#  OPT_OMPI_FLAGS	+=
  DEBUG_FLAGS		+= -Xa -sx6 -gv
#  DEBUG_OMP_FLAGS	+=
#  DEBUG_OMPI_FLAGS	+=
  PROF_FLAGS		+= -ftrace -DFTRACE
  PROF_LIBS		+= -lftrace
endif

###########################################################################
#
#  flags for IBM SP
#
###########################################################################

# Power4 regatta, 64 bit
ifeq (Power4-cc-64,${IMDSYS})
  CC_SERIAL	= xlc
  CC_MPI	= mpcc_r
  CC_OMP	= xlc_r
  CC_OMPI	= mpcc_r
  BIN_DIR	= ${HOME}/bin/powerpc
  OPT_FLAGS	+= -O4 -w -q64 -qunroll -qalign # -qtbtable=full
  OMP_FLAGS	+= -DUSE_WALLTIME -qsmp=omp
  OMPI_FLAGS	+= -qsmp=omp
  DEBUG_FLAGS	+= -g -q64
  PROF_FLAGS	+= -p
  PROF_LIBS	+= -lpl
  PAPI_LIBS     =  -lpapi64 -lpmapi
  LIBS          += # -lmpiprof
  LFLAGS        += -q64
endif

# Power4 regatta, 32bit
ifeq (Power4-cc-32,${IMDSYS})
  CC_SERIAL	= xlc
  CC_MPI	= mpcc
  CC_OMP	= xlc_r
  CC_OMPI	= mpcc_r
  BIN_DIR	= ${HOME}/bin/powerpc
  OPT_FLAGS	+= -O4 -w -q32 -qunroll -qalign # -qtbtable=full
  OMP_FLAGS	+= -DUSE_WALLTIME -qsmp=omp
  OMPI_FLAGS	+= -qsmp=omp
  DEBUG_FLAGS	+= -g -q64
  PROF_FLAGS	+= -p
  PROF_LIBS	+= -lpl
  PAPI_LIBS     =  -lpapi32 -lpmapi
  LIBS          += # -lmpiprof
  LFLAGS        += -q32
endif

# IBM GlueGene/L
ifeq (bgl-cc,${IMDSYS})
  BGLSYS 	= /bgl/BlueLight/ppcfloor/bglsys
  CC_SERIAL	= dummy
  CC_MPI	= blrts_xlc -DBG -DBGL -I${BGLSYS}/include
  BIN_DIR	= ${HOME}/bin
  OPT_FLAGS	+= -O4 -qarch=440 -qtune=440
  DEBUG_FLAGS	+= -g
  LIBS          += -L${BGLSYS}/lib -lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts
  LFLAGS        += -O4
endif

# IBM GlueGene/P
ifeq (bgp-cc,${IMDSYS})
  BGPSYS 	= /bgsys/drivers/ppcfloor/arch
  CC_SERIAL	= dummy
  CC_MPI	= mpixlc -DBG -DBGP -I${BGPSYS}/include
  BIN_DIR	= ${HOME}/bin
  OPT_FLAGS	+= -O4 -qarch=450 -qtune=450
  DEBUG_FLAGS	+= -g
  LIBS          += # -L${BGPSYS}/lib
  LFLAGS        += -O4
endif


###########################################################################
#
#  Flags for Cell Broad Band Engine
#
###########################################################################

# Cell Broad Band Engine with IBM XLC compiler
ifeq (cbe-xlc,${IMDSYS})
  # Path to Cell SDK
  CELLSDK=/opt/ibm/cell-sdk/prototype

  # Headers for libsync
  CELLSDKLIB=${CELLSDK}/src/lib

  # The spu-timing utility provided by the Cell SDK
  SPUTIMING=${CELLSDK}/bin/spu_timing -running-count

  SPUCC         = spuxlc
  CC_SERIAL     = ppuxlc -q32
  CC_MPI        = mpicc  -q32

  OPT_FLAGS     += -O3 -qstrict -qcpluscmt -qflag=w:w -qproto -qro -qwarn64  -qhot=level=1:simd  -qreport -qinline -qunroll=yes -qlibansi -qipa=level=1
  SPUCFLAGS     += -O3 -qstrict -qcpluscmt -qflag=w:w -qhot=level=1 -qproto -qreport -qinline  -qinline-dma64:mdma64:wait_dma -I${CELLSDKLIB}/sync -I${CELLSDKLIB}
  PP_FLAGS      += -I${CELLSDKLIB}/sync
  EMBEDSPU      = ppu-embedspu -m32

  BIN_DIR	= ${HOME}/bin/cbe

  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -pg
  PROF_LIBS	+= -pg

  SPULDFLAGS    +=
  LFLAGS        += -O3 -qlibansi -qipa=level=1
  LIBS          += -lspe2 -lpthread
endif


# Cell Broad Band Engine with GCC
ifeq (cbe-gcc,${IMDSYS})
  # Path to Cell SDK
  CELLSDK=/opt/ibm/cell-sdk/prototype

  # Headers for libsync
  CELLSDKLIB=${CELLSDK}/src/lib

  # The spu-timing utility provided by the Cell SDK
  SPUTIMING=${CELLSDK}/bin/spu_timing -running-count

  # Use the PPU compiler ppu-gcc, also for MPI!
  # as main compiler as it defines the __PPU__ macro
  # The ordinary gcc does not set this macro

  CC_SERIAL     = ppu-gcc -m32
  CC_MPI        = mpicc   -m32
  SPUCC         = spu-gcc
  OPT_FLAGS     += -pipe  -std=gnu89 -pedantic -Wall -Wno-unused -Wcast-qual -Wno-long-long   -O3 -funroll-loops
  PP_FLAGS      += -I${CELLSDKLIB}/sync
  SPUCFLAGS     += -pipe  -std=gnu89 -pedantic -Wall -Wunused -Wshadow -Wcast-qual -Wno-long-long -O2 -funroll-loops  -I${CELLSDKLIB}/sync/spu -I${CELLSDKLIB}
  EMBEDSPU      = ppu-embedspu -m32

  BIN_DIR	= ${HOME}/bin/cbe

  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
  PROF_LIBS	+= -g3 -pg

  SPULDFLAGS    +=
  LFLAGS        +=
  LIBS          += -lspe2
endif


###########################################################################
#
#  flags for generic version of IMD (slow)
#
###########################################################################

ifeq (generic-gcc,${IMDSYS})
  CC_SERIAL     = gcc
#  BIN_DIR       = ../Run
  BIN_DIR       = ${HOME}/bin/${HOSTTYPE}
  OPT_FLAGS     += -O
  DEBUG_FLAGS   += -g
endif

ifeq (generic-cc,${IMDSYS})
  CC_SERIAL     = cc
#  BIN_DIR       = ../Run
  BIN_DIR       = ${HOME}/bin/${HOSTTYPE}
  OPT_FLAGS     += -O
  DEBUG_FLAGS   += -g
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
PP_FLAGS += -DMPI
endif
# OpenMP
ifneq (,$(strip $(findstring omp,${MAKETARGET})))
PARALLEL = OMP
PP_FLAGS += -DOMP
endif
# MPI + OpenMP
ifneq (,$(strip $(findstring ompi,${MAKETARGET})))
PARALLEL = OMPI
PP_FLAGS += -DMPI -DOMPI
endif
# PACX
ifneq (,$(strip $(findstring pacx,${MAKETARGET})))
PARALLEL = PACX
PP_FLAGS += -DMPI -DPACX
endif


###########################################################################
#
#  Compiler, flags, libraries
#
###########################################################################

# compiler; if empty, we issue an error later
CC = ${CC_${PARALLEL}}

# optimization flags
OPT_FLAGS   += ${${PARALLEL}_FLAGS} ${OPT_${PARALLEL}_FLAGS}
DEBUG_FLAGS += ${${PARALLEL}_FLAGS} ${DEBUG_${PARALLEL}_FLAGS}

# libraries
LIBS += ${${PARALLEL}_LIBS}

# optimization or debug
CFLAGS := ${FLAGS}
PP_FLAGS += ${FLAGS}
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += ${DEBUG_FLAGS}
else
CFLAGS += ${OPT_FLAGS}
endif

# profiling support
ifneq (,$(findstring prof,${MAKETARGET}))
CFLAGS += ${PROF_FLAGS}
LIBS   += ${PROF_LIBS}
endif

# MPE logging
ifneq (,$(findstring mpel,${MAKETARGET}))
PP_FLAGS += -DMPE
MPI_LIBS += -llmpe -lmpe
endif

# MPE tracing
ifneq (,$(findstring mpet,${MAKETARGET}))
PP_FLAGS += -DMPE
MPI_LIBS += -ltmpe -lmpe
endif

ifneq (,$(findstring debug,${MAKETARGET}))
PP_FLAGS += -DDEBUG
endif

###########################################################################
#
# IMD sources
#
###########################################################################

IMDHEADERS      = config.h globals.h imd.h makros.h potaccess.h \
                  prototypes.h types.h

SOURCES         = imd_maxwell.c imd_integrate.c imd_misc.c \
	          imd_param.c imd_alloc.c imd_io.c imd_io_3d.c \
                  imd_potential.c imd_time.c imd_generate.c \
                  imd_distrib.c imd_main_3d.c
SOURCES2D       = ${SOURCES} imd_geom_2d.c imd_pictures_2d.c
SOURCES3D       = ${SOURCES} imd_geom_3d.c imd_pictures_3d.c
RISCSOURCES2D   = imd_main_risc_2d.c
RISCSOURCES3D   = imd_main_risc_3d.c

MPISOURCES      = imd_mpi_util.c
MPISOURCES2D    = ${MPISOURCES} imd_main_mpi_2d.c imd_geom_mpi_2d.c \
		  imd_comm_force_2d.c
MPISOURCES3D    = ${MPISOURCES} imd_main_mpi_3d.c imd_geom_mpi_3d.c \
		  imd_comm_force_3d.c imd_fix_cells_3d.c

NBLSOURCES3D    = imd_geom_mpi_3d.c imd_comm_force_3d.c imd_fix_cells_3d.c

PAIRSOURCES     = imd_forces.c
EAM2SOURCES     = imd_forces_eam2.c
MEAMSOURCES     = imd_forces_meam.c
CGSOURCES	= imd_cg.c
COVALENTSOURCES = imd_forces_covalent.c
UNIAXSOURCES    = imd_forces_uniax.c imd_gay_berne.c
EWALDSOURCES    = imd_forces_ewald.c

CNASOURCES      = imd_cna.c

EPITAXSOURCES   = imd_epitax.c

FEFLSOURCES	= imd_fefl.c

DEFORMSOURCES   = imd_deform.c

SOCKHEADERS     = sockets.h sockutil.h socket_io.h
SOCKSOURCES     = socket_io.c sockutil.c

QUASISOURCES    = imd_qc.c

CORRSOURCES     = imd_correl.c

TRANSSOURCES    = imd_transport.c

LASERSOURCES	= imd_laser.c

TTMSOURCES	= imd_ttm.c

SMSOURCES	= imd_sm.c

FCSSOURCES      = imd_forces_fcs.c

#BBOOSTSOURCES	= imd_bboost.c imd_bb_core1.c imd_bb_core2.c
BBOOSTSOURCES	= imd_bboost.c

#########################################################
#
# IMD Configuration rules
#
#########################################################

HEADERS := ${IMDHEADERS}

# twod or not twod
ifneq (,$(strip $(findstring 2d,${MAKETARGET})))

  # 2d, serial or mpi
  ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
    SOURCES	:= ${SOURCES2D} ${MPISOURCES2D}
    else
    SOURCES	:= ${SOURCES2D} ${RISCSOURCES2D}
  endif
  PP_FLAGS  += -DTWOD

else

  # 3d, vec, serial, mpi
  ifneq (,$(strip $(findstring vec,${MAKETARGET})))

    SOURCES	:= ${SOURCES3D} ${NBLSOURCES3D}
    ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
      SOURCES	+= ${MPISOURCES}
    endif

  else

    ifneq (,$(strip $(findstring cbe,${MAKETARGET})))

      SOURCES	:= ${SOURCES3D} ${NBLSOURCES3D}
      ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
        SOURCES	+= ${MPISOURCES}
      endif

    else

      ifneq (,$(strip $(findstring nbl,${MAKETARGET})))

        SOURCES	:= ${SOURCES3D} ${NBLSOURCES3D}
        ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
          SOURCES += ${MPISOURCES}
        endif

      else

	ifneq (,$(strip $(findstring kim, ${MAKETARGET})))
	  SOURCES := ${SOURCES3D} ${NBLSOURCES3D}
	  ifneq (,$(strip $(findstring mpi, ${MAKETARGET})))
	    SOURCES += ${MPISOURCES}
          endif
	else
          ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
            SOURCES := ${SOURCES3D} ${MPISOURCES3D}
          else
            SOURCES := ${SOURCES3D} ${RISCSOURCES3D}
          endif
        endif
      endif
    endif
  endif
endif

SRCMAIN = imd.c
ifneq (,$(strip $(findstring py,${MAKETARGET})))
  SRCMAIN  = py_imd.c
  CFLAGS  += -fPIC
endif
ifneq (,$(strip $(findstring jvis,${MAKETARGET})))
  SRCMAIN   = jvis_imd.c
  PP_FLAGS += -DNVE -DNVT -DNPT -DNPT_iso -DREFPOS
endif
SOURCES += ${SRCMAIN}

###  INTERACTION  #######################################

# pair interaction is the default
ifneq (,$(strip $(findstring vec,${MAKETARGET})))
  FORCESOURCES = imd_main_vec_3d.c
else
  ifneq (,$(strip $(findstring cbe,${MAKETARGET})))
    FORCESOURCES = imd_forces_cbe.c
    PP_FLAGS += -DCBE
  else
    ifneq (,$(strip $(findstring nbl,${MAKETARGET})))
      FORCESOURCES = imd_forces_nbl.c
      PP_FLAGS += -DNBL
      ifneq (,$(strip $(findstring kim, ${MAKETARGET})))
	ERROR = "KIM is not compatible with any NBL force routine"
      endif
    else
      FORCESOURCES = ${PAIRSOURCES}
    endif
  endif
endif

# PAIR
ifneq (,$(strip $(findstring pair,${MAKETARGET})))
  ifneq (,$(strip $(findstring kim,${MAKETARGET})))
    ERROR = "KIM is not compatible with PAIR force routines"
  endif
  PP_FLAGS  += -DPAIR
endif

# EAM2 or EAM  -  this is now the same
ifneq (,$(strip $(findstring eam,${MAKETARGET})))
  ifneq (,$(strip $(findstring kim,${MAKETARGET})))
    ERROR = "KIM is not compatible with EAM force routines"
  endif
  # EEAM
  ifneq (,$(strip $(findstring eeam,${MAKETARGET})))
    PP_FLAGS  += -DEEAM
  endif
  # MEAM
  ifneq (,$(strip $(findstring meam,${MAKETARGET})))
    PP_FLAGS  += -DMEAM
    FORCESOURCES += ${MEAMSOURCES}  ${COVALENTSOURCES}
  else
    PP_FLAGS  += -DEAM2
    ifneq (,$(strip $(findstring asympot,${MAKETARGET})))
      PP_FLAGS += -DASYMPOT
      FORCESOURCES += ${EAM2SOURCES}
    else
      ifeq (,$(strip $(findstring nbl,${MAKETARGET})))
        ifeq (,$(strip $(findstring vec,${MAKETARGET})))
          FORCESOURCES += ${EAM2SOURCES}
        endif
      endif
    endif
  endif
endif

# ADP
ifneq (,$(strip $(findstring adp,${MAKETARGET})))
  ifneq (,$(strip $(findstring kim,${MAKETARGET})))
    ERROR = "KIM is not compatible with ADP force routines"
  endif
PP_FLAGS  += -DADP
endif

# TTBP
ifneq (,$(strip $(findstring ttbp,${MAKETARGET})))
PP_FLAGS  += -DTTBP
FORCESOURCES += ${COVALENTSOURCES}
endif

# TTBPXT
ifneq (,$(strip $(findstring ttbpxt,${MAKETARGET})))
PP_FLAGS  += -DTTBP -DXT
FORCESOURCES += ${COVALENTSOURCES}
endif

# STIWEB
ifneq (,$(strip $(findstring stiweb,${MAKETARGET})))
  ifneq (,$(strip $(findstring kim,${MAKETARGET})))
    ERROR = "KIM is not compatible with STIWEB force routines"
  endif
PP_FLAGS  += -DSTIWEB
FORCESOURCES += ${COVALENTSOURCES}
endif

# TERNBCC
ifneq (,$(strip $(findstring ternbcc,${MAKETARGET})))
PP_FLAGS  += -DSTIWEB -DTERNBCC
FORCESOURCES += ${COVALENTSOURCES}
endif

# TERSOFF
ifneq (,$(strip $(findstring tersoffmod,${MAKETARGET})))
  ifneq (,$(strip $(findstring kim,${MAKETARGET})))
    ERROR = "KIM is not compatible with TERSOFFMOD force routines"
  endif
PP_FLAGS  += -DTERSOFFMOD
FORCESOURCES += ${COVALENTSOURCES}
ifneq (,$(strip $(findstring tersoffmod2,${MAKETARGET})))
PP_FLAGS  += -DTERSOFFMOD2
endif
else
ifneq (,$(strip $(findstring tersoff,${MAKETARGET})))
  ifneq (,$(strip $(findstring kim,${MAKETARGET})))
    ERROR = "KIM is not compatible with TERSOFF force routines"
  endif
PP_FLAGS  += -DTERSOFF
FORCESOURCES += ${COVALENTSOURCES}
ifneq (,$(strip $(findstring tersoff2,${MAKETARGET})))
PP_FLAGS  += -DTERSOFF2
endif
endif
endif

# BRENNER
ifneq (,$(strip $(findstring brenner,${MAKETARGET})))
PP_FLAGS  += -DBRENNER
FORCESOURCES += ${COVALENTSOURCES}
ifneq (,$(strip $(findstring brenner2,${MAKETARGET})))
PP_FLAGS  += -DBRENNER2
endif
endif

# KEATING
ifneq (,$(strip $(findstring keating,${MAKETARGET})))
PP_FLAGS  += -DKEATING
FORCESOURCES += ${COVALENTSOURCES}
endif

# LASER & LASERYZ
ifneq (,$(strip $(findstring laseryz,${MAKETARGET})))
        ifneq (,$(strip $(findstring laser,${MAKETARGET})))
		CFLAGS += -DLASERYZ -DLASER
		SOURCES += ${LASERSOURCES}
	else
	       	CFLAGS += -DLASERYZ -DLASER
		SOURCES += ${LASERSOURCES}
	endif
else
	ifneq (,$(strip $(findstring laser,${MAKETARGET})))
		PP_FLAGS += -DLASER
		SOURCES += ${LASERSOURCES}
	endif
endif

# PDECAY
ifneq (,$(findstring pdecay,${MAKETARGET}))
PP_FLAGS += -DPDECAY
endif

# TWO TEMPERATURE MODEL TTM
ifneq (,$(strip $(findstring ttm,${MAKETARGET})))
PP_FLAGS += -DTTM
SOURCES += ${TTMSOURCES}
endif

# STREITZ MINTMIRE MODEL SM
ifneq (,$(strip $(findstring sm,${MAKETARGET})))
  PP_FLAGS += -DSM
  SOURCES += ${SMSOURCES}
  ifeq (,$(strip $(findstring nbl,${MAKETARGET})))
    SOURCES += ${EWALDSOURCES}
  endif
endif

# UNIAX
ifneq (,$(strip $(findstring uniax,${MAKETARGET})))
PP_FLAGS  += -DUNIAX
FORCESOURCES = ${UNIAXSOURCES}
endif

# EWALD
ifneq (,$(strip $(findstring ewald,${MAKETARGET})))
PP_FLAGS  += -DEWALD
FORCESOURCES  += ${EWALDSOURCES}
endif

# FCS
ifneq (,$(strip $(findstring fcs,${MAKETARGET})))
PP_FLAGS      += -DUSEFCS
FORCESOURCES  += ${FCSSOURCES}
LIBS          += $(shell pkg-config --libs scafacos-fcs)
CFLAGS        += $(shell pkg-config --cflags scafacos-fcs)
endif

# VARCHG
ifneq (,$(strip $(findstring varchg,${MAKETARGET})))
CFLAGS  += -DVARCHG
endif

# DIPOLE
ifneq (,$(strip $(findstring dipole,${MAKETARGET})))
CFLAGS  += -DDIPOLE
endif

# BUCK
ifneq (,$(strip $(findstring buck,${MAKETARGET})))
CFLAGS  += -DBUCK
endif

# MORSE
ifneq (,$(strip $(findstring morse,${MAKETARGET})))
CFLAGS  += -DMORSE
endif

# EXTF
ifneq (,$(strip $(findstring extf,${MAKETARGET})))
CFLAGS  += -DEXTF
endif

# LJ - computed Lennard-Jones (for vector versions only)
ifneq (,$(strip $(findstring lj,${MAKETARGET})))
PP_FLAGS  += -DLJ
endif

ifneq (,$(strip $(findstring kim,${MAKETARGET})))
FORCESOURCES = imd_forces_kim.c
PP_FLAGS += -DKIM
LIBS          += -L${KIM_DIR}/KIM_API -lkim
CFLAGS        += -I${KIM_DIR}/KIM_API
endif

SOURCES += ${FORCESOURCES}

###  ENSEMBLES  #########################################

ifneq (,$(findstring nve,${MAKETARGET}))
PP_FLAGS += -DNVE
endif

ifneq (,$(findstring mik,${MAKETARGET}))
PP_FLAGS += -DMIK
endif


# CG
ifneq (,$(strip $(findstring cg,${MAKETARGET})))
SOURCES += ${CGSOURCES}
IMDHEADERS += cg_util.h
PP_FLAGS  += -DCG
endif
ifneq (,$(strip $(findstring acg,${MAKETARGET})))
PP_FLAGS  += -DACG
endif

ifneq (,$(findstring nvt,${MAKETARGET}))
PP_FLAGS += -DNVT
endif

ifneq (,$(findstring npt_iso,${MAKETARGET}))
PP_FLAGS += -DNPT -DNPT_iso
endif

ifneq (,$(findstring npt_axial,${MAKETARGET}))
PP_FLAGS += -DNPT -DNPT_axial
endif

ifneq (,$(findstring frac,${MAKETARGET}))
PP_FLAGS += -DFRAC
# SOURCES += ${FRACSOURCES}
endif

ifneq (,$(findstring damp,${MAKETARGET}))
PP_FLAGS += -DDAMP
endif

ifneq (,$(findstring ftg,${MAKETARGET}))
PP_FLAGS += -DFTG
endif

ifneq (,$(findstring finnis,${MAKETARGET}))
PP_FLAGS += -DFINNIS
endif

ifneq (,$(findstring stm,${MAKETARGET}))
PP_FLAGS += -DSTM
endif

###  OPTIONS  ############################################

# vector mode
ifneq (,$(findstring vec,${MAKETARGET}))
PP_FLAGS += -DVEC

ifneq (,$(findstring vec2,${MAKETARGET}))
PP_FLAGS += -DVEC2
else
ifneq (,$(findstring vec3,${MAKETARGET}))
PP_FLAGS += -DVEC3
endif
endif

endif

# nudged elastic band method
ifneq (,$(findstring neb,${MAKETARGET}))
SOURCES  += imd_neb.c
PP_FLAGS += -DNEB
CC        = ${CC_MPI}
endif

# common-neighbour analysis
ifneq (,$(findstring cna,${MAKETARGET}))
SOURCES += ${CNASOURCES}
SOURCES += ${COVALENTSOURCES}
PP_FLAGS += -DCNA
endif

# socket interface
ifneq (,$(strip $(findstring sock,${MAKETARGET})))
HEADERS += ${SOCKHEADERS}
SOURCES += ${SOCKSOURCES}
PP_FLAGS  += -DSOCKET_IO
endif

# timing
ifneq (,$(findstring timing,${MAKETARGET}))
PP_FLAGS += -DTIMING
endif

ifneq (,$(findstring and,${MAKETARGET}))
PP_FLAGS += -DAND
endif

ifneq (,$(findstring ber,${MAKETARGET}))
PP_FLAGS += -DBER
endif

ifneq (,$(findstring fbc,${MAKETARGET}))
PP_FLAGS += -DFBC
endif

ifneq (,$(findstring bend,${MAKETARGET}))
PP_FLAGS += -DBEND
PP_FLAGS += -DFBC
endif

ifneq (,$(findstring bboost,${MAKETARGET}))
PP_FLAGS += -DBBOOST
SOURCES += ${BBOOSTSOURCES}
endif

ifneq (,$(findstring debugLo,${MAKETARGET}))
PP_FLAGS += -DdebugLo
endif


ifneq (,$(findstring zapp,${MAKETARGET}))
PP_FLAGS += -DZAPP
endif

ifneq (,$(findstring flagatoms,${MAKETARGET}))
PP_FLAGS += -DFLAGEDATOMS
endif


ifneq (,$(findstring rigid,${MAKETARGET}))
PP_FLAGS += -DRIGID
endif

ifneq (,$(findstring sendrec,${MAKETARGET}))
PP_FLAGS += -DSR
endif

ifneq (,$(findstring einstein,${MAKETARGET}))
PP_FLAGS += -DEINSTEIN
endif

ifneq (,$(findstring fnorm,${MAKETARGET}))
PP_FLAGS += -DFNORM
endif

ifneq (,$(findstring getsaddle,${MAKETARGET}))
PP_FLAGS += -DGETSADDLE
endif
ifneq (,$(findstring getmin,${MAKETARGET}))
PP_FLAGS += -DGETMIN
endif


ifneq (,$(findstring relaxinfo,${MAKETARGET}))
PP_FLAGS += -DRELAXINFO
endif

ifneq (,$(findstring norhoh,${MAKETARGET}))
PP_FLAGS += -DNORHOH
endif

ifneq (,$(findstring glok,${MAKETARGET}))
PP_FLAGS += -DGLOK
endif

ifneq (,$(findstring adaptglok,${MAKETARGET}))
PP_FLAGS += -DADAPTGLOK
endif

ifneq (,$(findstring mix,${MAKETARGET}))
PP_FLAGS += -DMIX
endif

ifneq (,$(findstring fire,${MAKETARGET}))
PP_FLAGS += -DMIK
PP_FLAGS += -DMIX
PP_FLAGS += -DADAPTGLOK
PP_FLAGS += -DGLOK
PP_FLAGS += -DRELAXINFO
PP_FLAGS += -DFNORM
endif

ifneq (,$(findstring efilter,${MAKETARGET}))
PP_FLAGS += -DEFILTER
endif

ifneq (,$(findstring clone,${MAKETARGET}))
PP_FLAGS += -DCLONE
endif

ifneq (,$(findstring nnbr,${MAKETARGET}))
PP_FLAGS += -DNNBR
endif

ifneq (,$(findstring writef,${MAKETARGET}))
PP_FLAGS += -DWRITEF
endif

ifneq (,$(findstring deform,${MAKETARGET}))
PP_FLAGS += -DDEFORM
ifeq (,$(findstring homdef,${MAKETARGET}))
SOURCES += ${DEFORMSOURCES}
endif
endif

ifneq (,$(findstring cycle,${MAKETARGET}))
PP_FLAGS += -DCYCLE
PP_FLAGS += -DHOMDEF
SOURCES += ${DEFORMSOURCES}
endif
ifneq (,$(findstring homdef,${MAKETARGET}))
PP_FLAGS += -DHOMDEF
SOURCES += ${DEFORMSOURCES}
endif

ifneq (,$(findstring shock,${MAKETARGET}))
PP_FLAGS += -DSHOCK
endif

ifneq (,$(findstring stress,${MAKETARGET}))
PP_FLAGS += -DSTRESS_TENS
endif

ifneq (,$(findstring quasi,${MAKETARGET}))
PP_FLAGS += -DQUASI
SOURCES += ${QUASISOURCES}
endif

ifneq (,$(findstring disloc,${MAKETARGET}))
PP_FLAGS += -DDISLOC
endif

ifneq (,$(findstring sllod,${MAKETARGET}))
PP_FLAGS += -DSLLOD -DNVT
endif

ifneq (,$(findstring avpos,${MAKETARGET}))
PP_FLAGS += -DAVPOS
endif

ifneq (,$(findstring force,${MAKETARGET}))
PP_FLAGS += -DFORCE
endif

ifneq (,$(findstring nmoldyn,${MAKETARGET}))
PP_FLAGS += -DNMOLDYN
endif

ifneq (,$(findstring dsf,${MAKETARGET}))
PP_FLAGS += -DDSF
LFLAGS += -lfftw3
endif

ifneq (,$(findstring atdist,${MAKETARGET}))
PP_FLAGS += -DATDIST
endif

ifneq (,$(findstring diffpat,${MAKETARGET}))
PP_FLAGS += -DDIFFPAT -I ${FFTW_DIR}/include
ifneq (,$(findstring omp,${MAKETARGET}))
LIBS   += -L ${FFTW_DIR}/lib -lfftw3f_threads -lfftw3f -lgomp
else
LIBS   += -L ${FFTW_DIR}/lib -lfftw3f
endif
endif

ifneq (,$(findstring ordpar,${MAKETARGET}))
PP_FLAGS += -DORDPAR
endif

# EPITAX
ifneq (,$(strip $(findstring epitax,${MAKETARGET})))
PP_FLAGS  += -DEPITAX
SOURCES += ${EPITAXSOURCES}
endif

# FEFL
ifneq (,$(findstring fefl,${MAKETARGET}))
PP_FLAGS += -DFEFL
SOURCES += ${FEFLSOURCES}
endif

# Correlation
ifneq (,$(strip $(findstring corr,${MAKETARGET})))
PP_FLAGS  += -DCORRELATE
SOURCES += ${CORRSOURCES}
endif

# heat transport
ifneq (,$(findstring nvx,${MAKETARGET}))
  PP_FLAGS += -DNVX
  SOURCES  += ${TRANSSOURCES}
else
  ifneq (,$(strip $(findstring hc,${MAKETARGET})))
    PP_FLAGS += -DHC
    SOURCES  += ${TRANSSOURCES}
  endif
endif

# mean square displacement
ifneq (,$(strip $(findstring msqd,${MAKETARGET})))
PP_FLAGS  += -DMSQD
SOURCES += ${CORRSOURCES}
endif

# MONOLJ Case
ifneq (,$(findstring monolj,${MAKETARGET}))
PP_FLAGS += -DMONOLJ
endif

# Single precision
ifneq (,$(findstring single,${MAKETARGET}))
PP_FLAGS += -DSINGLE
endif

# monoatomic system (performance tweak)
ifneq (,$(findstring mono,${MAKETARGET}))
PP_FLAGS += -DMONO
endif

# high precision output (checkpoints)
ifneq (,$(findstring hpo,${MAKETARGET}))
PP_FLAGS += -DHPO
endif

# use 4-point-interpolation
ifneq (,$(findstring 4point,${MAKETARGET}))
PP_FLAGS += -DFOURPOINT
endif

# use spline interpolation
ifneq (,$(findstring spline,${MAKETARGET}))
PP_FLAGS += -DSPLINE
endif

# use papi
ifneq (,$(findstring papi,${MAKETARGET}))
PP_FLAGS += -DPAPI ${PAPI_INC}
LIBS   += ${PAPI_LIBS}
endif

# use external potential for indenters
ifneq (,$(findstring extpot,${MAKETARGET}))
CFLAGS += -DEXTPOT
SOURCES += imd_extpot.c
endif

# directory for python modules
ifeq (,${IMDPYDIR})
  PYDIR = ${HOME}/python/imd
else
  PYDIR = ${IMDPYDIR}
endif

# Added by Frank Pister
ifneq (,$(findstring cbe,${MAKETARGET}))
SOURCES += spu.c imd_cbe_calc_ppu.c imd_cbe_util.c
endif

ifneq (,$(findstring cbe2,${MAKETARGET}))
PP_FLAGS  += -DCBE2
endif

ifneq (,$(findstring ppu,${MAKETARGET}))
PP_FLAGS += -DON_PPU
endif

#Angle distribution & slip vector analysis
ifneq (,$(findstring ada,${MAKETARGET}))
PP_FLAGS += -DADA
SOURCES += imd_ada.c
ifeq (,$(findstring ${COVALENTSOURCES},${SOURCES}))
	SOURCES += ${COVALENTSOURCES}
endif
endif

ifneq (,$(findstring nye,${MAKETARGET}))
PP_FLAGS += -DNYETENSOR
SOURCES += imd_nyeTensorAnalysis_3d.c
ifeq (,$(findstring ${COVALENTSOURCES},${SOURCES}))
	SOURCES += ${COVALENTSOURCES}
endif
ifeq (,$(findstring imd_ada.c,${SOURCES}))
	SOURCES += imd_ada.c
	PP_FLAGS += -DADA
endif
endif

# processor affinity
ifneq (,$(findstring aff,${MAKETARGET}))
PP_FLAGS += -DAFF
endif


# Substitute .o for .c to get the names of the object files
OBJECTS := $(subst .c,.o,${SOURCES})


###########################################################################
#
#	 Rules
#
###########################################################################

# all objects depend on headers
${OBJECTS}: ${HEADERS}

# How to compile *.c files
# special rules for force computation
imd_forces.o: imd_forces.c
	${CC} ${CFLAGS} ${PP_FLAGS} ${RCD_FLAGS} -c imd_forces.c

imd_forces_nbl.o: imd_forces_nbl.c
	${CC} ${CFLAGS} ${PP_FLAGS} ${RCD_FLAGS} ${NOALIAS} -c imd_forces_nbl.c

# Uncommented by Frank Pister
imd_forces_cbe.o: imd_forces_cbe.c
	${CC}   ${CFLAGS} ${PP_FLAGS} -c imd_forces_cbe.c

# rule to make imd_forces_cbe.o

imd_forces_eam2.o: imd_forces_eam2.c
	${CC} ${CFLAGS} ${PP_FLAGS} ${RCD_FLAGS} -c imd_forces_eam2.c

imd_forces_covalent.o: imd_forces_covalent.c
	${CC} ${CFLAGS} ${PP_FLAGS} ${RCD_FLAGS} -c imd_forces_covalent.c

# generic compilation rule
.c.o:
	${CC} ${CFLAGS} ${PP_FLAGS} -c $<

# How to link
ifneq (,$(strip $(findstring py,${MAKETARGET})))
${MAKETARGET}: ${OBJECTS}
	@echo ${PYDIR} ${IMDPYDIR}
	swig ${PP_FLAGS} -python -module $@ -outdir ${PYDIR} imd.i
	${CC} ${CFLAGS} ${PP_FLAGS} -I/usr/include/python -c imd_wrap.c
	${CC} ${LFLAGS} -fPIC -shared -o _$@.so ${OBJECTS} imd_wrap.o ${LIBS}
	${MV} _$@.so ${PYDIR}; rm -f _$@.so
else
${MAKETARGET}: ${OBJECTS}
	${CC} ${LFLAGS} -o $@ ${OBJECTS} ${LIBS}
	${MV} $@ ${BIN_DIR}; rm -f $@
endif

# First recursion only set the MAKETARGET Variable
.DEFAULT:
ifneq (,${CC})
	./version.sh
	@echo \#define COMPILE_TARGET \"$@\" >> version.h
	${MAKE} MAKETARGET='$@' STAGE2
else
ifneq (,${IMDSYS})
	@echo "IMDSYS variable ${IMDSYS} is not recognized"
else
	@echo "IMDSYS variable is not set"
endif
endif

# Second recursion sets MAKETARGET variable and compiles
# An empty MAKETARGET variable would create an infinite recursion, so we check
STAGE2:
ifneq (,${ERROR})
	@echo -e "\nError: ${ERROR}\n"
else
ifneq (,${MAKETARGET})
	${MAKE} MAKETARGET='${MAKETARGET}' ${MAKETARGET}
else
	@echo 'No TARGET specified.'
endif
endif

###########################################################################
#
#	 Misc. TARGETs
#
###########################################################################

clean: clean_spu
	rm -f *.o *.u *~ \#* *.V *.T *.O *.il

# Remove SPU-related temp. files (e.g. assembler/timing files)
clean_spu:
	rm -f  spu spu1 spu2   spu.map   spu.s  imd_cbe_util_spu.s imd_cbe_calc_spu.s     spu.s.timing imd_cbe_util_spu.s.timing imd_cbe_calc_spu.s.timing *.lst



help:
	@echo "Usage: gmake imd[_<parallel>][_<option>[_<option>...]]"

socktest:
	gcc -o ${BINDIR}/socktest sockutil.c socktest.c









# SPU specific targets for CBE
# The variables SPUCXX, SPUCC must be set to the SPU C++ and C compilers
# e.g. to spu-c++ and spu-cc respectively.
# SPUCXXFLAGS and SPUCFLAGS are the corresponding flags passed to the
# compiler which may contain optimization or warning flags, for instance.
# SPUCXX and SPUCC may also be set to use the IBM  compilers.
# EMBEDSPU must contain the name of the ppu-embedspu executable which
# generates object files containing code for the SPU and which may
# be linked with the main program.



# The SPU binary containing the main program as well as some tools
# and an spu version of calc_wp
# The SPU object file which is linked to the main programm
spu.o: spu1
	${EMBEDSPU} hndle_cbe_calc spu1 spu.o



# 2 SPUlets created in a different manner

# Compile first translate C to assembler, then compile assembler sources
spu1:  spu.s  imd_cbe_util_spu.s  imd_cbe_calc_spu.s
	${SPUCC} -o spu1   ${SPUCFLAGS} ${SPULDFLAGS}     spu.s  imd_cbe_calc_spu.s  imd_cbe_util_spu.s

# Compile the C sources without generating assembler
spu2: spu.c imd_cbe_util.c imd_cbe_calc_spu.c imd_cbe.h  config.h
	${SPUCC} -o spu2 ${SPUCFLAGS} ${SPULDFLAGS} ${PP_FLAGS}   spu.c imd_cbe_util.c imd_cbe_calc_spu.c


# Generate sssembler code for some SPU sources
# This assembler code may be analyzed using the spu-timing tool
imd_cbe_calc_spu.s: imd_cbe_calc_spu.c imd_cbe.h  config.h
	${SPUCC} -o imd_cbe_calc_spu.s  ${SPUCFLAGS} ${PP_FLAGS}   -S imd_cbe_calc_spu.c

imd_cbe_util_spu.s: imd_cbe_util.c imd_cbe.h  config.h
	${SPUCC} -o imd_cbe_util_spu.s  ${SPUCFLAGS} ${PP_FLAGS}    -S imd_cbe_util.c

spu.s: spu.c imd_cbe.h  config.h
	${SPUCC} -o spu.s ${SPUCFLAGS} ${PP_FLAGS}  -S  spu.c


# Timing of SPU assembly code (for profiling purposes)
spu.s.timing: spu.s
	${SPUTIMING} spu.s

imd_cbe_calc_spu.s.timing: imd_cbe_calc_spu.s
	${SPUTIMING} imd_cbe_calc_spu.s

imd_cbe_util_spu.s.timing: imd_cbe_util_spu.s
	${SPUTIMING} imd_cbe_util_spu.s


spu_timing: spu.s.timing  imd_cbe_calc_spu.s.timing  imd_cbe_util_spu.s.timing
