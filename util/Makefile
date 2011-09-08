############################################################################
#
# IMD -- The ITAP Molecular Dynamics Program
#
# Copyright 1996-2007 Institute for Theoretical and Applied Physics,
# University of Stuttgart, D-70550 Stuttgart
#
#
# IMD Utilities Makefile - works only with GNU make
#
# $Revision$
# $Date$
#
############################################################################

# Compilers
CFLAGS = ${FLAGS}
CC     = gcc
FLAGS  = -O

ifneq (,$(strip $(findstring alpha,${HOSTTYPE})))
CC      = cc
FLAGS  += -DALPHA -O3 -float -fp_reorder -arch ev56 -tune host -w
endif

ifneq (,$(strip $(findstring sparc,${HOSTTYPE})))
CC      = cc
FLAGS  += -fast -xtarget=ultra3 -xarch=v9b
endif

# where to put executable
BINDIR = ${HOME}/bin/${HOSTTYPE}

############################################################################

# sample constructor
make_config: make_config.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ make_config.c -lm

# sample converter
atcon: atcon.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ atcon.c -lm

# sample converter
procdata: procdata.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ procdata.c -lm

# analyse/convert atom density distribution files
atdist: atdist.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ atdist.c -lm

# analyse/convert diffraction pattern files
diffpat: diffpat.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ diffpat.c -lm

# convert VASP simulation to atom density distribution file
vasp2atdist: vasp2atdist.c param.c param.h
	${CC} ${CFLAGS} -o ${BINDIR}/$@ param.c vasp2atdist.c -lm

# convert VASP simulation to sequence of IMD checkpoints
vasp2chkpt: vasp2chkpt.c param.c param.h
	${CC} ${CFLAGS} -o ${BINDIR}/$@ param.c vasp2chkpt.c -lm

# convert VASP simulation to msqd fqile
vasp2msqd: vasp2msqd.c param.c param.h
	${CC} ${CFLAGS} -o ${BINDIR}/$@ param.c vasp2msqd.c -lm

# convert VASP charge density to volume file
chg2volume: chg2volume.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ chg2volume.c -lm

# histograms of square displacement files
sqd2hist: sqd2hist.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ sqd2hist.c -lm

# smoothly cut potential table
cutpot: cutpot.c param.c param.h
	${CC} ${CFLAGS} -o ${BINDIR}/$@ param.c cutpot.c -lm

# compute power spectrum in .eng file
imd_power: fft.c imd_power.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ imd_power.c -lm

#ppm_on_ppm
ppm_on_ppm: ppm_on_ppm.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ ppm_on_ppm.c -lm

#imd2ppm
imd2ppm: imd2ppm.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ imd2ppm.c -lm

#dynsf
dynsf: dynsf.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ dynsf.c -lfftw3 -lm

#conf2conf
conf2conf: conf2conf.c conf_tools.c conf_tools.h misc.c misc.h
	${CC} ${CFLAGS} -DC2C -o ${BINDIR}/$@ misc.c conf_tools.c conf2conf.c -lm

#hccorr
hccorr: hccorr.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ hccorr.c -lfftw3 -lm

#hcprof
hcprof: hcprof.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ hcprof.c

############################################################################
#
# distribution tools
#
###########################################################################

DISTSOURCE = dist_tools.c dist_tools.h

dist2volume: ${DISTSOURCE} dist2volume.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ dist2volume.c -lm

dist2gplt: ${DISTSOURCE} dist2gplt.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ dist2gplt.c -lm

dist2ppm: ${DISTSOURCE} dist2ppm.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ dist2ppm.c -lm

dist2ascii: ${DISTSOURCE} dist2ascii.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ dist2ascii.c -lm

dist2hist: ${DISTSOURCE} dist2hist.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ dist2hist.c -lm

dist_extract: ${DISTSOURCE} dist_extract.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ dist_extract.c -lm

# quick hack; superseeded by dist2volume
dist2rvf: dist2rvf.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ dist2rvf.c -lm

# make distribution from checkpoint file; quick hack
sample2dist: sample2dist.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ sample2dist.c -lm


###########################################################################
#
# imd_pair and its decendants
#
###########################################################################

SOURCES = misc.c param.c util_param.c util_common.c conf_tools.c
HEADERS = misc.h util.h param.h conf_tools.h
LIBS    = -lm

ifneq (,$(strip $(findstring 2d,${MAKETARGET})))
FLAGS += -DTWOD
endif

ifneq (,$(strip $(findstring angle,${MAKETARGET})))
FLAGS += -DANGLE
SOURCES += imd_angle.c
endif

ifneq (,$(strip $(findstring cna,${MAKETARGET})))
FLAGS += -DCNA -DCOVALENT
SOURCES += imd_cna.c
endif

ifneq (,$(strip $(findstring conn,${MAKETARGET})))
FLAGS += -DCONN
SOURCES += imd_conn.c
endif

ifneq (,$(strip $(findstring coord,${MAKETARGET})))
FLAGS += -DCOORD
SOURCES += imd_coord.c
ifneq (,$(strip $(findstring tersoff,${MAKETARGET})))
FLAGS += -DTERSOFF
endif
endif

ifneq (,$(strip $(findstring elco,${MAKETARGET})))
FLAGS += -DELCO
SOURCES += imd_elco.c imd_stress.c

ifneq (,$(strip $(findstring eam,${MAKETARGET})))
FLAGS += -DEAM -DPAIR_POT
SOURCES += imd_elco_eam.c imd_elco_pair.c
HEADERS += potaccess.h

else 
ifneq (,$(strip $(findstring tersoff,${MAKETARGET})))
FLAGS += -DTERSOFF -DCOVALENT 
SOURCES += imd_elco_tersoff.c

else 
ifneq (,$(strip $(findstring stiweb,${MAKETARGET})))
FLAGS += -DSTIWEB -DCOVALENT 
SOURCES += imd_elco_stiweb.c

else 
ifneq (,$(strip $(findstring keating,${MAKETARGET})))
FLAGS += -DKEATING -DCOVALENT 
SOURCES += imd_elco_keating.c

else 
FLAGS += -DPAIR_POT
SOURCES += imd_elco_pair.c
HEADERS += potaccess.h

endif
endif
endif
endif

endif
ifneq (,$(strip $(findstring pair,${MAKETARGET})))
FLAGS += -DPAIR
SOURCES += imd_pair.c
endif

ifneq (,$(strip $(findstring remove,${MAKETARGET})))
FLAGS += -DREMAT
SOURCES += imd_remove_atoms.c
endif

ifneq (,$(strip $(findstring ps,${MAKETARGET})))
FLAGS += -DPS -DCOVALENT
SOURCES += imd_ps.c imd_ps_main.c
endif

ifneq (,$(strip $(findstring ring,${MAKETARGET})))
FLAGS += -DRING -DCOVALENT
SOURCES += imd_ring.c
ifneq (,$(strip $(findstring tersoff,${MAKETARGET})))
FLAGS += -DTERSOFF
endif
endif

ifneq (,$(strip $(findstring strain,${MAKETARGET})))
FLAGS += -DSTRAIN
SOURCES += imd_strain.c
endif

ifneq (,$(strip $(findstring stress,${MAKETARGET})))
FLAGS += -DSTRESS
SOURCES += imd_stress.c
endif

ifneq (,$(strip $(findstring torsion,${MAKETARGET})))
FLAGS += -DTORSION
SOURCES += imd_torsion.c
endif

###########################################################################
#
#	 Rules
#
###########################################################################

# Substitute .o for .c to get the names of the object files
OBJECTS := $(subst .c,.o,${SOURCES})

# all objects depend on headers
${OBJECTS}: ${HEADERS}

.SUFFIXES:

# generic compilation rule
%.o : %.c
	${CC} ${CFLAGS} -c $<

# How to link
${MAKETARGET}: ${OBJECTS} 
	${CC} ${FLAGS} -o $@ ${OBJECTS} ${LIBS}
	mv $@ ${BINDIR}; rm -f $@

# First recursion only set the MAKETARGET Variable
.DEFAULT:
	${MAKE} MAKETARGET='$@' STAGE2

# Second recursion sets MAKETARGET variable and compiles
# An empty MAKETARGET variable would create an infinite recursion, so we check
STAGE2:
ifneq (,${MAKETARGET})
	${MAKE} MAKETARGET='${MAKETARGET}' ${MAKETARGET}
else
	@echo 'No TARGET specified.'
endif

###########################################################################
#
#	 Misc. TARGETs
#
###########################################################################

# clean
clean:
	rm -f *.o *~ 

count:
	wc -l Makefile *.[ch] *.csh *.sh *_uv *.perl rav
