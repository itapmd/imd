
CC     = gcc
CFLAGS = -O -g3 ${FLAGS}
BINDIR = ${HOME}/bin/${HOSTTYPE}

# sample constructor/converter
atcon: atcon.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ atcon.c -lm

# sample converter
procdata: procdata.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ procdata.c -lm

# analyse/convert atom density distribution files
atdist: atdist.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ atdist.c -lm

# histograms of square displacement files
sqd2hist: sqd2hist.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ sqd2hist.c -lm

# smoothly cut potential table
cutpot: cutpot.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ cutpot.c -lm

# compute power spectrum in .eng file
imd_power: fft.c imd_power.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ imd_power.c -lm


##########################################################
#
# distribution tools
#
DISTSOURCE = dist_tools.c dist_tools.h

dist2volume: ${DISTSOURCE} dist2volume.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ dist2volume.c -lm

dist2ppm: ${DISTSOURCE} dist2ppm.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ dist2ppm.c -lm

dist2ascii: ${DISTSOURCE} dist2ascii.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ dist2ascii.c -lm

dist2hist: ${DISTSOURCE} dist2hist.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ dist2hist.c -lm

dist_extract: ${DISTSOURCE} dist_extract.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ dist_extract.c -lm

# quick hack; needs to be replaced by better version
dist2rvf: dist2rvf.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ dist2rvf.c -lm

# make distribution from checkpoint file; quick hack
sample2dist: sample2dist.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ sample2dist.c -lm


###########################################################
#
# imd_pair and its decendands
#
imd_pair.c: util.h imd_pair.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ imd_pair.c -lm

imd_angle.c: util.h imd_angle.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ imd_angle.c -lm

imd_conn.c: util.h imd_conn.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ imd_conn.c -lm

imd_coord.c: util.h imd_coord.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ imd_coord.c -lm

imd_elco.c: util.h imd_elco.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ imd_elco.c -lm

imd_ps.c: util.h imd_ps_main.c imd_ps.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ imd_ps.c -lm

imd_ring.c: util.h imd_ring.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ imd_ring.c -lm

imd_strain.c: util.h imd_strain.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ imd_strain.c -lm

imd_stress.c: util.h imd_stress.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ imd_stress.c -lm

imd_torsion.c: util.h imd_torsion.c
	${CC} ${CFLAGS} -o ${BINDIR}/$@ imd_torsion.c -lm

# clean
clean:
	rm -f *.o *~ 

