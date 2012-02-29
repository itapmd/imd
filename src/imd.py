#!/usr/bin/python

##############################################################################
#
# IMD -- The ITAP Molecular Dynamics Program
#
# Copyright 1996-2011 Institute for Theoretical and Applied Physics,
# University of Stuttgart, D-70550 Stuttgart
#
##############################################################################

##############################################################################
# $Revision$
# $Date$
##############################################################################

import sys, imp, getopt
from time import time

# usage message
def usage():
   print '\n Usage:\n'
   print '   ', sys.argv[0], '<module> [-r <n>] -p <paramfile>\n'
   print '        <module>          IMD python module, e.g. py_imd_nbl_nve'
   print '        -r <n>            restart from checkpoint number <n>'
   print '        -p <paramfile>    IMD parameter file (required)\n'
   print ' Parallel execution on <np> processors:\n'
   print '    mpirun -np <np> bwpython', sys.argv[0], \
                     '<module> [-r <n>] -p <paramfile>\n'
   print '    <module> must then be compiled with MPI enabled\n'
   return

# parse command line
restart = 0
if len(sys.argv) > 3:
    mod  = sys.argv[1]
    args = sys.argv[2:]
    try:
        opts, args = getopt.getopt(args,'r:p:')
    except getopt.GetoptError:
        usage()
        sys.exit()
    for o, a in opts:
        if o == '-r':
            if a.isdigit():
                restart = int(a)
            else:
                usage()
                sys.exit()
        if o == '-p':
            paramfile = a
else:
    usage()
    sys.exit()

# load module
try:
    lst = imp.find_module(mod)
except ImportError:
    print 'module ' + mod + ' not found'
    usage()
    sys.exit()
try:
    IMD = imp.load_module(mod,lst[0],lst[1],lst[2])
except ImportError:
    print 'module ' + mod + ' not found'
    usage()
    sys.exit()
finally:
    lst[0].close()

# determine available options
funcs = dir(IMD)
def_EFILTER     = 'write_config_ef'       in funcs
def_NNBR        = 'write_config_nb'       in funcs
def_ATDIST      = 'write_atdist'          in funcs
def_DIFFPAT     = 'write_diffpat'         in funcs
def_DISLOC      = 'write_config_dem'      in funcs
def_AVPOS       = 'write_config_avpos'    in funcs
def_STRESS_TENS = 'write_config_press'    in funcs
def_MSQD        = 'write_msqd'            in funcs
def_NMOLDYN     = 'write_nmoldyn'         in funcs
def_HOMDEF      = 'lin_deform'            in funcs
def_DEFORM      = 'deform_sample'         in funcs
def_FBC         = 'init_fbc'              in funcs
def_SOCKET_IO   = 'init_socket'           in funcs
def_EWALD       = 'init_ewald'            in funcs
def_REFPOS      = 'init_refpos'           in funcs
def_CG          = 'cg_step'               in funcs
def_ACG         = 'acg_step'              in funcs
def_GLOK        = 'update_glok'           in funcs
def_NBLIST      = 'check_nblist'          in funcs
def_CORRELATE   = 'init_correl'           in funcs
def_RELAX       = 'check_relaxed'         in funcs
def_RIGID       = 'calc_superforces'      in funcs
def_FORCE       = 'write_config_force'    in funcs
def_WRITEF      = 'write_config_wf'       in funcs
def_TEMPCONTROL = 'increment_temperature' in funcs
def_MPI         = 'init_mpi'              in funcs

# initialize MPI
if def_MPI:
    IMD.init_mpi()

# all public IMD global variables
G = IMD.cvar

# print IMD version
if G.myid == 0:
    print IMD.COMPILE_TARGET + ' of ' + IMD.DATE

# set restart and simulation phase
G.imdrestart = restart
simulation   = 1

# read parameter file
finished = IMD.read_parameters(paramfile, simulation)

# setup potentials
IMD.setup_potentials()

# read or generate atoms
if '_' == G.infilename[0]:
    IMD.generate_atoms(G.infilename)
else:
    IMD.read_atoms(G.infilename)

# initializations done once
# CBE omitted
# EPITAX omitted
if def_EWALD:
    IMD.init_ewald()

if def_SOCKET_IO:
    if G.myid == 0 and G.socket_int > 0:
        IMD.init_socket()

if G.imdrestart == 0 and G.eng_int > 0:
    IMD.write_eng_file_header()

if def_REFPOS:
    IMD.init_refpos()

# CNA omitted
# LASER omitted
# TTM omitted

# loop over simulation phases
while simulation == 1 or not finished:

    if G.myid == 0:
        print 'starting simulation', simulation
    
    # read parameters for new simulation phase
    if simulation > 1:
        finished = IMD.read_parameters(paramfile, simulation)
        IMD.make_box()
    
    # initializations for new phase
    # SHOCK init omitted
    # TTM init omitted
    # LASER init omitted
    # FRAC or FTG omitted
    if def_FBC:
        IMD.init_fbc()
    if def_CORRELATE or def_MSQD:
        IMD.init_correl(G.ncorr_rmax, G.ncorr_tmax)
    if def_NMOLDYN:
        if G.nmoldyn_int > 0:
            IMD.init_nmoldyn()
    if def_ATDIST:
        if G.atdist_int > 0:
            init_atdist()
    if def_DIFFPAT:
        if G.diffpat_int > 0:
            init_diffpat()
    if def_CG:
        if G.ensemble == IMD.ENS_CG:
            IMD.reset_cg()
    if def_ACG:
        G.acg_alpha = G.acg_init_alpha
    if def_DEFORM:
        G.deform_int = 0
    if def_RELAX:
        G.is_relaxed = 0
    
    # we should use the IMD variables steps_min, steps_max and steps
    # for the minimal, maximal, and current step number, as these
    # are accessed (and sometimes even changed) directly by some IMD routines
    start_time = time()
    G.steps = G.steps_min
    while G.steps <= G.steps_max:
        # SHOCK omitted
        # determine whether pressure tensor shall be computed in this step
        if def_STRESS_TENS:
            flag = G.eng_int  > 0 and 0 == G.steps % G.eng_int  or \
                   G.dist_int > 0 and 0 == G.steps % G.dist_int or G.relax_rate > 0.0
            if flag:
                G.do_press_calc = 1
            else:
                G.do_press_calc = 0
        # EPITAX omitted
        # deform sample with HOMDEF option
        if def_HOMDEF:
            if G.lindef_int > 0 and 0 == G.steps % G.lindef_int:
                if IMD.DIM == 2:
                    IMD.lin_deform(G.lindef_x, G.lindef_y,             G.lindef_size)
                else:
                    IMD.lin_deform(G.lindef_x, G.lindef_y, G.lindef_z, G.lindef_size)
        # deform sample with DEFORM option
        if def_DEFORM:
            if G.max_deform_int > 0:
                if G.is_relaxed or G.deform_int == G.max_deform_int:
                    IMD.deform_sample()
                    G.deform_int = 0
                    if G.ensemble == IMD.ENS_CG:
                        IMD.reset_cg()
                G.deform_int += 1
        # CNA omitted
        # update average positions
        if def_AVPOS:
            if G.steps == G.steps_min or G.steps == G.avpos_start:
                IMD.update_avpos()
        # update atom distribution
        if def_ATDIST:
            if G.atdist_int > 0 and G.steps >= G.atdist_start and G.steps <= G.atdist_end:
                G.update_atdist()
        # update diffraction pattern
        if def_DIFFPAT:
            if G.diffpat_int > 0 and G.steps >= G.diffpat_start and G.steps <= G.diffpat_end:
                G.update_diffpat(G.steps)
        # compute updated boundary forces; these are then applied in integrator
        if def_FBC:
            IMD.update_fbc()
        
        # calculate forces
        # beware: FBC forces and mobility restrictions are applied in integrator
        # for CG-like integrators, forces are computed from integrator
        if G.ensemble == IMD.ENS_CG:
            if def_CG:
                IMD.cg_step(G.steps)
            elif def_ACG:
                IMD.acg_step(G.steps)
        else:
            IMD.calc_forces(G.steps)
        # calculate forces on superparticles
        if def_RIGID:
            IMD.calc_superforces()
        
        # write forces for potfit before moving atoms
        if def_FORCE:
            if G.force_int > 0 and 0 == G.steps % G.force_int:
                IMD.write_config_force( G.steps / G.force_int)
        # write forces on boundary particles (only for postprocessing)
        # this does not yet include FBC forces!
        if def_WRITEF: 
            IMD.write_config_wf(G.steps)
        # EPITAX omitted
        if def_DISLOC:
            if G.steps == G.reset_Epot_step and 1 == G.calc_Epot_ref:
                IMD.reset_Epot_ref()
        # CNA omitted
        # update correlation functions
        if def_CORRELATE or def_MSQD:
            if G.steps >= G.correl_start and G.steps < G.correl_end or 0 == G.correl_end:
                istep = G.steps - G.correl_start
                if 0 == istep % G.correl_ts:
                    IMD.correlate(G.steps, G.correl_refstep, istep / G.correl_ts)
                if 0 != G.correl_int and G.steps - G.correl_refstep + 1 >= G.correl_int: 
                    G.correl_refstep += G.correl_int
        # update state of (adaptive) GLOK integrator
        if G.ensemble == IMD.ENS_GLOK:
            IMD.update_glok()
        # LASER omitted
        # TTM omitted
        
        # move atoms
        # for CG-like ensembles, atoms are moved together with force computation
        if G.ensemble != IMD.ENS_CG:
            IMD.move_atoms()
        
        # EPITAX omitted
        # increment temperature
        if def_TEMPCONTROL:
            IMD.increment_temperature()
        
        # Periodic I/O
        # write chechpoint
        if G.checkpt_int  > 0 and 0 == G.steps % G.checkpt_int:
            IMD.write_config(G.steps / G.checkpt_int, G.steps)
        # write entry to .eng file
        if G.eng_int  > 0 and 0 == G.steps % G.eng_int:
            IMD.write_eng_file(G.steps)
        # write distribution file
        if G.dist_int > 0 and 0 == G.steps % G.dist_int:
            IMD.write_distrib(G.steps)
        # write picture file
        if G.pic_int > 0 and 0 == G.steps % G.pic_int:
            IMD.write_pictures(G.steps)
        # TTM omitted
        # write energy-filtered atoms
        if def_EFILTER:
            if G.ef_checkpt_int > 0 and 0 == G.steps % G.ef_checkpt_int:
                IMD.write_config_ef(G.steps / G.ef_checkpt_int)
        # write coordination-filtered atoms
        if def_NNBR:
            if G.nb_checkpt_int > 0 and 0 == G.steps % G.nb_checkpt_int:
                IMD.write_config_nb(G.steps / G.nb_checkpt_int)
        # write positions for atom distributions
        if def_ATDIST:
            if G.atdist_pos_int > 0 and 0 == G.steps % G.atdist_pos_int:
                IMD.write_config_atdist_pos(G.steps / G.atdist_pos_int)
        # write differential energy and displacements; update reference positions
        if def_DISLOC:
            if G.steps == G.up_ort_ref:
                IMD.update_ort_ref()
            if G.dem_int > 0 and 0 == G.steps % G.dem_int:
                IMD.write_config_dem(G.steps)
            if G.dsp_int > 0 and G.steps > G.up_ort_ref and 0 == G.steps % G.dsp_int:
                IMD.write_config_dsp(G.steps)
        # write and update average positions
        if def_AVPOS:
            if G.steps > G.avpos_start and G.steps <= G.avpos_end:
                stp = G.steps - G.avpos_start
                if G.avpos_res > 0 and 0 == stp % G.avpos_res:
                    IMD.add_positions()
                if G.avpos_int > 0 and 0 == stp % G.avpos_int:
                    IMD.write_config_avpos(stp / G.avpos_int)
        # NVX omitted
        # write atoms with per-atom stress tensor
        if def_STRESS_TENS:
            if G.press_int > 0 and 0 == G.steps % G.press_int:
                IMD.write_config_press(G.steps / G.press_int)
        # write entry to nMoldyn trajectory
        if def_NMOLDYN:
            if G.nmoldyn_int > 0 and 0 == G.steps % G.nmoldyn_int:
                IMD.write_nmoldyn(G.steps)
        # DSF omitted
        # check for request on socket
        if def_SOCKET_IO:
            if G.socket_int > 0 and 0 == G.steps % G.socket_int:
                IMD.check_socket()
        # change box to relax pressure
        if def_HOMDEF:
            if G.relax_rate > 0.0:
                IMD.relax_pressure()
        # check if sample is relaxed
        if def_RELAX:
            IMD.check_relaxed()
        # check neighbor list, apply PBC, fix cell distribution
        if def_NBLIST:
            IMD.check_nblist()
        else:
            IMD.fix_cells()
        # write atom distribution at end of simulation
        if def_ATDIST:
            if G.atdist_int > 0 and G.steps == G.atdist_end:
                IMD.write_atdist()
        # write diffraction pattern at end of simulation
        if def_DIFFPAT:
            if G.diffpat_int > 0 and G.steps == G.diffpat_end:
                IMD.write_diffpat()
        # write square displacements at end of simulation
        if def_MSQD:
            if G.correl_end > 0 and G.steps == G.correl_end or \
            G.correl_end == 0 and G.steps == G.steps_max:
                IMD.write_config_sqd(0)
        # check whether a checkpoint shall be written
        if G.watch_int > 0 and 0 == G.steps % G.watch_int:
            IMD.check_write();
        # check whether we shall stop simulation
        if G.stop_int > 0 and 0 == G.steps % G.stop_int:
            finished = IMD.check_stop()
            if finished:
                break
        # check whether we have run out of time
        if G.maxwalltime > 0:
            finished = IMD.check_walltime()
            if finished:
                break
        # increment time step
        G.steps += 1

    # write timing information at end of each simulation phase
    if G.myid == 0:
        print 'finished simulation', simulation
        atsteps = ((G.steps_max - G.steps_min) * G.natoms)
        diff_time = 1000000.0 * (time() - start_time) / atsteps
        print G.num_cpus * diff_time, "microseconds per step and atom"
    simulation += 1

IMD.close_files()
