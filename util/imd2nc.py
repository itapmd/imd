#!/usr/bin/python

###############################################################################
#
# Convert an IMD trajectory file to a netCDF file for nMolDyn
#
###############################################################################

# import the modules we need
import sys, string, time              # basic services
from Numeric import *                 # easy access to Numeric
from Scientific.IO.NetCDF import *    # easy access to NetCDF
from Scientific.IO.TextFile import TextFile
from MMTK import *
from MMTK.Trajectory import Trajectory, TrajectoryOutput
from MMTK.Trajectory import SnapshotGenerator 
from MMTK.Trajectory import LogOutput, TrajectoryVariable
from MMTK.ParticleProperties import ParticleVector

# usage message
def usage():
   print 'Usage:', sys.argv[0], '<types> <trajectory>'
   print
   print '   <types>:       comma-separated list of chemical symbols'
   print
   print '   <trajectory>:  IMD trajectory file (*.nmoldyn)'
   print
   print 'Example:  imd2nc.py Al,Ni trajectory.nmoldyn'
   print
   return

# check number of arguments
if len(sys.argv) != 3:
   print 'Wrong number of arguments.\n'
   usage()
   sys.exit()

# defaults for units - customize if necessary
tu = 0.01018         # time unit of input in ps
lu = 0.1             # length unit of input in nm
vu = lu / tu         # velocity unit of input in nm/ps

# get command line arguments
atom_types      = sys.argv[1].split(',')
trajectory_file = sys.argv[2]
output_file     = trajectory_file + '.nc'

# number of atoms and types
try:
   traj_file = TextFile(trajectory_file, 'r')
except:
   print '\nCould not open file', trajectory_file,'\n'
   usage()
   sys.exit()
line      = traj_file.readline()
numtypes  = map( string.atoi, line.split())
ntypes    = len(numtypes)
natoms    = sum(numtypes)
line      = traj_file.readline()
box       = map( string.atof, line.split())
have_box  = len(box) == 3
if have_box:
   do_pbc = True
   line   = traj_file.readline()
line      = traj_file.readline()
nitems    = len( line.split() )
traj_file.close()
if len(atom_types) != len(numtypes):
   print 'Wrong number of atoms types.\n'
   usage()
   sys.exit()

# build universe (periodic, of possible)
if do_pbc:
   box = ( box[0]*lu, box[1]*lu, box[2]*lu )
   universe = OrthorhombicPeriodicUniverse(box,None)
   print 'Building periodic universe...'
else:
   universe = InfiniteUniverse(None)
   print 'Building infinite universe...'

# create atoms, with positions initialized
for t in range(ntypes):    
   for i in range(numtypes[t]):
      universe.addObject(Atom(atom_types[t],position=Vector(0.0,0.0,0.0)))

trajectory = Trajectory(universe, output_file,'w')
atoms = universe.atomList()

# if we have velocities, add them to the universe
if nitems == 6:
   vel = ParticleVector(universe)
   for i in range(natoms):
      vel[i] = Vector( 0.0, 0.0, 0.0 )
   universe.setVelocities(vel)
   snapshot = SnapshotGenerator(universe,
              actions=[TrajectoryOutput(trajectory,
              ['configuration','velocities','time'], 0, None, 1)])
else:
   snapshot = SnapshotGenerator(universe,
              actions=[TrajectoryOutput(trajectory,
              ['configuration','time'], 0, None, 1)])

# open trajectory again
traj_file = TextFile (trajectory_file, 'r')
traj_file.readline()
if have_box:
   traj_file.readline()

# convert trajectory
finished = False
while finished == False:
   try:
      line = traj_file.readline().split()
   except IOError:
      finished = True
      break
   if len(line) == 0:
      finished = True
      break
   if nitems == 6:   # with velocities
      for i in range(natoms):
         lst = map(string.atof, traj_file.readline().split())
         x, y, z, vx, vy, vz = lst
         atoms[i].setPosition( Vector( x*lu,   y*lu,  z*lu ) )
         vel[i] = Vector( vx*vu, vy*vu, vz*vu )
      universe.setVelocities(vel)
   else:             # without velocities
      for i in range(natoms):
         x,y,z = map(string.atof, traj_file.readline().split())
         atoms[i].setPosition( Vector( x*lu,   y*lu,  z*lu ) )
   time = string.atof(line[0]) * tu
   snapshot(data={'time':time})

# close files
traj_file.close()    
trajectory.close()

