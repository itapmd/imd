#!/usr/bin/perl -w
#
# script to convert IMD atom files to AtomEye extended cfg format
#
# Authors: Erik Bitzek, Franz G?hler
#
# currently supports only 3D ASCII atom files having a correct header
#
# $Revision$
# $Date$
#
# Usage:  imd2cfg.perl -i <infile> -t <types> <options>
#
# <types> is a space separated list of atom types, like: -t Ni Al
# Sufficiently many atom types must be specified. 
#
# Supported options are:
#
#  -v           include also velocities, if available (default: no velocities)
#  -s           sort according to atom types (default: no sort)
#  -Lx <x y z>  first  box vector (default: read from header)
#  -Ly <x y z>  second box vector (default: read from header)
#  -Lz <x y z>  third  box vector (default: read from header)

# The system of units used in the input file must be specified, 
# like in the following examples.
# Otherwise, AtomEye will give wrong information on units, but atoms
# are displayed anyway.

# Angstrom-eV-amu units
$mass_unit   = 1.0;         # mass unit in amu
$length_unit = 1.0;         # Length unit in Agstrom
$time_unit   = 0.01018;     # time unit in ps
$energy_unit = "eV";        # name of energy unit

# Angstrom-eV-ps units
# $mass_unit   = 0.00010365;  # mass unit in amu
# $length_unit = 1.0;         # Length unit in Agstrom
# $time_unit   = 1.0;         # time unit in ps
# $energy_unit = "eV";        # name of energy unit

# some variables we will use...
$sort=0;
$ntypes=0;
$natoms=0;
$incveloc=0;
@H0=0.0;

# conversion factor of velocities to A/fs
$veloc_unit  = "A/fs";                     # name of velocity unit
$v0 = 1000.0 * $time_unit / $length_unit;  # conversion factor of to A/fs

# usage message
$usage = "Usage: $0 -i <infile> -t <types> [-v] [-s] [-Lx <x y z>] [-Ly <x y z>] [-Lz <x y z>]";

# first we parse the command line
for ($i = 0; $i <= $#ARGV; $i++) {
    if ($ARGV[$i] =~ /-t/) {
        while (($i < $#ARGV) && ($ARGV[$i+1] !~ /-/)) {
            $atomtype[$ntypes++]=$ARGV[++$i];
        }
    }
    elsif ($ARGV[$i] =~ /-i/) {
        $infile = $ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /-v/) {
        $incveloc = 1;
    }
    elsif ($ARGV[$i] =~ /-s/) {
        $sort = 1;
    }
    elsif ($ARGV[$i] =~ /-Lx/) {
        for($j=1; $j<=3; $j++) {
            $H0[1][$j] = $ARGV[$i+$j];
            $i++;
        }
    }
    elsif ($ARGV[$i] =~ /-Ly/) {
        for($j=1; $j<=3; $j++) {
            $H0[2][$j]=$ARGV[$i+$j];
            $i++;
        }
    }
    elsif ($ARGV[$i]=~ /-Lz/) {
        for($j=1; $j<=3; $j++) {
            $H0[3][$j] = $ARGV[$i+$j];
            $i++;
        }
    }
    else {
        die "Unknown argument $ARGV[$i]\n" . $usage;
    }
}

# check input
if (!($ntypes)) {
    die "No atom types\n" . $usage;
}
if (!($infile)) {
    die "No input file\n" . $usage;
}

$outfile=$infile.".cfg";

open(IMD,$infile) or die "can't open iputfile $infile\n";

# check for information in the header of the imdfile
while ($line = <IMD>) {

    # parse header lines
    if ($line =~ /\#/) {

        @zeile=split(/[\ \t\n]+/,$line);

        # format line:
        if($zeile[0] eq "#F") {
            if($zeile[1] ne "A") {
                die "Sorry - binary files not (yet) supported\n";
            }
            $nnumber = $zeile[2];
            $nsorte  = $zeile[3];
            $nmass   = $zeile[4];
            $npos    = $zeile[5];
            $nveloc  = $zeile[6];
            $ndata   = $zeile[7];
            if ($npos ne 3) {
                die "3 atom coordinates are required!";
            }
        }

        # comments:
        if ($zeile[0] eq "#C") {
            $ncontents=0;
            while ($zeile[$ncontents+1]) {
                $contents[$ncontents] = $zeile[$ncontents+1];
                $ncontents++;
            }
            if (!($nnumber+$nsorte+$nmass+$npos+$nveloc+$ndata eq $ncontents)){
                die "Wrong number of items in contents line\n";
            }
        }

        # box vectors:
        if ($zeile[0] eq "#X") {
            $H0[1][1] = $zeile[1];
            $H0[1][2] = $zeile[2];
            $H0[1][3] = $zeile[3];
        }
        elsif ($zeile[0] eq "#Y") {
            $H0[2][1] = $zeile[1];
            $H0[2][2] = $zeile[2];
            $H0[2][3] = $zeile[3];
        }
        elsif ($zeile[0] eq "#Z") {
            $H0[3][1] = $zeile[1];
            $H0[3][2] = $zeile[2];
            $H0[3][3] = $zeile[3];
        }

    }

    # count atoms
    else {
        $natoms++;
    }
}
close IMD;

if (($H0[1][1]==0.0) || ($H0[2][2]==0.0) || ($H0[3][3]==0.0)) {
    die "information about box size is lacking";
}

# dual basis, first unnormalized
$B[1][1] = $H0[2][2] * $H0[3][3] - $H0[2][3] * $H0[3][2];
$B[1][2] = $H0[2][3] * $H0[3][1] - $H0[2][1] * $H0[3][3];
$B[1][3] = $H0[2][1] * $H0[3][2] - $H0[2][2] * $H0[3][1];

$B[2][1] = $H0[3][2] * $H0[1][3] - $H0[3][3] * $H0[1][2];
$B[2][2] = $H0[3][3] * $H0[1][1] - $H0[3][1] * $H0[1][3];
$B[2][3] = $H0[3][1] * $H0[1][2] - $H0[3][2] * $H0[1][1];

$B[3][1] = $H0[1][2] * $H0[2][3] - $H0[1][3] * $H0[2][2];
$B[3][2] = $H0[1][3] * $H0[2][1] - $H0[1][1] * $H0[2][3];
$B[3][3] = $H0[1][1] * $H0[2][2] - $H0[1][2] * $H0[2][1];

# box volume
$det = $B[1][1] * $H0[1][1] + $B[1][2] * $H0[1][2] + $B[1][3] * $H0[1][3];

# normalize dual basis
$B[1][1] = $B[1][1] / $det;
$B[1][2] = $B[1][2] / $det;
$B[1][3] = $B[1][3] / $det;

$B[2][1] = $B[2][1] / $det;
$B[2][2] = $B[2][2] / $det;
$B[2][3] = $B[2][3] / $det;

$B[3][1] = $B[3][1] / $det;
$B[3][2] = $B[3][2] / $det;
$B[3][3] = $B[3][3] / $det;

open(OUT,">$outfile");

# write header of extended cfg file
print OUT "Number of particles = $natoms\n";
print OUT "A = $length_unit Angstrom (basic length-scale)\n";
for ($i=1; $i<=3; $i++) {
    for($j=1; $j<=3; $j++) {
        printf(OUT "H0(%d,%d) = %f A\n", $i, $j, $H0[$i][$j]);
    }
}
# printf(OUT "R = %f [ps^-1]\n", 1.0/$time_unit);

# write property information to header of extended cfg file
$off  = $nnumber + $nsorte + $nmass + $npos;
$naux = 0;
print OUT ".NO_VELOCITY.\n";
if (($nveloc) && ($incveloc)) {
    $ec  = $nnumber + $npos + $nveloc + 1 + $ndata;
    print  OUT "entry_count = $ec\n";
    printf(OUT "auxiliary[%d] = %s [%s]\n", $naux++, "vx",   $veloc_unit);
    printf(OUT "auxiliary[%d] = %s [%s]\n", $naux++, "vy",   $veloc_unit);
    printf(OUT "auxiliary[%d] = %s [%s]\n", $naux++, "vz",   $veloc_unit);
    printf(OUT "auxiliary[%d] = %s [%s]\n", $naux++, "Ekin", $energy_unit);
    $off = $off + $nveloc;
}
else {
    $off = $off + $nveloc;
    if ($nveloc) {
        $ec  = $nnumber + $npos + 1 + $ndata;
        print  OUT "entry_count = $ec\n";
        printf(OUT "auxiliary[%d] = %s [%s]\n", $naux++, "Ekin", $energy_unit);
    }
    else {
        $ec  = $nnumber + $npos + $ndata;
        print  OUT "entry_count = $ec\n";
    }
}
for ($i=0; $i<$ndata; $i++) {
    if ($contents[$off] =~ /Epot/) {
        $unit = $energy_unit;
    }
    else {
        $unit = "reduced units";
    }
    printf(OUT "auxiliary[%d] = %s [%s]\n", $naux++, $contents[$off++], $unit);
}
print OUT "auxiliary[$naux] = number []\n";

# sorting of file to reduce data size
if ($sort) {
    $infiles = $infile.".s";
    system "sort -n -k 2 $infile -o $infiles";
    open(IMD2,$infiles) or die "can't open iputfile  $infiles\n";
}
else {
    open(IMD2,$infile) or die "can't open iputfile  $infile\n";
}
$lasttyp=-1;

while (($line = <IMD2>)) {

    if ($line !~ /\#/) {

        # default values
        $mass = 1.0; $typ = 0;
        @zeile=split(/[\ \t\n]+/,$line);

        # imd format:
        $i  = 0;
        if ($nnumber) { 
            $nr   = $zeile[$i++]; 
        }
        if ($nsorte) { 
            $typ  = $zeile[$i++];
            if ($typ >= $ntypes) {
                die "not enough atom types!";
            }
        } 
        if ($nmass) { 
            $mass = $zeile[$i++];
        }
        $x = $zeile[$i++];
        $y = $zeile[$i++];
        $z = $zeile[$i++];
        if ($nveloc) {
            $vx = $zeile[$i++];
            $vy = $zeile[$i++];
            $vz = $zeile[$i++];
            $ekin = 0.5 * $mass * ($vx * $vx + $vy * $vy + $vz * $vz);
        }
        for ($j=0; $j<$ndata; $j++) {
            $aux[$j] = $zeile[$i++];
        }

        # positions in reduced coordinates
        $s1 = $x  * $B[1][1] + $y  * $B[1][2] + $z  * $B[1][3]; 
        $s2 = $x  * $B[2][1] + $y  * $B[2][2] + $z  * $B[2][3]; 
        $s3 = $x  * $B[3][1] + $y  * $B[3][2] + $z  * $B[3][3]; 

        if ($typ != $lasttyp) {
            printf (OUT "%f \n%s\n", $mass/$mass_unit, $atomtype[$typ]);
            $lasttyp=$typ;
        }
        printf (OUT "%f %f %f ", $s1, $s2, $s3);
        if ($nveloc) {
            if ($incveloc) {
                printf (OUT "%f %f %f ", $vx*$v0, $vy*$v0, $vz*$v0);
            }
            printf(OUT "%f ", $ekin);
        }
        for ($i=0; $i<$ndata; $i++) {
            printf (OUT "%f ", $aux[$i]);
        }
        if ($nnumber) {
             printf(OUT "%d ", $nr);
        }
        print OUT "\n";
    }
}

if ($sort) {
    system "rm $infiles";
}
