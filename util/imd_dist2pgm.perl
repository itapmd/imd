#!/client/bin/perl
# 
# imd_dist2pgm.perl konviertiert eine 2d Verteilungsfunktion 
# in ein Graymap
# 
# J. Stadler, Sommer 96
#             
# liest von STDIN, schreibt auf STDOUT

# Header einlese
$_     = <STDIN>;
@zeile = split;
$xres = $zeile[2];
$yres = $zeile[3];
	

# Werte in Array einlesen

for( $i=0; $i<$xres; $i++ ) {
    for( $j=0; $j<$yres; $j++ ) {
	$_ = <STDIN>;	

	s/\s*//;             
	@zeile = split;
	if ("-p" eq $ARGV[0] ) {
          $bitmap[$i][$j] = $zeile[0] > 0 ? $zeile[0] : - $zeile[0]; 
	} else {
	  $bitmap[$i][$j] =  $zeile[1]; 
        };				
    };
};
# Min/Max ausrechnen

$xmin=$bitmap[0][0]; 
$xmax=$bitmap[0][0]; 

if ("-p" eq $ARGV[0] ) {
    $factor = 10;		
} else {
    $factor = 1024;
};				

				
for( $i=0; $i<$xres; $i++ ) {
    for( $j=0; $j<$yres; $j++ ) {
    if ($bitmap[$i][$j] > $xmax) { $xmax = $bitmap[$i][$j]; };
    if ($bitmap[$i][$j] < $xmin) { $xmin = $bitmap[$i][$j]; };
};				
};

# Write pgm output

print("P2 $xres $yres \n");
printf("%d\n",$xmax * $factor);

for( $j=0; $j < $yres; $j++ ) {
    for( $i=0; $i < $xres; $i++ ) { printf("%d ",$bitmap[$i][$j] * $factor ) };	
    print("\n");
};

