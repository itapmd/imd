#!/usr/local/bin/tcsh
# concatenate per-CPU output files
foreach file (*.0)
  # strip off .0 at the end of $file
  set new=`echo $file | sed -e 's/\.0$//'`
  touch $new.head
  cat $new.head $new.[0-9]* > $new
  rm -f $new.*
end
