#!/bin/sh
# concatenate per-CPU output files
for file in *.0; do
  # strip off .0 at the end of $file
  new=`echo $file | sed -e 's/\.0$//'`
  touch $new.head
  cat $new.head $new.[0-9]* > $new
  rm -f $new.*
done
