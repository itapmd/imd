#!/bin/sh
# concatenate per-CPU output files
for file in *.0; do
  # strip off .0 at the end of $file
  new=`echo $file | sed -e 's/\.0$//'`
  head=`echo $new | sed -e 's/\.chkpt$/\.head/'`
  touch $head
  cat $head $new.* > $new
  rm -f $new.*
  rm -f $head
done
