#!/bin/sh

rm -f version.h
grep=`which grep`
date=`$grep '\$Date' *.[ch] Makefile | awk '{print $3}' | sort -r | head -1`
echo \#define DATE \"$date\" > version.h
