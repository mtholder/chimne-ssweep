#!/bin/sh
f=$1
echo '#NEXUS'
echo 'BEGIN CHIMNeSSweep;'
tail -n1 $f
echo 'END;'
