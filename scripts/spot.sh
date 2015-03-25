#!/bin/bash

tags=$1
hotspots=$2

n=`bedops -e -1 $tags $hotspots | awk '{ t += $5; } END { print t; }'`
total=`bedops -u $tags | awk '{ t += $5; } END { print t; }'`

echo "SPOT: $(echo "scale=4; $n/$total" | bc)"