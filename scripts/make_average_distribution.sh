
#!/usr/bin/env bash

PREFIX=$1
ARGS=$@
INPUTBGS=$(echo $@ | cut -d' ' -f2-)

sort -m -k1,1 -k2,2n -o ${PREFIX}.bg $INPUTBGS
