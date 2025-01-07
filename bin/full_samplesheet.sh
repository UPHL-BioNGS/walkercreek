#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
script=$SCRIPT_DIR/create_samplesheet.pl
echo "sample,r1,r2,r3,r4"

for VAR in "$@"
do
    $script -i $VAR -f
    for DIR in `find $VAR -type d`; do
        $script -i $DIR -f
    done
done | sort | uniq

# Taken from https://github.com/CDCgov/mycosnp-nf/blob/master/bin/mycosnp_full_samplesheet.sh
