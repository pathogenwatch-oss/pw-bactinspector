#!/bin/bash

set -e -o pipefail

export PYTHONPATH=/:/bactinspector/

echoerr() { printf "%s\n" "$*" >&2; }

cat - > /tmp/input.fasta

reference=$( mash dist /tmp/input.fasta /bactinspector/libraries/filter.k21s1000.msh -d 0.075 | sort -gk3 | head -n 1 | cut -f 2 )

echoerr ${reference}

next_params=$( python3 /bactinspector/get_genus.py ${reference} )

next_library=$( echo "${next_params}" | cut -f 1 )
next_threshold=$( echo "${next_params}" | cut -f 2 )

echoerr ${next_library}

python3 /bactinspector/run_bactinspector.py check_species \
        -f input.fasta \
        -i /tmp/ \
        -s \
        -l /bactinspector/libraries/${next_library}.k21s1000 \
        -d ${next_threshold} \
        > /tmp/result.out

echoerr $( cat /tmp/result.out )

cat /tmp/result.out | python3 /bactinspector/map_to_pw.py