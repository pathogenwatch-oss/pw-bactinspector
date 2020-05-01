#!/bin/bash

export PYTHONPATH=/:/bactinspector/

cat - > /tmp/input.fasta

python3 /bactinspector/run_bactinspector.py check_species \
        -f input.fasta \
        -i /tmp/ \
        -s \
        -l /bactinspector/data/all_complete_refseq.k21s1000 \
        -d 0.075 \
        | python3 /bactinspector/map_to_pw.py