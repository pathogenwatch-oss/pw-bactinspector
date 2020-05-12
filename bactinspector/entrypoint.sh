#!/bin/bash

set -e -o pipefail

export PYTHONPATH=/:/bactinspector/

echoerr() { printf "%s\n" "$*" >&2; }

cat - > /tmp/input.fasta

python3 speciator.py /tmp/input.fasta /bactinspector/libraries
