#!/bin/bash

#Script to generate permutations for all subjects (searchlight):

set -e #stop script if an error occurs.

python3 wholebrain_rawdist_perm.py 001
python3 wholebrain_rawdist_perm.py 004
python3 wholebrain_rawdist_perm.py 005
python3 wholebrain_rawdist_perm.py 006
python3 wholebrain_rawdist_perm.py 007
python3 wholebrain_rawdist_perm.py 008
python3 wholebrain_rawdist_perm.py 009
python3 wholebrain_rawdist_perm.py 010
python3 wholebrain_rawdist_perm.py 012
python3 wholebrain_rawdist_perm.py 013

echo "All scripts have been executed."
