#!/bin/bash

# Read the list of subjects from subjects.txt and process them line by line
while IFS= read -r subject; do
    echo "Running script for subject: $subject"
    #python rsa_plots_crossnobis_pseudorun.py "$subject"
    python rsa_modelcomparisonV2.py
done < subjects2.txt

