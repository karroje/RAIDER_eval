#!/bin/bash
############################################
# Nathan Figueroa - Miami University
# April 01 2013
#
# Simple script for automating all of 
# RepeatScout's procedures
#
# Takes one parameter: Sequence file name
#
############################################

SEQUENCE=$1

../build_lmer_table -sequence $SEQUENCE -freq freq_table
../RepeatScout -sequence $SEQUENCE -freq freq_table -output repscout.fa
cat repscout.fa | perl ../filter-stage-1.prl > filtered_repscout.fa

