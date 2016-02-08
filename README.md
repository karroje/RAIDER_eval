RAIDER_eval
===========

This repository contains various developmental versions of the RAIDER tool for de novo sequence identification, as well as as automated tests for the analysis RAIDER results.

**RAIDER code**:
* ~~RAIDER_src: Source code for original RAIDER prototype.~~ 
* ~~RAIDER_prescan: Source cod for a low-memory version of RAIDER.~~
* ~~RAIDER_srv_v2: Adaptation of RAIDER to spaced seed searches.~~
* phRAIDER: Current invocation of phRAIDER.

Simulation code:
* chromosome_simulator.py: 
  * Generates simulated chromosomes based on an order-k Markov chain estimated from a template organism genome.
  * Leaves transposable element sequences intact, while generating inter-TE sequences from the Markov chain.
* markov_gen.py: Libraries for working with k-th order Markov chains:

Pipeline code:
* redhawk.py: Python library for interfacing with a pbs job management system
* ~~Raider_eval.py: Pipeline for automated comparison of de novo TE identification tools (deprecated)~~
* testing_pipeline.py: New pipeline for initiating large batches of TE identification tool runs in parallel on a PBS cluster
* collect_stats.py: Analysis scripts -- collect results from testing_pipeline.py output and create result tables.

Analysis code:
* analysis*.R: R code for generating plots and statistical analysis from the pipeline output


**Authors**:

* Nathanial Figueroa 
  * Co-author: RAIDER
* Xiaolin Liu
* Carly Schaeffer 
  * Co-author: phRAIDER
  * Contributed to developement of pipeline code.
* **Dr. John Karro: Primary investigator**
  * co-author: RAIDER, RAIDER_eval
  * Primary developer of simulation code, pipeline code, analysis code
