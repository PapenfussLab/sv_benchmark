# Structural Variant Benchmarking

This project contains a collection of scripts used for SV benchmarking and visualisation of the results. The project structure is as follows:

## Scripts for creating benchmarking VCFs

* input.*: reference data for a given sample
* data.*: directory containing the sequencing data and variant calls for a given sample
* data.*/*.metadata: generated file continaing bash variables containing the metadata about the file (e.g. library fragment size, aligner, variant caller, ...)
* alignbam.sh: script for running aligners (fastq to BAM)
* common.sh: helper utility script. Contains generic functions such as loading/saving metadata and executing jobs.
* call_*.sh: bash shell scripts for running the relevant variant caller (fastq/BAM to VCF)
* bamtofastq.sh: conversion script so de novo assembly based variant callers can be run on samples with input BAM files
* *2vcf.py: helper scripts to convert the results of a variant caller using a custom format to VCF
* clean*.sh: helper scripts to remove incomplete data
* setting.*: settings file for a given sample.
* xcall_*.sh: variant caller script that aren't working/don't work with the latest version of the variant caller

## Scripts for processing benchmarking VCFs

* R/global.R: contains constants and common data
* R/contig.R: contains instance-dependent configurations data such as directory locations
* R/install.R: helper script to install require packages in a new R environment
* R/lib*.R: collections of functions used
* R/manuscript_figures.R: script for generating non-simulation benchmarking paper figures
* R/shiny.Rproj: shiny app.
* R/server.R: shiny app
* R/ui.R: shiny app
* R/shinyCache2.R: caching infrastructure. We run out of memory if we load 5000+ VCFs at the same time so we heavily cache intermediate results
* R/precache.R: script for pre-caching shiny app results. Usability is poor when it takes 4 hours respond to a change in a drop-down so by precomputing all possible combinations of input elements on the cluster, we reduce both deployed app CPU usage, and disk usage (since we only need to deploy the final cached results, not all raw and intermediate results as well). Use --plot to generate simulation benchmarking paper figures
