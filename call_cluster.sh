#!/bin/bash
# Invokes all indel callers
#
PPN=4
PPN_MULTI=20
MEMGB=64
echo ./call_breakdancer.sh $1 -l nodes=1:ppn=${PPN},mem=${MEMGB}gb
#./call_clever.sh $1  -l nodes=1:ppn=${PPN},mem=${MEMGB}gb
echo ./call_cortex.sh $1  -l nodes=1:ppn=${PPN},mem=750gb
echo ./call_crest.sh $1  -l nodes=1:ppn=${PPN},mem=${MEMGB}gb
#./call_defuse.sh $1  -l nodes=1:ppn=${PPN},mem=${MEMGB}gb
echo ./call_delly.sh $1  -l nodes=1:ppn=${PPN},mem=${MEMGB}gb
#./call_gasv.sh $1  -l nodes=1:ppn=${PPN},mem=${MEMGB}gb
echo ./call_gridss.sh $1  -l nodes=1:ppn=${PPN_MULTI},mem=${MEMGB}gb
echo ./call_hydra.sh $1  -l nodes=1:ppn=${PPN_MULTI},mem=${MEMGB}gb
echo ./call_lumpy.sh $1  -l nodes=1:ppn=${PPN},mem=${MEMGB}gb
echo ./call_manta.sh $1  -l nodes=1:ppn=${PPN_MULTI},mem=${MEMGB}gb
#./call_meerkat.sh $1  -l nodes=1:ppn=${PPN},mem=${MEMGB}gb
echo ./call_pindel.sh $1  -l nodes=1:ppn=${PPN_MULTI},mem=${MEMGB}gb
#./call_samtools.sh $1  -l nodes=1:ppn=${PPN},mem=${MEMGB}gb
echo ./call_socrates.sh $1  -l nodes=1:ppn=${PPN},mem=${MEMGB}gb
#./call_tigra.sh $1  -l nodes=1:ppn=${PPN},mem=${MEMGB}gb
#./call_variationhunter.sh $1  -l nodes=1:ppn=${PPN},mem=${MEMGB}gb
#./call_whamg.sh $1  -l nodes=1:ppn=${PPN},mem=${MEMGB}gb
echo ./call_samtools.sh $1  -l nodes=1:ppn=${PPN},mem=${MEMGB}gb
