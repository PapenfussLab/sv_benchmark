#!/bin/bash
#
# quername sorts files already sorted by coordinate
. common.sh

for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	# So we don't overwrite the output from the alignment itself
	XC_SUFFIX=.sc2sq
	XC_OUTPUT=$CX.sq.bam
	XC_MULTICORE=1
	XC_SCRIPT="
	mkdir $CX 2>/dev/null
	module add novoalign
	novosort -n -3 -t $CX -c \$(nproc) -o $XC_OUTPUT.tmp.bam $BAM && \
	mv $XC_OUTPUT.tmp.bam $XC_OUTPUT
	"
	xc_exec
done
