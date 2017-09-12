#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
PATH=$PATH

for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=bcftools/1.3.1
	CX_CALLER_OPTIONS=-m
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	module remove samtools bcftools
	module add samtools/1.5 bcftools/1.3.1
	samtools mpileup -R -u -f $CX_REFERENCE $BAM | bcftools call $CX_CALLER_OPTIONS -v -O v --threads \$(nproc) -o $CX/calls.vcf && \
	mv $CX/calls.vcf $XC_OUTPUT"
	xc_exec
done

