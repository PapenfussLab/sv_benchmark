#!/bin/bash
#
# runs whhamg against bams
#
. common.sh
CALLER=whamg/1.8.0
WHAM_ROOT=$BASE_DIR/tools/$CALLER
export PATH=$WHAM_ROOT/bin:$PATH
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	CX_BAM2CFG_FLAGS=-m
	CX_CALLER_FLAGS=""
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	module add openmpi perl
	ln -s $BAM input.bam
	ln -s $BAM.bai input.bam.bai
	whamg -a $CX_REFERENCE -f input.bam | perl $WHAM_ROOT/utils/filtWhamG.pl > wham.vcf && \
	cp wham.vcf $XC_OUTPUT
	"
	xc_exec
	exit
done

