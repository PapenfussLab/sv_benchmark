#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
CALLER=breakdancer/1.4.5
# Need Statistics::Descriptive Math::CDF perl packages
export PATH=$BASE_DIR/tools/$CALLER/bin:$PATH
export PATH=$BASE_DIR/tools/$CALLER/perl:$PATH
# Need old version of GD
export PERL5LIB=$BASE_DIR/tools/breakdancer/perllib:$PERL5LIB
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	CX_BAM2CFG_FLAGS=-m
	CX_CALLER_FLAGS=""
	cx_save
	XC_OUTPUT=$CX.vcf
	#echo "readgroup:1	platform:illumina	map:$CX_BAM	readlen:$CX_READ_LENGTH	mean:$CX_READ_FRAGMENT_LENGTH	std:$CX_READ_FRAGMENT_STDDEV	exe:samtools view" > $CX.cfg
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	module add samtools perl
	ln -s $BAM input.bam
	bam2cfg.pl $CX_BAM2CFG_FLAGS $BAM > input.cfg
	breakdancer-max $CX_CALLER_FLAGS input.cfg > ouput.ctx
	$BASE_DIR/breakdancer2vcf.py \$(cut -f 10 < $CX/input.cfg | cut -b 5- ) < ouput.ctx > $XC_OUTPUT
	"
	xc_exec
done

