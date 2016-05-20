#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
CALLER=breakpointer/1.0

# https://www.broadinstitute.org/cancer/cga/dranger
# dRanger is not yet available to external users.

for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	# create split reference in required format
	if [[ ! -f  $CX_REFERENCE.split/$(head -1 $CX_REFERENCE.fai | cut -f 1).txt ]] ; then
		if [[ ! -f  $CX_REFERENCE.split/$(head -1 $CX_REFERENCE.fai | cut -f 1).fa ]] ; then
			echo "Missing reference split by chr in $CX_REFERENCE.split"
			exit
		fi
		for FA in $CX_REFERENCE.split/*.fa ; do
			tail -n +2 $FA | tr -d '\r\n' > $(dirname $FA)/$(basename $FA .fa).txt
		done
	fi
	exit
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	CX_BAM2CFG_FLAGS=-m
	CX_CALLER_FLAGS=""
	cx_save
	XC_OUTPUT=$CX.vcf
	echo "readgroup:1	platform:illumina	map:$CX_BAM	readlen:$CX_READ_LENGTH	mean:$CX_READ_FRAGMENT_LENGTH	std:$CX_READ_FRAGMENT_STDDEV	exe:samtools view" > $CX.cfg
	XC_SCRIPT="module add samtools $CALLER ; rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	ln -s $BAM input.bam
	bam2cfg.pl $CX_BAM2CFG_FLAGS $BAM > input.cfg
	breakdancer-max $CX_CALLER_FLAGS input.cfg > ouput.ctx
	$BASE_DIR/breakdancer2vcf.py $CX_READ_FRAGMENT_LENGTH < ouput.ctx > $XC_OUTPUT
	"
	xc_exec
done

