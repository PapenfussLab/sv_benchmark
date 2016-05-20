#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
CALLER=manta/0.29.6
export PATH=$BASE_DIR/tools/$CALLER/bin:$PATH
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	ln -s $BAM input.bam
	ln -s $BAM.bai input.bam.bai
	configManta.py \
		--bam input.bam \
		--referenceFasta $CX_REFERENCE \
		--runDir $CX && \
	$CX/runWorkflow.py -m local -j \$(nproc)
	"
	xc_exec
done

