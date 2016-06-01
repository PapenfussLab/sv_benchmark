#!/bin/bash
#
# runs manta
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
	# TODO: should we merge diploidSV.vcf and candidateSV.vcf ?
	# This either stripping the genotyping from diploid, or
	# putting fake genotyping into candidate.
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	ln -s $BAM input.bam
	ln -s $BAM.bai input.bam.bai
	configManta.py \
		--bam input.bam \
		--referenceFasta $CX_REFERENCE \
		--runDir $CX && \
	$CX/runWorkflow.py -m local -j \$(nproc) && \
	$BASE_DIR/manta_combine.sh $CX/results > $CX.vcf
	"
	xc_exec
done

