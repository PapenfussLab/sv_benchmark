#!/bin/bash
#
# runs socrates against bams
#
. common.sh
VERSION=1.13.1
CALLER=socrates/$VERSION
SOCRATES_JAR=$BASE_DIR/tools/socrates-$VERSION-jar-with-dependencies.jar
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	if [[ -f $BAM.bt2.bam ]] ; then
		echo "Using $BAM.bt2.bam as proxy"
		CX_BAM=$BAM.bt2.bam
	fi
	CX_CALLER=$CALLER
	CX_CALLER_ARGS=
	if [[ "$CX_ALIGNER_SOFTCLIP" == 0 ]] ; then
		echo "Socrates: skipping end-to-end aligned $BAM"
		continue
	fi
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	ln -s $CX_BAM $CX/input.bam
	ln -s $CX_BAM.bai $CX/input.bam.bai
	module remove bowtie2
	module add bowtie2/2.3.2
	java -Xmx30g -jar $SOCRATES_JAR -t \$(nproc) $CX_REFERENCE $CX/input.bam && \
	$BASE_DIR/socrates2vcf.py $CX/results_Socrates_paired_*.txt $CX/results_Socrates_unpaired_*.txt > $CX.vcf
	"
	xc_exec
done

