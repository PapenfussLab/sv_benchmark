#!/bin/bash
#
# converts a BAM back to FASTQ
# This approach allows downsampling to use the same set of reads as BWA alignment
# so all callers get the same set of input reads
#
. common.sh

for BAM in $DATA_DIR/*.sq.bam ; do
	echo Processing $BAM
	cx_load $BAM
	if [[ "$CX_ALIGNER" != "$2" ]] ; then
		echo -n #continue
	fi
	# unset context added by aligner in alignbam.sh
	unset CX_FQ1
	unset CX_FQ2
	unset CX_ALIGNER
	unset CX_ALIGNER_MODE
	unset CX_ALIGNER_SOFTCLIP
	unset CX_ALIGNER_VERSION
	cx_save
	XC_OUTPUT=$CX.1.fq.gz
	XC_SCRIPT="
	# FU=$CX.u.fq.tmp # shouldn't have any unpaired, and if so, we can't handle downstream
	#SamToFastq VALIDATION_STRINGENCY=SILENT INPUT=$BAM F=$CX.1.fq.tmp F2=$CX.2.fq.tmp
	bedtools bamtofastq -i $BAM -fq $CX.1.fq.tmp -fq2 $CX.2.fq.tmp && 
	mv $CX.1.fq.tmp $CX.1.fq &&
	mv $CX.2.fq.tmp $CX.2.fq
	"
	xc_exec
done
