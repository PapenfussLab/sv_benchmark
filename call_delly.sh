#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
#CALLER=delly/0.6.8
CALLER=delly/0.7.6
export PATH=$BASE_DIR/tools/$CALLER/bin:$PATH
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	cx_save
	XC_OUTPUT=$CX.vcf
	if [[ $CALLER =~ "0.6.8" ]] ; then
		XC_SCRIPT="module add picard-tools; rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
		delly -t DEL -o del.vcf -g $CX_REFERENCE $CX_BAM && \
		delly -t DUP -o dup.vcf -g $CX_REFERENCE $CX_BAM && \
		delly -t INV -o inv.vcf -g $CX_REFERENCE $CX_BAM && \
		delly -t TRA -o tra.vcf -g $CX_REFERENCE $CX_BAM && \
		SortVcf I=del.vcf I=dup.vcf I=inv.vcf I=tra.vcf O=merged.vcf SEQUENCE_DICTIONARY=${CX_REFERENCE/.fa/.dict} && \
		cp merged.vcf $XC_OUTPUT
		"
	else
		XC_SCRIPT="module add picard-tools bcftools; rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
		delly call -t DEL -o del.bcf -g $CX_REFERENCE $CX_BAM && \
		delly call -t DUP -o dup.bcf -g $CX_REFERENCE $CX_BAM && \
		delly call -t INV -o inv.bcf -g $CX_REFERENCE $CX_BAM && \
		delly call -t TRA -o tra.bcf -g $CX_REFERENCE $CX_BAM && \
		bcftools merge --force-samples del.bcf dup.bcf inv.bcf tra.bcf > merged.vcf && \
		cp merged.vcf $XC_OUTPUT
		"
	fi
	xc_exec
done

