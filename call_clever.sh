#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
CALLER=clever/270dee03
PATH=$BASE_DIR/tools/clever-toolkit/scripts:$PATH
for BAM in $DATA_DIR/*.sq.bam ; do
	cx_load $BAM
	if [ "$CX_READ_FRAGMENT_LENGTH" == "" ] ; then
		echo "$CX missing read fragment length information: skipping"
		continue
	fi
	if [ $(($CX_READ_FRAGMENT_LENGTH - (2 * $CX_READ_LENGTH) )) -le 5 ] ; then
		# % postprocess-predictions  --vcf --stddev 50.1498 7d1f1aa985583461bcb1d006eb3700bf/predictions.raw.txt -8.13184 > predictions.vcf
		# % postprocess-predictions: error: no such option: -8
		# postprocess-predictions [options] <predictions(.gz)> <insert-size-mean>
		# crashes when a negative insert-size-mean is passed as it tries to parse
		# to parameter as an option since it has a leading minus sign
		echo "clever requires a positive average of unread sequence within fragments. Skipping Read Length ${CX_READ_LENGTH}bp for ${CX_READ_FRAGMENT_LENGTH}bp fragments."
		continue
	fi
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	# TODO: use --use_xa --sorted flags if SAM XA field is written (ie: BWA)
	cx_save
	XC_OUTPUT=$CX.vcf
	# Clever does not write the VCF SV headers - we need to add them ourselves
	XC_SCRIPT="module add bwa ; rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	ln -s $CX_BAM input.bam
	clever input.bam $CX_REFERENCE results && \
	cp indel.vcf $CX.vcf && \
	tail -n +2 $CX/predictions.vcf >> $CX.vcf
	# need to add missing VCF headers
	"
	xc_exec
	#mkdir -p results/work
	#cd results/work
	#add-score-tags-to-bam -s $CX_REFERENCE < ../../raw.bam > input.bam && \
	#insert-length-histogram -m lengths.mean-and-sd < input.bam > lengths.histogram && \
	#bam-to-alignment-priors lengths.histogram < input.bam | split-priors-by-chromosome -z priors
done

