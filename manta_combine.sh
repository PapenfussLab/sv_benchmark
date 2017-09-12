#!/bin/bash
# combines the three VCFs output by Manta
#
# $1 manta results directory
DIR="$1/variants"
if [[ ! -d "$DIR" ]] ; then
	echo "Could not find $DIR"
	exit 1
fi

gunzip - < $DIR/diploidSV.vcf.gz > $DIR/diploidSV.vcf && \
gunzip - < $DIR/candidateSV.vcf.gz > $DIR/candidateSV.vcf && \
gunzip - < $DIR/candidateSmallIndels.vcf.gz > $DIR/candidateSmallIndels.vcf && \
grep -E  '^##' $DIR/diploidSV.vcf | grep -vE '^((##FILTER)|(##INFO)|(##FORMAT)|(##ALT))' && \
echo '##FILTER=<ID=candidateSV,Description="Variant output to candidateSV call set">' && \
echo '##FILTER=<ID=candidateSmallIndels,Description="Variant output to candidateSmallIndels call set">' && \
grep -hE '^((##FILTER)|(##INFO)|(##FORMAT)|(##ALT))' $DIR/diploidSV.vcf $DIR/candidateSV.vcf $DIR/candidateSmallIndels.vcf | sort | uniq && \
grep -vE '^##' $DIR/diploidSV.vcf && \
grep -vE '^#' $DIR/candidateSV.vcf | awk 'BEGIN { FS="\t"; OFS="\t" } { $7="candidateSV" ; $10="" ; print }' && \
grep -vE '^#' $DIR/candidateSmallIndels.vcf  | awk 'BEGIN { FS="\t"; OFS="\t" } { $7="candidateSmallIndels" ; $10="" ; print }'