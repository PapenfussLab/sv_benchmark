#!/bin/bash
#
# runs variantionhunter
#
. common.sh

CALLER=variationhunter/0.04
VH_DIR=$BASE_DIR/tools/variationhunter/CommonLawRelease
HG19_DIR=$VH_DIR/Hg19_NecessaryFiles
export PATH=$VH_DIR/selection:$VH_DIR/clustering:$PATH

for DIVET in $DATA_DIR/*_DIVET.vh ; do
	cx_load $DIVET
	CX_DIVET=$DIVET
	CX_CALLER=$CALLER
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="mkdir -p $CX $CX/align $CX 2>/dev/null;
# run VariationHunter
cd $CX
echo 1 > $CX/sample.lib
echo lib sample $DIVET $(( CX_READ_FRAGMENT_LENGTH - 3 * CX_READ_FRAGMENT_STDDEV )) $(( CX_READ_FRAGMENT_LENGTH + 3 * CX_READ_FRAGMENT_STDDEV )) $CX_READ_LENGTH >> $CX/sample.lib

awk -f $BASE_DIR/variationhunter_divit_filter.awk $DIVET > $CX/filtered.DIVET.vh

VH \
	--chro $HG19_DIR/AllChro \
	--init $HG19_DIR/initInfo \
	--lib $CX/sample.lib \
	--repeat $HG19_DIR/Hg19.Satellite \
	--gap $HG19_DIR/hg19_Gap.Table.USCS.Clean \
	--output $CX/sample.cluster.out \
	--outputRead $CX/sample.name \
	--maxmapping 500 \
	--prunprob 0.001 && \
multiInd_SetCover \
	-l $CX/sample.lib \
	-r $CX/sample.name \
	-c $CX/sample.cluster.out \
	-t 50000 \
	-o $CX/out.sv && \
$BASE_DIR/variationhunter2vcf.py < out.sv > $XC_OUTPUT
	"
	xc_exec
done

