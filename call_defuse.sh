#!/bin/bash
#
#
. common.sh
CALLER=defuse/0.8.0
DEFUSE_DIR=$BASE_DIR/tools/defuse/0.8.0
DEFUSE_DATADIR=$BASE_DIR/tools/defuse/dataset_directory
export PATH=$DEFUSE_DIR/scripts:$PATH

# Reference preprocessing
# defuse_create_ref.pl -d $DEFUSE_DATADIR -c $DEFUSE_DIR/scripts/config.txt

for FQ1 in $DATA_DIR/*.1.fq ; do
	cx_load $FQ1
	CX_CALLER=$CALLER
	CX_FQ1=$FQ1
	CX_FQ2=${FQ1/.1./.2.}
	cx_save
	XC_OUTPUT=$CX/results.tsv
	XC_SCRIPT="
#rm -rf $CX;
mkdir -p $CX 2>/dev/null; cd $CX;
module add samtools bowtie ucsc-tools R gmap-gsnap 

defuse_run.pl \
	-c $DEFUSE_DIR/scripts/config.txt \
	-d $DEFUSE_DATADIR \
	-1 $CX_FQ1 \
	-2 $CX_FQ2 \
	-o $CX \
	-p $MAX_CORES

"
	xc_exec
done

