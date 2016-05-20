#!/bin/bash
#
#
. common.sh
CALLER=soapsv/1.02
PATH=$BASE_DIR/tools/$CALLER/final:$PATH
PATH=/usr/local/bioinfsoftware/SOAPdenovo/SOAPdenovo2-bin-LINUX-generic-r240:$PATH
for FQ1 in $DATA_DIR/*.1.fq ; do
	cx_load $FQ1
	CX_FQ1=$FQ1
	CX_FQ2=${FQ1/.1.fq/.2.fq}
	CX_CALLER=$CALLER
	
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="
rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX;
mkdir ref

	"
	xc_exec
done

