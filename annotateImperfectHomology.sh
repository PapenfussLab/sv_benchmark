#!/bin/bash
#
#
. common.sh
FULL_JAR=~/bin/gridss-1.4.2-SNAPSHOT-jar-with-dependencies.jar

for VCF in $DATA_DIR/*.vcf ; do
	cx_load $VCF
	XC_OUTPUT=$VCF
	XC_OUTPUT=$CX.annimphom.bedpe
	CX_INPUT_VCF=$VCF
	CX_IHOMANN=$XC_OUTPUT
	cx_save
	XC_MULTICORE=1
	XC_SCRIPT="
module add java R
rm -rf $CX;
mkdir $CX 2>/dev/null; cd $CX;
Rscript $BASE_DIR/tobedpe.R $VCF $CX/in.bedpe && \
java \\
	-ea \\
	-Xmx3g \\
	-cp $FULL_JAR \\
	gridss.AnnotateInexactHomologyBedpe \\
	INPUT=$CX/in.bedpe \\
	OUTPUT=$CX/out.bedpe \\
	REFERENCE_SEQUENCE=$CX_REFERENCE \\
	VERBOSITY=DEBUG \\
	&& \\
cp $CX/out.bedpe $XC_OUTPUT
		"
	xc_exec
done

