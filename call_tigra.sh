#!/bin/bash
#
# runs TIGRA off breakdancer, pindel, and delly data
#
. common.sh

CALLER=tigra/0.3.7
# bwa mem already on path
export PATH=/usr/local/bioinfsoftware/samtools/samtools-0.1.6:$PATH # 0.3.7
export LD_LIBRARY_PATH=/usr/local/bioinfsoftware/samtools/samtools-0.1.6:$LD_LIBRARY_PATH # 0.3.7
#export PATH=/usr/local/bioinfsoftware/samtools/samtools-0.1.19/bin:$PATH # 0.4.0
export PATH=$BASE_DIR/tools/$CALLER:$PATH
export PATH=$BASE_DIR/tools/tigra/tigra-ext:$PATH

CHUNK_SIZE=1000
# for UPSTREAM_METADATA in $(grep -E "(breakdancer)|(pindel)|(delly)" $DATA_DIR/*.metadata) ; do
for UPSTREAM_METADATA in $(grep -L tigra $(grep -lE "(breakdancer)|(pindel)|(delly)" $DATA_DIR/*.metadata)) ; do
	cx_load $UPSTREAM_METADATA
	CX_UPSTREAM_VCF=$CX.vcf
	CX_UPSTREAM=$CX
	if [[ ! -f $CX_UPSTREAM_VCF ]] ; then
		echo Missing $CX_CALLER output $CX.vcf
		continue
	fi
	CX_CALLER=$CALLER/$CX_CALLER
	cx_save
	XC_OUTPUT=$CX.vcf
	if [[ "$CX_CALLER" =~ "breakdancer" ]] ; then
		ARGS="-b "
		CX_UPSTREAM_FILE="$CX_UPSTREAM/ouput.ctx"
	else
		ARGS=""
		CX_UPSTREAM_FILE="$CX_UPSTREAM_VCF"
	fi
	# parallel version
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
		mkdir -p $CX/split
		grep -E \"^##\" $CX_UPSTREAM_FILE > $CX/header.vcf
		split -d -a 8 -l $CHUNK_SIZE $CX_UPSTREAM_FILE $CX/split/input
		for file in $CX/split/input* ; do cat $CX/header.vcf \$file > \$file.vcf ; done
		find $CX/split -name 'input*.vcf' | parallel TIGRA-ext.pl \
			-d $CX \
			-r $CX_REFERENCE \
			-o $CX/sv{#}.vcf \
			-I /usr/local/bioinfsoftware/BEDTools/BEDTools-Version-2.11.2/bin/intersectBed \
			$ARGS \
			{} \
			$CX_BAM && \
		cat $BASE_DIR/tigra.vcf $CX/sv*.vcf > $CX.vcf
	"
	# serial version
	XC_SCRIPT="#rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
		if [[ ! -s b.bed ]] ; then
			TIGRA-ext.pl \
				-d $CX \
				-r $CX_REFERENCE \
				-o $CX/sv.vcf \
				-I /usr/local/bioinfsoftware/BEDTools/BEDTools-Version-2.11.2/bin/intersectBed \
				$ARGS \
				$CX_UPSTREAM_FILE \
				$CX_BAM
		fi
		if [[ ! -s $CX/intersect.bed ]] ; then
			if [[ ! -f $CX/a.bed.original ]] ; then
				mv $CX/a.bed $CX/a.bed.original
			fi
			if [[ ! -f $CX/b.bed.original ]] ; then
				mv $CX/b.bed $CX/b.bed.original
			fi
			awk '((\$2 <= \$3) && (\$2 > 0))' $CX/a.bed.original > $CX/a.bed
			awk '((\$2 <= \$3) && (\$2 > 0))' $CX/b.bed.original > $CX/b.bed
			intersectBed -wo -a $CX/a.bed -b $CX/b.bed > $CX/intersect.bed
			TIGRA-ext_no_tigra_bwa.pl \
			-d $CX \
			-r $CX_REFERENCE \
			-o $CX/sv.vcf \
			-I /usr/local/bioinfsoftware/BEDTools/BEDTools-Version-2.11.2/bin/intersectBed \
			$ARGS \
			$CX_UPSTREAM_FILE \
			$CX_BAM
		fi
		cat $BASE_DIR/tigra.vcf $CX/sv.vcf > $CX.vcf
	"
	xc_exec
	exit
done

