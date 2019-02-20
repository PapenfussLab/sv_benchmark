#!/bin/bash
#
# Aligns paired fastq reads
#
. common.sh

export PATH=$BASE_DIR/tools/variationhunter/mrfast-2.6.0.1:$PATH

function align_subread {
	if [ "$3" != "local" ] ; then
		return
	fi
	cx_load $1
	CX_FQ1=$1
	CX_FQ2=$2
	CX_ALIGNER=subread
	CX_ALIGNER_MODE=$3
	cx_save
	INDEX=$CX_REFERENCE.subread
	if [ ! -f $INDEX.00.b.array ] ; then
		echo "Unable to find subread index for $CX_REFERENCE"
		echo subread-buildindex -o $CX_REFERENCE.subread $CX_REFERENCE 
		exit 1
	fi
	if [[ "$1" == *.gz ]] ; then
		echo "subread does not handle .gz compressed input, skipping $1"
		return
	fi 
	XC_MULTICORE=1
	XC_OUTPUT=$CX.su.bam
	XC_SCRIPT="module add subread
	subread-align -T \$(nproc) -i $INDEX -r $1 -R $2 |
		AddOrReplaceReadGroups I=/dev/stdin O=$CX.tmp.bam $READ_GROUP_PARAMS &&
	mv $CX.tmp.bam $XC_OUTPUT
	"
	xc_exec
}
function align_bwa {
	if [[ "$3" =~ "local" ]] ; then
		cx_load $1
	else 
		return
	fi
	CX_FQ1=$1
	CX_FQ2=$2
	CX_ALIGNER=bwamem
	CX_ALIGNER_MODE=$3
	ARGS=""
	if [ "$3" == "local,multimapping" ] ; then
		ARGS="-a"
	fi
	cx_save
	INDEX=${CX_REFERENCE}.bwa
	if [ ! -f "$INDEX.bwt" ] ; then
		INDEX=${CX_REFERENCE}
		if [ ! -f "$INDEX.bwt" ] ; then
			echo "Unable to find bwa index for $CX_REFERENCE"
			echo bwa index $CX_REFERENCE
			exit 1
		fi
	fi
	XC_MULTICORE=1
	XC_OUTPUT=$CX.su.bam
	XC_SCRIPT="
	bwa mem -t \$(nproc) $ARGS $INDEX $1 $2 > $CX.tmp.sam && \
		AddOrReplaceReadGroups I=$CX.tmp.sam O=$CX.tmp.bam $READ_GROUP_PARAMS &&
	mv $CX.tmp.bam $XC_OUTPUT &&
	rm $CX.tmp.sam
	"
	xc_exec
}
function align_bowtie2 {
	cx_load $1
	CX_FQ1=$1
	CX_FQ2=$2
	CX_ALIGNER=bowtie2
	CX_ALIGNER_MODE=$3
	if [[ "$CX_ALIGNER_MODE" =~ "global" ]] ; then
		ALIGNER_FLAGS="--end-to-end"
	elif [[ "$CX_ALIGNER_MODE" =~ "local" ]] ; then
		ALIGNER_FLAGS="--local"
	else
		return
	fi
	if [[ "$CX_ALIGNER_MODE" =~ "multimapping" ]] ; then
		ALIGNER_FLAGS="$ALIGNER_FLAGS -k $MULTIMAPPING_LOCATIONS"
		CX_MULTIMAPPING_LOCATIONS=$MULTIMAPPING_LOCATIONS
	fi
	cx_save
	#INDEX=$(ls -1 ${CX_REFERENCE/.fa/}*.rev.1.bt2)
	#INDEX=${INDEX/.rev.1.bt2/}
	INDEX=$CX_REFERENCE
	if [[ ! -f $INDEX.rev.1.bt2 ]] ; then
		echo "Unable to find bowtie2 index for $CX_REFERENCE"
		echo bowtie2-build $CX_REFERENCE $CX_REFERENCE
		exit 1
	fi
	XC_MULTICORE=1
	XC_OUTPUT=$CX.su.bam
	XC_SCRIPT="
	module add bowtie2
	bowtie2 --threads \$(nproc) --mm $ALIGNER_FLAGS -x $INDEX -1 $1 -2 $2 |
		AddOrReplaceReadGroups I=/dev/stdin O=$CX.tmp.bam $READ_GROUP_PARAMS &&
	mv $CX.tmp.bam $XC_OUTPUT
	"
	xc_exec
}
function align_novoalign {
	cx_load $1
	CX_FQ1=$1
	CX_FQ2=$2
	CX_ALIGNER=novoalign
	CX_ALIGNER_MODE=$3
	if [[ "$CX_ALIGNER_MODE" =~ "global" ]] ; then
		ALIGNER_FLAGS="-o FullNW"
	elif [[ "$CX_ALIGNER_MODE" =~ "local" ]] ; then
		ALIGNER_FLAGS="-o Softclip"
	else
		return
	fi
	if [[ "$CX_ALIGNER_MODE" =~ "multimapping" ]] ; then
		ALIGNER_FLAGS="$ALIGNER_FLAGS -r All -e $MULTIMAPPING_LOCATIONS"
		CX_MULTIMAPPING_LOCATIONS=$MULTIMAPPING_LOCATIONS
	fi
	cx_save
	INDEX=$(ls -1 ${CX_REFERENCE/.fa/}*.novoindex)
	if [ ! -f "$INDEX" ] ; then
		echo "Unable to find novoalign index for $CX_REFERENCE"
		novoindex $CX_REFERENCE.novoindex $CX_REFERENCE 
		exit 1
	fi
	if [ $(( $(stat -c%s $INDEX) / 1000000 )) -gt $(( $MAX_MEMORY - 1024 )) ] ; then
		echo "$(basename $CX): Not scheduling novoalign as index $INDEX requires more memory than is available: $(($(stat -c%s $INDEX) / 1000000))Mb required, $(( $MAX_MEMORY  - 1024 ))Mb free"
		return
	fi
	XC_MULTICORE=1
	XC_OUTPUT=$CX.su.bam
	XC_SCRIPT="
	novoalign -c \$(nproc) -o SAM $ALIGNER_FLAGS -F STDFQ -i PE $CX_READ_FRAGMENT_LENGTH,$CX_READ_FRAGMENT_STDDEV -d $INDEX -f $1 $2 |
		AddOrReplaceReadGroups I=/dev/stdin O=$CX.tmp.bam $READ_GROUP_PARAMS &&
	mv $CX.tmp.bam $XC_OUTPUT
	"
	xc_exec
}

function align_mrfast {
	cx_load $1
	INDEX=$(ls -1 $CX_REFERENCE.index)
	if [ ! -f "$INDEX" ] ; then
		echo "Fatal: Unable to find mrfast index for $CX_REFERENCE"
		echo mrfast --index $CX_REFERENCE
		exit 1
	fi
	CX_ALIGNER_MIN_CONCORD=$((CX_READ_FRAGMENT_LENGTH - 4 * CX_READ_FRAGMENT_STDDEV ))
	CX_ALIGNER_MAX_CONCORD=$((CX_READ_FRAGMENT_LENGTH + 4 * CX_READ_FRAGMENT_STDDEV ))
	CX_FQ1=$1
	CX_FQ2=$2
	CX_ALIGNER=mrfast
	CX_ALIGNER_MODE=global
	cx_save
	XC_OUTPUT=$CX.su.bam
	# 3 stddev = concordant
	MIN=$((CX_READ_FRAGMENT_LENGTH - 3 * CX_READ_FRAGMENT_STDDEV ))
	MAX=$((CX_READ_FRAGMENT_LENGTH + 3 * CX_READ_FRAGMENT_STDDEV ))
	# http://mrfast.sourceforge.net/manual.html
	# Paired End Mode: The number of reads in each file should not exceed 1 million (500,000 pairs),
	# however chunk size of 500,000 reads (250,000 pairs) is recommended.
	CHUNK_READS=250000
	CHUNK_LINES=$((4 * $CHUNK_READS))
	XC_MULTICORE=1
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	split -a 6 -l $CHUNK_LINES $CX_FQ1 split1
	split -a 6 -l $CHUNK_LINES $CX_FQ2 split2
	for CHUNK in split1* ; do
		CHUNK=\${CHUNK/split1/}
		FQ1=\$CHUNK.1.fq
		FQ2=\$CHUNK.2.fq
		mv split1\$CHUNK \$FQ1
		mv split2\$CHUNK \$FQ2
	done
	rm -f $CX/mrfast.sh
	for FQ1 in $CX/*.1.fq ; do
		CHUNK=\${FQ1/.1.fq/}
		FQ2=\${FQ1/.1./.2.}
		echo mrfast --search $CX_REFERENCE --pe --seq1 \$FQ1 --seq2 \$FQ2 --discordant-vh --min $MIN --max $MAX --maxdis $MULTIMAPPING_LOCATIONS -o \$CHUNK >> $CX/mrfast.sh
	done
	parallel -j \$(nproc) < $CX/mrfast.sh
	cat *_DIVET.vh > output_DIVET.vh && 
	mv output_DIVET.vh ${CX}_DIVET.vh
	"
	xc_exec
}

function align_razers3 {
	echo "TODO: implement aligner"
	# http://bioinformatics.oxfordjournals.org/content/early/2012/08/23/bioinformatics.bts505
}
function align_segemehl {
	echo "TODO: implement aligner"
}
function align_masai {
	echo "TODO: implement aligner"
}
function align_soap2 {
	echo "TODO: implement aligner"
}
function align_hisat {
	echo "TODO: implement aligner"
}
function align_star {
	echo "TODO: implement aligner"
}
function align_tophat2 {
	echo "TODO: implement aligner"
}
function align_bbmap {
	echo "TODO: implement aligner"
}


# Align reads
for FQ1 in $(ls -1 $DATA_DIR/*.1.fastq.gz $DATA_DIR/*.1.fq 2>/dev/null) ; do
	FQ2=${FQ1/.1./.2.}
	cx_load $FQ1
	if [ ! -z "$CX_DOWNSAMPLE_FROM" ] ; then
		# skip downsampled bams as they can be efficiently generated
		# by downsamplebam.sh
		continue
	fi
	#FQ1="<( bamToFastq -i $CX_BAM -fq /dev/stdout -fq2 /dev/null )"
	#FQ2="<( bamToFastq -i $CX_BAM -fq /dev/null -fq2 /dev/stdout )"
	RGSM=$FQ1
	if [[ "$CX_REFERENCE_VCF" != "" ]] ; then
		RGSM=$(basename $CX_REFERENCE_VCF)
	fi
	RGLB=$FQ1
	if [[ "$CX_REFERENCE_VCF" != "" ]] ; then
		RGLB="$(basename $CX_REFERENCE_VCF)_rl${CX_READ_LENGTH}_f${CX_READ_FRAGMENT_LENGTH}"
	fi
	READ_GROUP_PARAMS="RGPL=illumina RGPU=sim RGSM=$RGSM RGLB=$RGLB RGCN=sim${CX_READ_SIM}"
	if [[ "$2" == "" ]] ; then
		ALIGNERS_TO_PROCESS="$ALIGNERS"
	else
		ALIGNERS_TO_PROCESS="$2"
	fi
	# $2 now used by qsub
	ALIGNERS_TO_PROCESS="$ALIGNERS"
	for ALIGNER in $ALIGNERS_TO_PROCESS ; do
		echo "Processing $ALIGNER"
		for ALIGNER_MODE in $ALIGNER_MODES ; do
			echo "Processing $ALIGNER_MODE"
			align_$ALIGNER $FQ1 $FQ2 $ALIGNER_MODE
		done
	done
done

