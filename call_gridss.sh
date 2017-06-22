#!/bin/bash
#
#
. common.sh
VERSION=1.4.1
CALLER=gridss/$VERSION
FULL_JAR=~/bin/gridss-$VERSION-jar-with-dependencies.jar

for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	if [[ "$CX_ALIGNER_SOFTCLIP" == 0 ]] ; then
		echo "GRIDSS: skipping end-to-end aligned $BAM"
		continue
	fi
	BLACKLIST=$CX_BLACKLIST
	if [[ "$BLACKLIST" == "" ]] ; then
		BLACKLIST="null"
	fi
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	for CX_ASSEMBLY_METHOD in $GRIDSS_METHODS ; do
		for CX_K in $GRIDSS_KMER ; do
			for CX_MODEL in $GRIDSS_MODEL ; do
				for CX_EXCLUSION in $GRIDSS_EXCLUSION ; do
					for CX_REALIGNER in $GRIDSS_ALIGNER ; do
						CX_CALLER_ARGS="$CX_K,$CX_ASSEMBLY_METHOD,$CX_MODEL,$CX_EXCLUSION,$CX_REALIGNER"
						cx_save
						XC_OUTPUT=$CX.vcf
						XC_MULTICORE=1
						FSARGS=""
						if [[ "$CX_READ_FRAGMENT_LENGTH" != "" ]] ; then
							FSARGS="INPUT_MIN_FRAGMENT_SIZE=$(( CX_READ_FRAGMENT_LENGTH - 3 * CX_READ_FRAGMENT_STDDEV ))"
							FSARGS="$FSARGS INPUT_MAX_FRAGMENT_SIZE=$(( CX_READ_FRAGMENT_LENGTH + 3 * CX_READ_FRAGMENT_STDDEV ))"
						fi
						# TODO: fix realignment for GRIDSS 1.0 command line
						# realignment.aligner=$CX_REALIGNER
						XC_SCRIPT="
#rm -rf $CX;
mkdir $CX 2>/dev/null; cd $CX;
# locking already taken care of
rm -r $CX/gridss.lock.breakend.vcf
cat > $CX/gridss_custom.properties << EOF
#visualisation.assemblyProgress = true
#visualisation.buffers = true
#visualisation.bufferTrackingItervalInSeconds = 10
assembly.method=$CX_ASSEMBLY_METHOD
assembly.k=$CX_K
scoring.model=$CX_MODEL
variantcalling.minSize=$MIN_EVENT_SIZE
EOF
if [[ \"$CX_EXCLUSION\" == \"SC\" ]] ; then
	echo scoring.exclude=SplitRead >> $CX/gridss_custom.properties
	echo scoring.exclude=Indel >> $CX/gridss_custom.properties
	echo scoring.exclude=SoftClip >> $CX/gridss_custom.properties
fi
if [[ \"$CX_EXCLUSION\" == \"RP\" ]] ; then
	echo scoring.exclude=UnmappedMate >> $CX/gridss_custom.properties
	echo scoring.exclude=DiscordantPair >> $CX/gridss_custom.properties
	echo scoring.exclude=SoftClip >> $CX/gridss_custom.properties
fi
if [[ \"$CX_MODEL\" == \"ReadCount\" ]] ; then
	echo variantcalling.minScore=2 >> $CX/gridss_custom.properties
	echo variantcalling.lowQuality=13 >> $CX/gridss_custom.properties
fi
java \\
	-Dsamjdk.buffer_size=524288 \\
	-Dsamjdk.create_index=true \\
	-Dsamjdk.use_async_io_read_samtools=true \\
	-Dsamjdk.use_async_io_write_samtools=true \\
	-Dsamjdk.use_async_io_write_tribble=true \\
	-Dsamjdk.compression_level=1 \\
	-ea \\
	-Xmx31g \\
	-cp $FULL_JAR \\
	gridss.CallVariants \\
	TMP_DIR=$CX \\
	WORKING_DIR=$CX \\
	INPUT=$BAM \\
	$FSARGS \\
	ASSEMBLY=$CX/assembly.bam \\
	OUTPUT=$CX/breakend.vcf \\
	REFERENCE_SEQUENCE=$CX_REFERENCE \\
	CONFIGURATION_FILE=$CX/gridss_custom.properties \\
	VERBOSITY=DEBUG \\
	BLACKLIST=$BLACKLIST \\
	&& \\
cp $CX/breakend.vcf $CX.vcf
		"
						xc_exec
					done
				done
			done
		done
	done
done

