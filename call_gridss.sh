#!/bin/bash
#
#
. common.sh
VERSION=0.11.5
CALLER=gridss/$VERSION
#FULL_JAR=~/bin/${CALLER/\//-}-jar-with-dependencies.jar
FULL_JAR=~/bin/gridss-$VERSION-jar-with-dependencies.jar

for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	if [ "$CX_ALIGNER_SOFTCLIP" == 0 ] ; then
		echo "GRIDSS: skipping end-to-end aligned $BAM"
		continue
	fi
	BLACKLIST=$CX_BLACKLIST
	if [[ "$BLACKLIST" == "" ]] ; then
		BLACKLIST="none"
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
						FSARGS=""
						if [[ "$CX_READ_FRAGMENT_LENGTH" != "" ]] ; then
							FSARGS="INPUT_MIN_FRAGMENT_SIZE=$(( CX_READ_FRAGMENT_LENGTH - 3 * CX_READ_FRAGMENT_STDDEV ))"
							FSARGS="$FSARGS INPUT_MAX_FRAGMENT_SIZE=$(( CX_READ_FRAGMENT_LENGTH + 3 * CX_READ_FRAGMENT_STDDEV ))"
						fi
						XC_SCRIPT="
rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX;
cat > $CX/gridss_custom.properties << EOF
visualisation.assemblyProgress = true
visualisation.buffers = true
visualisation.bufferTrackingItervalInSeconds = 10
assembly.method=$CX_ASSEMBLY_METHOD
assembly.k=$CX_K
scoring.model=$CX_MODEL
variantcalling.minSize=$MIN_EVENT_SIZE
realignment.aligner=$CX_REALIGNER
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
	-ea \\
	-Xmx32g \\
	-cp $FULL_JAR \\
	au.edu.wehi.idsv.Idsv \\
	TMP_DIR=$CX \\
	WORKING_DIR=$CX \\
	INPUT=$BAM \\
	$FSARGS \\
	OUTPUT=$CX/breakend.vcf \\
	REFERENCE=$CX_REFERENCE \\
	CONFIGURATION_FILE=$CX/gridss_custom.properties \\
	PER_CHR=$GRIDSS_PER_CHR \\
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

