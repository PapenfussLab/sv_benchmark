#!/bin/bash
#
# runs GASVPro against bams
#
. common.sh
CALLER=gasvpro/20140228
BAMTOOLS=$BASE_DIR/tools/bamtools/2.4.0/bin/bamtools
JAVA="java -Xmx31g"
for BAM in $DATA_DIR/*.sq.bam ; do
	cx_load $BAM
	if [ $(($CX_READ_FRAGMENT_LENGTH - (2 * $CX_READ_LENGTH + 50 ) )) -le 0 ] ; then
		echo "GASV requires at least 50bp unread within fragments. Skipping Read Length ${CX_READ_LENGTH}bp for ${CX_READ_FRAGMENT_LENGTH}bp fragments."
		continue
	fi
	CHR_PREFIX=""
	if [[ "$(head -c 4 $CX_REFERENCE)" == ">chr" ]] ; then
		CHR_PREFIX="chr"
	fi
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	cx_save
	# use (offset) reference index of chromosomes as custom ordering
	XC_OUTPUT=$CX.vcf
mkdir -p $CX
# TODO: CHROMOSOME_NAMING_FILE
cat > $CX/GASVPro.sh << EOF
#!/bin/sh
BAMFILEHQ=$CX/unique.bam  ##BAMFILEHQ
BAMFILELQ=$CX/multi.bam  ##BAMFILELQ
MCMCTMPDIR=$CX/mcmc ##MCMCTMPDIR #GIVE FULL PATH!
GASVDIR=$BASE_DIR/tools/$CALLER    ##GASVDIRECTORY
EOF
	tail -n +35 $BASE_DIR/tools/$CALLER/bin/GASVPro.sh >> $CX/GASVPro.sh
	chmod +x $CX/GASVPro.sh
	XC_SCRIPT="module add samtools; mkdir $CX $CX/mcmc 2>/dev/null; cd $CX
	
	if [[ ! -f $CX/unique.bam ]] ; then
		if [[ \"$CX_MULTIMAPPING_LOCATIONS\" != \"\" ]] ; then
			echo Using bamtools to extract uniquely mapping reads
			$BAMTOOLS filter -in $CX_BAM -tag 'IH:>1' -out $CX/multi.bam &
			$BAMTOOLS filter -in $CX_BAM -tag 'IH:<=1' -out $CX/unique.bam &
			wait
		else
			samtools view -q 10 -u -o $CX/unique.bam -U $CX/multi.bam $CX_BAM
		fi
	fi
	$CX/GASVPro.sh && \
	$BASE_DIR/gasv2vcf.py < BAMToGASV.gasvpro.in.ALL.MCMCThreshold.clusters.pruned.clusters $BASE_DIR/gasv.CHROMOSOME_NAMING_FILE.txt > $CX.vcf"
	xc_exec
done

