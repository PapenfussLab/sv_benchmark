#!/bin/bash
#
# runs CREST against bams
#
. common.sh
CALLER=crest/0.0.1
export PATH=$BASE_DIR/tools/blat/34x12/bin:$PATH 
export PATH=$BASE_DIR/tools/cap3/2010827:$PATH 
export PATH=$BASE_DIR/tools/$CALLER:$PATH 
export PERL5LIB=$BASE_DIR/tools/$CALLER
# TODO: missing BioPerl
# cpan
# install Bio::DB::Sam

BLAT_PORT=24513
BLAT_POLL_INTERVAL=10
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=crest
	CX_CALLER_SETTINGS=multithreaded
	cx_save
	if [ ! -f $CX_REFERENCE.2bit ] ; then
		echo "Missing $CX_REFERENCE.2bit - unable run CREST"
		exit 1
	fi
	BLAT_TMP_FILE=/tmp/crest.blat.$BLAT_PORT.$(basename $CX).lock
	XC_OUTPUT=$CX.vcf
	XC_TRAP="rm $BLAT_TMP_FILE ; if ls /tmp/crest.blat.$BLAT_PORT.*.lock 2>/dev/null ; then echo 'Additional instance using blat server' ; else gfServer stop localhost $BLAT_PORT 2>/dev/null; pkill -s 0 gfServer; echo 'Blat server terminated' ; fi "
	XC_SCRIPT="
	#rm -rf $CX;
	mkdir $CX $CX/tmp 2>/dev/null; cd $CX
	export TMP=$CX/tmp
	export TEMP=$CX/tmp
	export TMPDIR=$CX/tmp
	module add samtools ucsc-tools perl parallel
	ln -s $BAM $CX/in.bam
	ln -s $BAM.bai $CX/in.bam.bai
	echo > $BLAT_TMP_FILE
	gfServer start localhost $BLAT_PORT $CX_REFERENCE.2bit & 
	while
		echo Waiting for ${BLAT_POLL_INTERVAL}s for BLAT server to accept requests 1>&2
		sleep $BLAT_POLL_INTERVAL
		gfClient localhost $BLAT_PORT $CX_REFERENCE.2bit <( echo -e '>\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' ) /dev/null
		[[ \$? != 0 ]]
	do
		:
	done
	echo Success connecting to BLAT server 1>&2
	# output file is in working directory so we need to move to the data directory before running CREST
	if [[ ! -f in.bam.cover ]] ; then
		cut -f 1 $CX_REFERENCE.fai | parallel extractSClip.pl -i in.bam --ref_genome $CX_REFERENCE -r {}
		cat in.bam.*.cover > in.bam.cover
		echo extractSClip.pl complete. Runnin CREST.pl
	else 
		echo Skipping extractSClip.pl - in.bam.cover already exists
	fi
	cut -f 1 $CX_REFERENCE.fai | parallel CREST.pl $CX_CALLER_ARGS -l $CX_READ_LENGTH --tmpdir $CX/tmp -f in.bam.cover -d in.bam --ref_genome $CX_REFERENCE -t $CX_REFERENCE.2bit --blatserver localhost --blatport $BLAT_PORT -r {} -p in.bam.{}
	cat in.bam.*.predSV.txt > in.bam.predSV.txt
	$BASE_DIR/crest2vcf.py < in.bam.predSV.txt > $CX.vcf
	"
	xc_exec
done

