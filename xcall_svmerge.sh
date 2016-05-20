#!/bin/bash
#
# runs svmerge against bams
#

# TODO: were can I find old versions of the required software?

. common.sh
CALLER=svmerge/1.2r37
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	PINDEL_DIR=$(grep -l CX_CALLER=pindel $(grep -l "CX_BAM=$CX_BAM" $DATA_DIR/*.metadata) | cut -b 1-32)
	BREAKDANCER_DIR=$(grep -l CX_CALLER=pindel $(grep -l "CX_BAM=$CX_BAM" $DATA_DIR/*.metadata) | cut -b 1-32)
	if [[ "$PINDEL_DIR" ==  "" ]] ; then
		echo "Missing pindel output for $BAM - skipping"
		continue
	fi
	if [[ "$BREAKDANCER_DIR" ==  "" ]] ; then
		echo "Missing pindel output for $BAM - skipping"
		continue
	fi
	cx_save
	XC_MULTICORE=1
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
cat > project.config << EOF
## General parameters

project=project
name=project
version=1
svdir=sv_calls
projdir=$CX/project
chrRange=1-22
chrOther=X Y
gender=female
exedir=$BASE_DIR/tools/$CALLER/
species=homo_sapiens

defaultQueue=bioinformatics

## Location of reference files

bam=$CX_BAM
bai=$CX_BAM.bai
#bamdir=/path/to/directory/of/chromosome/bams/
# specify either chrrefdir or reffile
#chrrefdir=/path/to/chomosome/dir/
reffile=$CX_REFERENCE

## Location of SV calling software

# Samtools Version 0.1.12-10 and older (must have samtools pileup)
samtools=/usr/local/bioinfsoftware/samtools/samtools-0.1.10/samtools 

cnddir=/path/to/cnD/

bam2conf=/path/to/breakdancer/bam2cfg.pl
bdexe=/path/to/breakdancer/BreakDancer_Max

pinexe=/path/to/pindel_0.2.3/pindel

secexe=/path/to/SECluster.pl

#rdxdir=/path/to/rdxplorer/

## List SV callers 

# breakdancer
# pindel
# sec
# rdx
# cnd

callerlist=pindel breakdancer sec cnd 

## SV caller parameters

defaultQueue=normal

## BreakDancerMax
BDconf=1
BDconfParams=-c 7 -n 10000
# do not include -o below in BDparams
BDparams=-c 7 -m 1000000 -q 20  
BDcopynum=2
BDmem=3000
BDqueue=bioinformatics



# Pindel
PDconf=/path/to/pindel_config
PDoptParams=-x 5 -v 1000
PDfiltermem=3000
PDmem=2000
PDqueue=normal

# SEC - single end clusters
SECfilter=1
SECqual=20
SECmin=5
SECminCluster=5
SECmax=500
SECfilterMem=2000
SECfilterQueue=normal
SECmem=4000
SECqueue=normal


# cnD
CNDnohet=1
CNDpileup=1
CNDgccorrect=1
CNDsnprate=0.001
# do not include --prefix below
CNDparams=--smooth=100 --repeat-cutoff=0.35 
CNDpileupMem=3000
CNDpileupQueue=normal
CNDmem=2000
CNDqueue=normal

# RDXplorer

#RDXSplitBam=1
#RDXqueue=long
#RDXmem=2000


## Filter and Merge Steps

bedexe=/path/to/BEDTools/bin/intersectBed
gaps=/path/to/NCBI37/hg19_gap.txt
gapsBuffer=600
centel=/path/to/NCBI37/hg19_cen_tel.txt
centelBuffer=1000000
#filterOther=/path/to/other/regions.txt
#filterOtherBuffer=100

# fraction of the SV call that must overlap with gap/cen/tel
gapOverlap=0.25

# breakdancer
BDscore=25
BDrs=2

# pindel
PDscore=30
PDsupports=10

# rdx
#RDXscore=5
#RDXout=NA18506.rdx.txt

otherFilterGaps=1
otherDelCalls=newcaller/del.tab
otherInsCalls=newcaller/ins.tab
otherInvCalls=newcaller/inv.tab
otherLossCalls
otherGainCalls


## Local assembly

# Assemble from specific config files only

#assemMin=17
#assemMax=17

outdir=.
makeConfig=1
joblimit=75

submatrix=/path/to/SVMerge/submat.txt
exonerateExe=exonerate-2.2.0

# Assemble only if the exonerate output doesn't exist
checkdone=1

# bsub queue and memory
assemQueue=normal
#assemMem=2000

# subseq: alignments to part of a chr rather than whole chr.
subseq=1

# Velvet parameters

velvet=1
velveth=velveth
velvetg=velvetg
hashlen=29
ins_len=220
exp_cov=35
cov_cutoff=2

# ABySS parameters
#abyss=1
#abyss-pe=/path/to/abyss-pe
kmer=25
npairs=10


# Make config file for running exonerate only 
#exonerate=1


## Alignment parsing

bamcheck=1
meanCov=42
zyg=het
offset=1
parseSplit=1
parseSplitLines=250
parseJoblimit=10
parseQueue=normal
parseMem=2000
EOF

$BASE_DIR/tools/$CALLER/
	
	"
	xc_exec
done

