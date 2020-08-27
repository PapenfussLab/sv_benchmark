#!/bin/bash
if [[ ! -f hg38.fa.out.gz ]] ; then
	wget http://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz
fi
if [[ ! -f hg19.fa.out.gz ]] ; then
	wget http://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm405-db20140131/hg19.fa.out.gz
fi
