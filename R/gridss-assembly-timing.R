library(GRanges)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggbio)
library(StructuralVariantAnnotation)
#library(BSgenome.Hsapiens.UCSC.hg38)

# load hg38 contigs used
hg38fai <- read.table("~/Papenfuss_lab/resources/gatk/hg38/v0/Homo_sapiens_assembly38.fasta.fai", col.names=c("seqname", "seqlength", "offset", "linebases", "linebytes"))
row.names(hg38fai) <- hg38fai$seqname

# load hg38 repeatmasker annotations
# pre-processed using tail -n +4 hg38.sorted.fa.out | sed -e 's/^/ /' | tr -s " " | cut -f 6,7,8,11,12 -d " " | tr " " "," > repeatmasker.csv
rmdt <- read.csv("~/Papenfuss_lab/projects/reference_genomes/human/hg38/repeatmasker/repeatmasker.csv", col.names=c("chr", "start", "end", "repeat", "repeatClass"))
grrm <- GRanges(seqnames=rmdt$chr,
                ranges=IRanges(start=rmdt$start, end=rmdt$end),
                repeatType=rmdt$repeat.,
                repeatClass=rmdt$repeatClass)
grrm$repeatSummary <- str_extract(grrm$repeatClass, "^[^?/]*")
grrm$repeatSummary[str_detect(grrm$repeatSummary, "RNA")] <- "*RNA"
remove(rmdt)


log <- read.csv("~/hg38.assembly.bam.events.csv", header=FALSE, stringsAsFactors=FALSE)

loadlog <- log %>% filter(V3 == "load") %>%
  mutate(V4=as.character(V4)) %>%
  mutate(V6=pmin(hg38fai[V4,]$seqlength, V6)) %>% # adjust to end of chr
  filter(V6>V5) %>% # TODO: why is this needed at all?!?
  arrange(match(V4, hg38fai$seqname), V5)
gr <- GRanges(
  seqnames=loadlog$V4,
  ranges=IRanges(start=loadlog$V5, end=loadlog$V6),
  strand="*",
  filtered=as.logical(loadlog$V8),
  density=loadlog$V7/(loadlog$V6-loadlog$V5),
  time=loadlog$V9)
genome(gr) <- "hg38"
seqlengths(gr) <- hg38fai[names(seqlengths(gr)),]$seqlength
gr <- gr[order(seqnames(gr))]
# Add RepeatMasker annotation
gr$repeatSummary <- ""
gr$repeatSummary <- grrm$repeatSummary[findOverlaps(gr, grrm, select="first")]
gr$repeatSummary[is.na(gr$repeatSummary)] <- ""

vcf <- readVcf("~/hg38.gridss.vcf", "hg38")
hcgr <- rowRanges(vcf)[rowRanges(vcf)$FILTER %in% c(".", "PASS"),]
gr$nearHighConfidenceCall <- overlapsAny(gr, hcgr, maxgap=500)


ggplot(data.frame(gr)) +
  aes(x=density, fill=nearHighConfidenceCall) +
  geom_histogram(bins=100) +
  facet_wrap(~ repeatSummary) +
  scale_x_log10(limits=c(1e-4, NA)) +
  scale_y_log10()

ggplot(data.frame(gr)[gr$nearHighConfidenceCall,]) +
  aes(x=density, fill=nearHighConfidenceCall) +
  geom_histogram(bins=100)

ggplot(data.frame(gr)) +
  aes(x=density, colour=repeatSummary) +
  geom_density(size=2) +
  scale_x_log10() +
  scale_colour_brewer(palette="Paired") + 
  coord_cartesian(xlim=c(1e-3, max(data.frame(gr)$density)))

plotGrandLinear(gr[gr$density>=1], aes(y=log10(density)),
                space.skip=0) +
  theme(axis.text.x=element_text(angle=-90, hjust=0))

ggplot(data.frame(gr)) +
  aes(x=time, y=density) +
  geom_hex(bins=100) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ repeatSummary)



ggplot(data.frame(gr)) +
  aes(x=density, fill=nearHighConfidenceCall) +
  geom_histogram(bins=100) +
  scale_x_log10(limits=c(1e-4, NA)) +
  scale_y_log10()









# Convert timing information into per-base estimage
# Note: it is actually an interval in which contigs are emitted out of order

# histogram of GRIDSS per base timing vs repeat-masker annotation

# histogram of maximum interval graph size vs repeat-masker annotation
