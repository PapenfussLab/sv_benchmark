#!/usr/bin/env Rscript
library("optparse")
parser <- OptionParser(option_list=list(
	make_option(c("-p", "--plot"), action="store_true", default=FALSE, help="Generate plots")
))
arguments = parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args
write(commandArgs(trailingOnly=TRUE), stderr())
write(names(opt), stderr())
write(opt$plot, stderr())

library(ggplot2)
library(scales)

calldf <- .CachedLoadCallsForId(
	datadir="./data.na12878",
	metadata=LoadCachedMetadata("./data.na12878"),
	id="676904146cce653bb2e31a77b4d3e037",
	maxgap=200,
	sizemargin=0.25,
	ignore.strand=TRUE,
	requiredHits=1,
	truthgr=truthgr,
	truthgrName=NULL,
	grtransform=simoptions$grtransform[["DAC"]],
	grtransformName="DAC",
	nominalPosition=FALSE)
calldf <- calldf %>%
	filter(abs(svLen) > 50)%>%
	mutate(Classification = ifelse(tp, "True Positive", ifelse(fp, "False Positive", "False Negative")))


# event size density
ggplot(calldf) +
	aes(x=svLen, y = ..ncount.., colour=Classification) +
	geom_histogram() +
	scale_y_continuous(labels = percent_format()) +
	facet_wrap( ~ Classification) +
	scale_x_log10(limits=c(50, max(calldf$svLen)))

ggplot(calldf %>% filter(!fn)) +
	aes(x=svLen, fill=tp) +
	geom_histogram() +
	scale_x_log10() +
	labs(title="Called Event Size Distribution", fill="Correctly Called")

ggplot(calldf %>% filter(!fp)) +
	aes(x=svLen, fill=tp) +
	geom_histogram() +
		scale_x_log10() +
	labs(title="Actual Event Size Distribution", fill="Correctly Called")

ggplot(calldf %>% filter(!fp)) +
	aes(x=svLen, fill=tp) +
	geom_histogram() +
	scale_x_log10() +
	facet_wrap(~ repeatClass) +
	labs(title="Actual Event Size Distribution", fill="Correctly Called")


# ROC by event size
calldf$svLenBin <- cut(abs(calldf$svLen), breaks=c(0, 50, 100, 200, 300, 400, 500, 1000, 10000,1000000000), labels=c("<=50", "51-100", "101-200", "201-300", "301-400", "401-500", "501-1,000", "1,001-10,000","10,000+"))
rocby <- function(df, ...) {
	groupingCols <- quos(...)
	calldf %>%
			select(Id, CallSet, !!!groupingCols, QUAL, tp, fp, fn) %>%
			rbind(calldf %>% select(Id, CallSet, !!!groupingCols) %>% distinct(Id, CallSet, !!!groupingCols) %>% mutate(QUAL=max(calldf$QUAL) + 1, tp=0, fp=0, fn=0)) %>%
			#filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
			group_by(Id, CallSet, !!!groupingCols) %>%
			arrange(desc(QUAL)) %>%
			mutate(events=sum(tp) + sum(fn)) %>%
			mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn)) %>%
			# each QUAL score is a point on the ROC plott
			group_by(Id, CallSet, !!!groupingCols, QUAL, events) %>%
			summarise(tp=max(tp), fp=max(fp), fn=max(fn)) %>%
			# QUAL scores with the same number of tp calls can be merged on the ROC plot
			group_by(Id, CallSet, !!!groupingCols, tp, events) %>%
			summarise(fp=max(fp), fn=max(fn), QUAL=min(QUAL))
}
ggplot(rocby(calldf)) +
	aes(y=tp, x=fp, colour=Id) +
	geom_line() +
	labs(title="ROC")

ggplot(rocby(calldf, repeatClass)) +
	aes(y=tp, x=fp, colour=Id) +
	geom_line() +
	facet_wrap(~ repeatClass) +
	labs(title="ROC")

ggplot(rocby(calldf, svLenBin)) +
	aes(y=tp, x=fp, colour=Id) +
	geom_line() +
	facet_wrap(~ svLenBin) +
	labs(title="ROC")

ggplot(rocby(calldf, svLenBin, repeatClass)) +
	aes(y=tp, x=fp, colour=Id) +
	geom_line() +
	facet_grid(svLenBin ~ repeatClass) +
	labs(title="ROC")

# what is the repeat context of the FP calls?

#############################
# sequencce homology context
calldf$ihomlenBin <- cut(abs(calldf$ihomlen), breaks=c(0, 5, 10, 25, 50, 100, 300, 1000))

ggplot(calldf) +
	aes(x=ihomlen, fill=Classification) +
	geom_histogram() +
	scale_x_log10()

calldf %>%
	group_by(ihomlenBin) %>%
	summarise(falsePositiveRate=sum(fp)/sum(fp+tp))
ggplot(calldf %>%
	group_by(ihomlenBin) %>%
	summarise(falseNegativeRate=sum(fn)/sum(tp+fn))) +
	aes(x=ihomlenBin, y=falseNegativeRate) +
	geom_point() +
	labs("False Negative Rate by size of breakpoint sequence homology")





# TODO:
 * SNP/INDEL context of SV calls
 * actual coverage vs mean coverage
 * mappabiliyumappability of surrounding bases
 * event size distribution vs real event size distribution































