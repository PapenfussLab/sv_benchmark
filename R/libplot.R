library(ggplot2)
library(rtracklayer)

theme_set(theme_bw())

scale_y_power4 <- scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)^4, labels=c("0", "", "", "", "", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"))
scale_y_power5 <- scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)^5, labels=c("0", "", "", "", "", "0.5", "", "0.7", "0.8", "0.9", "1.0"))

scale_y01_small <- scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, .5, .75, 1), labels=c("0", "", "0.5", "", "1"))

scale_x_svlen <- scale_x_continuous(
	breaks=2**(0:16),
	labels=c("1", "2", "4", "8", "16", "32", "64", "128", "256", "512", "1k", "2k", "4k", "8k", "16k", "32k", "64k"),
	minor_breaks=c(1,2,3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 24, 28, 32, 48, 64, 80, 96, 112, 128, 160, 192, 224, 256, 288, 320, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536),
	trans="log10",
	expand = c(0,0))
scale_x_svlen_short <- scale_x_continuous(
	breaks=2**(0:16),
	labels=c("1", "", "", "", "16", "", "", "", "256", "", "", "", "4k", "", "", "", "64k"),
	minor_breaks=c(1,2,3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 24, 28, 32, 48, 64, 80, 96, 112, 128, 160, 192, 224, 256, 288, 320, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536),
	trans="log10",
	expand = c(0,0))
scale_x_svlen_medium <- scale_x_continuous(
	breaks=2**(0:16),
	labels=c("1", "2", "4", "8", "16", "32", "64", "128", "", "512", "", "2k", "", "8k", "", "32k", ""),
	minor_breaks=c(1,2,3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 24, 28, 32, 48, 64, 80, 96, 112, 128, 160, 192, 224, 256, 288, 320, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536),
	trans="log10",
	expand = c(0,0))
scale_x_log_fp <- scale_x_continuous(breaks=c(0, 10, 100, 1000, 10000, 100000)+1,
	labels=c("0", "10", "100", "1k", "10k", "100k"),
	minor_breaks=1+c(1*c(2,4,6,8),10*c(2,4,6,8),100*c(2,4,6,8),1000*c(2,4,6,8),10000*c(2,4,6,8)),
	trans="log10")
scale_x_log_fp_short <- scale_x_continuous(breaks=c(1, 11, 101, 1001, 10001, 100001),
	labels=c("0", "", "100", "", "10k", ""),
	minor_breaks=c(),
	trans="log10")


saveplot <- function(file=file, ...) {
	ggsave(paste0("png/", file, ".png"), ...)
	ggsave(paste0("eps/", file, ".eps"), ...)
}


plotPrecRecall <- function(plotdf) {
	if (nrow(plotdf) == 0) {
		write("plotPrecRecall: no data!", stderr())
		return(NULL)
	}
	# display breakpoint counts instead of breakend counts
	plotdf$tp <- plotdf$tp/2
	plotdf$fp <- plotdf$fp/2
	p <- ggplot(plotdf %>% arrange(desc(QUAL))) +
		aes(group = paste(Id, CallSet), y = precision, x = tp, colour=caller, linetype=CallSet) +
		geom_line(size=1) +
		scale_colour_brewer(palette = "Paired") +
		scale_y_continuous(limits=c(0,1)) +
		labs(title = "", y = "Precision", x = "Recall (true positive count)")
	return(p)
}

plotPrecRecallRepeat <- function(plotdf) {
	if (nrow(plotdf) == 0) {
		write("plotPrecRecallRepeat: no data!", stderr())
		return(NULL)
	}
	# display breakpoint counts instead of breakend counts
	plotdf$tp <- plotdf$tp/2
	plotdf$fp <- plotdf$fp/2
	if (nrow(plotdf) == 0) return(NULL)
	p <- ggplot(plotdf %>% arrange(desc(QUAL))) +
		aes(group = paste(Id, CallSet), y = precision, x = tp, colour=caller, linetype=CallSet) +
		geom_line(size=1) +
		facet_wrap(~ repeatClass, scales="free") +
		scale_colour_brewer(palette = "Paired") +
		scale_y_continuous(limits=c(0,1)) +
		labs(title = "Precision-Recall by RepeatMasker repeat class", y = "Precision", x = "Recall (true positive count)")
	return(p)
}
plotRocLinear <- function(plotdf) {
	if (nrow(plotdf) == 0) {
		write("plotRoc: no data!", stderr())
		return(NULL)
	}
	# display breakpoint counts instead of breakend counts
	plotdf$tp <- plotdf$tp/2
	plotdf$fp <- plotdf$fp/2
	p <- ggplot(plotdf %>% arrange(desc(QUAL))) +
		aes(group = paste(Id, CallSet), y = tp, x = fp, colour=caller, linetype=CallSet) +
		geom_line(size=1) +
		coord_cartesian(ylim=c(0, max(plotdf$tp)), xlim=c(0, max(plotdf$tp))) +
		scale_colour_brewer(palette = "Paired") +
		labs(title = "", y = "True Positives", x = "False positives")
	return(p)
}
plotFacetedRocLogFP <- function(plotdf, linetype, colour, facet_eval_string=NULL) {
	if (nrow(plotdf) == 0) {
		write("plotFacetedRoc: no data!", stderr())
		return(NULL)
	}
	# display breakpoint counts instead of breakend counts
	plotdf$tp <- plotdf$tp/2
	plotdf$fp <- plotdf$fp/2
	p <- ggplot(plotdf %>% arrange(desc(QUAL))) +
		aes(group=paste(Id, CallSet), y=sens, x=fp+1) +
		aes_string(linetype=linetype, colour=colour) +
		geom_line() +
		scale_x_log_fp +
		labs(title="", y="Sensitivity", x="False Positives")
	if (!is.null(facet_eval_string)) {
		p <- p + facet_grid(facet_eval_string)#eval(parse(text=paste("caller ~ ", paste(simfacets[!(simfacets %in% c(input$simlinetype, input$simcolour))], collapse=" + ")))))
	}
	return(p)
}
plotFacettedSensByEventSize <- function(plotdf, linetype, colour, facet_eval_string=NULL) {
	if (nrow(plotdf) == 0) {
		write("plotPrecRecall: no data!", stderr())
		return(NULL)
	}
	p <- ggplot(plotdf) +
		aes(group=paste(Id, CallSet), x=abs(svLen), y=sens) +
		aes_string(linetype=linetype, colour=colour) +
		geom_line(size=0.5) +
		scale_x_svlen +
		labs(title="", y="Sensitivity", x="Event size")
	if (!is.null(facet_eval_string)) {
		p <- p + facet_grid(facet_eval_string)#eval(parse(text=paste("caller ~ ", paste(simfacets[!(simfacets %in% c(input$simlinetype, input$simcolour))], collapse=" + ")))))
	}
	return(p)
}


# TODO
# visibility of
# ..., Eichler 2016

