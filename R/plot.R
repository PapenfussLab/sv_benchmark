#load("C:\\dev\\ws\\indelappraisal\\details.RData")

#
# Data Transforms
#

# add calculated fields:
detail$tp$svLenError <- abs(detail$tp$svLen - detail$tp$svLenCalled)

# Group SVs into buckets
#svGroupLabels <- c("Tiny (1-10]", "Small (10-50]", "Medium (50-200]", "Large (200, 10000]", "Huge 10000+")
#svGroupBreaks <- c(0, 10, 50, 200, 10000, 1e15)
#detail$tp$svGroup <- cut(abs(detail$tp$svLen), breaks=svGroupBreaks, labels=svGroupLabels, ordered_result=TRUE)
#detail$fp$svGroup <- cut(abs(detail$fp$svLen), breaks=svGroupBreaks, labels=svGroupLabels, ordered_result=TRUE)
#detail$fn$svGroup <- cut(abs(detail$fn$svLen), breaks=svGroupBreaks, labels=svGroupLabels, ordered_result=TRUE)
# Aggregate detailed data into counts
#summaryExpected <- ApplyDetailAggregation(cast, detail, caller + aligner + softclipped + fragLength + readLength + readDepth + svType + svGroup ~ ., fun.aggregate=length, value="svLen")

#
# Data Sanity Checks
#
# Shouldn't be able to match SV without matching both ends within tolerance
detail$tp[detail$tp$svLenError > 300,]
# Counts should match underlying reference VCF
#refVcf <-	vcfs[["78927fc52fe33cf2f61898d028296e8e"]]
#refVcfdf <- data.frame(svLen=sapply(info(refVcf)$SVLEN, "[", 1), svType=info(refVcf)$SVTYPE)
#refVcfdf$svGroup <- cut(abs(refVcfdf$svLen), breaks=svGroupBreaks, labels=svGroupLabels, ordered_result=TRUE)
#summaryExpected <- cast(refVcfdf, svType + svLen ~ ., fun.aggregate=length, value="svLen")
#colnames(summaryExpected)[colnames(summaryExpected) == "(all)"] <- "expectedRefTotal"
#summaryExpected <- summaryExpected[summaryExpected$reftotal - summaryExpected$expectedRefTotal != 0,]

#
# Plots
#
if (Sys.info()[["sysname"]] == "Linux") {
	setwd("~/plots")
} else {
	setwd("C:/dev//plots")
}
#install.packages("ggplot2")
library(ggplot2)

# graph common elements
log2SvLenAxis <- scale_x_continuous(breaks = c(-17, -13, -9, -5, -1, 1, 5, 9, 13, 17), 
																		labels = c("-64k", "-4k", "-256", "-16", "-1", "1", "16", "256", "4k", "64k"))

ext <- ".jpg"
fragSize <- 400
dt <- ApplyDetailAggregation(count, detail, vars=c("caller", "aligner", "softclipped", "fragLength", "readLength", "readDepth", "svLen"))
#for (rd in unique(dt$readDepth)) { #rd <- 100
rd <- 100
# Caller sensitivity
dt$log2SvLen <- sign(dt$svLen) * (log2(abs(dt$svLen)) + 1) # +1 so 1 and -1 aren't both 0
dt$Alignment <- factor(ifelse(dt$softclipped, "Local", "Global"), levels=c("Local", "Global"))
dtSens <- dt[!is.nan(dt$sensitivity) & dt$fragLength==fragSize & dt$readDepth==rd, ]
#dtSens$softclipped <- factor(dtSens$softclipped, levels=c(TRUE, FALSE), ordered=TRUE) # flip order so TRUE has solid line
ggplot(dtSens[dtSens$Alignment=="Local", ], aes(x=log2SvLen, y=sensitivity)) +
	geom_line(aes(group=caller, colour=caller)) +
	facet_grid(aligner ~ readLength) +
	theme_bw() + ggtitle("Sensitivity per Aligner")
ggsave(paste0("rd", rd, "_", "sensitivity_aligners", ext), width=10, height=7.5)
ggplot(dtSens[dtSens$Alignment=="Local", ], aes(x=log2SvLen, y=sensitivity)) +
	geom_line(aes(group=interaction(aligner, readLength), colour=factor(readLength), linetype=aligner)) +
	facet_wrap( ~ caller) +
	theme_bw() + ggtitle("Caller sensitivity (soft-clipped alignment)")
ggsave(paste0("rd", rd, "_", "sensitivity_callers", ext), width=10, height=7.5)
ggplot(dtSens[dtSens$Alignment=="Global", ], aes(x=log2SvLen, y=sensitivity)) +
	geom_line(aes(group=interaction(aligner, readLength), colour=factor(readLength), linetype=aligner)) +
	facet_wrap( ~ caller) +
	theme_bw() + ggtitle("Caller sensitivity (end-to-end alignment)")
ggsave(paste0("rd", rd, "_", "sensitivity_callers_nosoftclip", ext), width=10, height=7.5)
for (aligner in unique(dtSens$aligner)) {
	ggplot(dtSens[dtSens$aligner==aligner, ], aes(x=log2SvLen, y=sensitivity)) +
		geom_line(aes(group=interaction(Alignment, readLength), colour=factor(readLength), linetype=Alignment)) +
		facet_wrap( ~ caller) +
		theme_bw() + 
		log2SvLenAxis + 
		labs(colour="Read Length", x="Variant Size", y="Sensitivity") + 
		ggtitle(paste("Caller sensitivity:", aligner))
	ggsave(paste0("rd", rd, "_", "sensitivity_aligner_", aligner, "", ext), width=10, height=7.5)
}
for (caller in unique(dtSens$caller)) {
	ggplot(dtSens[dtSens$caller==caller, ], aes(x=log2SvLen, y=sensitivity)) +
		geom_line(aes(group=interaction(Alignment, readLength), colour=factor(readLength), linetype=Alignment)) +
		facet_wrap( ~ aligner) +
		log2SvLenAxis + 
		labs(colour="Read Length", x="Variant Size", y="Sensitivity") + 
		theme_bw() + ggtitle(paste("Caller sensitivity:", caller))
	ggsave(paste0("rd", rd, "_", "sensitivity_caller_", caller, "", ext), width=10, height=7.5)
}
# Effect of soft-clipping on indel calling
dtscvs <- merge(dtSens[dtSens$Alignment=="Local", ], dtSens[dtSens$Alignment=="Global", ], by=c("caller", "aligner", "fragLength", "readLength", "readDepth", "svLen", "log2SvLen"))
dtscvs$softclipSensitivityImprovement <- dtscvs$sensitivity.x - dtscvs$sensitivity.y
ggplot(dtscvs, aes(x=log2SvLen, y=softclipSensitivityImprovement)) +
	geom_line(aes(group=readLength, colour=factor(readLength))) +
	facet_grid(aligner ~ caller) +
	theme_bw() +
	log2SvLenAxis + 
	labs(colour="Read Length", x="Variant Size", y="Local alignment sensitivity improvement") + 
	ggtitle("Effect Global vs Local alignment on called indel")
ggsave(paste0("rd", rd, "_", "sensitivity_softclipping", ext), width=10, height=7.5)

# Calling rates
for (aligner in unique(dt$aligner)) {
	for (softclip in unique(dt$Alignment)) {
		dtCallingRate <- dt[dt$aligner == aligner & dt$Alignment == softclip & dt$fragLength==fragSize & dt$readDepth==rd,]
		if (nrow(dtCallingRate) > 0) {
			# Calling rates
			dtcr <- ddply(dtCallingRate, ~ caller + readLength + svLen + log2SvLen, summarise, allcalls=sum(tpDup) + sum(fp), correctcalls=sum(tp))
			ggplot(dtcr, aes(x=log2SvLen)) +
				geom_line(aes(y=correctcalls, group=readLength, colour=interaction(factor(readLength)))) +
				geom_line(aes(y=allcalls, group=readLength, colour=interaction(factor(readLength))), linetype="dashed") +
				scale_y_log10() + 
				facet_wrap( ~ caller) +
				log2SvLenAxis + 
				labs(x="Variant Size") + 
				theme_bw() + ggtitle(paste0("Calling rate, ", rd, "x ", aligner, " ", softclip, " alignment"))
			ggsave(paste0("rd", rd, "_", "calling_rate_sc", softclip, "_", aligner, ext), width=10, height=7.5)
		}
	}
}
for (caller in unique(dt$caller)) {
	dtCallingRate <- dt[dt$caller == caller & dt$fragLength==fragSize & dt$readDepth==rd,]
	if (nrow(dtCallingRate) > 0) {
		dtcr <- ddply(dtCallingRate, ~ aligner + readLength + svLen + log2SvLen + Alignment, summarise, calls=sum(tp), fpCalls=sum(tpDup) + sum(fp))
		ggplot(dtcr, aes(x=log2SvLen)) +
			geom_line(aes(y=calls, group=readLength, colour=interaction(factor(readLength)))) +
			geom_line(aes(y=fpCalls, group=readLength, colour=interaction(factor(readLength))), linetype="dashed") +
			scale_y_log10() + 
			facet_grid(Alignment ~ aligner) +
			log2SvLenAxis + 
			labs(colour="Read Length", x="Variant Size") + 
			theme_bw() + ggtitle(paste0(caller, " calling rate, ", rd, "x "))
		ggsave(paste0("rd", rd, "_", "calling_rate_by_aligner_", caller, ext), width=10, height=7.5)
	}
}
# Sensitivity vs Read Depth 
for (aligner in unique(dt$aligner)) {
	for (softclip in unique(dt$Alignment)) {
		dtCurrentAligner <- dt[dt$aligner == aligner & dt$Alignment == softclip & dt$fragLength==fragSize, ]
		if (nrow(dtCurrentAligner) > 0) {
			# Calling rates
			dtrdcalls <- ddply(dtCurrentAligner, ~ readDepth + caller + readLength, summarise, correctcalls=sum(tp), falsepositivecalls=sum(tpDup) + sum(fp))
			ggplot(dtrdcalls, aes(x=readDepth)) +
				geom_line(aes(group=readLength, colour=interaction(factor(readLength)), y=correctcalls)) +
				geom_line(aes(group=readLength, colour=interaction(factor(readLength)), y=falsepositivecalls), linetype="dashed") +
				facet_wrap( ~ caller) +
				labs(colour="Read Length", x="Read Depth", y="Variant calls") + 
				scale_x_continuous(breaks = c(4, 8, 15, 30, 60, 100)) + 
				theme_bw() + ggtitle(paste0("Read Depth ", aligner, " ", softclip, " alignment"))
			ggsave(paste0("read_depth_sc", softclip, "_", aligner, ext), width=10, height=7.5)
		}
	}
}
dtrdcalls <- ddply(dt[dt$fragLength==fragSize, ], ~ readLength + readDepth + caller + aligner + Alignment, summarise, correctcalls=sum(tp), calls=sum(tpDup) + sum(fp), tpDupCalls=sum(tpDup))
ggplot(dtrdcalls, aes(x=readDepth, group=factor(readLength), colour=factor(readLength))) +
	geom_line(aes(y=calls), linetype="dashed") +
	geom_line(aes(y=tpDupCalls), linetype="dotted") +
	geom_line(aes(y=correctcalls)) +
	facet_grid(Alignment + aligner ~ caller) +
	theme_bw() + ggtitle(paste0("Read Depth"))
ggsave(paste0("read_depth", ext), width=10, height=7.5)

# calculate best aligner for each caller
# TODO!
ddply(dt[dt$fragLength==fragSize & dt$readDepth==rd,], ~ caller + aligner + Alignment + readLength, summarise, tp=sum(tp), fp=sum(fp), fn=sum(fn))



# Sensitivity vs False Discovery Rate
library(reshape)
#dtroc <- c("caller", "aligner", "softclipped", "fragLength", "readLength", "readDepth", "svLen", "sensitivity", "fdr")
rocfn <- detail$fn
rocfp <- detail$fp
roctp <- detail$tp
# force numeric quality scores field
rocfp$qual[!is.numeric(rocfp$qual) | is.na(rocfp$qual)] <- 0
roctp$qual[!is.numeric(roctp$qual) | is.na(roctp$qual)] <- 0
rocfn$qual <- -1 # add a sentinal qual so we have a non-zero fixed count of fn when we cumsum
# strip unwanted columns
groupingColumns <- c("caller", "aligner", "softclipped", "fragLength", "readLength", "readDepth", "svLen", "qual") # qual must be last
rocfn <- rocfn[groupingColumns]
rocfp <- rocfp[groupingColumns]
roctp <- roctp[c(groupingColumns, "isDuplicate")]
# TP duplicates are considered false positives
rocfp <- rbind(rocfp, roctp[roctp$isDuplicate, ][groupingColumns])
roctp <- roctp[!roctp$isDuplicate, ]
# Add counts and merge
rocfn <- ddply(rocfn, groupingColumns, nrow)
rocfp <- ddply(rocfp, groupingColumns, nrow)
roctp <- ddply(roctp, groupingColumns, nrow)
rocfn <- rename(rocfn, c("V1"="fn"))
rocfp <- rename(rocfp, c("V1"="fp"))
roctp <- rename(roctp, c("V1"="tp"))
dtroc <- merge(rocfn, merge(rocfp, roctp, all=TRUE), all=TRUE)
dtroc$tp[is.na(dtroc$tp)] <- 0
dtroc$fn[is.na(dtroc$fn)] <- 0
dtroc$fp[is.na(dtroc$fp)] <- 0

dtroc <- dtroc[order(-dtroc$qual), ]
dtroc <- ddply(dtroc, head(groupingColumns, -1), summarise,
							 sensitivity=cumsum(tp) / (sum(tp) + sum(fn)),
							 fdr=cumsum(fp) / (cumsum(tp) + cumsum(fp)),
							 precision=cumsum(tp)/(cumsum(tp)+cumsum(fp)))

dtroc <- dtroc[!(dtroc$sensitivity==0 & dtroc$fdr==0),]
# ROC curves

aligner <- "bwamem"
ggplot(dtroc[dtroc$aligner==aligner & dtroc$fragLength==fragSize & dtroc$readDepth==rd, ],
	aes(x=fdr, y=sensitivity, group=caller, colour=caller, shape=caller)) +
	geom_line() +
	geom_point() + 
	facet_grid(readLength ~ svLen) +
	theme_bw() + ggtitle(paste("ROC ", aligner))
ggsave(paste0("rd", rd, "_", "roc_aligner_bwamem", ext), width=10, height=7.5)


ggplot(dtroc,
			 aes(x=precision, y=sensitivity, colour=factor(svLen), shape=factor(readDepth))) +
	geom_line() +
	geom_point() + 
	facet_grid(readLength ~ aligner) +
	scale_x_reverse() + 
	theme_bw() + ggtitle(paste("ROC ", aligner))


head(dtroc[dtroc$aligner=="bwamem" & dtroc$fdr>0.4& dtroc$fdr<1,])
head(dtroc[dtroc$aligner=="bwamem" & dtroc$svLen==-20 & dtroc$readLength==36 & dtroc$readDepth==100,])


