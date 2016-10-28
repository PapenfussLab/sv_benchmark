#install.packages("xlsx")
library(ggplot2)
library(xlsx) # make sure jre\bin\server in on path: http://stackoverflow.com/questions/2399027/cannot-load-rjava-because-cannot-load-a-shared-library
# once-off vertor graphics conversion
#PostScriptTrace("../assets/iconmonstr-download-icon.ps", "../assets/iconmonstr-download-icon.rgml.xml")
#PostScriptTrace("../assets/iconmonstr-text-file-4-icon.ps", "../assets/iconmonstr-text-file-4-icon.rgml.xml")


setwd("C:/dev/ws/indelappraisal")
sw <- read.xlsx("../docs/IndelSoftware.xlsx", sheetName="SV", startRow=2, stringsAsFactors=FALSE)
# strip out all rows after the first blank in the Software column
sw <- sw[1:(min(which(is.na(sw$Software)))-1), ]
sw$Family <- ifelse(is.na(sw$Family), sw$Software, sw$Family)

# Group by Family
swMap <- data.frame(Software=sw$Software, Family=sw$Family, stringsAsFactors=FALSE)
swOrder <- data.frame(Family=swMap$Family, ordinal=rep(0, length(swMap$Family)))
swOrder <- swOrder[!duplicated(swOrder),]
swOrder$ordinal <- nrow(swOrder) - seq_len(nrow(swOrder)) + 1

pubs <- data.frame(Software=sw$Software, PublicationDate=sw$Date.of.Publication, Evidence=sw$Evidence, stringsAsFactors=FALSE)
pubs <- merge(pubs, swMap)

versions <- read.xlsx("../docs/IndelSoftware.xlsx", sheetName="SoftwareReleases")
versions <- versions[order(versions$ReleaseDate),]
versions <- versions[versions$Software!="MAQ",]

earliest <- rbind(data.frame(Family=versions$Family, Date=versions$ReleaseDate),
									data.frame(Family=pubs$Family, Date=pubs$PublicationDate))
earliest <- earliest[order(earliest$Date), ]
earliest <- earliest[!duplicated(earliest$Family), ]
earliest$Software <- earliest$Family
earliest <- merge(earliest, pubs)

events <- rbind(data.frame(Family=versions$Family, Date=versions$ReleaseDate, Type=rep("Software Release", nrow(versions))),
								data.frame(Family=pubs$Family, Date=pubs$PublicationDate, Type=rep("Publication", nrow(pubs))))



# add ordinal to output tables
versions <- merge(versions, swOrder)
pubs <- merge(pubs, swOrder)
earliest <- merge(earliest, swOrder)
events <- merge(events, swOrder)

# don't print the family publication withat publication date since we special case that
pubs$Software <- ifelse(pubs$Software == pubs$Family, rep("", nrow(pubs)), pubs$Software)
# remove ugly version numbers
versions$Version[versions$Family=="samtools" & versions$Version %in% c("0.1.2", "0.1.4", "0.1.10", "0.1.11", "0.1.12", "0.1.14", "0.1.15") ] <- NA
versions$Version[versions$Family=="VariationHunter" & versions$Version %in% c("0.3") ] <- NA
versions$Version[versions$Family=="SVMerge" & versions$Version %in% c("1.0r6", "1.0r10", "1.0r19", "1.1r36") ] <- NA
versions$Version[versions$Family=="VarScan2" & versions$Version %in% c("2.2.3", "2.2.7", "2.2.11", "2.2.12", "2.3.2") ] <- NA
versions$Version[versions$Family=="GASV" & versions$Version %in% c("0.9", "1.5.2","2.0.1") ] <- NA
versions$Version[versions$Family=="SOAPsplice" & versions$Version %in% c("1.2", "1.3","1.4") ] <- NA
versions$Version[versions$Family=="SNVer" & versions$Version %in% c("0.0.2", "0.1.1", "0.1.2", "0.2.1", "0.2.2", "0.3.1", "0.4.1") ] <- NA
versions$Version[versions$Family=="DELLY" & versions$Version %in% c("0.0.3", "0.0.4", "0.0.5", "0.0.6", "0.0.7", "0.0.8", "0.0.10") ] <- NA
versions$Version[versions$Family=="CLEVER" & versions$Version %in% c("1.1") ] <- NA
versions$Version[versions$Family=="RetroSeq" & versions$Version %in% c("1.31", "1.32", "1.33", "1.34", "1.41") ] <- NA


ggplot(data=NULL, aes(y=ordinal)) +	
	scale_x_date(limits=c(as.Date("2008-11-01"), as.Date("2013-10-01"))) +
	geom_text(data=versions,aes(label=Version, x=ReleaseDate, vjust=1.1), size=3) + 
	geom_text(data=pubs,aes(label=Software, x=PublicationDate, colour=Evidence, vjust=-0.4)) + 
	geom_text(data=earliest, aes(label=Software, x=Date, colour=Evidence, hjust=1.2)) + 
	geom_segment(data=earliest, aes(x=Date, y=ordinal, xend=as.Date("2013-10-01"), yend=ordinal), linetype=3) + 
	geom_point(data=events, aes(x=Date, shape=Type), colour="black", fill="black") + 
	scale_shape_manual(values=c(20, 25)) +
	scale_y_continuous(breaks=NULL) + 
	geom_hline(yintercept=7.5) + 
	annotate("text", label="Specialised", x=as.Date("2008-10-15"), y=3) + 
	annotate("text", label="Small indel", x=as.Date("2008-10-15"), y=38) + 
	geom_hline(yintercept=33.5) + 
	theme_bw() + theme(
		axis.text.y=element_blank(),
		axis.title.y=element_blank(),
		axis.ticks.y=element_blank(),
		panel.border = element_blank()) +
	#theme(legend.position="none") +
#	ggtitle("INDEL calling software") +
	xlab("Release / Publication Date")

	# TODO: add line from min(pub, version) to max(all(pub, version))
ggsave("publication.jpg", width=18, height=18)
ggsave("publication.pdf", width=18, height=18)

# http://stackoverflow.com/questions/2181902/how-to-use-an-image-as-a-point-in-ggplot
# https://github.com/hadley/ggplot2/blob/master/R/geom-point-.r
# http://cran.r-project.org/web/packages/grImport/grImport.pdf
#library(proto)
#library(grImport) #install.packages("grImport")
# GeomPointImage <- proto(ggplot2:::Geom, {
#	 objname <- "point_image"
#	 desc <- "Custom grob"
#	 draw_groups <- function(., ...) .$draw(...)
#	 draw <- function(., data, scales, coordinates, just=c("centre", "centre"), ...) {
#		 picture <- readPicture("../assets/iconmonstr-download-icon.rgml.xml")
#		 pictureGrob(picture)
#		 #with(coordinates$transform(data, scales), pictureGrob(picture, x, y))
#	 }
#	 required_aes <- c("x", "y")
#	 default_aes <- function(.) aes(x=0.5, y=0.5)
#	 default_stat <- function(.) StatIdentity
#	 guide_geom <- function(.) "custom"
#	 draw_legend <- function(., data, ...) {	}
#	 icon <- function(.) {	}
#	 examples <- function() { }		
# })
# geom_point_image <- function (mapping = NULL, data = NULL, stat = "identity", position = "identity", na.rm = FALSE, ...) {
#	 GeomPointImage$new(mapping = mapping, data = data, stat = stat, position = position, na.rm = na.rm, ...)
# }
#