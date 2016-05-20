dataDir <- "data.hg19chr12"
if (Sys.info()[["sysname"]] == "Linux") {
	setwd(paste0("~/dev/ws/indelappraisal/", dataDir))
} else {
	setwd(paste0("W:/i/", dataDir))
}
source("C:/dev/ws/indelappraisal/libindelappraisal.R")
metadata <- LoadMetadata()
vcfs <- LoadVcfs(metadata)
detail <- GenerateDetailedCallInformation(metadata, vcfs)
#save(metadata, detail, file=paste0("../", dataDir, ".RData"))
# load("../details.RData")
# TODO: gridss - remove duplicate identical breakend calls

# only IDSV
#caller <- "idsv"
#vcfList <- metadata[metadata$CX_CALLER==caller | metadata$Id %in% GetMetadataId(metadata[metadata$CX_CALLER==caller,]$CX_REFERENCE_VCF),]$File
#vcfList <- vcfList[!is.na(vcfList)]
#vcfList <- c(gsub(".metadata", ".vcf", vcfList), gsub(".metadata", ".reference.vcf", vcfList))
#vcfList <- vcfList[file.exists(vcfList)]
#vcfs <- LoadVcfs(metadata, vcfList)



