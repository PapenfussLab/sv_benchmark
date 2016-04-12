GetId <- function(filenames) {
	cf <- as.character(filenames)
	if (length(cf) == 0) {
		return(character(0))
	} else if (all(is.na(cf))) {
		return(cf)
	} else {
		return(str_match(basename(as.character(filenames)), "([^.]+)(\\..*)*$")[,2])
	}
}
test_that("GetId", {
	expect_equal(c("a", "m", "bin"),
		GetId(c("a", "C:\\directory\\m.as.ds.dsa.ss.vcf", "/usr/local/bin")))
	expect_equal(character(0), GetId(c()))
	expect_equal(character(0), GetId(NULL))
	expect_equal(NA_character_, GetId(NA))
})
# Load metadata into a dataframe
LoadMetadata <- function(directory) {
	write("Loading metadata", stderr())
	filenames <- list.files(directory, pattern="*.metadata$", full.names=TRUE)
	zeroSizeFiles = file.info(filenames)$size == 0
	if (any(zeroSizeFiles)) {
		warning(paste("Skipping files", filenames[zeroSizeFiles], "as they have 0 size."))
		filenames <- filenames[!zeroSizeFiles]
	}
	#metadata <- foreach (filename=filenames, .export=c("GetId"), .combine=rbind) %dopar% {
	metadata <- lapply(filenames, function(filename) {
		md <- read.csv(filename, header=FALSE, sep="=", quote = "\"'", col.names=c("CX", "V"))
		md$File <- filename
		md$Id <- GetId(filename)
		md
	})
	metadata <- do.call(rbind, metadata)
	metadata <- data.frame(lapply(metadata, as.character), stringsAsFactors=FALSE)
	metadata <- cast(metadata, File + Id ~ CX, value="V")  # pivot on context name
	rownames(metadata) <- metadata$Id
	# transform known numeric data to expected type
	metadata$CX_READ_FRAGMENT_LENGTH <- as.numeric(as.character(metadata$CX_READ_FRAGMENT_LENGTH))
	metadata$CX_READ_LENGTH <- as.numeric(as.character(metadata$CX_READ_LENGTH))
	metadata$CX_READ_DEPTH <- as.numeric(as.character(metadata$CX_READ_DEPTH))
	write(paste(nrow(metadata), "metadata files loaded"), stderr())
	return(metadata)
}
# infer proxy quality scores for ROC purposes based on strength of support
withqual <- function(vcf, caller) {
  if (is.null(rowRanges(vcf)$QUAL)) {
    rowRanges(vcf)$QUAL <- NA_real_
  }
  if (!is.na(caller) && !is.null(caller) && all(is.na(rowRanges(vcf)$QUAL))) {
    caller <- str_extract(caller, "^[^/]+") # strip version
    # use total read support as a qual proxy
    if (caller %in% c("delly")) {
      rowRanges(vcf)$QUAL <- ifelse(is.na(info(vcf)$PE), 0, info(vcf)$PE) + ifelse(is.na(info(vcf)$SR), 0, info(vcf)$SR)
    } else if (caller %in% c("crest")) {
      rowRanges(vcf)$QUAL <- ifelse(is.na(info(vcf)$right_softclipped_read_count), 0, info(vcf)$right_softclipped_read_count) + ifelse(is.na(info(vcf)$left_softclipped_read_count), 0, info(vcf)$left_softclipped_read_count)
    } else if (caller %in% c("pindel")) {
      rowRanges(vcf)$QUAL <- geno(vcf)$AD[,1,2]
    } else if (caller %in% c("lumpy")) {
      rowRanges(vcf)$QUAL <- unlist(info(vcf)$SU)
    }
  }
  if (any(is.na(rowRanges(vcf)$QUAL))) {
  	#if (is.null(caller) && is.na(caller)) {
	warning(paste("Missing QUAL scores for", caller))
  }
  return(vcf)
}
#' Loads a minimal structural variant GRanges from the VCF
LoadMinimalSVs <- function(filename, caller) {
	vcf <- readVcf(filename, "hg19")
	vcf <- withqual(vcf, caller)
	gr <- breakpointRanges(vcf)
	gr$paramRangeID <- NULL
	gr$REF <- NULL
	gr$ALT <- NULL
	gr$svtype <- NULL
	#gr$svLen <- NULL
	gr$insSeq <- NULL
	#gr$insLen <- NULL
	return(gr)
}
#vcf <- readVcf("C:/dev/sv_benchmark/data.aligner/5afa7ffdf2cc32602476526d5b477c5c.vcf", "hg19")
#' Loads structural variant GRanges from the VCFs in the given directory
LoadVCFs <- function(directory, pattern="*.vcf$", metadata=NULL, existingList=NULL) {
	write("Loading VCFs", stderr())
	filenames <- list.files(directory, pattern=pattern, full.names=TRUE)
	zeroSizeFiles = file.info(filenames)$size == 0
	if (any(zeroSizeFiles)) {
		write(paste("Skipping file", filenames[zeroSizeFiles], "due to 0 size.\n"))
		warning(paste("Skipping files", paste(filenames[zeroSizeFiles]), "due to 0 size.\n"))
		filenames <- filenames[!zeroSizeFiles]
	}
	# exclude already loaded VCFs
	filenames <- filenames[!(GetId(filenames) %in% names(existingList))]
	# only load VCFS that have metadata
	#if (!is.null(metadata)) {
	#	filenames <- filenames[GetId(filenames) %in% metadata$Id]
	#}
	#vcfs <- foreach (filename=filenames, .packages="VariantAnnotation") %dopar% { # Parallel load of VCFs
	grlist <- lapply(filenames, function(filename) {
		write(paste0("Loading ", filename), stderr())
		caller <- NULL
		if (!is.null(metadata)) {
			caller <- metadata$CX_CALLER[metadata$Id == GetId(filename)]
		}
		gr <- LoadMinimalSVs(filename, caller)
		return (gr)
	})
	names(grlist) <- GetId(filenames)
	grlist[sapply(grlist, is.null)] <- NULL # Remove NULL VCFs list
	write(paste("Loaded", length(grlist), "VCFs"), stderr())
	return(c(existingList, grlist))
}
