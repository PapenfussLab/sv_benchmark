# setwd("R")
source("libmanuscript_figures.R")

#generate_figures(datadir, sample_name, ids, truth_id, truth_name, grtransformName, allow_missing_callers=TRUE)

ignore.interchromosomal <- TRUE
mineventsize <- 51
maxeventsize <- NULL #TODO: what max event size should we use? Need one for chm since truth only has calls under a certain size
maxgap <- 200
sizemargin <- 0.25
ignore.strand <- TRUE
nominalPosition <- FALSE
longread_requiredHits <-5
longread_minMapq <- 0

# NA12878
na12878_truth <- c(
	"Parikh et al"="00000000000000000000000000000003",
	# TODO: which truth sets are these?
	"0001"="00000000000000000000000000000001", # Mills 2010?
	"0000"="00000000000000000000000000000000", # Kidd 2008?
	"0002"="00000000000000000000000000000002" # Sudmunt 2015?
		)
datadir <- "../data.na12878"
sample_name <- "NA12878"
ids <- c(
	#"2f4d15f5f6e428fb91e12dac971571b9", #bcftools/1.3.1
	"26511afa5055a939fdf99a2ac24938cc", #manta/1.1.1
	"6ae03359fcf0a39aa236a0aaea7ea915", #breakdancer/1.4.5
	"80dd0c2aa34964f330702ee3d0bfa53f", #hydra/master20160129
	"8aaf2886ffe782e666661d6b890b330a", #gridss/1.4.1
	"215312bf8994d7370151e8f9e6995ba5", #pindel/0.2.5b6 # broken"8d6a1ff660ab1a5374891f760903bdfb", #pindel/0.2.5b6
	"a6cbc4fc76337d871ef899d51fbae0d9", #socrates/1.13.1
	"b6aeda0d81ff8c839e381b57898dc3d8", #crest
	"e57987b7dec305f8be3a91fa95108c61", #cortex/1.0.5.14
	"f849b9493f2dd5dc20b6b7e49e1c89d7", #delly/0.7.6
	"fa5f9a52d3b67a206cb6e581053a9269") #lumpy/0.2.11
grtransformName <- "DAC"
longreadbedpedir <- paste0(datadir, "/../", "input.na12878/longread/")
for (i in seq_along(na12878_truth)) {
	truth_id <- na12878_truth[i]
	truth_name <- names(na12878_truth)[i]
	generate_figures(
		datadir, sample_name, ids, truth_id, truth_name,
		grtransformName,
		longreadbedpedir = NULL,
		#### LONG READS NOT WORKING ###
		# longreadbedpedir = longreadbedpedir,
		#### CORTEX NOT FINISHED ###
		allow_missing_callers = TRUE,
		eventtype = "DEL")
}

## chm datasets

longreadbedpedir <- NULL #TODO use Eichler PacBio reads directly

#chm13
datadir <- "../data.chm13"
sample_name <- "chm13"
ids <- c(
    "1bf3f864802712ffbe88f5bf29c7035a", # breakdancer/1.4.5
    "3b4917f2b3fc05710961bd98fe795863", # hydra/master20160129
    "48321dad78b7fdd801f6778c438799b3", # delly/0.7.6
    "4c7e0c2ae616fe15db24b2caf3786b93", # socrates/1.13.1
    "5ececd73193ffa02446ec5e2ea269fab", # pindel/0.2.5b6
    "6ab531d72c9b112e4a798e4285205cf3", # lumpy/0.2.11
    "907ff8ccd8657b8554ec31feb09f461b", # gridss/1.4.1
    "a6034953f1524610656f002b74faa1c4", # crest
    "b264ac30faf22d80ef51d8c4a25b1275"  # manta/1.1.1
    # "c639766990fcfbca2eb45f7806362fe6" # bcftools/1.3.1
	)
truth_id <- "00000000000000000000000000000013"
truth_name <- "Huddleston et al"
grtransformName <- "None"
generate_figures(
	datadir, sample_name, ids, truth_id, truth_name,
	grtransformName, allow_missing_callers = TRUE, eventtype = "DEL")

#chm1
datadir <- "../data.chm1"
sample_name <- "chm1"
ids <- c(
    "16c58fbcc5633564b10ebe8f78d87883", # manta/1.1.1
    "43a13d07730deb934e9fc01e3b3cd26f", # hydra/master20160129
    "74a311b09582df074b2bd839b06f8775", # pindel/0.2.5b6
    "8dcad8fe04f4ebc0ad3254ab4420cdc8", # breakdancer/1.4.5
    # "94bb6deef9f1bf1f9027a47e8488ae4f", # bcftools/1.3.1
    "9d134f160ac68c0445002fbb78db4a5e", # delly/0.7.6
    "acd889cc16741fb0fba62faa4f7005f3", # lumpy/0.2.11
    "b1112f1c3cbd28c464f58fc5c5c02f9b", # gridss/1.4.1
    "b142875fd36e402e4aff5f9e31f257f7")  # socrates/1.13.1

truth_id <- "00000000000000000000000000000001"
truth_name <- "Huddleston et al"
grtransformName <- "None"
generate_figures(
	datadir, sample_name, ids, truth_id, truth_name,
	grtransformName, allow_missing_callers = TRUE, eventtype = "DEL")

#chmboth
datadir <- "../data.chmboth"
sample_name <- "chmboth"
ids <- c(
	"21e805013e3d43cdb2123469a284154b", # breakdancer/1.4.5
	"c9e4d5d92fadf3165d86aaf28a7d57e7", # delly/0.7.6
	"97e18047b004eeafaa1f13465f2da01d", # gridss/1.4.1
	"267c637130c88d0ed28369aa25f4fc2b", # hydra/master20160129
	"b39d9b583e9bbb76e60e59182f25baa3", # lumpy/0.2.11
	"38e090223d3ce8084593742f95fe3571", # manta/1.1.1
	"0519837b3526f1bb8a92bd141bdc1d8b", # pindel/0.2.5b6
	"50c375a9aa5b71f643effb7c78d01504", # socrates/1.13.1
	"6a7cabe5c50fb7256eaf75f30a5a7fa1", # cortex/1.0.5.14
	"0fcd5cd6a5844ed4210db524ea841ca5"  # crest
	)

truth_id <- "00000000000000000000000000000014"
truth_name <- "Huddleston et al"
grtransformName <- "None"
generate_figures(
	datadir, sample_name, ids, truth_id, truth_name,
	grtransformName, allow_missing_callers = TRUE, eventtype = "DEL")

#chmboth - manually merged truth set
datadir <- "../data.chmboth"
sample_name <- "chmboth"
ids <- c(
	"21e805013e3d43cdb2123469a284154b", # breakdancer/1.4.5
	"c9e4d5d92fadf3165d86aaf28a7d57e7", # delly/0.7.6
	"97e18047b004eeafaa1f13465f2da01d", # gridss/1.4.1
	"267c637130c88d0ed28369aa25f4fc2b", # hydra/master20160129
	"b39d9b583e9bbb76e60e59182f25baa3", # lumpy/0.2.11
	"38e090223d3ce8084593742f95fe3571", # manta/1.1.1
	"0519837b3526f1bb8a92bd141bdc1d8b", # pindel/0.2.5b6
	"50c375a9aa5b71f643effb7c78d01504", # socrates/1.13.1
	"6a7cabe5c50fb7256eaf75f30a5a7fa1", # cortex/1.0.5.14
	"0fcd5cd6a5844ed4210db524ea841ca5"  # crest
)

truth_id <- "00000000000000000000000000000003"
truth_name <- "Manually merged truth set"
grtransformName <- "None"
generate_figures(
	datadir, sample_name, ids, truth_id, truth_name,
	grtransformName, allow_missing_callers = TRUE, eventtype = "DEL")

#HG002
datadir <- "../data.HG002"
sample_name <- "HG002"
ids <- c(
	"274d2a1693dde71a4d1817f6418c004f", # breakdancer/1.4.5
	"27fb72e39c6ba30b8e728e63deae4640", # hydra/master20160129
	"06b54896d6319ad0a3a25891d6d9bb03", # lumpy/0.2.11 # vcf empty
	"465e01c683f87eeab9f899d4ad7fe65f", # manta/1.1.1
	"476bb92e736e6842977cf0f6bfca72cc", # delly/0.7.6
	"4d3c8e2d0f784d80c7f73e534e0a7d22", # crest
	"5a52bfca501ce56c15771a96bac99fa0", # cortex/1.0.5.14
	"81bfafd1a366f2288882b93e3e0fd56e", # gridss/1.6.1-SNAPSHOT
	# "8e7ed64cfb252d7a39bfdb5b406b8d4f", # crest # broken - alt contigs only
	# "99ffc54b5253251c802fc3db2adf0e5c", # bcftools/1.3.1 # Presumably don't want/need this
	"eded64e93b9fa957f80bf040baec9da4"  # pindel/0.2.5b6
)

truth_id <- "00000000000000000000000000000002"
truth_name <- "GiaB NIST Tier 1 v0.6"
grtransformName <- "HG002_NIST_T1"
generate_figures(
	datadir, sample_name, ids, truth_id, truth_name,
	grtransformName, allow_missing_callers = TRUE, eventtype = "DEL")

