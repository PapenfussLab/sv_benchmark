# Importing metadata

metadata_hg002_files <- dir("../data.HG002/") %>% str_subset(".*metadata") %>% str_c("../data.HG002/", .)

read_metadata_file <- function(f) {
	read_delim(f, delim = "=", col_names = c("variable", "value")) %>%
	mutate(filename = f,
		   Id = str_remove(basename(f), ".metadata"),
		   corresponding_vcf = str_c(str_remove(f, ".metadata"), ".vcf"),
		   corresponding_vcf_exists = file.exists(corresponding_vcf)) %>%
	spread(key = variable, value = value)
	}

metadata_hg002 <- map(metadata_hg002_files, read_metadata_file) %>% bind_rows()

caller_versions_na12878 <- c(
	"manta/1.1.1",
	"breakdancer/1.4.5",
	"hydra/master20160129",
	"gridss/1.4.1",
	"pindel/0.2.5b6",
	"socrates/1.13.1",
	"crest",
	"cortex/1.0.5.14",
	"delly/0.7.6",
	"lumpy/0.2.11"
)

caller_versions_chm13 <- c(
	"breakdancer/1.4.5",
	"hydra/master20160129",
	"delly/0.7.6",
	"socrates/1.13.1",
	"pindel/0.2.5b6",
	"lumpy/0.2.11",
	"gridss/1.4.1",
	"crest",
	"manta/1.1.1"
)

caller_versions_chm1 <- c(
	"manta/1.1.1",
	"hydra/master20160129",
	"pindel/0.2.5b6",
	"breakdancer/1.4.5",
	"delly/0.7.6",
	"lumpy/0.2.11",
	"gridss/1.4.1",
	"socrates/1.13.1"
)

caller_versions_chmboth <- c(
	"breakdancer/1.4.5",
	"delly/0.7.6",
	"gridss/1.4.1",
	"hydra/master20160129",
	"lumpy/0.2.11",
	"manta/1.1.1",
	"pindel/0.2.5b6",
	"socrates/1.13.1"
)

all(caller_versions_chm13 %in% caller_versions_na12878)
all(caller_versions_chm1 %in% caller_versions_na12878)
all(caller_versions_chmboth %in% caller_versions_na12878)
# so can take na12878 as the master list

caller_versions_na12878[!(caller_versions_na12878 %in% (metadata_hg002 %>% filter(corresponding_vcf_exists))$CX_CALLER)]

caller_names_na12878 <- str_extract(caller_versions_na12878, "^[^/]+")

caller_versions_hg002

caller_names_na12878[!(caller_names_na12878 %in% (
	metadata_hg002 %>%
		filter(corresponding_vcf_exists) %>%
		mutate(caller=str_extract(CX_CALLER, "^[^/]+")) %>%
		pull(caller)))]

caller_versions_na12878
caller_versions_hg002 <-
	unique((
		metadata_hg002 %>%
			filter(corresponding_vcf_exists)$CX_CALLER[
				str_extract((metadata_hg002 %>%
							 	filter(corresponding_vcf_exists))$CX_CALLER, "^[^/]+")
				%in% caller_names_na12878]))

metadata_hg002 %>%
	filter(corresponding_vcf_exists) %>%
	mutate(caller = str_extract(CX_CALLER, "^[^/]+")) %>%
	select(caller, CX_READ_LENGTH, CX_REFERENCE, corresponding_vcf_exists) %>%
	arrange(caller) %>% View()

read_lengths <- unique(metadata_hg002$CX_READ_LENGTH)

cross_df(
	list(caller = caller_names_na12878,
		 CX_READ_LENGTH = read_lengths)) %>%
	drop_na() %>% arrange(caller) %>%
	left_join(
		metadata_hg002 %>%
			mutate(caller = str_extract(CX_CALLER, "^[^/]+"))) %>%
	select(caller, CX_READ_LENGTH, CX_READ_DEPTH, corresponding_vcf_exists, CX_REFERENCE, Id) %>%
	arrange(CX_READ_LENGTH, CX_READ_DEPTH)

# Using this for HG002:

metadata_hg002 %>%
	filter(CX_READ_LENGTH == 150,
		   CX_READ_DEPTH == 60,
		   corresponding_vcf_exists,
		   str_detect(CX_BAM, "cameron.d")) %>%
	select(Id, CX_CALLER)
