library(ggplot2)
library(rtracklayer)

theme_set(theme_bw())

scale_y_power4 <- scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)^4, labels=c("0", "", "", "", "", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"))
scale_y_power5 <- scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)^5, labels=c("0", "", "", "", "", "0.5", "", "0.7", "0.8", "0.9", "1.0"))

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
scale_x_log_fp <- scale_x_continuous(breaks=c(1, 11, 101, 1001, 10001, 100001),
	labels=c("0", "10", "100", "1k", "10k", "100k"),
	minor_breaks=c(),
	trans="log10")

saveplot <- function(file=file, ...) {
  ggsave(paste0("png/", file, ".png"), ...)
  ggsave(paste0("eps/", file, ".eps"), ...)
}
