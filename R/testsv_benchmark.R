source("sv_benchmark.R")

test_that("GetId", {
  expect_equal(c("a", "m", "bin"),
               GetId(c("a", "/dir/m.as.ds.dsa.ss.vcf", "/usr/local/bin")))
  expect_equal(character(0), GetId(c()))
  expect_equal(character(0), GetId(NULL))
  expect_equal(NA_character_, GetId(NA))
  if (as.character(Sys.info())[1] == "Windows") {
    expect_equal(c("m"), GetId(c("C:\\directory\\m.as.ds.dsa.ss.vcf")))
  }
})
test_that(".distance", {
  expect_equal(data.frame(
    min=c(1, 1, 0),
    max=c(1, 20, 5)
  ), .distance(IRanges(start=c(1, 1, 1), end=c(1, 10, 10)),
               IRanges(start=c(2, 11, 5), end=c(2, 21, 5))))
  
  expect_equal(data.frame(
    min=c(1, 1, 0),
    max=c(1, 20, 5)
  ), .distance(IRanges(start=c(2, 11, 5), end=c(2, 21, 5)),
               IRanges(start=c(1, 1, 1), end=c(1, 10, 10))))
})


test_that("findSpanningBreakpoints", {
  query <- GRanges(seqnames=c("chr1", "chr1", "chr1", "chr1"),
                   ranges=IRanges(start=c(1000, 2000, 2200, 3000), width=1),
                   strand=c("-", "+", "-", "+"),
                   partner=c("a2", "a1", "b2", "b1"))
  names(query) <- c("a1", "a2", "b1", "b2")
  subject <- GRanges(seqnames=c("chr1", "chr1"),
                     ranges=IRanges(start=c(1000, 3000), width=1),
                     strand=c("-", "+"),
                     partner=c("c2", "c1"))
  names(subject) <- c("c1", "c2")
  spanning <- findSpanningBreakpoints(query, subject, maxSpanningFragmentSize=300)
  expect_equal(spanning$subjectHits, c(1, 2))
  expect_equal(spanning$localBreakend, c("a1", "b2"))
  expect_equal(spanning$remoteBreakend, c("b2", "a1"))
  
  # intermediate fragment is too big
  expect_equal(nrow(findSpanningBreakpoints(query, subject, maxSpanningFragmentSize=198)), 0)
})



