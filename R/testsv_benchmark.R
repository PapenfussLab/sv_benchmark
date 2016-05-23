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

