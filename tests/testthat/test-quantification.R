test_that("FeatureMatrix works", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  fragments <- CreateFragmentObject(fpath)
  mat <- FeatureMatrix(
    fragments = fragments,
    features = granges(atac_small)
  )
  expect_equal(dim(mat), c(323, 50))
})
test_that("FeatureMatrix works on grange on diff seqnames", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  fragments <- CreateFragmentObject(fpath)
  features <- suppressWarnings(c(
    granges(atac_small),
    GenomicRanges::GRanges(
      seqnames = "fake_chr",
      ranges = IRanges::IRanges(start = 1, end = 1000)
    )
  ))
  mat <- FeatureMatrix(
    fragments = fragments,
    features = features
  )
  expect_equal(dim(mat), c(323, 50))
  mat <- FeatureMatrix(
    fragments = fragments,
    features = features,
    keep_all_features = TRUE
  )
  expect_equal(dim(mat), c(324, 50))
})