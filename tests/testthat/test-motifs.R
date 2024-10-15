pwm <- readRDS("../testdata/pwm_2motifs.rds")
genome.fasta <- system.file("extdata", "chr1_start.fa", package = "Signac")
genome <- Rsamtools::FaFile(genome.fasta)

test_that("AddMotifs works", {
  motif <- AddMotifs(
    atac_small[["peaks"]],
    genome,
    pwm,
    verbose = FALSE
  )
  expect_equal(dim(Motifs(motif)), c(323, 2))
})

test_that("AddMotifs works with fakechr", {
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
    features = features,
    keep_all_features = TRUE,
    verbose = FALSE
  )
  object <- CreateChromatinAssay(
    counts = mat,
    sep = c(":", "-"),
    fragments = fragments,
    verbose = FALSE
  )
  motif <- AddMotifs.default(
    granges(object)[100:nrow(object)],
    genome,
    pwm,
    verbose = FALSE
  )
  expect_equal(dim(motif), c(225, 2))
  test <- SetAssayData(
    object = object,
    layer = "motifs",
    new.data = motif,
    verbose = FALSE
  )
  expect_equal(dim(Motifs(test)), c(324, 2))
  # If some features are not in the ChromatinAssay
  # They should be removed
  test <- SetAssayData(
    object = atac_small,
    layer = "motifs",
    new.data = motif,
    verbose = FALSE
  )
  expect_equal(dim(Motifs(test)), c(323, 2))
})
