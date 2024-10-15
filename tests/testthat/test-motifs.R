test_that("AddMotifs works", {
  library(JASPAR2022)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg19)

  pwm <- getMatrixSet(
    x = JASPAR2022,
    opts = list(species = 9606, all_versions = FALSE)
  )
  motif <- AddMotifs(
    atac_small[["peaks"]],
    BSgenome.Hsapiens.UCSC.hg19,
    pwm
  )
  expect_equal(dim(Motifs(motif)), c(323, 692))
})

test_that("AddMotifs works with fakechr", {
  library(JASPAR2022)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg19)
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
    keep_all_features = TRUE
  )
  object <- CreateChromatinAssay(
    counts = mat,
    sep = c(":", "-"),
    fragments = fragments
  )

  pwm <- getMatrixSet(
    x = JASPAR2022,
    opts = list(species = 9606, all_versions = FALSE)
  )
  motif <- AddMotifs.default(
    granges(object)[100:nrow(object)],
    BSgenome.Hsapiens.UCSC.hg19,
    pwm
  )
  expect_equal(dim(motif), c(225, 692))
  test <- SetAssayData(
    object = object,
    layer = 'motifs',
    new.data = motif
  )
  expect_equal(dim(Motifs(test)), c(324, 692))
  # If some features are not in the ChromatinAssay
  # They should be removed
  test <- SetAssayData(
    object = atac_small,
    layer = "motifs",
    new.data = motif
  )
  expect_equal(dim(Motifs(test)), c(323, 692))
})
