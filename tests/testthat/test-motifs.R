pwm <- TFBSTools::PFMatrixList(
  "MA0030.1" = TFBSTools::PFMatrix(
    ID = "MA0030.1",
    name = "FOXF2",
    matrixClass = "Fork head/winged helix factors", strand = "+",
    bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
    tags = list(
      family = "FOX", species = "9606",
      tax_group = "vertebrates", medline = "7957066",
      tfbs_shape_id = "37", remap_tf_name = "FOXF2"
    ),
    profileMatrix = matrix(
      c(
        1L, 10L, 17L, 13L, 3L, 7L, 0L, 27L, 27L, 27L, 0L, 27L, 16L, 7L,
        10L, 7L, 4L, 5L, 11L, 0L, 0L, 0L, 0L, 0L, 25L, 0L, 4L, 4L,
        7L, 5L, 2L, 5L, 8L, 20L, 0L, 0L, 0L, 0L, 0L, 0L, 2L, 6L,
        9L, 5L, 4L, 5L, 0L, 0L, 27L, 0L, 0L, 0L, 2L, 0L, 5L, 10L
      ),
      byrow = TRUE, nrow = 4,
      dimnames = list(c("A", "C", "G", "T"))
    )
  ),
  "MA0031.1" = TFBSTools::PFMatrix(
    ID = "MA0031.1",
    name = "FOXD1",
    matrixClass = "Fork head/winged helix factors", strand = "+",
    bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
    tags = list(
      family = "FOX", species = "9606",
      tax_group = "vertebrates", medline = "7957066",
      tfbs_shape_id = "38", remap_tf_name = "FOXD1"
    ),
    profileMatrix = matrix(
      c(
        1L, 0L, 19L, 20L, 18L, 1L, 20L, 7L,
        1L, 0L, 1L, 0L, 1L, 18L, 0L, 2L,
        17L, 0L, 0L, 0L, 1L, 0L, 0L, 3L,
        1L, 20L, 0L, 0L, 0L, 1L, 0L, 8L
      ),
      byrow = TRUE, nrow = 4,
      dimnames = list(c("A", "C", "G", "T"))
    )
  )
)

test_that("AddMotifs works", {
  library(BSgenome.Hsapiens.UCSC.hg19)
  motif <- AddMotifs(
    atac_small[["peaks"]],
    BSgenome.Hsapiens.UCSC.hg19,
    pwm,
    verbose = FALSE
  )
  expect_equal(dim(Motifs(motif)), c(323, 2))
})

test_that("AddMotifs works with fakechr", {
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
    BSgenome.Hsapiens.UCSC.hg19,
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
