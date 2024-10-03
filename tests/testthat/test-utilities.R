library(GenomicRanges)

test_that("AverageCounts works", {
  expect_equal(
    object = as.vector(x = AverageCounts(object = atac_small)),
    expected = c(162.80, 118.74),
    tolerance = 1 / 1000
  )
})

test_that("CellsPerGroup works", {
  expect_equal(
    object = as.vector(x = CellsPerGroup(object = atac_small)),
    expected = c(50, 50, 0)
  )
})

test_that("SortIdents works",{
  set.seed(1)
  atac_small$test <- sample(1:10, ncol(atac_small), replace = TRUE)
  atac_small <- SortIdents(object = atac_small, label = "test")
  expect_equal(
    object = levels(atac_small$test),
    expected = c("10", "6", "3", "7", "9", "5", "4", "8", "1", "2")
  )
  Idents(atac_small) <- sample(1:10, ncol(atac_small), replace = TRUE)
  atac_small <- SortIdents(object = atac_small)
  expect_equal(
    object = levels(Idents(atac_small)),
    expected = c("1", "9", "10", "2", "3", "5", "4", "7", "6", "8")
  )
})

test_that("GRanges conversion works", {
  correct_string <- c("chr1-1-10", "chr2-12-3121")
  correct_ranges <- GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges(start = c(1, 12), end = c(10, 3121))
  )
  granges <- StringToGRanges(regions = correct_string)
  string_ranges <- GRangesToString(grange = correct_ranges)
  expect_equal(object = granges, expected = correct_ranges)
  expect_equal(object = string_ranges, expected = correct_string)
})

test_that("ChunkGRanges works", {
  granges <- GRanges(
    seqnames = c("chr1"),
    ranges = IRanges(start = seq(1, 10), end = seq(1, 10) + 1)
  )
  split_ranges <- ChunkGRanges(granges = granges, nchunk = 3)
  correct_split <- list(
    GRanges(seqnames = "chr1", ranges = IRanges(start = 1:3, end = (1:3) + 1)),
    GRanges(seqnames = "chr1", ranges = IRanges(start = 4:6, end = (4:6) + 1)),
    GRanges(seqnames = "chr1", ranges = IRanges(start = 7:10, end = (7:10) + 1))
  )
  expect_equal(object = split_ranges, expected = correct_split)
})

test_that("Extend works", {
  granges <- GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges(
      start = c(1000, 1200),
      end = c(10000, 312100)
    ),
    strand = c("+", "-")
  )
  correct_extended <- GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges(
      start = c(500, 200),
      end = c(11000, 312600)
    ),
    strand = c("+", "-")
  )
  extended_granges <- Extend(x = granges, upstream = 500, downstream = 1000)
  expect_equal(object = extended_granges, expected = correct_extended)
})

test_that("ExtractCell works", {
  expect_equal(
    object = ExtractCell(x = "chr1\t1\t300\tTGCA\t1"), expected = "TGCA"
  )
})

test_that("CreateBWGroup works with single tile", {
  outdir <- file.path(tempdir(), "createBW")
  dir.create(outdir, showWarnings = FALSE)
  fake.bed.data <- data.frame(
    seqnames = rep("chr1", 5),
    start = c(0, 10, 100, 110, 300),
    end = c(100, 150, 200, 250, 500),
    cell_name = rep("fake_cell", 5),
    nb = 1:5
  )
  write.table(
    fake.bed.data, file.path(outdir, "0.bed"),
    col.names = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE
  )
  CreateBWGroup(
    groupNamei = "0",
    availableChr = "chr1",
    chromLengths = c("chr1" = 249250621),
    tiles = GRanges(seqnames = "chr1", ranges = IRanges(start = 1, end = 249250621)),
    normBy = NULL,
    tileSize = 249250621,
    normMethod = 'RC',
    cutoff = NULL,
    outdir = outdir
  )
  expect_equal(object = length(list.files(outdir)), expected = 2)
  expect(file.exists(file.path(outdir, "0-TileSize-249250621-normMethod-rc.bw")), "File does not exist.")
  bw <- rtracklayer::import.bw(file.path(outdir, "0-TileSize-249250621-normMethod-rc.bw"))
  expect_equal(object = bw$score, 20000)
})

test_that("CreateBWGroup works with 100bp tile", {
  outdir <- file.path(tempdir(), "createBW2")
  dir.create(outdir, showWarnings = FALSE)
  fake.bed.data <- data.frame(
    seqnames = rep("chr1", 5),
    start = c(0, 10, 100, 110, 300),
    end = c(100, 150, 200, 250, 500),
    cell_name = rep("fake_cell", 5),
    nb = 1:5
  )
  write.table(
    fake.bed.data, file.path(outdir, "0.bed"),
    col.names = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE
  )
  CreateBWGroup(
    groupNamei = "0",
    availableChr = "chr1",
    chromLengths = c("chr1" = 249250621),
    tiles = GRanges(seqnames = "chr1", ranges = IRanges(start = seq(1, 249250621, 100), end = c(seq(100, 249250621, 100), 249250621))),
    normBy = NULL,
    tileSize = 100,
    normMethod = "RC",
    cutoff = NULL,
    outdir = outdir
  )
  expect_equal(object = length(list.files(outdir)), expected = 2)
  expect(file.exists(file.path(outdir, "0-TileSize-100-normMethod-rc.bw")), "File does not exist.")
  bw <- rtracklayer::import.bw(file.path(outdir, "0-TileSize-100-normMethod-rc.bw"))
  expect_equal(object = bw$score, c(6000, 8000, 2000, 0))
})

test_that("CreateBWGroup works with seqlength equal to final pos", {
  outdir <- file.path(tempdir(), "createBW2")
  dir.create(outdir, showWarnings = FALSE)
  fake.bed.data <- data.frame(
    seqnames = rep("chr1", 5),
    start = c(0, 10, 100, 110, 300),
    end = c(100, 150, 200, 250, 500),
    cell_name = rep("fake_cell", 5),
    nb = 1:5
  )
  write.table(
    fake.bed.data, file.path(outdir, "0.bed"),
    col.names = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE
  )
  CreateBWGroup(
    groupNamei = "0",
    availableChr = "chr1",
    chromLengths = c("chr1" = 500),
    tiles = GRanges(seqnames = "chr1", ranges = IRanges(start = seq(1, 499, 100), end = c(seq(100, 500, 100)))),
    normBy = NULL,
    tileSize = 100,
    normMethod = "RC",
    cutoff = NULL,
    outdir = outdir
  )
  expect_equal(object = length(list.files(outdir)), expected = 2)
  expect(file.exists(file.path(outdir, "0-TileSize-100-normMethod-rc.bw")), "File does not exist.")
  bw <- rtracklayer::import.bw(file.path(outdir, "0-TileSize-100-normMethod-rc.bw"))
  expect_equal(object = bw$score, c(6000, 8000, 2000))
})

test_that("ExportGroupBW works", {
  outdir <- file.path(tempdir(), "ExportGroupBW")
  dir.create(outdir, showWarnings = FALSE)
  fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
  cells <- colnames(x = atac_small)
  names(x = cells) <- cells
  frags <- CreateFragmentObject(
    path = fpath,
    cells = cells,
    verbose = FALSE,
    validate = FALSE
  )
  Fragments(atac_small) <- frags
  ExportGroupBW(
    object = atac_small,
    assay = NULL,
    group.by = NULL,
    idents = NULL,
    normMethod = "RC",
    tileSize = 100,
    minCells = 5,
    cutoff = NULL,
    chromosome = NULL,
    outdir = outdir,
    verbose = TRUE
  )
  expect_equal(object = length(list.files(outdir)), expected = 4)
  expect(
    file.exists(file.path(outdir, "0-TileSize-100-normMethod-rc.bw")),
    "File does not exist."
  )
  bw <- rtracklayer::import.bw(file.path(outdir, "0-TileSize-100-normMethod-rc.bw"))
  expect_equal(object = length(seqlengths(bw)), 298)
})