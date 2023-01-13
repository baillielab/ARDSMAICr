## overlap_test()

### test that subsetting n_genes works, universe size passes correctly, and a valid GeneOverlap obj
### is returned

test_that("overlap_test works correctly", {
  set.seed(1)
  data_genes <- data.frame(gene = sample(1:19220, 1000))
  overlap_genes <- data.frame(gene = sample(1:19220, 100))
  go_obj_res <- suppressWarnings(overlap_test(data_genes, overlap_genes, n_data_genes = 500, universe = "HGNC"))
  go <- attributes(go_obj_res)
  expect_true(nrow(data_genes) == 1000)
  expect_true(nrow(overlap_genes) == 100)
  expect_true(length(go$intersection) >= 1 && length(go$intersection) <= 100)
  expect_true(go$pval >= 0 && go$pval <= 1)
})

### test n_genes error

test_that("n_genes error correct", {
  expect_error(suppressWarnings(overlap_test(data_genes, data_biolitmine, n_data_genes = 10000)))
})

### test "FANTOM_L" works

test_that("FANTOM_L set correctly", {
  test <- suppressWarnings(overlap_test(data_genes, data_biolitmine, n_data_genes = "100", universe = "FANTOM_L"))
  test_res <- attributes(test)
  genome_size <- test_res$genome.size
  expect_equal(genome_size, 12606)
})

### test "FANTOM_S" works

test_that("FANTOM_L set correctly", {
  test <- suppressWarnings(overlap_test(data_genes, data_biolitmine, n_data_genes = "100", universe = "FANTOM_S"))
  test_res <- attributes(test)
  genome_size <- test_res$genome.size
  expect_equal(genome_size, 11512)
})

### test "HGNC" works

test_that("FANTOM_L set correctly", {
  test <- suppressWarnings(overlap_test(data_genes, data_biolitmine, n_data_genes = "100", universe = "HGNC"))
  test_res <- attributes(test)
  genome_size <- test_res$genome.size
  expect_equal(genome_size, 19220)
})