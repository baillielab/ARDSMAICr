## count_lists()

### test unavailable lists are filtered out

test_that("available lists are filtered", {
  test_data <- readRDS(test_path("testdata", "count_lists_available_equal.rds"))
  expect_equal(suppressWarnings(count_lists(test_data)), 2)
})

### test behaviour when all lists are unavailable

test_that("available lists are filtered", {
  test_data <- readRDS(test_path("testdata", "count_lists_all_false.rds"))
  expect_equal(suppressWarnings(count_lists(test_data)), 0)
})

### test that error is thrown when duplicate lists are found

test_that("uIDs are unique", {
  test_data <- readRDS(test_path("testdata", "count_lists_duplicate_uID.rds"))
  expect_error(suppressWarnings(count_lists(test_data)), "Duplicate lists detected...")
})

## lists_per_method()

### test unavailable lists are filtered out

test_that("available lists are filtered", {
  test_data <- readRDS(test_path("testdata", "count_lists_available_equal.rds"))
  test_data_2 <- readRDS(test_path("testdata", "lists_per_method_tibble.rds"))
  expect_identical(suppressWarnings(lists_per_method(test_data)), test_data_2)
})

### test error is thrown when a method entry is NA

test_that("method NAs", {
  test_data <- readRDS(test_path("testdata", "missing_methods.rds"))
  expect_error(suppressWarnings(lists_per_method(test_data)))
})

## count_genes()

### test behaviour when there is a single unique gene

test_that("single gene", {
  test_data <- readRDS(test_path("testdata", "same_genes.rds"))
  expect_equal(suppressWarnings(count_genes(test_data)), 1)
})

### test case sensitivity

test_that("case sensitivity", {
  test_data <- readRDS(test_path("testdata", "lowercase_genes.rds"))
  expect_equal(suppressWarnings(count_genes(test_data)), 1)
})

## genes_per_list()

### test list genes sum to total genes

test_that("genes sum correctly", {
  test_data <- readRDS(test_path("testdata", "sample_genes.rds"))
  test_data <- as.data.frame(test_data)
  list_genes <- suppressWarnings(genes_per_list(test_data))
  list_genes_sum <- sum(list_genes$n_genes)
  num_cols <- unlist(lapply(test_data, is.numeric))
  num_cols_df <- test_data[ , num_cols]
  num_cols_df <- num_cols_df[ , !names(num_cols_df) %in% c("maic_score")]
  col_sums <- colSums(num_cols_df  != 0)
  test_total <- sum(col_sums)
  expect_equal(list_genes_sum, test_total)
})

## genes_per_method()

### test list genes sum to total genes

test_that("genes sum correctly", {
  test_data_genes <- readRDS(test_path("testdata", "per_method_genes.rds"))
  test_data_study<- readRDS(test_path("testdata", "per_method_study.rds"))
  test_data_genes <- as.data.frame(test_data_genes)
  list_genes <- suppressWarnings(genes_per_method(test_data_study, test_data_genes))
  list_genes_sum <- sum(list_genes$n_genes)
  num_cols <- unlist(lapply(test_data_genes, is.numeric))
  num_cols_df <- test_data_genes[ , num_cols]
  num_cols_df <- num_cols_df[ , !names(num_cols_df) %in% c("maic_score")]
  col_sums <- colSums(num_cols_df  != 0)
  test_total <- sum(col_sums)
  expect_equal(list_genes_sum, test_total)
})

## lists_per_gene()

### test lists sum to total lists

test_that("lists sum correctly", {
  test_data <- readRDS(test_path("testdata", "sample_genes.rds"))
  test_data <- as.data.frame(test_data)
  lists <- suppressWarnings(lists_per_gene(test_data))
  lists_sum <- sum(lists$n_lists)
  num_cols <- unlist(lapply(test_data, is.numeric))
  num_cols_df <- test_data[ , num_cols]
  num_cols_df <- num_cols_df[ , !names(num_cols_df) %in% c("maic_score")]
  col_sums <- colSums(num_cols_df  != 0)
  test_total <- sum(col_sums)
  expect_equal(lists_sum, test_total)
})

## lists_per_method_per_gene()

### test lists sum to total lists

test_that("lists sum correctly", {
  test_data_genes <- readRDS(test_path("testdata", "per_method_genes.rds"))
  test_data_study<- readRDS(test_path("testdata", "per_method_study.rds"))
  test_data_genes <- as.data.frame(test_data_genes)
  lists <- suppressWarnings(lists_per_method_per_gene(test_data_study, test_data_genes))
  lists_sum <- sum(lists$n_lists)
  num_cols <- unlist(lapply(test_data_genes, is.numeric))
  num_cols_df <- test_data_genes[ , num_cols]
  num_cols_df <- num_cols_df[ , !names(num_cols_df) %in% c("maic_score")]
  col_sums <- colSums(num_cols_df  != 0)
  test_total <- sum(col_sums)
  expect_equal(lists_sum, test_total)
})

## methods_per_gene()

### test methods sum to total methods

test_that("methods sum correctly", {
  test_data <- readRDS(test_path("testdata", "sample_genes.rds"))
  test_data <- as.data.frame(test_data)
  methods <- suppressWarnings(methods_per_gene(test_data))
  methods_sum <- sum(methods$n_methods)
  expect_equal(methods_sum, 43)
})




