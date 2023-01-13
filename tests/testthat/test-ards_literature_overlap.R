## overlap_table()

### test biolitmine = TRUE behaviour if not correct format

test_that("check biolitmine = TRUE is correct", {
  data <- readRDS(test_path("testdata", "biolitmine_fail.rds"))
  expect_error(suppressWarnings(overlap_table(data_genes, data, biolitlimine = TRUE)))
})

### test that columns are correct biolitmine = TRUE

test_that("check that gene column is returned", {
  data <- data_biolitmine[1:50,]
  table <- suppressWarnings(overlap_table(data_genes, data, biolitmine = TRUE))
  col_list <- table$x$tag$attribs$columns
  first_author <- col_list[[c(1,1)]]
  expect_equal(first_author, "gene")
})

### ### test that columns are correct biolitmine = FALSE

test_that("check that gene column is returned", {
  data <- readRDS(test_path("testdata", "biolitmine_fail.rds"))
  table <- suppressWarnings(overlap_table(data_genes, data, biolitmine = FALSE))
  col_list <- table$x$tag$attribs$columns
  first_author <- col_list[[c(1,1)]]
  expect_equal(first_author, "gene")
})

## overlap_unidentified_table()

### test biolitmine = TRUE behaviour

test_that("check biolitmine = TRUE is correct", {
  data <- readRDS(test_path("testdata", "biolitmine_fail.rds"))
  expect_error(suppressWarnings(overlap_unidentified_table(data_genes, data, biolitlimine = TRUE)))
})

### test that columns are correct biolitmine = TRUE

test_that("check that gene column is returned", {
  data <- data_biolitmine[1:50,]
  table <- suppressWarnings(overlap_unidentified_table(data_genes, data, biolitmine = TRUE))
  col_list <- table$x$tag$attribs$columns
  first_author <- col_list[[c(1,1)]]
  expect_equal(first_author, "Gene")
})

### test that if biolitmine = FALSE a tibble is returned

test_that("tibble if biolitmine = FALSE", {
  data <- data_biolitmine[1:50,]
  table <- suppressWarnings(overlap_unidentified_table(data_genes, data, biolitmine = FALSE))
  expect_type(table, "list")
})

## overlap_venn()

### test if biolitmine = TRUE

test_that("creates venn as expected biolitmine = TRUE", {
  data <- data_biolitmine[1:50,]
  vdiffr::expect_doppelganger(
    title = "create upset TRUE ",
    fig = suppressWarnings(overlap_venn(data_genes, data, biolitmine = TRUE)),
  )
})

### test if biolitmine = FALSE

test_that("creates venn as expected biolitmine = FALSE", {
  data <- readRDS(test_path("testdata", "biolitmine_fail.rds"))
  vdiffr::expect_doppelganger(
    title = "create upset FALSE",
    fig = suppressWarnings(overlap_venn(data_genes, data, biolitmine = FALSE)),
  )
})

### test biolitmine = TRUE behaviour

test_that("check biolitmine = TRUE is correct", {
  data <- readRDS(test_path("testdata", "biolitmine_fail.rds"))
  expect_error(suppressWarnings(overlap_venn(data_genes, data, biolitlimine = TRUE)))
})

## overlap_count()

### test if the function returns the correct values for a given input

test_that("overlap_count returns the correct values", {
  data_genes <- tibble::tibble(gene = c("A", "B", "C", "D"))
  data_alternative <- tibble::tibble(gene = c("A", "B", "E", "F"))
  result <- suppressWarnings(overlap_count(data_genes, data_alternative))
  expect_equal(result, list(alternative_n = 4, maic_overlap_n = 2, maic_overlap_percentage = 50))
})

### test if the function returns the correct values when there are no common genes

test_that("overlap_count returns the correct values when no common genes", {
  data_genes <- tibble::tibble(gene = c("A", "B", "C", "D"))
  data_alternative <- tibble::tibble(gene = c("E", "F", "G", "H"))
  result <- overlap_count(data_genes, data_alternative)
  expect_equal(result, list(alternative_n = 4, maic_overlap_n = 0, maic_overlap_percentage = 0))
})

