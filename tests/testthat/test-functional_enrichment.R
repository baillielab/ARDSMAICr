## gp_enrichment()

### test that P values are limited < 1

test_that("correct P value", {
  expect_error(suppressWarnings(gp_enrichment(data_genes, p_threshold = 1, source = "GO")))
})

### test that sources are correct

test_that("correct sources", {
  expect_error(suppressWarnings(gp_enrichment(data_genes, p_threshold = 0.1, source = "ABC")))
})

### check args are passed to gost()

test_that("correct sources", {
  result <- suppressWarnings(gp_enrichment(data_genes, p_threshold = 0.0001, source = "REAC"))
  pval <- result$meta$query_metadata$user_threshold
  sources <- result$meta$query_metadata$sources
  expect_equal(pval, 0.0001, tolerance = 1e-5)
  expect_equal(sources, "REAC")
})

### check func returns a list

test_that("correct sources", {
  result <- suppressWarnings(gp_enrichment(data_genes, p_threshold = 0.0001, source = "REAC"))
  expect_type(result, "list")
})

## gp_plot()

test_that("creates manhatten plot as expected", {
  result <- readRDS(test_path("testdata", "gp_result.rds"))
  vdiffr::expect_doppelganger(
    title = "create manhatten",
    fig = suppressWarnings(gp_plot(result)),
  )
})

## gp_tidy_results()

### test that func returns expected result

test_that("gp_tidy_results returns expected output", {
  gp_result <- list(result = data.frame(
    p_value = c(0.001, 0.05), term_size = c(10, 20),
    query_size = c(5, 15), intersection_size = c(3, 8),
    precision = c(0.6, 0.53), recall = c(0.3, 0.53),
    term_id = c("GO:0000001", "GO:0000002"),
    source = c("GO", "GO"),
    term_name = c("term 1", "term 2")
  ))
  expected_output <- data.frame(
    Rank = c("1", "2"), term_name = c("Term 2", "Term 1"),
    source = c("GO", "GO"), term_id = c("GO:0000002", "GO:0000001"),
    intersection_size = c(8, 3), term_size = c(20, 10),
    recall = c(0.53, 0.3), precision = c(0.53, 0.6),
    query_size = c(15, 5), p_value = c("5.00e-02", "1.00e-03")
  )
  expected_output$source <- as.factor(expected_output$source)
  expect_equal(suppressWarnings(gp_tidy_results(gp_result)), expected_output)
})

## gp_table()

### test that columns are correct

test_that("first author column is returned", {
  result <- suppressWarnings(gp_enrichment(data_genes, p_threshold = 0.00001, source = "GO"))
  tidy_result <- suppressWarnings(gp_tidy_results(result))
  table <- suppressWarnings(gp_table(tidy_result))
  col_list <- table$x$tag$attribs$columns
  rank <- col_list[[c(1,1)]]
  expect_equal(rank, "Rank")
})

## gp_count_pathways()

### check func returns correct type when by_pathway = TRUE

test_that("correct type wehn by_pathway = TRUE", {
gp_res <- readRDS(test_path("testdata", "gp_result.rds"))
count_gp <- suppressWarnings(gp_count_pathways(gp_res, by_pathway = TRUE))
expect_type(count_gp, "list")
})

### check func returns correct type when by_pathway = FALSE

test_that("correct type wehn by_pathway = TRUE", {
  gp_res <- readRDS(test_path("testdata", "gp_result.rds"))
  count_gp <- suppressWarnings(gp_count_pathways(gp_res))
  expect_type(count_gp, "integer")
})