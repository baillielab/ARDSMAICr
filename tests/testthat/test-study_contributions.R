## contributions_calculation()

### test that the output has the correct number of studies

test_that("test correct number of studies", {
  data <- readRDS(test_path("testdata", "sample_genes.rds"))
  data <- as.data.frame(data)
  output <- suppressWarnings(contributions_calculation(data))
  data_rows <- unlist(lapply(data, is.numeric))
  data_rows_df <- data[ , data_rows]
  data_rows_df <- data_rows_df[ , !names(data_rows_df) %in% c("maic_score")]
  expected_rows <- ncol(data_rows_df)
  expect_equal(nrow(output), expected_rows)
})

### test that the output has the correct number of columns

test_that("test correct number of columns", {
  data <- readRDS(test_path("testdata", "sample_genes.rds"))
  output <- suppressWarnings(contributions_calculation(data))
  expect_equal(ncol(output), 5)
})


