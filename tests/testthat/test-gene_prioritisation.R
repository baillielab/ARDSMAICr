## inflection_point()

### test that correct inflection point is returned

test_that("Data set with non-numeric rownames returns correct inflection point", {
  data <- readRDS(test_path("testdata", "inflection_test.rds"))
  expected_output <- list(maic_score = 1.84258, gene_number = 1308)
  result <- suppressWarnings(inflection_point(data_genes))
  expect_equal(result, expected_output, tolerance = 1e-5)
})

## inflection_point_plot

test_that("creates inflection plot as expected", {
  sample_genes <- readRDS(test_path("testdata", "large_sample_genes.rds"))
  vdiffr::expect_doppelganger(
    title = "create plot",
    fig = suppressWarnings(inflection_point_plot(sample_genes, first_break = 100, increment = 50)),
  )
})