## methods_upset()

test_that("creates upset plot as expected", {
  sample_genes <- readRDS(test_path("testdata", "sample_genes.rds"))
  vdiffr::expect_doppelganger(
    title = "create upset",
    fig = suppressWarnings(methods_upset(sample_genes)),
  )
})
