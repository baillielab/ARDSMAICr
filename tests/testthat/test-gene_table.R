## gene_table()

### check unsupported methods error

test_that("error for unsupported methods", {
  data <- readRDS(test_path("testdata", "gene_table_methods.rds"))
  expect_error(suppressWarnings(gene_table(data)))
})

### test that columns are correct

test_that("check that rank column is returned", {
  data <- data_genes[1:50,]
  table <- suppressWarnings(gene_table(data))
  col_list <- table$x$tag$attribs$columns
  first_author <- col_list[[c(1,1)]]
  expect_equal(first_author, "rank")
})
