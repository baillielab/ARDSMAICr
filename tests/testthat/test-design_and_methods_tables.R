## design_table()

### test that duplicate PMIDS are removed

test_that("the function subsets the table correctly", {
  data <- data.frame(
    First_author = c("Author 1", "Author 2", "Author 3"),
    PMID = c(1234567, 7654321, 1234567),
    Year = c(2020, 2019, 2018),
    Focus = c("ARDS", "ARDS", "ARDS"),
    ARDS_definition = c("Berlin definition", "AECC definition", "Berlin definition"),
    ARDS_pts = c(100, 200, 50)
  )
  table <- suppressWarnings(design_table(data))
  auth <- table$x$tag$attribs$data
  dup_removed <- stringr::str_detect(auth, c("Author 3"))
  expect_false(dup_removed)
})

## methods_table()

### test for duplicates error

test_that("entries are unique", {
  test_data <- readRDS(test_path("testdata", "duplicate_methods_study.rds"))
  expect_error(suppressWarnings(methods_table(test_data)), "Duplicates detected...if there are none please qualify Technologies")
})

### test that columns are correct

test_that("first author column is returned", {
  table <- suppressWarnings(methods_table(data_study))
  col_list <- table$x$tag$attribs$columns
  first_author <- col_list[[c(1,1)]]
  expect_equal(first_author, "First_author")
})