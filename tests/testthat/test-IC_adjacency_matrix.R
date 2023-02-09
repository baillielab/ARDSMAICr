## inf_adjacency_matrix()

### test that default function returns expected values

test_that("Correct shared information content returned", {
  data <- readRDS(test_path("testdata", "adj_check.rds"))
  result <- suppressWarnings(inf_adj_matrix(data))
  result_A <- result[1, 1]
  result_B <- result[3, 2]
  result_C <- result[1, 3]
  result_D <- result[2, 4]
  expect_equal(result_A, 0)
  expect_equal(result_B, 2.5)
  expect_equal(result_C, 2)
  expect_equal(result_D, 2.5)
})

### test that expected contributions are returned when contribution == TRUE

test_that("Correct shared information contribution returned", {
  data <- readRDS(test_path("testdata", "adj_check.rds"))
  result <- suppressWarnings(inf_adj_matrix(data, contribution = TRUE))
  result_C <- result[4, 3]
  expect_equal(result_C, 4)
})

### test that expected contents are returned when unique == FALSE

test_that("Correct shared information content without unique genes returned", {
  data <- readRDS(test_path("testdata", "adj_check.rds"))
  result <- suppressWarnings(inf_adj_matrix(data, unique = FALSE))
  result_D <- result[4, 4]
  expect_equal(result_D, 0)
})

### test that correct type returned when as_list = FALSE

test_that("Correct shared information content without unique genes returned", {
  data <- readRDS(test_path("testdata", "adj_check.rds"))
  result <- suppressWarnings(inf_adj_matrix(data))
  expect_type(result, "double")
})

### test that correct type returned when as_list = TRUE

test_that("Correct shared information content without unique genes returned", {
  data <- readRDS(test_path("testdata", "adj_check.rds"))
  result <- suppressWarnings(inf_adj_matrix(data, as_list = TRUE))
  expect_type(result, "list")
})


## normalise_gene_scores()

test_that("Correct normalise values returned returned", {
  data <- readRDS(test_path("testdata", "adj_check.rds"))
  result <- suppressWarnings(normalise_gene_scores(data))
  result_D1 <- result |> purrr::pluck(5,1)
  result_D2 <- result |> purrr::pluck(5,2)
  expect_equal(result_D1, 0.2)
  expect_equal(result_D2, 0.2)
})