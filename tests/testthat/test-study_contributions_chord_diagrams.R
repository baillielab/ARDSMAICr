## contributions_chord_bystudy()

### tests output consistent

test_that("Output is consistent", {
  vdiffr::expect_doppelganger(
    title = "plot",
    fig = suppressWarnings(contributions_chord_bystudy(data_contributions)),
  )
})

## contributions_chord_bymethod()

### tests output consistent

test_that("Output is consistent", {
  vdiffr::expect_doppelganger(
    title = "plot",
    fig = suppressWarnings(contributions_chord_bymethod(data_contributions, data_study)),
  )
})
