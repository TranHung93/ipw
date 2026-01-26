library(testthat)
library(ipw)

test_that("tstartfun correctly identifies starting times and handles -1 default", {
  # Create a simple test dataset
  mydata <- data.frame(
    patient = c(1, 1, 1, 2, 2),
    time.days = c(10, 20, 30, 5, 15)
  )

  # Run function
  tstart <- tstartfun(patient, time.days, mydata)

  # 1. Check for the -1 default on first records
  # For Patient 1, time 10 should start at -1
  # For Patient 2, time 5 should start at -1
  expect_equal(tstart[1], -1)
  expect_equal(tstart[4], -1)

  # 2. Check that subsequent records use the previous timevar value
  # Patient 1: [-1, 10], [10, 20], [20, 30]
  expect_equal(tstart[2], 10)
  expect_equal(tstart[3], 20)

  # Patient 2: [-1, 5], [5, 15]
  expect_equal(tstart[5], 5)
})

test_that("tstartfun maintains original data order", {
  # Create unsorted data
  unsorted_data <- data.frame(
    id = c(1, 2, 1),
    time = c(20, 5, 10)
  )

  # The function should return tstart in the order: [10, -1, -1]
  # Because:
  # Row 1 (id 1, time 20) is the 2nd record for id 1 -> tstart = 10
  # Row 2 (id 2, time 5) is the 1st record for id 2 -> tstart = -1
  # Row 3 (id 1, time 10) is the 1st record for id 1 -> tstart = -1

  res <- tstartfun(id, time, unsorted_data)

  expect_equal(res, c(10, -1, -1))
})

test_that("tstartfun handles single-record individuals", {
  single_dat <- data.frame(id = 101, time = 50)
  res <- tstartfun(id, time, single_dat)

  expect_equal(res, -1)
})
