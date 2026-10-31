test_that("tstartfun computes correct interval start times for multi-record subjects", {
  mydata <- data.frame(
    patient = c(1, 1, 1, 1, 2, 2, 2),
    time.days = c(14, 34, 41, 56, 0, 11, 28)
  )

  tstart <- tstartfun(patient, time.days, mydata)

  # First record for each patient must be -1
  expect_equal(tstart[1], -1)
  expect_equal(tstart[5], -1)

  # Subsequent records must equal previous time.days value
  expect_equal(tstart[2:4], c(14, 34, 41))
  expect_equal(tstart[6:7], c(0, 11))
})

test_that("tstartfun maintains original dataset row ordering when unsorted", {
  unsorted_data <- data.frame(
    id = c(1, 2, 1, 2, 1),
    time = c(20, 15, 10, 5, 30)
  )

  res <- tstartfun(id, time, unsorted_data)

  # Expected tstart for id=1 (times 10, 20, 30 -> tstarts -1, 10, 20)
  # and id=2 (times 5, 15 -> tstarts -1, 5)
  # Row 1: id=1, time=20 -> 2nd record -> tstart = 10
  # Row 2: id=2, time=15 -> 2nd record -> tstart = 5
  # Row 3: id=1, time=10 -> 1st record -> tstart = -1
  # Row 4: id=2, time=5  -> 1st record -> tstart = -1
  # Row 5: id=1, time=30 -> 3rd record -> tstart = 20
  expect_equal(res, c(10, 5, -1, -1, 20))
})

test_that("tstartfun handles single-record individuals", {
  single_dat <- data.frame(id = c(101, 102), time = c(50, 12))
  res <- tstartfun(id, time, single_dat)

  expect_equal(res, c(-1, -1))
})

test_that("tstartfun supports both unquoted symbols and character string column names", {
  mydata <- data.frame(patient = c(1, 1, 2), time.days = c(10, 20, 15))

  res_sym <- tstartfun(patient, time.days, mydata)
  res_str <- tstartfun("patient", "time.days", mydata)

  expect_equal(res_sym, res_str)
})

test_that("tstartfun throws helpful errors on invalid inputs", {
  mydata <- data.frame(patient = c(1, 1), time.days = c(10, 20))

  # Missing data argument
  expect_error(tstartfun(patient, time.days), "data")

  # Non-dataframe data
  expect_error(tstartfun(patient, time.days, data = list(patient = 1, time.days = 10)))

  # Missing values in id or timevar
  data_na_id <- data.frame(patient = c(1, NA), time.days = c(10, 20))
  expect_error(tstartfun(patient, time.days, data_na_id), "missing values")

  data_na_time <- data.frame(patient = c(1, 1), time.days = c(10, NA))
  expect_error(tstartfun(patient, time.days, data_na_time), "missing values")
})
