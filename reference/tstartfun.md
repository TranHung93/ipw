# Compute Starting Time For Counting Process Notation

Computes starting time for follow-up intervals in longitudinal data
using counting process notation. Within each subject, starting time
equals -1 for the initial record and the previous record's follow-up
time for subsequent records.

## Usage

``` r
tstartfun(id, timevar, data)
```

## Arguments

- id:

  Patient or unit identifier (unquoted symbol, string, or vector).

- timevar:

  Follow-up time variable (unquoted symbol, string, or vector).

- data:

  Data frame containing `id` and `timevar`.

## Value

Numeric vector of starting times in the same row order as `data`.

## Author

Willem M. van der Wal <willem@vanderwalresearch.com>, Ronald B. Geskus
<rgeskus@oucru.org>

## Examples

``` r
mydata1 <- data.frame(
  patient = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2),
  time.days = c(14, 34, 41, 56, 72, 98, 0, 11, 28, 35)
)
mydata1$tstart <- tstartfun(patient, time.days, mydata1)
mydata1
#>    patient time.days tstart
#> 1        1        14     -1
#> 2        1        34     14
#> 3        1        41     34
#> 4        1        56     41
#> 5        1        72     56
#> 6        1        98     72
#> 7        2         0     -1
#> 8        2        11      0
#> 9        2        28     11
#> 10       2        35     28
```
