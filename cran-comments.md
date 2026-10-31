## CRAN Resubmission / Update Notes for ipw 1.3.1

This is a patch release updating the `ipw` package from version 1.3.0 to 1.3.1.

### Summary of Changes
* Refactored core functions (`ipwtm`, `ipwpoint`, `ipwplot`, `tstartfun`) to enhance code maintainability, error handling, input validation, and formula construction safety.
* Fixed variable name deparsing bugs when arguments are supplied as symbols/variables.
* Expanded unit test suite using `testthat` (3rd edition) to cover all supported model families and edge cases.

### Test Environments
* Local macOS / Linux (R 4.x.x)
* Windows Server 2022 (win-builder R-devel and R-release)
* R-hub multi-platform builder

### R CMD check results
There were 0 ERRORS, 0 WARNINGS, 0 NOTES.

### Reverse dependencies
No reverse dependency issues detected.
