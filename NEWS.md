# ipw 1.3.1

* **Core Refactoring (`ipwtm`, `ipwpoint`, `ipwplot`, `tstartfun`)**:
  * Improved defensive input validation, parameter checking, and error messages across all core functions.
  * Modernized formula construction (`as.formula()`) to eliminate fragile `eval(parse(text = ...))` string manipulations.
  * Fixed variable symbol evaluation issues when passing variables as function arguments rather than raw names.
  * Modularized `ipwtm` into dedicated family modeling engines (`binomial`, `survival`, `multinomial`, `ordinal`, `gaussian`).
  * Replaced manual survival hazard step-differencing with continuous baseline hazard interpolation using `survival::basehaz` and `stats::approx`.
  * Refactored `tstartfun` to use safe column evaluation and direct vector lag logic (`c(-1, t[-length(t)])`).

* **Test Suite Expansion**:
  * Upgraded `testthat` suite to 3rd edition specs, expanding coverage for `multinomial`, `ordinal`, `gaussian`, unstabilized models, weight truncation, and invalid parameter rejection.
