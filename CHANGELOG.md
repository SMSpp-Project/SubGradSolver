# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- the module layout of SMS++ (`include/`, `src/`, `test/`, the CMake and
  makefile builds, the CI), with a unit test solving three small problems
  with every rule and scheme
- the convex combination of linearizations kept in the global pools for the
  primal recovery, either that of the direction or the ergodic average
  (`intPrimalRec`), given by `get_dual_solution()` as the important
  linearization of each component
- a valid lower bound from the direction on the box of the variables
- the handling of the Modification of the functions, of the variables and
  of their bounds
- the addition of variables to the components (strongly quasi-additive,
  with 0 in their box) is followed rather than starting again: the values,
  the points, the subgradient and the direction are kept and extended to
  the new variables, whose entries come from the linearizations in the
  global pools, the direction starting again only if it is not there, and
  the Stepsize is formatted; the other changes of the structure still read
  it again and forget what is known about the functions
- `dblNZEps`: when the norm of the direction projected on the tangent cone
  of the box is at most `dblNZEps`, the value at the stability center minus
  the linearization error is a lower bound, and that norm is also the one of
  the stopping condition, whatever `intProjection` is

- the C++20 module `smspp.subgradsolver` (`src/SubGradSolver.cppm`), which
  exports every name that the headers declare at namespace scope and
  re-exports `smspp`, so that a source file can write `import
  smspp.subgradsolver;` in place of the `#include` of the headers; it is
  built with the CMake option `SMSpp_BUILD_MODULE`, as the target
  `SMS++::SubGradSolver_module`

### Changed

- the algorithm re-ported from the `SubGrad` of NDO onto the interface of
  `C05Function`: the subgradient is the sum of those of all the components,
  the 0-th one included, and the linearization errors are those at the
  stability center
- the parameters have descriptive names (`intProjection`, `intScheme`,
  `dbltStar`, `strStepsize`, `strDeflection`, ...), which replace the old
  ones, and those of the rules are parameters of `SubGradSolver`
- `Stepsize` and `Deflection` have their own factories, in which
  `PrimalDual` is registered as both

### Fixed

- `intFumeroEta1` and `intFumeroEta2` of `FumeroTV` were swapped with
  respect to the paper: `intFumeroEta2` is the count of failures of the
  first phase, `intFumeroEta1` that of the second

- the subgradient was never computed, the projections were never written
  back, the bounds of the variables had the wrong sign and a concave
  function was not maximized

## [0.3.0] - 2021-02-05

### Changed

- Maintenance release.

## [0.2.0] - 2020-03-06

### Added

- Changelog.

## [0.1.0] - 2019-11-05

### Added

- First test release.

[Unreleased]: https://gitlab.com/smspp/subgradsolver/-/compare/0.3.0...develop
[0.3.0]: https://gitlab.com/smspp/subgradsolver/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/subgradsolver/-/compare/0.1.0...0.2.0
[0.1.0]: https://gitlab.com/smspp/subgradsolver/-/tags/0.1.0
