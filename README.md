# SubGradSolver

[![Pipeline](https://gitlab.com/smspp/subgradsolver/badges/develop/pipeline.svg)](https://gitlab.com/smspp/subgradsolver/-/pipelines?ref=develop)
[![Release](https://gitlab.com/smspp/subgradsolver/-/badges/release.svg)](https://gitlab.com/smspp/subgradsolver/-/releases)
[![Documentation](https://img.shields.io/badge/docs-Doxygen-blue)](https://smspp.gitlab.io/smspp-project/)
[![C++ Standard](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

`SubGradSolver` is a SMS++ `CDASolver` for the minimization of a convex
function, or the maximization of a concave one, given as the sum of
`C05Function` (the components) and of a linear or separable quadratic 0-th
component, over a box of the variables; the Block is read as `BundleSolver`
reads it, so that `SubGradSolver` can be used wherever `BundleSolver` is,
in particular as the inner `Solver` of a `LagrangianDualSolver`. The method
is the unified subgradient method of

> A. Frangioni, E. Gorgone, B. Gendron. "On the computational efficiency of
> subgradient methods: a case study with Lagrangian bounds", Mathematical
> Programming Computation 9(4), 573-604, 2017

i.e., a (deflected, projected, incremental) subgradient method whose
stepsize and deflection rules are objects of the classes `Stepsize` and
`Deflection`, chosen by name with the parameters `strStepsize` and
`strDeflection` among those registered in their factories:

- the stepsize rules `Polyak`, `ColorTV` (that of the Volume algorithm) and
  `FumeroTV`, which are target value rules, the first and the last needing
  a finite lower bound on the function (from the Block or from the
  components);

- the deflection rule `Volume` (that of the modified Volume algorithm, the
  default with `ColorTV`), or none;

- `PrimalDual`, the Primal-Dual subgradient method with simple or weighted
  averages, which chooses both the stepsize and the deflection.

The best value found is a valid bound whenever `compute()` stops; the other
bound is only known when the Block gives one, or when the direction, a
convex combination of linearizations, has a finite minimum on the box of the
variables, which is always the case when all of them have finite bounds.
For each component with a global pool of linearizations, `SubGradSolver`
keeps there a convex combination of them, either that of the direction or
the ergodic average weighted by the stepsizes (`intPrimalRec`), and gives it
as the important linearization of the component in `get_dual_solution()`,
which is how a `LagrangianDualSolver` recovers its primal solution. When
the components only get new variables, e.g., the multipliers of new
constraints in a Lagrangian dual, `SubGradSolver` keeps what it knows and
extends it to them, the new entries of the subgradient and of the direction
coming from the global pools, rather than starting again. All the
parameters, those of the rules included, are parameters of `SubGradSolver`,
and therefore can be given in its `ComputeConfig`; the header
`SubGradSolver.h` describes them.


## Getting started

These instructions will let you build the `SubGradSolver` module on
your system.

### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and its
  requirements.

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
cmake --install .
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(SubGradSolver)
target_link_libraries(<my_target> SMS++::SubGradSolver)
```

### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `SubGradSolver` has to include a "main makefile" of
the module, which typically is either [makefile-c](makefile-c) including all
necessary libraries comprised the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). These in turn recursively include all the
required other makefiles, hence one should only need to edit the "main
makefile" for compilation type (C++ compiler and its options) and it all
should be good to go. In case some of the external libraries are not at their
default location, it should only be necessary to create the
`../extlib/makefile-paths` out of the `extlib/makefile-default-paths-*` for
your OS `*` and edit the relevant bits (commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.

### Tests

The unit test in `test/` solves three small problems with known optimum,
one of which has the structure of a Lagrangian dual and is changed between
the solutions, with one `SubGradSolver` for each rule and scheme, all
registered by the `BlockSolverConfig` `test/SGBSCfg.txt`.


## Getting help

If you need support, you want to submit bugs or propose a new feature, you
can [open a new issue](https://gitlab.com/smspp/subgradsolver/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

- **Enrico Gorgone**  
  Dipartimento di Matematica ed Informatica  
  Università di Cagliari

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors

- **Laura Fiorella**  
  Dipartimento di Matematica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.


## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
