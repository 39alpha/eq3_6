# eq3_6

[![Eq3/6](https://github.com/39alpha/eq3_6/actions/workflows/build.yml/badge.svg)](https://github.com/39alpha/eq3_6/actions/workflows/build.yml)

A Software Package for Geochemical Modeling

This is an updated version of [EQ3/6 Version 8.0a](https://www-gs.llnl.gov/energy-homeland-security/geochemistry)
developed at Lawrence Livermore National Laboratories. The changes include:

1. **Transitioning to Fortran 2008** though legacy `go to` and `COMMON` still linger, for example.
2. **Modernizing the build system and directory structure** using [`meson`](https://mesonbuild.com) for builds.
3. **Adding integration tests** against legacy outputs and CI-driven snapshots
4. **Minor bug fixes,** e.g. [b84f80c](https://github.com/39alpha/eq3_6/commit/b84f80c4b329145238336d0a536a6ec97811eefe)
5. **Usability features.** For example, the primary executables now take their inputs at the command line rather than
   relying on specifically named input files, e.g.

```bash
$ eqpt test/data/eqpt/cmp.d0
$ eq3nr cmp.d1 test/data/eq3nr/cmp/cristoba.3i
$ eq6 cmp.d1 test/data/eq6/cmp/crisqtz.6i
```

> **NOTE**: We will continue to apply bug fixes to the FORTRAN 77 sources, but those changes will be maintained on the
[v9.x.x](https://github.com/39alpha/eq3_6/tree/v9.x.x) branch and continue to use Makefile-based builds.

## Installation

```bash
$ git clone https://github.com/39alpha/eq3_6
$ just install
```
This will build all of the executables and libraries and install the following
directories:
- `bin` - contains the executable files `eq3nr`, `eq6`, `eqpt`, `xcon3` and `xcon6`.
- `lib` - contains both shared and static libraries, e.g. `libeq.a` and `libeq.so`.
- `include` - contains the various header files you may need to use the libraries.

By default `just install` will place all installation files in `/usr/local`, e.g. `/usr/local/bin/eq3nr`. If you wish to
change that, you can specify an installation prefix:
```bash
$ just install --prefix=/opt/eq3_6 -C build
```
which will place the files in `/opt/eq3_6`, e.g. `/opt/eq3_6/bin/eq3nr`.

## Build Dependencies

- [`meson`](https://mesonbuild.com/) `pip install meson`
- [`ninja`](https://ninja-build.org/) `pip install ninja`
- `gfortran` (pretty much any version as far as we can tell)
- `coreutils` (MacOS, testing only - `brew install coreutils`)
- `parallel` (MacOS, testing only - `brew install parallel`)

## Licensing

This repo includes the source files for [EQ3/6 v8.0a](https://www-gs.llnl.gov/energy-homeland-security/geochemistry),
a software package for geochemical modeling developed by Thomas Wolery at Lawrence Livermore National Laboratory under
the [BSD 3-Clause License](https://github.com/39alpha/eq3_6/blob/main/EQ36_LICENSE).

## References

- Wolery, T. J., and USDOE. EQ3/6 A Software Package for Geochemical Modeling. Computer software. December 13, 2010.
  https://www.osti.gov//servlets/purl/1231666. doi:https://doi.org/10.11578/dc.20210416.44.
- Wolery, T. J. and R. L. Jarek. Software User's Manual EQ36, Version 8.0. U.S. Tech. Rep. 2003. Department of Energy,
  Office of Civilian Radioactive Waste Management, Office of Repository Development. 10813-UM-8.0-00.
