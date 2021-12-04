# eq3_6
A Software Package for Geochemical Modeling

This is a slightly updated version of
[EQ3/6 Version 8.0a](https://www-gs.llnl.gov/energy-homeland-security/geochemistry)
developed at Lawrence Livermore National Laboratories. At this point, the only changes made have
been to directory structure and adding a Make-based build system for easy building and installation.
Further changes will be made to address bug fixes and enhancements.

## Installation
```bash
$ make
$ make test
$ make install
```
This will build all of the executables and libraries and install the following directories:
- `bin` - contains the executable files `eq3nr`, `eq6`, `eqpt`, `xcon3` and `xcon6`
- `lib` - contains the library files `libeqlib.a`, `libeqlibu.a` and `libeqlibg.a`
- `include` - contains the various header files you may need to use the libraries
- `share` - contains various data files and the like

By default `make install` will place all installation files in `/usr/local`, e.g.
`/usr/local/bin/eq3nr`. If you wish to change that, you may specify an installation prefix:
```bash
$ make PREFIX=/opt/eq3_6 install
```
which will place the files in `/opt/eq3_6`, e.g. `/opt/eq3_6/bin/eq3nr`.

Of course, if you don't have write access to the install location, you will need to `sudo`.

By default, `make` will compile all of the targets with optimimization level-3 `-O3`. If you'd like
to change the compilation flags, e.g. to enable debugging `-g` or change the level `-O0`, you can do
that with the `FFLAGS` option:
```bash
$ make FFLAGS=-g

## Dependencies

- `gfortran` (pretty much any version as far as we can tell)
- `coreutils` (MacOS, testing only - `brew install coreutils`)

## References

- Wolery, T. J., and USDOE. EQ3/6 A Software Package for Geochemical Modeling. Computer software. December 13, 2010. https://www.osti.gov//servlets/purl/1231666. doi:https://doi.org/10.11578/dc.20210416.44.
- Wolery, T. J. and R. L. Jarek. Software User's Manual EQ36, Version 8.0. U.S. Tech. Rep. 2003. Department of Energy, Office of Civilian Radioactive Waste Management, Office of Repository Development. 10813-UM-8.0-00.
