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

## Dependencies

- `gfortran` (pretty much any version as far as we can tell)

## References

[NEED REFERENCE FOR ORIGINAL EQ3/6]
