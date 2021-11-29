# eq3_6
A Software Package for Geochemical Modeling

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

By default `make install` will place all installation files in `/usr`, e.g. `/usr/bin/eq3nr`. If you
wish to change that, you may specify an installation prefix:
```bash
$ make PREFIX=/opt/eq3_6 install
```
which will place the files in `/opt/eq3_6`, e.g. `/opt/eq3_6/bin/eq3nr`.

Of course, if you don't have write access to the install location, you will need to `sudo`.
