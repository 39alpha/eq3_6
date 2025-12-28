TARGET=../../bin/eq3nr
SOURCES=$(wildcard src/*.f90)
OBJECTS=$(subst src/,obj/,$(SOURCES:.f90=.o))
LIBRARIES=../../lib/libeqlib.a ../../lib/libeqlibu.a ../../lib/libeqlibg.a
INCLUDES=-I../../include

include ../Makefile.exe.in
