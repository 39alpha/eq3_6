TARGET=../../bin/eqpt
SOURCES=$(wildcard src/*.f)
OBJECTS=$(subst src/,obj/,$(SOURCES:.f=.o))
LIBRARIES=../../lib/libeqlibu.a
FFLAGS:=-O3

include ../Makefile.exe.in
