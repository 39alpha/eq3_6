TARGET=../../bin/xcon6
SOURCES=$(wildcard src/*.f90)
OBJECTS=$(subst src/,obj/,$(SOURCES:.f90=.o))
LIBRARIES=../../lib/libeqlibu.a
INCLUDES=-I../../include

include ../Makefile.exe.in
