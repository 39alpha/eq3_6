TARGET=../../bin/eq6
MODSOURCES=src/mod6pt.f90 src/mod6xf.f90
SOURCES=$(wildcard src/*.f90)
OBJECTS=$(subst src/,obj/,$(SOURCES:.f90=.o))
MODOBJECTS=$(subst src/,obj/,$(MODSOURCES:.f90=.o))
LIBRARIES=../../lib/libeqlib.a ../../lib/libeqlibg.a ../../lib/libeqlibu.a 
INCLUDES=-I../../include

include ../Makefile.exe.in
