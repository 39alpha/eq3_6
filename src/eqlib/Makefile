TARGET=../../lib/libeqlib.a
SOURCES=$(wildcard src/*.f)
OBJECTS=$(subst src/,obj/,$(SOURCES:.f=.o))
INCLUDES=-I../../include
FFLAGS:=-O3

include ../Makefile.lib.in
