TARGET=../../bin/xcon6
SOURCES=$(wildcard src/*.f)
OBJECTS=$(subst src/,obj/,$(SOURCES:.f=.o))
LIBRARIES=../../lib/libeqlibu.a
INCLUDES=-I../../include
FFLAGS:=-O3

include ../Makefile.exe.in
