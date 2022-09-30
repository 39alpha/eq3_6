STATIC_TARGET=../../lib/libeqlib.a
SHARED_TARGET=../../lib/libeqlib.so
SOURCES=$(wildcard src/*.f)
OBJECTS=$(subst src/,obj/,$(SOURCES:.f=.o))
INCLUDES=-I../../include

include ../Makefile.lib.in
