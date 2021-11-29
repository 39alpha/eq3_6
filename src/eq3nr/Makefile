TARGET=../../bin/eq3nr
SOURCES=$(wildcard src/*.f)
OBJECTS=$(subst src/,obj/,$(SOURCES:.f=.o))
LIBRARIES=../../lib/libeqlib.a ../../lib/libeqlibu.a ../../lib/libeqlibg.a

all: $(TARGET)

$(TARGET): $(OBJECTS) $(LIBRARIES)
	gfortran -o $(TARGET) $^

obj/%.o: src/%.f
	@mkdir -p obj
	gfortran -c -o $@ $^

../../lib/libeqlib.a:
	+make -C ../eqlib

../../lib/libeqlibu.a:
	+make -C ../eqlibu

../../lib/libeqlibg.a:
	+make -C ../eqlibg

clean:
	@rm -rf obj
