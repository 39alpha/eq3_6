TARGET=../../lib/libeqlibg.a
SOURCES=$(wildcard src/*.f90)
OBJECTS=$(subst src/,obj/,$(SOURCES:.f90=.o))
INCLUDES=-I../../include

all: $(TARGET)

$(TARGET): $(OBJECTS)
	ar rsc $(TARGET) $^

obj/%.o: src/%.f90
	@mkdir -p obj
	gfortran $(FFLAGS) $(INCLUDES) -c -o $@ $^

clean:
	@rm -rf obj

.PHONY: clean
