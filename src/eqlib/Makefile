TARGET=../../lib/libeqlib.a
SOURCES=$(wildcard src/*.f)
OBJECTS=$(subst src/,obj/,$(SOURCES:.f=.o))

all: $(TARGET)

$(TARGET): $(OBJECTS)
	ar rsc $(TARGET) $^

obj/%.o: src/%.f
	@mkdir -p obj
	gfortran -c -o $@ $^

clean:
	@rm -rf obj
