PREFIX:=/usr

all:
	@mkdir -p bin lib
	+make -C src/eq3nr
	+make -C src/eq6
	+make -C src/eqpt
	+make -C src/xcon3
	+make -C src/xcon6

clean:
	@rm -rf bin/* lib/*
	+make -C src/eqlib clean
	+make -C src/eqlibu clean
	+make -C src/eqlibg clean
	+make -C src/eq3nr clean
	+make -C src/eq6 clean
	+make -C src/eqpt clean
	+make -C src/xcon3 clean
	+make -C src/xcon6 clean

install:
	mkdir -p $(PREFIX)/bin $(PREFIX)/lib $(PREFIX)/include
	cp -R bin/* $(PREFIX)/bin
	cp -R lib/* $(PREFIX)/lib
	cp -R include/* $(PREFIX)/include

.PHONY: clean install
