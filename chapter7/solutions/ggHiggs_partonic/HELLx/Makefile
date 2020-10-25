CC = g++

PKG = hell-x

CFLAGS = -O2 -Wall -I./include

obj = objects

OBJECTS = $(obj)/hell-x.o $(obj)/special_functions.o $(obj)/expansionSFs.o

examples: examples/test_hell-x

all: examples
lib: lib$(PKG).a

examples/%: examples/%.cc lib$(PKG).a
	$(CC) $(CFLAGS) -o $@ $< lib$(PKG).a

lib$(PKG).a: $(OBJECTS)
	ar -r lib$(PKG).a $^

$(obj)/%.o: src/%.cc include/%.hh include/math/matrix.hh include/expansionSFs.hh include/version.hh
	$(CC) $(CFLAGS) -o $@ -c $<

$(obj)/%.o: src/math/%.cc include/math/%.hh
	$(CC) $(CFLAGS) -o $@ -c $<




clean: clean-tilde
	rm -f $(examples)
	rm -f lib$(PKG).a
	rm -f $(obj)/*.o

clean-tilde:
	rm -f *~
	rm -f src/*~
	rm -f include/*~
	rm -f include/math/*~


.PHONY: clean-tilde clean lib examples all
