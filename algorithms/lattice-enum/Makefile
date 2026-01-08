.PHONY: all clean

all: lattice-enum lattice-enum_debug

lattice-enum: lattice-enum.c Makefile
	gcc -O2 lattice-enum.c -o lattice-enum

lattice-enum_debug: lattice-enum.c Makefile
	gcc -O0 -g lattice-enum.c -o lattice-enum_debug

clean:
	rm -rf lattice-enum lattice-enum_debug
