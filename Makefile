#Compiler options
GCC = gcc
CFLAGS = -Wall -Wshadow=global -fopenmp -march=native -O4

all:
	$(GCC) src/random.c -c -o random.o $(CFLAGS)
	$(GCC) src/anyrng.c -o anyrng random.o -lm $(CFLAGS)

example:
	$(GCC) src/example.c -o example $(CFLAGS)

clean:
	rm -f random.o
	rm -f anyrng
	rm -f example
