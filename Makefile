CC=gcc
CFLAGS=-I. -Wall -lz -O2

.PHONY: clean

2d_automaton: 2d_automaton.c compress.c
	$(CC) $(CFLAGS) -o 2d_automaton 2d_automaton.c compress.c

wolfram_automaton: wolfram_automaton.c compress.c
	$(CC) $(CFLAGS) -o wolfram_automaton wolfram_automaton.c compress.c

clean:
	  $(RM) 2d_automaton wolfram_automaton
