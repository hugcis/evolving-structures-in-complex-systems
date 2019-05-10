CC=gcc
LD=gcc
CFLAGS=-I. -Wall -O2 -march=native -funroll-loops -ffast-math
LDFLAGS=-I. -Wall -lz
SOURCES:=$(wildcard src/*.c)
OBJS:=$(SOURCES:src/%.c=build/%.o)
TARGET:=bin/automaton

all: $(TARGET)

$(TARGET): $(OBJS)
	@mkdir -p bin
	$(LD) $(LDFLAGS) $+ -o $@

build/%.o: src/%.c
	@mkdir -p build
	$(CC) $(CFLAGS) $+ -c -o $@

.PHONY: clean
clean:
	rm -rf $(OBJS)
	rm -rf $(TARGET)
	rmdir bin
	rmdir build

.PHONY: rebuild
rebuild:
	$(MAKE) clean
	$(MAKE) all
