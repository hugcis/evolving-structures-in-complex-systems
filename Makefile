CC=gcc
LD=gcc
CFLAGS=-Isrc -Wall -O3 -march=native -funroll-loops -ffast-math -flto=thin
LDFLAGS=-Wall -lz  -framework Accelerate -lgsl -O3 -flto=thin

SRCDIR:=src
BUILDDIR:=build
BINDIR:=bin

SRCEXT:=c
MKDIR:=@mkdir -p
SOURCES:=$(shell find $(SRCDIR) -type f -name "*.$(SRCEXT)")
OBJS:=$(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

DIRS:= automaton utils nn search
SOURCEDIRS:=$(foreach dir,$(DIRS),$(addprefix $(SRCEDIR)/, $(dir)))
TARGETDIRS:=$(foreach dir,$(DIRS),$(addprefix $(BUILDDIR)/, $(dir)))

TARGET:=$(BINDIR)/automaton

all: directories $(TARGET)
	$(MAKE) tools/viz/step_to_ppm

$(TARGET): $(OBJS)
	$(MKDIR) $(BINDIR)
	$(LD) $(LDFLAGS) $+ -o $@

.PHONY: directories
directories:
	$(MKDIR) $(TARGETDIRS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	$(MKDIR) $(BUILDDIR)
	$(CC) $(CFLAGS) $(XCFLAGS) $+ -c -o $@

.PHONY: clean
clean:
	rm -rf $(OBJS)
	rm -rf $(TARGET)
	rmdir $(TARGETDIRS)
	rmdir $(BINDIR)
	rmdir $(BUILDDIR)

tools/viz/step_to_ppm:
	gcc ./tools/viz/step_to_ppm.c -o ./tools/viz/step_to_ppm -lnetpbm

.PHONY: rebuild
rebuild:
	$(MAKE) clean
	$(MAKE) all
