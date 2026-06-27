# fcube — Makefile
#
# Targets:
#   make            → build library (default)
#   make test       → build and run all tests
#   make clean      → remove build artefacts
#   make distclean  → clean + remove the library archive
#
# Compiler selection (default: gfortran):
#   make FC=ifort   (Intel classic)
#   make FC=ifx     (Intel oneAPI)
#   make FC=nagfor  (NAG)
#
# Build type (default: release):
#   make BUILD=debug    → -O0 -g -fcheck=all -fbacktrace
#   make BUILD=release  → -O2

# ------------------------------------------------------------------ #
#  Directories — declared first so FFLAGS can reference BUILDDIR
# ------------------------------------------------------------------ #

SRCDIR   := src
TESTDIR  := test
BUILDDIR := build

# ------------------------------------------------------------------ #
#  Compiler
# ------------------------------------------------------------------ #

# Guard against the make built-in FC=f77 default and environment noise.
ifeq ($(origin FC),default)
  FC := gfortran
endif
ifeq ($(FC),f77)
  FC := gfortran
endif

# ------------------------------------------------------------------ #
#  Flags
# ------------------------------------------------------------------ #

BUILD ?= release

FFLAGS_COMMON := -std=f2018 -Wall -Wextra -Wno-maybe-uninitialized -J$(BUILDDIR)

ifeq ($(BUILD),debug)
  FFLAGS := $(FFLAGS_COMMON) -O0 -g -fcheck=all -fbacktrace
else
  FFLAGS := $(FFLAGS_COMMON) -O2
endif

# ------------------------------------------------------------------ #
#  Sources — listed in dependency order (each module before its users)
# ------------------------------------------------------------------ #

LIB_SRCS := \
  $(SRCDIR)/cube_kinds.f90 \
  $(SRCDIR)/cube_data.f90  \
  $(SRCDIR)/cube_io.f90    \
  $(SRCDIR)/cube_arith.f90 \
  $(SRCDIR)/cube_fft.f90   \
  $(SRCDIR)/cube_diff.f90

LIB_OBJS := $(patsubst $(SRCDIR)/%.f90, $(BUILDDIR)/%.o, $(LIB_SRCS))

# ------------------------------------------------------------------ #
#  Library archive
# ------------------------------------------------------------------ #

LIB := $(BUILDDIR)/libfcube.a

# ------------------------------------------------------------------ #
#  Test sources and binaries
# ------------------------------------------------------------------ #

TEST_SRCS := \
  $(TESTDIR)/test_kinds.f90 \
  $(TESTDIR)/test_data.f90  \
  $(TESTDIR)/test_io.f90    \
  $(TESTDIR)/test_arith.f90 \
  $(TESTDIR)/test_diff.f90

TEST_BINS := $(patsubst $(TESTDIR)/%.f90, $(BUILDDIR)/%, $(TEST_SRCS))

# ------------------------------------------------------------------ #
#  Default target
# ------------------------------------------------------------------ #

# ------------------------------------------------------------------ #
#  Application
# ------------------------------------------------------------------ #

APP_SRC := app/main.f90
APP_BIN := $(BUILDDIR)/fcube_example

$(APP_BIN): $(APP_SRC) $(LIB) | $(BUILDDIR)
	$(FC) $(FFLAGS) $< -o $@ -L$(BUILDDIR) -lfcube

.PHONY: all test example clean distclean

all: $(LIB)

example: $(APP_BIN)
	@echo ""
	./$(APP_BIN)

# ------------------------------------------------------------------ #
#  Library build
# ------------------------------------------------------------------ #

$(LIB): $(LIB_OBJS)
	ar rcs $@ $^
	@echo "Built $@"

# Compile a library source into build/
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# ------------------------------------------------------------------ #
#  Module dependency chain (USE order)
# ------------------------------------------------------------------ #

$(BUILDDIR)/cube_data.o:  $(BUILDDIR)/cube_kinds.o
$(BUILDDIR)/cube_io.o:    $(BUILDDIR)/cube_data.o
$(BUILDDIR)/cube_arith.o: $(BUILDDIR)/cube_data.o
$(BUILDDIR)/cube_fft.o:   $(BUILDDIR)/cube_kinds.o
$(BUILDDIR)/cube_diff.o:  $(BUILDDIR)/cube_fft.o   $(BUILDDIR)/cube_data.o  $(BUILDDIR)/cube_kinds.o

# ------------------------------------------------------------------ #
#  Test binaries
# ------------------------------------------------------------------ #

$(BUILDDIR)/test_%: $(TESTDIR)/test_%.f90 $(LIB) | $(BUILDDIR)
	$(FC) $(FFLAGS) $< -o $@ -L$(BUILDDIR) -lfcube

# Per-test module dependencies (so make knows when to relink)
$(BUILDDIR)/test_kinds: $(BUILDDIR)/cube_kinds.o
$(BUILDDIR)/test_data:  $(BUILDDIR)/cube_data.o  $(BUILDDIR)/cube_kinds.o
$(BUILDDIR)/test_io:    $(BUILDDIR)/cube_io.o    $(BUILDDIR)/cube_data.o  \
                        $(BUILDDIR)/cube_kinds.o
$(BUILDDIR)/test_arith: $(BUILDDIR)/cube_arith.o $(BUILDDIR)/cube_data.o  \
                        $(BUILDDIR)/cube_kinds.o
$(BUILDDIR)/test_diff:  $(BUILDDIR)/cube_diff.o  $(BUILDDIR)/cube_fft.o   \
                        $(BUILDDIR)/cube_data.o   $(BUILDDIR)/cube_kinds.o

# ------------------------------------------------------------------ #
#  Run tests
# ------------------------------------------------------------------ #

# Tests need to find the .cube fixture relative to the project root,
# so we run each binary from the project root directory.
test: $(TEST_BINS)
	@echo ""
	@echo "=== Running tests ==="
	@fail=0; \
	for bin in $(TEST_BINS); do \
	  echo "--- $$bin ---"; \
	  if ./$$bin; then :; else fail=$$((fail+1)); fi; \
	done; \
	echo ""; \
	if [ $$fail -eq 0 ]; then \
	  echo "All test suites passed."; \
	else \
	  echo "$$fail suite(s) FAILED."; \
	  exit 1; \
	fi

# ------------------------------------------------------------------ #
#  Housekeeping
# ------------------------------------------------------------------ #

clean:
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.mod $(TEST_BINS)

distclean: clean
	rm -f $(LIB)
	rmdir $(BUILDDIR) 2>/dev/null || true
