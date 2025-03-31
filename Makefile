.SUFFIXES:

FC := mpif90
CC := mpicc
INCDIRS := p4est/build/local/include
LIBDIRS := p4est/build/local/lib
LIBS := p4est sc z m
CFLAGS := -Wall -O2 -g
TARGET := test_refinement test_advection test_xdmf_writer benchmark_ghostcell

.PHONY: all
all: $(TARGET)

# Determine compiler brand
compiler_version = $(shell $(FC) --version)
compiler_brand = $(word 1, $(compiler_version))

ifeq ($(compiler_brand), GNU)
	FFLAGS ?= -Wall -Wextra -Wrealloc-lhs -O0 -g -fcheck=all -Jsrc
	# FFLAGS ?= -Wall -O2 -g -Jsrc
else ifeq ($(compiler_brand), nvfortran)
	FFLAGS ?= -Wall -acc=gpu -fast -Mpreprocess -static-nvidia -g -module src
endif

$(TARGET): src/m_foap4.o src/p4est_wrapper.o src/m_xdmf_writer.o

# Dependencies
src/m_foap4.o: src/m_xdmf_writer.mod
src/test_refinement.o: src/m_foap4.mod
src/benchmark_ghostcell.o: src/m_foap4.mod
src/test_advection.o: src/m_foap4.mod
src/test_xdmf_writer.o: src/m_foap4.mod

.PHONY: clean
clean:
	$(RM) $(TARGET) src/*.o src/*.mod src/*.smod

# How to get .o object files from .c source files
src/%.o: src/%.c
	$(CC) -c -o $@ $(OBJ) $(CFLAGS) $(addprefix -I,$(INCDIRS)) -o $@ $<

# How to get .o object files from .f90 source files
src/%.o: src/%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get .mod files from .f90 source files (remake only if they have been
# removed, otherwise assume they are up to date)
src/%.mod: src/%.f90 src/%.o
	@test -f $@ || $(FC) -c -o $(@:.mod=.o) $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get Fortran executables from .o object files
%: src/%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
