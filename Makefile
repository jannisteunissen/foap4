.SUFFIXES:

FC := mpif90
CC := mpicc
FYPPFLAGS := -n
INCDIRS := p4est/build/local/include
LIBDIRS := p4est/build/local/lib
LIBS := p4est sc z m
CFLAGS := -Wall -O2 -g
TARGETS := test_refinement test_advection test_xdmf_writer test_euler	\
	benchmark_ghostcell

.PHONY: all
all: $(TARGETS)

# Determine compiler brand
compiler_version = $(shell $(FC) --version)
compiler_brand = $(word 1, $(compiler_version))

ifeq ($(compiler_brand), GNU)
	FFLAGS ?= -Wall -O2 -g -Jsrc -cpp $(FFLAGS_USER)
	ifeq ($(DEBUG), 1)
		FFLAGS += -O0 -fcheck=all
	endif
else ifeq ($(compiler_brand), nvfortran)
	FFLAGS ?= -Wall -acc=gpu -gpu=ccall -fast -Mpreprocess -static-nvidia	\
	-g -module src $(FFLAGS_USER)
else ifeq ($(compiler_brand), pgfortran)
	FFLAGS ?= -Wall -acc=gpu -fast -gpu=ccall -Mpreprocess -static-nvidia	\
	-g -module src $(FFLAGS_USER)
endif

# Dependencies
$(TARGETS): src/m_foap4.o src/p4est_wrapper.o src/m_xdmf_writer.o src/m_config.o
$(addsuffix .o,$(addprefix src/,$(TARGETS))): src/m_foap4.mod
src/test_euler.o: src/m_euler.mod
test_euler: src/m_euler.o
src/m_foap4.o: src/m_xdmf_writer.mod

.PHONY: clean
clean:
	$(RM) $(TARGETS) src/*.o src/*.mod src/*.smod

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

.PRECIOUS: src/%.f90
src/%.f90: src/%.fpp
	fypp $(FYPPFLAGS) $< $@
