.SUFFIXES:

FC := mpif90
CC := mpicc
INCDIRS := $(HOME)/opt/include
LIBDIRS := $(HOME)/opt/lib
LIBS := p4est sc z m
CFLAGS := -Wall -O0 -g
# FFLAGS := -Wall -O0 -g -fcheck=all
TARGET := first_test

# Determine compiler brand
compiler_version = $(shell $(FC) --version)
compiler_brand = $(word 1, $(compiler_version))

ifeq ($(compiler_brand), GNU)
	FFLAGS ?= -Wall -O0 -g -fcheck=all
else ifeq ($(compiler_brand), nvfortran)
	FFLAGS ?= -Wall -acc=gpu -fast -Mpreprocess -static-nvidia -g
endif

first_test: m_foap4.o p4est_wrapper.o m_xdmf_writer.o

m_foap4.o: m_xdmf_writer.mod

.PHONY: all
all: $(TARGET)

.PHONY: clean
clean:
	$(RM) $(TARGET) *.o *.mod *.smod

# How to get .o object files from .c source files
%.o: %.c
	$(CC) -c -o $@ $(OBJ) $(CFLAGS) $(addprefix -I,$(INCDIRS)) -o $@ $<

# How to get .o object files from .f90 source files
%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get .mod files from .f90 source files (remake only if they have been
# removed, otherwise assume they are up to date)
%.mod: %.f90 %.o
	@test -f $@ || $(FC) -c -o $(@:.mod=.o) $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get Fortran executables from .o object files
%: %.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
