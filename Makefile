.SUFFIXES:

FC := mpif90
CC := mpicc
FYPPFLAGS := -n
INCDIRS := p4est/build/local/include
LIBDIRS := p4est/build/local/lib
LIBS := p4est sc z m
CFLAGS := -Wall -O2 -g
TARGETS_2D := build/test_refinement_2d build/test_advection_2d build/test_xdmf_writer_2d \
	build/test_euler_2d build/benchmark_ghostcell_2d
TARGETS_3D := build/test_refinement_3d build/test_advection_3d build/test_euler_3d \
	build/test_xdmf_writer_3d

OBJDIR := build
SRCDIR := src

.PHONY: all clean
all: $(TARGETS_2D) $(TARGETS_3D)

# Determine compiler brand
compiler_version = $(shell $(FC) --version)
compiler_brand = $(word 1, $(compiler_version))

ifeq ($(compiler_brand), GNU)
	FFLAGS ?= -Wall -O2 -g -J$(OBJDIR) -cpp $(FFLAGS_USER) -Wno-unused-dummy-argument -Wl,--no-warn-execstack
	ifeq ($(DEBUG), 1)
		FFLAGS += -O0 -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan
		CFLAGS += -O0
	endif
else ifeq ($(compiler_brand), nvfortran)
	FFLAGS ?= -Minform=warn -acc=gpu,strict -fast -gpu=ccnative -Mpreprocess \
	-static-nvidia -g -module $(OBJDIR) $(FFLAGS_USER)
else ifeq ($(compiler_brand), pgfortran)
	FFLAGS ?= -Minform=warn -acc=gpu,strict -fast -gpu=ccnative -Mpreprocess \
	-static-nvidia -g -module $(OBJDIR) $(FFLAGS_USER)
endif

# General dependencies
$(TARGETS_2D): $(OBJDIR)/m_foap4_2d.o $(OBJDIR)/p4est_wrapper_2d.o $(OBJDIR)/m_config.o
$(TARGETS_3D): $(OBJDIR)/m_foap4_3d.o $(OBJDIR)/p4est_wrapper_3d.o $(OBJDIR)/m_config.o
$(addsuffix .o,$(TARGETS_2D)): $(OBJDIR)/m_foap4_2d.o $(OBJDIR)/m_config.o
$(addsuffix .o,$(TARGETS_3D)): $(OBJDIR)/m_foap4_3d.o $(OBJDIR)/m_config.o

# Dependencies for 2D targets
$(OBJDIR)/test_advection_2d: $(OBJDIR)/m_physics_advection.o
$(OBJDIR)/test_euler_2d: $(OBJDIR)/m_physics_euler_2d.o

# Dependencies for 3D targets
$(OBJDIR)/test_advection_3d: $(OBJDIR)/m_physics_advection.o
$(OBJDIR)/test_euler_3d: $(OBJDIR)/m_physics_euler_3d.o

# Ensure build directory exists
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Compile Fortran source files to object files
$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# Compile Fortran executables from .o object files
$(OBJDIR)/%: $(OBJDIR)/%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

# Compile C source files to object files for 2D
$(OBJDIR)/%_2d.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) -c -o $@ $< $(CFLAGS) -DNDIM=2 $(addprefix -I,$(INCDIRS))

# Compile C source files to object files for 3D
$(OBJDIR)/%_3d.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) -c -o $@ $< $(CFLAGS) -DNDIM=3 $(addprefix -I,$(INCDIRS))

# Clean target
clean:
	$(RM) $(OBJDIR)/*

.PRECIOUS: $(SRCDIR)/%.f90 $(SRCDIR)/%_2d.f90 $(SRCDIR)/%_3d.f90
$(SRCDIR)/%.f90: $(SRCDIR)/%.fpp
	fypp $(FYPPFLAGS) $(FYPP_USER) $< $@
$(SRCDIR)/%_2d.f90: $(SRCDIR)/%.fpp
	fypp $(FYPPFLAGS) $(FYPP_USER) -D NDIM=2 $< $@
$(SRCDIR)/%_3d.f90: $(SRCDIR)/%.fpp
	fypp $(FYPPFLAGS) $(FYPP_USER) -D NDIM=3 $< $@
