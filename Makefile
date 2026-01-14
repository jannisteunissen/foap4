.SUFFIXES:

# ==============================================================================
# Compiler and tool definitions
# ==============================================================================
F90C ?= mpif90
CC := mpicc
FYPP := fypp

# ==============================================================================
# Directory structure
# ==============================================================================
BUILDDIR := build
OBJDIR := $(BUILDDIR)/obj
BINDIR := $(BUILDDIR)/bin
GENDIR := $(BUILDDIR)/generated

SRCDIR := src
COREDIR := $(SRCDIR)/core
PHYSICSDIR := $(SRCDIR)/physics
NUMERICSDIR := $(SRCDIR)/numerics
TESTSDIR := $(SRCDIR)/tests
UTILSDIR := $(SRCDIR)/utils

vpath %.fpp $(COREDIR) $(PHYSICSDIR) $(NUMERICSDIR) $(TESTSDIR) $(UTILSDIR)
vpath %.f90 $(COREDIR) $(PHYSICSDIR) $(NUMERICSDIR) $(TESTSDIR) $(UTILSDIR)

# ==============================================================================
# External libraries
# ==============================================================================
INCDIRS := p4est/build/local/include
INCDIRS += $(COREDIR) $(PHYSICSDIR) $(UTILSDIR) $(NUMERICSDIR)
INCFLAGS := $(addprefix -I,$(INCDIRS))
LIBDIRS := p4est/build/local/lib
LIBS := p4est sc z m

# ==============================================================================
# Preprocessor flags
# ==============================================================================
FYPPFLAGS := -n $(INCFLAGS)

# ==============================================================================
# C compiler flags
# ==============================================================================
CFLAGS := -Wall -O2 -g

# ==============================================================================
# Compiler detection and Fortran flags
# ==============================================================================
compiler_version := $(shell $(F90C) --version 2>/dev/null)
compiler_brand ?= $(word 1, $(compiler_version))

ifeq ($(compiler_brand),GNU)
    FFLAGS ?= -Wall -O2 -g -J$(OBJDIR) -cpp -Wno-unused-dummy-argument -Wl,--no-warn-execstack
    ifeq ($(DEBUG),1)
        FFLAGS += -O0 -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan
        CFLAGS += -O0
    endif
else ifneq (,$(filter $(compiler_brand),nvfortran pgfortran))
    # Will be used for nvfortran and pgfortran
    FFLAGS ?= -Minform=warn -acc=gpu,strict -fast -gpu=ccnative -Mpreprocess \
        -static-nvidia -g -module $(OBJDIR)
else ifeq ($(compiler_brand),Cray)
    FFLAGS ?= -M878 -O2 -eT -ea -ef -ffree -h acc \
        -h acc_model=auto_async_none:fast_addr:no_deep_copy -J$(OBJDIR)
else
    $(warning Unknown compiler "$(compiler_brand)", using default flags)
    FFLAGS ?= -O2 -g -J$(OBJDIR)
endif

# Add user flags
FFLAGS += $(FFLAGS_USER)

# ==============================================================================
# Target definitions
# ==============================================================================
TARGETS_2D := $(addprefix $(BINDIR)/,\
    test_refinement_2d \
    test_advection_2d \
    test_xdmf_writer_2d \
    test_euler_2d \
    test_benchmark_ghostcell_2d)

TARGETS_3D := $(addprefix $(BINDIR)/,\
    test_refinement_3d \
    test_advection_3d \
    test_euler_3d \
    test_xdmf_writer_3d)

# ==============================================================================
# Common object files
# ==============================================================================
COMMON_OBJS_2D := $(OBJDIR)/m_foap4_2d.o $(OBJDIR)/p4est_wrapper_2d.o $(OBJDIR)/m_config.o
COMMON_OBJS_3D := $(OBJDIR)/m_foap4_3d.o $(OBJDIR)/p4est_wrapper_3d.o $(OBJDIR)/m_config.o

# ==============================================================================
# Phony targets
# ==============================================================================
.PHONY: all clean 2d 3d help

all: $(TARGETS_2D) $(TARGETS_3D)

2d: $(TARGETS_2D)

3d: $(TARGETS_3D)

help:
	@echo "Usage: make [target] [options]"
	@echo ""
	@echo "Targets:"
	@echo "  all      - Build all 2D and 3D targets (default)"
	@echo "  2d       - Build only 2D targets"
	@echo "  3d       - Build only 3D targets"
	@echo "  clean    - Remove all build artifacts"
	@echo "  help     - Show this help message"
	@echo ""
	@echo "Options:"
	@echo "  DEBUG=1      - Enable debug flags"
	@echo "  F90C=<comp>  - Set Fortran compiler (default: mpif90)"
	@echo "  FFLAGS_USER= - Additional Fortran flags"
	@echo "  FYPP_USER=   - Additional fypp flags"
	@echo ""
	@echo "Directory structure:"
	@echo "  $(SRCDIR)/core/     - Core framework files"
	@echo "  $(SRCDIR)/physics/  - Physics modules"
	@echo "  $(SRCDIR)/tests/    - Test programs"
	@echo "  $(SRCDIR)/numerics/ - Numerical methods"
	@echo "  $(BUILDDIR)/        - All build artifacts"
	@echo ""
	@echo "Detected compiler: $(compiler_brand)"

clean:
	$(RM) -r $(BUILDDIR)

# ==============================================================================
# Directory creation
# ==============================================================================
$(OBJDIR) $(BINDIR) $(GENDIR):
	mkdir -p $@

# ==============================================================================
# Fypp preprocessing rules
# ==============================================================================
.PRECIOUS: $(GENDIR)/%.f90 $(GENDIR)/%_2d.f90 $(GENDIR)/%_3d.f90

# Generic fypp rule (no dimension)
$(GENDIR)/%.f90: %.fpp | $(GENDIR)
	$(FYPP) $(FYPPFLAGS) $(FYPP_USER) $< $@

# 2D fypp rule
$(GENDIR)/%_2d.f90: %.fpp | $(GENDIR)
	$(FYPP) $(FYPPFLAGS) $(FYPP_USER) -D NDIM=2 $< $@

# 3D fypp rule
$(GENDIR)/%_3d.f90: %.fpp | $(GENDIR)
	$(FYPP) $(FYPPFLAGS) $(FYPP_USER) -D NDIM=3 $< $@

# ==============================================================================
# Fortran compilation rules
# ==============================================================================
$(OBJDIR)/%.o: $(GENDIR)/%.f90 | $(OBJDIR)
	$(F90C) -c -o $@ $< $(FFLAGS) $(INCFLAGS)

$(OBJDIR)/%.o: %.f90 | $(OBJDIR)
	$(F90C) -c -o $@ $< $(FFLAGS) $(INCFLAGS)

# ==============================================================================
# C compilation rules
# ==============================================================================
$(OBJDIR)/%_2d.o: $(COREDIR)/%.c | $(OBJDIR)
	$(CC) -c -o $@ $< $(CFLAGS) -DNDIM=2 $(INCFLAGS)

$(OBJDIR)/%_3d.o: $(COREDIR)/%.c | $(OBJDIR)
	$(CC) -c -o $@ $< $(CFLAGS) -DNDIM=3 $(INCFLAGS)

# ==============================================================================
# Linking rules
# ==============================================================================
LDFLAGS := $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

# Generic linking rule for 2D targets
$(BINDIR)/%_2d: $(OBJDIR)/%_2d.o | $(BINDIR)
	$(F90C) -o $@ $^ $(FFLAGS) $(LDFLAGS)

# Generic linking rule for 3D targets
$(BINDIR)/%_3d: $(OBJDIR)/%_3d.o | $(BINDIR)
	$(F90C) -o $@ $^ $(FFLAGS) $(LDFLAGS)

# ==============================================================================
# Special dependencies for 2D targets
# ==============================================================================
$(OBJDIR)/test_advection_2d.o: $(COMMON_OBJS_2D) $(OBJDIR)/m_physics_advection.o
$(BINDIR)/test_advection_2d: $(COMMON_OBJS_2D) $(OBJDIR)/m_physics_advection.o

$(OBJDIR)/test_euler_2d.o: $(COMMON_OBJS_2D) $(OBJDIR)/m_physics_euler_2d.o
$(BINDIR)/test_euler_2d: $(COMMON_OBJS_2D) $(OBJDIR)/m_physics_euler_2d.o

$(OBJDIR)/test_refinement_2d.o: $(COMMON_OBJS_2D)
$(BINDIR)/test_refinement_2d: $(COMMON_OBJS_2D)

$(OBJDIR)/test_xdmf_writer_2d.o: $(COMMON_OBJS_2D)
$(BINDIR)/test_xdmf_writer_2d: $(COMMON_OBJS_2D)

$(OBJDIR)/test_benchmark_ghostcell_2d.o: $(COMMON_OBJS_2D)
$(BINDIR)/test_benchmark_ghostcell_2d: $(COMMON_OBJS_2D)

# ==============================================================================
# Special dependencies for 3D targets
# ==============================================================================
$(OBJDIR)/test_advection_3d.o: $(COMMON_OBJS_3D) $(OBJDIR)/m_physics_advection.o
$(BINDIR)/test_advection_3d: $(COMMON_OBJS_3D) $(OBJDIR)/m_physics_advection.o

$(OBJDIR)/test_euler_3d.o: $(COMMON_OBJS_3D) $(OBJDIR)/m_physics_euler_3d.o
$(BINDIR)/test_euler_3d: $(COMMON_OBJS_3D) $(OBJDIR)/m_physics_euler_3d.o

$(OBJDIR)/test_refinement_3d.o: $(COMMON_OBJS_3D)
$(BINDIR)/test_refinement_3d: $(COMMON_OBJS_3D)

$(OBJDIR)/test_xdmf_writer_3d.o: $(COMMON_OBJS_3D)
$(BINDIR)/test_xdmf_writer_3d: $(COMMON_OBJS_3D)
