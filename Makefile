.SUFFIXES:

# ==============================================================================
# Compiler and tool definitions
# ==============================================================================
F90C ?= mpif90
CC := mpicc
FYPP := fypp
AR := ar
ARFLAGS := rcs

# ==============================================================================
# External libraries (paths not depending on build dirs)
# ==============================================================================
LIBDIRS := p4est/build/local/lib
LIBS := p4est sc z m

# ==============================================================================
# Floating point type for block data
# ==============================================================================
FLOAT_BITS ?= 64

# ==============================================================================
# GPU / offloading configuration
# One canonical variable OFFLOAD_KIND drives everything: acc | omp | ompcpu | none
# By default, use OpenACC (unless OpenMP or OPENMP_CPU is explicitly requested)
# ==============================================================================
ifeq ($(OPENMP),1)
    OFFLOAD_KIND := omp
else ifeq ($(OPENMP_CPU),1)
    OFFLOAD_KIND := ompcpu
else ifeq ($(OPENACC),0)
    OFFLOAD_KIND := none
else
    OFFLOAD_KIND := acc
endif

ifeq ($(OFFLOAD_KIND),acc)
    OFFLOAD_FYPP  := -D USE_OPENACC
    OFFLOAD_MODEL := "OpenACC GPU"
    OFFLOAD_TAG   := -acc
else ifeq ($(OFFLOAD_KIND),omp)
    OFFLOAD_FYPP  := -D USE_OPENMP
    OFFLOAD_MODEL := "OpenMP GPU"
    OFFLOAD_TAG   := -omp
else ifeq ($(OFFLOAD_KIND),ompcpu)
    OFFLOAD_FYPP  := -D USE_OPENMP_CPU
    OFFLOAD_MODEL := "OpenMP CPU"
    OFFLOAD_TAG   := -ompcpu
else
    OFFLOAD_FYPP  :=
    OFFLOAD_MODEL := "None"
    OFFLOAD_TAG   := -cpu
endif

# ==============================================================================
# Compiler detection
# ==============================================================================
compiler_version := $(shell $(F90C) --version 2>/dev/null)
compiler_brand ?= $(word 1, $(compiler_version))

# ==============================================================================
# Build configuration tag and build directory
# ==============================================================================
CONFIG_TAG := $(compiler_brand)-f$(FLOAT_BITS)$(OFFLOAD_TAG)
ifeq ($(DEBUG),1)
  CONFIG_TAG := $(CONFIG_TAG)-dbg
endif

BUILDDIR ?= build-$(CONFIG_TAG)

# ==============================================================================
# Directory structure
# ==============================================================================
OBJDIR := $(BUILDDIR)/obj
BINDIR := $(BUILDDIR)/bin
LIBDIR := $(BUILDDIR)/lib
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
# Include directories and flags (depend on directory structure)
# ==============================================================================
INCDIRS := p4est/build/local/include
INCDIRS += $(COREDIR) $(PHYSICSDIR) $(UTILSDIR) $(NUMERICSDIR) $(GENDIR)
INCFLAGS := $(addprefix -I,$(INCDIRS))

# ==============================================================================
# Preprocessor flags
# ==============================================================================
FYPPFLAGS := -n $(INCFLAGS) -D FLOAT_BITS=$(FLOAT_BITS) $(OFFLOAD_FYPP)

# ==============================================================================
# C compiler flags
# ==============================================================================
CFLAGS := -Wall -O2 -g

# ==============================================================================
# Fortran flags (depend on compiler_brand and OBJDIR)
# ==============================================================================
ifeq ($(compiler_brand),GNU)
    FFLAGS ?= -Wall -g -J$(OBJDIR) -cpp -Wno-unused-dummy-argument -Wl,--no-warn-execstack
    ifeq ($(DEBUG),1)
        FFLAGS += -O0 -fcheck=all
        CFLAGS += -O0
    endif
    ifeq ($(SAFE),1)
        FFLAGS += -O2
    else
        FFLAGS += -Ofast -march=native
    endif
    ifeq ($(OFFLOAD_KIND),acc)
        FFLAGS += -fopenacc -foffload=nvptx-none
        CFLAGS += -fopenacc -foffload=nvptx-none
    else ifneq (,$(filter $(OFFLOAD_KIND),omp ompcpu))
        FFLAGS += -fopenmp
        CFLAGS += -fopenmp
    endif

else ifneq (,$(filter $(compiler_brand),nvfortran pgfortran))
    NV_COMMON := -Minform=warn -fast -Mpreprocess -static-nvidia -g -module $(OBJDIR)
    ifeq ($(OFFLOAD_KIND),acc)
        FFLAGS ?= $(NV_COMMON) -acc=gpu,strict -gpu=ccnative
    else ifeq ($(OFFLOAD_KIND),omp)
        FFLAGS ?= $(NV_COMMON) -mp=gpu -gpu=ccnative
    else ifeq ($(OFFLOAD_KIND),ompcpu)
        FFLAGS ?= $(NV_COMMON) -mp
    else
        FFLAGS ?= $(NV_COMMON)
    endif

else ifeq ($(compiler_brand),Cray)
    FFLAGS ?= -M878 -O2 -eT -ea -ef -ffree -J$(OBJDIR)
    ifeq ($(OFFLOAD_KIND),acc)
        FFLAGS += -h acc -h acc_model=auto_async_none:fast_addr:no_deep_copy
    else ifneq (,$(filter $(OFFLOAD_KIND),omp ompcpu))
        FFLAGS += -h omp
    endif

else
    $(warning Unknown compiler "$(compiler_brand)", using default flags)
    FFLAGS ?= -O2 -g -J$(OBJDIR)
    ifneq ($(OFFLOAD_KIND),none)
        $(warning $(OFFLOAD_MODEL) flags unknown for this compiler; add manually via FFLAGS_USER)
    endif
endif

FFLAGS += $(FFLAGS_USER)

# ==============================================================================
# Library definitions
# ==============================================================================
LIB_2D := $(LIBDIR)/libfoap4_2d.a
LIB_3D := $(LIBDIR)/libfoap4_3d.a

BASE_OBJS := m_foap4_types
BASE_OBJS_2D := $(addprefix $(OBJDIR)/,$(addsuffix _2d.o,$(BASE_OBJS)))
BASE_OBJS_3D := $(addprefix $(OBJDIR)/,$(addsuffix _3d.o,$(BASE_OBJS)))

COMMON_OBJS := \
    m_foap4 \
    p4est_wrapper \
    m_io \
    m_rk \
    m_amr_flags \
    m_physics_advection \
    m_physics_euler
COMMON_OBJS_2D := $(addprefix $(OBJDIR)/,$(addsuffix _2d.o,$(COMMON_OBJS)))
COMMON_OBJS_3D := $(addprefix $(OBJDIR)/,$(addsuffix _3d.o,$(COMMON_OBJS)))

LIB_OBJS_2D := $(BASE_OBJS_2D) $(COMMON_OBJS_2D) \
               $(OBJDIR)/m_physics_shallow_water_2d.o \
               $(OBJDIR)/m_config.o

LIB_OBJS_3D := $(BASE_OBJS_3D) $(COMMON_OBJS_3D) $(OBJDIR)/m_config.o

# ==============================================================================
# Target definitions
# ==============================================================================
TARGETS_2D := $(addprefix $(BINDIR)/,\
    test_refinement_2d \
    test_advection_2d \
    test_xdmf_writer_2d \
    test_euler_2d \
    test_benchmark_ghostcell_2d \
    test_shallow_water_2d \
    test_euler_dmr_2d \
)

TARGETS_3D := $(addprefix $(BINDIR)/,\
    test_refinement_3d \
    test_advection_3d \
    test_xdmf_writer_3d \
    test_euler_3d \
)

# ==============================================================================
# Ensure limiter and flux schemes in NUMERICSDIR are preprocessed
# ==============================================================================
NUMERICS_FPP := $(wildcard $(NUMERICSDIR)/limiter_*.fpp) \
	$(wildcard $(NUMERICSDIR)/flux_scheme_*.fpp)
NUMERICS_GEN := $(patsubst $(NUMERICSDIR)/%.fpp,$(GENDIR)/%.f90,$(NUMERICS_FPP))

# ==============================================================================
# Phony targets
# ==============================================================================
.PHONY: all clean 2d 3d libs lib2d lib3d help build-summary

all: libs $(TARGETS_2D) $(TARGETS_3D)
	@$(MAKE) --no-print-directory build-summary

2d: lib2d $(TARGETS_2D)
	@$(MAKE) --no-print-directory build-summary

3d: lib3d $(TARGETS_3D)
	@$(MAKE) --no-print-directory build-summary

build-summary:
	@echo "=============================================="
	@echo "Build complete."
	@echo "  Build directory : $(BUILDDIR)"
	@echo "  Compiler        : $(compiler_brand)"
	@echo "  Offload model   : $(OFFLOAD_MODEL)"
	@echo "  Float bits      : $(FLOAT_BITS)"
	@echo "  Binaries in     : $(BINDIR)"
	@echo "=============================================="

libs: lib2d lib3d

lib2d: $(LIB_2D)

lib3d: $(LIB_3D)

help:
	@echo "Usage: make [target] [options]"
	@echo ""
	@echo "Targets:"
	@echo "  all      - Build all libraries and targets (default)"
	@echo "  2d       - Build 2D library and targets"
	@echo "  3d       - Build 3D library and targets"
	@echo "  libs     - Build both static libraries"
	@echo "  lib2d    - Build 2D static library ($(LIB_2D))"
	@echo "  lib3d    - Build 3D static library ($(LIB_3D))"
	@echo "  clean    - Remove all build artifacts in BUILDDIR"
	@echo "  help     - Show this help message"
	@echo ""
	@echo "Options:"
	@echo "  BUILDDIR=<dir>  - Set build directory (default: build)"
	@echo "  DEBUG=1         - Enable debug flags"
	@echo "  SAFE=1          - Use -O2 instead of -Ofast"
	@echo "  OPENACC=1       - Enable OpenACC GPU offloading (default; disabled with OPENMP)"
	@echo "  OPENMP=1        - Enable OpenMP target GPU offloading (disables OPENACC)"
	@echo "  OPENMP_CPU=1    - Enable OpenMP CPU threading (disables OPENACC)"
	@echo "  FLOAT_BITS=N    - Floating point precision for block data; 32 or 64"
	@echo "  F90C=<comp>     - Set Fortran compiler (default: mpif90)"
	@echo "  FFLAGS_USER=... - Additional Fortran flags appended to FFLAGS"
	@echo ""
	@echo "Examples:"
	@echo "  make 2d FLOAT_BITS=32"
	@echo "  make all DEBUG=1"
	@echo "  make 3d F90C=nvfortran OPENACC=1"
	@echo "  make 3d F90C=nvfortran OPENMP=1"
	@echo "  make BUILDDIR=build_debug DEBUG=1"
	@echo ""
	@echo "Detected compiler: $(compiler_brand)"
	@echo "Current BUILDDIR: $(BUILDDIR)"
	@echo "Current FLOAT_BITS: $(FLOAT_BITS)"
	@echo "Offload model: $(OFFLOAD_MODEL)"
	@echo "Current FFLAGS: $(FFLAGS)"


clean:
	$(RM) -r $(BUILDDIR)

# ==============================================================================
# Directory creation
# ==============================================================================
$(OBJDIR) $(BINDIR) $(GENDIR) $(LIBDIR):
	mkdir -p $@

# ==============================================================================
# Static library rules
# ==============================================================================
$(LIB_2D): $(LIB_OBJS_2D) | $(LIBDIR)
	$(AR) $(ARFLAGS) $@ $^

$(LIB_3D): $(LIB_OBJS_3D) | $(LIBDIR)
	$(AR) $(ARFLAGS) $@ $^

# ==============================================================================
# Fypp preprocessing rules
# ==============================================================================
.PRECIOUS: $(GENDIR)/%.f90 $(GENDIR)/%_2d.f90 $(GENDIR)/%_3d.f90

$(GENDIR)/%.f90: %.fpp | $(GENDIR)
	$(FYPP) $(FYPPFLAGS) $(FYPP_USER) $< $@

$(GENDIR)/%_2d.f90: %.fpp | $(GENDIR)
	$(FYPP) $(FYPPFLAGS) $(FYPP_USER) -D NDIM=2 $< $@

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

$(BINDIR)/%_2d: $(OBJDIR)/%_2d.o $(LIB_2D) | $(BINDIR)
	$(F90C) -o $@ $< $(LIB_2D) $(FFLAGS) $(LDFLAGS)

$(BINDIR)/%_3d: $(OBJDIR)/%_3d.o $(LIB_3D) | $(BINDIR)
	$(F90C) -o $@ $< $(LIB_3D) $(FFLAGS) $(LDFLAGS)

# ==============================================================================
# Dependencies for the library
# ==============================================================================
$(OBJDIR)/m_foap4_2d.o: $(OBJDIR)/m_foap4_types_2d.o
$(OBJDIR)/m_foap4_3d.o: $(OBJDIR)/m_foap4_types_3d.o

$(OBJDIR)/m_rk_2d.o: $(OBJDIR)/m_foap4_types_2d.o
$(OBJDIR)/m_rk_3d.o: $(OBJDIR)/m_foap4_types_3d.o

$(OBJDIR)/m_amr_flags_2d.o: $(OBJDIR)/m_foap4_2d.o
$(OBJDIR)/m_amr_flags_3d.o: $(OBJDIR)/m_foap4_3d.o

$(OBJDIR)/m_physics_advection_2d.o: $(OBJDIR)/m_foap4_types_2d.o
$(OBJDIR)/m_physics_advection_3d.o: $(OBJDIR)/m_foap4_types_3d.o

$(OBJDIR)/m_physics_euler_2d.o: $(OBJDIR)/m_foap4_types_2d.o
$(OBJDIR)/m_physics_euler_3d.o: $(OBJDIR)/m_foap4_types_3d.o

$(OBJDIR)/m_physics_shallow_water_2d.o: $(OBJDIR)/m_foap4_types_2d.o
$(OBJDIR)/m_physics_shallow_water_3d.o: $(OBJDIR)/m_foap4_types_3d.o

$(COMMON_OBJS_2D): $(BASE_OBJS_2D)
$(COMMON_OBJS_3D): $(BASE_OBJS_3D)

# ==============================================================================
# Dependencies for the targets
# ==============================================================================

# All 2D target objects depend on 2D library
$(patsubst $(BINDIR)/%,$(OBJDIR)/%.o,$(TARGETS_2D)): $(LIB_2D) $(NUMERICS_GEN)

# All 3D target objects depend on 3D library
$(patsubst $(BINDIR)/%,$(OBJDIR)/%.o,$(TARGETS_3D)): $(LIB_3D) $(NUMERICS_GEN)
