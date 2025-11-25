.SUFFIXES:

FC := mpif90
CC := mpicc
FYPPFLAGS := -n
INCDIRS := p4est/build/local/include src/physics
LIBDIRS := p4est/build/local/lib
LIBS := p4est sc z m
CFLAGS := -Wall -O2 -g
TARGETS_2D := test_refinement_2d test_advection_2d test_xdmf_writer_2d test_euler_2d	\
	benchmark_ghostcell_2d
TARGETS_3D := test_refinement_3d test_advection_3d test_euler_3d

.PHONY: all
all: $(TARGETS_2D) $(TARGETS_3D)

# Determine compiler brand
compiler_version = $(shell $(FC) --version)
compiler_brand = $(word 1, $(compiler_version))

ifeq ($(compiler_brand), GNU)
	FFLAGS ?= -Wall -O2 -g -Jsrc -cpp $(FFLAGS_USER)	\
	-Wno-unused-dummy-argument -Wl,--no-warn-execstack
	ifeq ($(DEBUG), 1)
		FFLAGS += -O0 -fcheck=all -ffpe-trap=invalid,zero,overflow	\
		-finit-real=snan
		CFLAGS += -O0
	endif
else ifeq ($(compiler_brand), nvfortran)
	FFLAGS ?= -Minform=warn -acc=gpu -fast -gpu=ccnative -Mpreprocess	\
	-static-nvidia -g -module src $(FFLAGS_USER)
else ifeq ($(compiler_brand), pgfortran)
	FFLAGS ?= -Minform=warn -acc=gpu,verystrict -fast -gpu=ccall	\
	-Mpreprocess -static-nvidia -g -module src $(FFLAGS_USER)
endif

# Dependencies
$(TARGETS_2D): src/m_foap4_2d.o src/p4est_wrapper_2d.o src/m_xdmf_writer.o src/m_config.o
$(TARGETS_3D): src/m_foap4_3d.o src/p4est_wrapper_3d.o src/m_xdmf_writer.o src/m_config.o

$(addsuffix .o,$(addprefix src/,$(TARGETS_2D))): src/m_foap4_2d.mod src/m_config.mod
$(addsuffix .o,$(addprefix src/,$(TARGETS_3D))): src/m_foap4_3d.mod src/m_config.mod

src/test_euler_2d.o: src/physics/m_euler_2d.o src/physics/euler_2d.f90
src/test_euler_3d.o: src/physics/m_euler_3d.o
test_euler_2d: src/physics/m_euler_2d.o
test_euler_3d: src/physics/m_euler_3d.o

src/test_advection_2d.o: src/m_advection.mod
src/test_advection_3d.o: src/m_advection.mod
test_advection_2d: src/m_advection.o
test_advection_3d: src/m_advection.o

src/m_foap4_2d.o: src/m_xdmf_writer.mod src/p4est_wrapper_2d.o
src/m_foap4_3d.o: src/m_xdmf_writer.mod src/p4est_wrapper_3d.o

.PHONY: clean
clean:
	$(RM) $(TARGETS) src/*.o src/*.mod src/*.smod src/physics/*.o src/limiters/*.o

# How to get .o object files from .c source files
src/%_2d.o: src/%.c
	$(CC) -c -o $@ $(OBJ) $(CFLAGS) -DNDIM=2 $(addprefix -I,$(INCDIRS)) -o $@ $<
src/%_3d.o: src/%.c
	$(CC) -c -o $@ $(OBJ) $(CFLAGS) -DNDIM=3 $(addprefix -I,$(INCDIRS)) -o $@ $<

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

.PRECIOUS: src/%.f90 src/%_2d.f90 src/%_3d.f90
src/%.f90: src/%.fpp
	fypp $(FYPPFLAGS) $(FYPP_USER) $< $@
src/%_2d.f90: src/%.fpp
	fypp $(FYPPFLAGS) $(FYPP_USER) -D NDIM=2 $< $@
src/%_3d.f90: src/%.fpp
	fypp $(FYPPFLAGS) $(FYPP_USER) -D NDIM=3 $< $@
