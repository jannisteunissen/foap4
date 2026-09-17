#:mute
! Central macro definitions for GPU/CPU parallelism.
! Provides OpenACC, OpenMP target, and OpenMP CPU variants.
! Macro overview:
!   GPU_IFDEF / GPU_ENDIF           - Conditional compilation guards for GPU code.
!   DEFAULT_PRESENT                 - Sets default data mapping (present/defaultmap).
!   COPYIN(varlist)                 - Maps variables to device with copy-in semantics.
!   DECLARE_DEVICE(varlist)         - Declares device-resident variables.
!   ENTER_DATA_COPYIN(varlist)      - Explicit data region: copyin.
!   ENTER_DATA_CREATE(varlist)      - Explicit data region: allocate on device.
!   EXIT_DATA_DELETE(varlist)       - Remove device data.
!   UPDATE_DEVICE(varlist)          - Copy host -> device.
!   UPDATE_SELF(varlist)            - Copy device -> host.
!   PARALLEL(clauses='')            - Begin a parallel region.
!   END_PARALLEL()                  - End a parallel region.
!   LOOP_INNER(clauses='')          - Inner loop parallelism.
!   LOOP_FLAT(clauses='')           - Flat loop parallelism.
!   LOOP_OUTER(clauses='')          - Outer loop parallelism.
!   PARALLEL_LOOP_FLAT(clauses='')  - Combined parallel+flat loop.
!   PARALLEL_LOOP_OUTER(clauses='') - Combined parallel+outer loop.
!   ATOMIC()                        - Atomic operation directive.
!   ROUTINE_SEQ()                   - Marks routine as sequential on device.
!   HOST_DATA_USE_DEVICE(varlist)   - Host-data region using device memory.
!   END_HOST_DATA()                 - End of host-data region.
#:if defined('USE_OPENACC')

#:def GPU_IFDEF()
#ifdef _OPENACC
#:enddef
#:def GPU_ENDIF()
#endif
#:enddef

#:def DEFAULT_PRESENT()
default(present)
#:enddef

#:def COPYIN(varlist)
copyin(${varlist}$)
#:enddef

#:def DECLARE_DEVICE(varlist)
!$acc declare create(${varlist}$)
#:enddef

#:def ENTER_DATA_COPYIN(varlist)
!$acc enter data copyin(${varlist}$)
#:enddef

#:def ENTER_DATA_CREATE(varlist)
!$acc enter data create(${varlist}$)
#:enddef

#:def EXIT_DATA_DELETE(varlist)
!$acc exit data delete(${varlist}$)
#:enddef

#:def UPDATE_DEVICE(varlist)
!$acc update device(${varlist}$)
#:enddef

#:def UPDATE_SELF(varlist)
!$acc update self(${varlist}$)
#:enddef

#:def PARALLEL(clauses='')
!$acc parallel ${clauses}$
#:enddef

#:def END_PARALLEL()
!$acc end parallel
#:enddef

#:def LOOP_INNER(clauses='')
!$acc loop independent ${clauses}$
#:enddef

#:def LOOP_FLAT(clauses='')
!$acc loop independent ${clauses}$
#:enddef

#:def LOOP_OUTER(clauses='')
!$acc loop independent ${clauses}$
#:enddef

#:def PARALLEL_LOOP_FLAT(clauses='')
!$acc parallel loop independent ${clauses}$
#:enddef

#:def PARALLEL_LOOP_OUTER(clauses='')
!$acc parallel loop independent ${clauses}$
#:enddef

#:def ATOMIC()
!$acc atomic
#:enddef

#:def ROUTINE_SEQ()
!$acc routine seq
#:enddef

#:def HOST_DATA_USE_DEVICE(varlist)
!$acc host_data use_device(${varlist}$)
#:enddef

#:def END_HOST_DATA()
!$acc end host_data
#:enddef

#:elif defined('USE_OPENMP')
#:def GPU_IFDEF()
#ifdef _OPENMP
#:enddef
#:def GPU_ENDIF()
#endif
#:enddef

#:def DEFAULT_PRESENT()
defaultmap(present)
#:enddef

#:def COPYIN(varlist)
map(to: ${varlist}$)
#:enddef

#:def DECLARE_DEVICE(varlist)
!$omp declare target(${varlist}$)
#:enddef

#:def ENTER_DATA_COPYIN(varlist)
!$omp target enter data map(to: ${varlist}$)
#:enddef

#:def ENTER_DATA_CREATE(varlist)
!$omp target enter data map(alloc: ${varlist}$)
#:enddef

#:def EXIT_DATA_DELETE(varlist)
!$omp target exit data map(delete: ${varlist}$)
#:enddef

#:def UPDATE_DEVICE(varlist)
!$omp target update to(${varlist}$)
#:enddef

#:def UPDATE_SELF(varlist)
!$omp target update from(${varlist}$)
#:enddef

#:def PARALLEL(clauses='')
!$omp target teams ${clauses}$
#:enddef

#:def END_PARALLEL()
!$omp end target teams
#:enddef

#:def LOOP_INNER(clauses='')
!$omp parallel do ${clauses}$
#:enddef

#:def LOOP_FLAT(clauses='')
!$omp distribute parallel do ${clauses}$
#:enddef

#:def LOOP_OUTER(clauses='')
!$omp distribute ${clauses}$
#:enddef

#:def PARALLEL_LOOP_FLAT(clauses='')
!$omp target teams distribute parallel do ${clauses}$
#:enddef

#:def PARALLEL_LOOP_OUTER(clauses='')
!$omp target teams distribute ${clauses}$
#:enddef

#:def ATOMIC()
!$omp atomic
#:enddef

#:def ROUTINE_SEQ()
!$omp declare target
#:enddef

#:def HOST_DATA_USE_DEVICE(varlist)
!$omp target data use_device_addr(${varlist}$)
#:enddef

#:def END_HOST_DATA()
!$omp end target data
#:enddef

#:elif defined('USE_OPENMP_CPU')
! No GPU offloading - emit nothing
#:def GPU_IFDEF()
#if 0
#:enddef
#:def GPU_ENDIF()
#endif
#:enddef

#:def DEFAULT_PRESENT()
#:enddef

#:def COPYIN(varlist)
#:enddef

#:def DECLARE_DEVICE(varlist)
#:enddef

#:def ENTER_DATA_COPYIN(varlist)
#:enddef

#:def ENTER_DATA_CREATE(varlist)
#:enddef

#:def EXIT_DATA_DELETE(varlist)
#:enddef

#:def UPDATE_DEVICE(varlist)
#:enddef

#:def UPDATE_SELF(varlist)
#:enddef

#:def PARALLEL(clauses='')
!$omp parallel ${clauses}$
#:enddef

#:def END_PARALLEL()
!$omp end parallel
#:enddef

! Don't use nested parallelism on CPU
#:def LOOP_INNER(clauses='')
#:enddef

#:def LOOP_FLAT(clauses='')
!$omp do ${clauses}$
#:enddef

#:def LOOP_OUTER(clauses='')
!$omp do ${clauses}$
#:enddef

#:def PARALLEL_LOOP_FLAT(clauses='')
!$omp parallel do ${clauses}$
#:enddef

#:def PARALLEL_LOOP_OUTER(clauses='')
!$omp parallel do ${clauses}$
#:enddef

#:def ATOMIC()
!$omp atomic
#:enddef

#:def ROUTINE_SEQ()
#:enddef

#:def HOST_DATA_USE_DEVICE(varlist)
#:enddef

#:def END_HOST_DATA()
#:enddef

#:else

#:def GPU_IFDEF()
#if 0
#:enddef
#:def GPU_ENDIF()
#endif
#:enddef
#:def DEFAULT_PRESENT()
#:enddef
#:def COPYIN(varlist)
#:enddef
#:def DECLARE_DEVICE(varlist)
#:enddef
#:def ENTER_DATA_COPYIN(varlist)
#:enddef
#:def ENTER_DATA_CREATE(varlist)
#:enddef
#:def EXIT_DATA_DELETE(varlist)
#:enddef
#:def UPDATE_DEVICE(varlist)
#:enddef
#:def UPDATE_SELF(varlist)
#:enddef
#:def PARALLEL(clauses='')
#:enddef
#:def END_PARALLEL()
#:enddef
#:def LOOP_INNER(clauses='')
#:enddef
#:def LOOP_FLAT(clauses='')
#:enddef
#:def LOOP_OUTER(clauses='')
#:enddef
#:def PARALLEL_LOOP_FLAT(clauses='')
#:enddef
#:def PARALLEL_LOOP_OUTER(clauses='')
#:enddef
#:def ATOMIC()
#:enddef
#:def ROUTINE_SEQ()
#:enddef
#:def HOST_DATA_USE_DEVICE(varlist)
#:enddef
#:def END_HOST_DATA()
#:enddef
#:endif

#:endmute
