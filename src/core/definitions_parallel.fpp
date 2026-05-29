#:mute

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

#:def LOOP(clauses='')
!$acc loop independent ${clauses}$
#:enddef

#:def LOOP_OUTER(clauses='')
!$acc loop independent ${clauses}$
#:enddef

#:def PARALLEL_LOOP(clauses='')
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

#:def WAIT()
!$acc wait
#:enddef

#:def ASYNC()
async
#:enddef

#:elif defined('USE_OPENMP')
#:def GPU_IFDEF()
#ifdef _OPENMP
#:enddef
#:def GPU_ENDIF()
#endif
#:enddef

#:def DEFAULT_PRESENT()
#! TODO: check if defaultmap(present) works
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

#:def LOOP(clauses='')
!$omp parallel do ${clauses}$
#:enddef

#:def LOOP_OUTER(clauses='')
!$omp distribute ${clauses}$
#:enddef

#:def PARALLEL_LOOP(clauses='')
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

#:def WAIT()
!$omp taskwait
#:enddef

#:def ASYNC()
nowait
#:enddef

#:else
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
#:enddef
#:def END_PARALLEL()
#:enddef
#:def LOOP(clauses='')
#:enddef
#:def LOOP_OUTER(clauses='')
#:enddef
#:def PARALLEL_LOOP(clauses='')
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
#:def WAIT()
#:enddef
#:def ASYNC()
#:enddef
#:endif

#:endmute
