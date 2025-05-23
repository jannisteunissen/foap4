#:def GPU_DECLARE_COPYIN(*args)
#:if defined("USE_OPENACC")
!$acc declare copyin(${', '.join(args)}$)
#:elif defined("USE_OPENMP")
!$omp declare target(${', '.join(args)}$)
#:endif
#:enddef

#:def GPU_LOOP()
#:if defined("USE_OPENACC")
!$acc loop
#:elif defined("USE_OPENMP")
!$omp loop
#:else
! NO_GPU
#:endif
#:enddef

#:def GPU_PARALLEL_LOOP(nowait=False)
#:if defined("USE_OPENACC")
!$acc parallel loop#{if nowait}# async#{endif}#
#:elif defined("USE_OPENMP")
!$omp target teams loop#{if nowait}# nowait#{endif}#
#:else
! NO_GPU
#:endif
#:enddef

#:def GPU_WAIT_ASYNC()
#:if defined("USE_OPENACC")
!$acc wait
#:elif defined("USE_OPENMP")
!$omp taskwait
#:endif
#:enddef

#:def GPU_DATA_CREATE(*args)
#:if defined("USE_OPENACC")
!$acc enter data create(${', '.join(args)}$)
#:elif defined("USE_OPENMP")
!$omp target enter data map(alloc:${', '.join(args)}$)
#:endif
#:enddef

#:def GPU_DATA_COPYIN(*args)
#:if defined("USE_OPENACC")
!$acc enter data copyin(${', '.join(args)}$)
#:elif defined("USE_OPENMP")
!$omp target enter data map(to:${', '.join(args)}$)
#:endif
#:enddef

#:def GPU_UPDATE_DEVICE(*args)
#:if defined("USE_OPENACC")
!$acc update device(${', '.join(args)}$)
#:elif defined("USE_OPENMP")
!$omp target update to(${', '.join(args)}$)
#:endif
#:enddef

#:def GPU_UPDATE_HOST(*args)
#:if defined("USE_OPENACC")
!$acc update self(${', '.join(args)}$)
#:elif defined("USE_OPENMP")
!$omp target update from(${', '.join(args)}$)
#:endif
#:enddef

#:def GPU_DELETE(*args)
#:if defined("USE_OPENACC")
!$acc exit data delete(${', '.join(args)}$)
#:elif defined("USE_OPENMP")
!$omp target exit data map(delete:${', '.join(args)}$)
#:endif
#:enddef

#:def GPU_START_PARALLEL()
#:if defined("USE_OPENACC")
!$acc parallel
#:elif defined("USE_OPENMP")
!$omp target teams
#:endif
! NO_GPU
#:enddef

#:def GPU_END_PARALLEL()
#:if defined("USE_OPENACC")
!$acc end parallel
#:elif defined("USE_OPENMP")
!$omp end target teams
#:endif
#:enddef

#:def GPU_USE_DEVICE_PTR(*args)
#:if defined("USE_OPENACC")
!$acc host_data use_device(${', '.join(args)}$)
#:elif defined("USE_OPENMP")
!$omp target data use_device_ptr(${', '.join(args)}$)
#:endif
#:enddef

#:def GPU_END_USE_DEVICE_PTR()
#:if defined("USE_OPENACC")
!$acc end host_data
#:elif defined("USE_OPENMP")
!$omp end target data
#:endif
#:enddef

#:def GPU_ROUTINE_SEQ()
#:if defined("USE_OPENACC")
!$acc routine seq
#:endif
!$omp declare target
#:enddef
