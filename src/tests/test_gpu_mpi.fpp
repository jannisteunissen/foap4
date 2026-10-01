#:include 'definitions_parallel.fpp'
! Simple test for GPU-direct MPI communication. If the MPI library is not
! CUDA/GPU-aware, this will either produce wrong results, hang, or crash.
program test_gpu_mpi
  use m_foap4_${NDIM}$d
  use mpi_f08

  implicit none

  integer, parameter    :: n = 10000
  integer               :: mpirank, mpisize
  real(fp), allocatable :: send_buf(:), recv_buf(:)

  allocate(send_buf(n), recv_buf(n))
  ${ENTER_DATA_CREATE('send_buf, recv_buf')}$

  call MPI_Init()
  call MPI_Comm_rank(MPI_COMM_WORLD, mpirank)
  call MPI_Comm_size(MPI_COMM_WORLD, mpisize)

  call test_ring_with_wrapper(mpirank, mpisize, n, send_buf, recv_buf)
  call test_ring_no_wrapper(mpirank, mpisize, n, send_buf, recv_buf)

  call MPI_Finalize()

  ${EXIT_DATA_DELETE('send_buf, recv_buf')}$
  deallocate(send_buf, recv_buf)

contains

  subroutine test_ring_no_wrapper(mpirank, mpisize, n, send_buf, recv_buf)
    integer, intent(in)     :: mpirank, mpisize
    integer, intent(in)     :: n
    real(fp), intent(inout) :: send_buf(n), recv_buf(n)
    integer                 :: dest, source, tag
    integer                 :: i, n_errors, n_errors_global
    real(fp)                :: expected
    type(MPI_Request)       :: req_send, req_recv
    type(MPI_Request)       :: reqs(2)
    type(MPI_Status)        :: statuses(MPI_STATUS_SIZE)

    tag      = 123
    n_errors = 0

    if (mpisize < 2) then
       if (mpirank == 0) print *, "This test requires at least 2 MPI ranks"
       call MPI_Abort(MPI_COMM_WORLD, 1)
    end if

    if (mpirank == 0) print *, &
         "Testing GPU-direct MPI communication with", mpisize, "ranks"

    ! Fill the send buffer on the device
    ${PARALLEL_LOOP_FLAT('private(i)')}$ ${DEFAULT_PRESENT()}$
    do i = 1, n
       send_buf(i) = real(mpirank * n + i, fp)
    end do

    dest   = mod(mpirank + 1,            mpisize)
    source = mod(mpirank - 1 + mpisize,  mpisize)

    ${HOST_DATA_USE_DEVICE('recv_buf, send_buf')}$
    call MPI_Irecv(recv_buf, n, MPI_DOUBLE, source, tag, &
         MPI_COMM_WORLD, req_recv)
    call MPI_Isend(send_buf, n, MPI_DOUBLE, dest, tag, &
         MPI_COMM_WORLD, req_send)
    ${END_HOST_DATA()}$

    reqs(1) = req_send
    reqs(2) = req_recv
    call MPI_Waitall(2, reqs, statuses)

    ! Verify received data
    ${PARALLEL_LOOP_FLAT('private(expected) reduction(+:n_errors)')}$ ${DEFAULT_PRESENT()}$
    do i = 1, n
       expected = real(source * n + i, fp)
       if (abs(recv_buf(i) - expected) > 1.0e-12_fp * abs(expected)) then
          n_errors = n_errors + 1
       end if
    end do

    call MPI_Allreduce(n_errors, n_errors_global, 1, MPI_INTEGER, MPI_SUM, &
         MPI_COMM_WORLD)

    if (mpirank == 0) then
       if (n_errors_global == 0) then
          print *, "PASS"
       else
          print *, "FAIL: found", n_errors_global, "mismatching elements"
       end if
    end if

    if (n_errors_global > 0) then
       call MPI_Abort(MPI_COMM_WORLD, 1)
    end if
  end subroutine test_ring_no_wrapper

  subroutine test_ring_with_wrapper(mpirank, mpisize, n, send_buf, recv_buf)
    integer, intent(in)     :: mpirank, mpisize
    integer, intent(in)     :: n
    real(fp), intent(inout) :: send_buf(n), recv_buf(n)
    integer                 :: dest, source, tag
    integer                 :: i, n_errors, n_errors_global, nreqs
    real(fp)                :: expected
    type(MPI_Request)       :: requests(2)
    type(MPI_Status)        :: statuses(MPI_STATUS_SIZE)

    tag      = 456
    n_errors = 0
    nreqs    = 0

    if (mpisize < 2) then
       if (mpirank == 0) print *, "This test requires at least 2 MPI ranks"
       call MPI_Abort(MPI_COMM_WORLD, 1)
    end if

    if (mpirank == 0) print *, &
         "Testing GPU-direct MPI via wrappers with ", mpisize, "ranks"

    ${PARALLEL_LOOP_FLAT('private(i)')}$ ${DEFAULT_PRESENT()}$
    do i = 1, n
       send_buf(i) = real(mpirank * n + i, fp)
    end do

    dest   = mod(mpirank + 1,           mpisize)
    source = mod(mpirank - 1 + mpisize, mpisize)

    ! Wrappers append their request to requests and increase nreqs
    call f4_mpi_irecv_wrapper(recv_buf, int(n, MPI_COUNT_KIND), source, tag, &
         MPI_COMM_WORLD, requests, nreqs)
    call f4_mpi_isend_wrapper(send_buf, int(n, MPI_COUNT_KIND), dest, tag, &
         MPI_COMM_WORLD, requests, nreqs)

    call MPI_Waitall(nreqs, requests, statuses)

    ${PARALLEL_LOOP_FLAT('private(expected) reduction(+:n_errors)')}$ ${DEFAULT_PRESENT()}$
    do i = 1, n
       expected = real(source * n + i, fp)
       if (abs(recv_buf(i) - expected) > 1.0e-12_fp * abs(expected)) then
          n_errors = n_errors + 1
       end if
    end do

    call MPI_Allreduce(n_errors, n_errors_global, 1, MPI_INTEGER, MPI_SUM, &
         MPI_COMM_WORLD)

    if (mpirank == 0) then
       if (n_errors_global == 0) then
          print *, "PASS"
       else
          print *, "FAIL: found", n_errors_global, "mismatching elements"
       end if
    end if

    if (n_errors_global > 0) call MPI_Abort(MPI_COMM_WORLD, 1)

  end subroutine test_ring_with_wrapper

end program test_gpu_mpi
