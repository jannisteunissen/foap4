program test_xdmf_writer
  use m_xdmf_writer
  use mpi_f08

  implicit none
  integer, parameter :: dp = kind(0.0d0)

  integer :: ierr, mpisize

  call MPI_init(ierr)

  call MPI_comm_size(MPI_COMM_WORLD, mpisize, ierr);
  if (mpisize /= 1) then
     call MPI_finalize(ierr)
     print *, "Number of MPI tasks should be 1"
     stop
  end if

  call multi_block_test("output/xdmf_test_single", [1, 1], [16, 16])
  call multi_block_test("output/xdmf_test_multiple", [4, 2], [8, 8])

  call MPI_finalize(ierr)

contains

  subroutine multi_block_test(fname, n_blocks_dim, nx)
    character(len=*), intent(in) :: fname
    integer, intent(in)          :: n_blocks_dim(2)
    integer, intent(in)          :: nx(2)
    integer                      :: n_blocks
    integer, parameter           :: n_cc        = 2
    character(len=10)            :: cc_names(n_cc) = ["rho", "phi"]
    integer, parameter           :: n_gc        = 1
    real(dp), allocatable        :: origin(:, :), dr(:, :)
    real(dp), allocatable        :: cc_data(:, :, :, :)
    real(dp), parameter          :: time        = 1.0_dp
    integer                      :: i, j, ii, jj, i_block
    integer                      :: lo(2), hi(2)
    real(dp)                     :: rr(2)
    real(dp), parameter          :: pi = acos(-1.0_dp)
    real(dp), parameter          :: domain_length(2) = [1.0_dp, 1.0_dp]

    lo = 1 - n_gc
    hi = nx + n_gc

    n_blocks = product(n_blocks_dim)
    allocate(origin(2, n_blocks), dr(2, n_blocks))
    allocate(cc_data(lo(1):hi(1), lo(2):hi(2), n_cc, n_blocks))

    do jj = 1, n_blocks_dim(2)
       do ii = 1, n_blocks_dim(1)
          i_block = (jj - 1) * n_blocks_dim(1) + ii

          dr(:, i_block) = domain_length / (nx * n_blocks_dim)
          origin(:, i_block) = [ii-1, jj-1] * dr(:, i_block) * nx

          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rr = origin(:, i_block) + [i-0.5_dp, j-0.5_dp] * dr(:, i_block)
                cc_data(i, j, 1, i_block) = product(sin(rr * pi))
                cc_data(i, j, 2, i_block) = product(cos(rr * pi))
             end do
          end do
       end do
    end do

    call xdmf_write_blocks_2DCoRect(MPI_COMM_WORLD, trim(fname), n_blocks, &
         nx+2*n_gc, n_cc, cc_names, n_gc, origin, dr, cc_data)

    ! call xdmf_write_blocks_2DCoRect(mpicomm, trim(full_fname), &
    !      f4%n_blocks, f4%bx+2*f4%n_gc, f4%n_vars, &
    !      f4%var_names, f4%n_gc, f4%block_origin(:, 1:f4%n_blocks), dr, &
    !      cc_data=f4%uu(:, :, :, 1:f4%n_blocks), time=time, viewer=viewer)

  end subroutine multi_block_test

end program test_xdmf_writer
