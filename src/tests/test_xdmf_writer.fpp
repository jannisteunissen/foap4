program test_xdmf_writer_${NDIM}$d
  use mpi_f08
  use m_foap4_${NDIM}$d
  use m_io_${NDIM}$d

  implicit none
  integer :: ierr, mpisize

  ! Change this to paraview to view the output in paraview
  character(len=20) :: viewer = "visit"

  call MPI_init(ierr)

  call MPI_comm_size(MPI_COMM_WORLD, mpisize, ierr);
  if (mpisize /= 1) then
     call MPI_finalize(ierr)
     print *, "Number of MPI tasks should be 1"
     stop
  end if

#:if NDIM == 2
  call multi_block_test_2d("output/xdmf_test_single_2d_0gc", [1, 1], [16, 16], 0)
  call multi_block_test_2d("output/xdmf_test_single_2d_1gc", [1, 1], [16, 16], 1)
  call multi_block_test_2d("output/xdmf_test_multiple_2d_0gc", [4, 2], [8, 8], 0)
  call multi_block_test_2d("output/xdmf_test_multiple_2d_1gc", [4, 2], [8, 8], 1)
#:elif NDIM == 3
  call multi_block_test_3d("output/xdmf_test_single_3d_0gc", [1, 1, 1], [16, 16, 16], 0)
  call multi_block_test_3d("output/xdmf_test_single_3d_1gc", [1, 1, 1], [16, 16, 16], 1)
  call multi_block_test_3d("output/xdmf_test_multiple_3d_0gc", [4, 2, 2], [8, 8, 8], 0)
  call multi_block_test_3d("output/xdmf_test_multiple_3d_1gc", [4, 2, 2], [8, 8, 8], 1)
#:endif

  call MPI_finalize(ierr)

contains

#:if NDIM == 2
  subroutine multi_block_test_2d(fname, n_blocks_dim, nx, n_gc_out)
    character(len=*), intent(in) :: fname
    integer, intent(in)          :: n_blocks_dim(2)
    integer, intent(in)          :: nx(2)
    integer, intent(in)          :: n_gc_out
    integer                      :: n_blocks
    integer, parameter           :: n_cc        = 2
    character(len=10)            :: cc_names(n_cc) = ["rho", "phi"]
    integer, parameter           :: n_gc        = 1
    real(dp), allocatable        :: origin(:, :), dr(:, :)
    real(fp), allocatable        :: cc_data(:, :, :, :)
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
                cc_data(i, j, 1, i_block) = real(product(sin(rr * pi)), fp)
                cc_data(i, j, 2, i_block) = real(product(cos(rr * pi)), fp)
             end do
          end do
       end do
    end do

    call io_xdmf_write_blocks_2DCoRect(MPI_COMM_WORLD, trim(fname), n_blocks, &
         nx, n_cc, cc_names, n_gc, n_gc_out, origin, dr, cc_data, time=time, &
         viewer=viewer)

  end subroutine multi_block_test_2d
#:elif NDIM == 3
  subroutine multi_block_test_3d(fname, n_blocks_dim, nx, n_gc_out)
    character(len=*), intent(in) :: fname
    integer, intent(in)          :: n_blocks_dim(3)
    integer, intent(in)          :: nx(3)
    integer, intent(in)          :: n_gc_out
    integer                      :: n_blocks
    integer, parameter           :: n_cc        = 2
    character(len=10)            :: cc_names(n_cc) = ["rho", "phi"]
    integer, parameter           :: n_gc        = 1
    real(dp), allocatable        :: origin(:, :), dr(:, :)
    real(fp), allocatable        :: cc_data(:, :, :, :, :)
    real(dp), parameter          :: time        = 1.0_dp
    integer                      :: i, j, k, ii, jj, kk, i_block
    integer                      :: lo(3), hi(3)
    real(dp)                     :: rr(3)
    real(dp), parameter          :: pi = acos(-1.0_dp)
    real(dp), parameter          :: domain_length(3) = [1.0_dp, 1.0_dp, 1.0_dp]

    lo = 1 - n_gc
    hi = nx + n_gc

    n_blocks = product(n_blocks_dim)
    allocate(origin(3, n_blocks), dr(3, n_blocks))
    allocate(cc_data(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), n_cc, n_blocks))

    do kk = 1, n_blocks_dim(3)
       do jj = 1, n_blocks_dim(2)
          do ii = 1, n_blocks_dim(1)
             i_block = (kk - 1) * n_blocks_dim(2) * n_blocks_dim(1) + &
                  (jj - 1) * n_blocks_dim(1) + ii

             dr(:, i_block) = domain_length / (nx * n_blocks_dim)
             origin(:, i_block) = [ii-1, jj-1, kk-1] * dr(:, i_block) * nx

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      rr = origin(:, i_block) + &
                           [i-0.5_dp, j-0.5_dp, k-0.5_dp] * dr(:, i_block)
                      cc_data(i, j, k, 1, i_block) = real(product(sin(rr * pi)), fp)
                      cc_data(i, j, k, 2, i_block) = real(product(cos(rr * pi)), fp)
                   end do
                end do
             end do
          end do
       end do
    end do

    call io_xdmf_write_blocks_3DCoRect(MPI_COMM_WORLD, trim(fname), n_blocks, &
         nx, n_cc, cc_names, n_gc, n_gc_out, origin, dr, cc_data, time=time, &
         viewer=viewer)

  end subroutine multi_block_test_3d
#:endif

end program test_xdmf_writer_${NDIM}$d
