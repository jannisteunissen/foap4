!> Foap4 stands for "Fortran OpenACC p4est", and implements MPI-parallel
!> quadtree/octree adaptive mesh refinement with support for OpenACC. Together
!> with p4est_wrapper.c, this module implements the required data structures and
!> methods.
!>
!> Author(s): Jannis Teunissen
#:include 'definitions.fpp'
module m_foap4_${NDIM}$d
  use mpi_f08
  use, intrinsic :: iso_c_binding
  use iso_fortran_env, only: error_unit
  use m_foap4_types_${NDIM}$d

  implicit none
  private

  ! Public types exported from m_foap4_types
  public :: ndim
  public :: dp, fp, f4_mpi_float
  public :: foap4_t
  public :: f4_bc_dirichlet, f4_bc_neumann
  public :: f4_bc_linear_extrap, f4_bc_fixed_value
  public :: f4_face_xlo, f4_face_xhi
  public :: f4_face_ylo, f4_face_yhi
#:if NDIM == 3
  public :: f4_face_zlo, f4_face_zhi
#:endif

  ! Public methods
  public :: f4_initialize
  public :: f4_destroy
  public :: f4_finalize
  public :: f4_construct_brick
  public :: f4_set_bc_scalar
  public :: f4_reset_wtime
  public :: f4_print_wtime
  public :: f4_get_mesh_revision
  public :: f4_get_global_highest_level
  public :: f4_get_num_local_blocks
  public :: f4_get_num_global_blocks
  public :: f4_cell_coord
  public :: f4_block_face_coord
  public :: f4_exchange_buffers
  public :: f4_update_ghostcells
  public :: f4_adjust_refinement
  public :: f4_partition
  public :: f4_get_load_imbalance
  public :: f4_fix_c2f_flux
  public :: f4_compute_sum
  public :: f4_compute_max

contains

  !> Initialize p4est and MPI
  subroutine f4_initialize(f4, p4est_log_type)
    type(foap4_t), intent(inout) :: f4
    character(len=*), intent(in) :: p4est_log_type
    integer                      :: log_level, ierr

    ! Different log levels defined in sc.h
    select case (p4est_log_type)
    case ("default")
       log_level = (-1) ! Selects the SC default threshold.
    case ("always")
       log_level = 0 ! Log absolutely everything.
    case ("trace")
       log_level = 1 ! Prefix file and line number.
    case ("debug")
       log_level = 2 ! Any information on the internal state.
    case ("verbose")
       log_level = 3 ! Information on conditions, decisions.
    case ("info")
       log_level = 4 ! Most relevant things a function is doing.
    case ("statistics")
       log_level = 5 ! Important for consistency/performance.
    case ("production")
       log_level = 6 ! A few lines at most for a major api function.
    case ("essential")
       log_level = 7 ! Log a few lines max (version info) per program.
    case ("error")
       log_level = 8 ! Log errors only.
    case default
       error stop "Unknow value for p4est_log_type (try default, error, ...)"
    end select

    call pw_initialize(f4%pw, f4%mpicomm%MPI_VAL, log_level)

    ! Set MPI rank and size
    call MPI_Comm_rank(f4%mpicomm, f4%mpirank, ierr)
    call MPI_Comm_size(f4%mpicomm, f4%mpisize, ierr)

    ! Set maximum number of MPI requests for MPI communication. There is at
    ! most one message between each pair of ranks, but messages can be split
    ! if they exceed 16 GB.
    f4%max_requests = 2 * f4%mpisize + 10

#ifdef _OPENACC
    call set_openacc_device(f4)
#endif

  end subroutine f4_initialize

#ifdef _OPENACC
  !> Set the device to be used by OpenACC
  !> Inspired by https://enccs.github.io/gpu-programming/10-multiple_gpu/
  subroutine set_openacc_device(f4)
    use openacc
    type(foap4_t), intent(in) :: f4
    type(MPI_Comm)            :: host_comm
    integer                   :: host_rank, host_size, ierr
    integer                   :: my_device, num_devices, tasks_per_device
    integer(acc_device_kind)  :: dev_type

    ! Split communicator into subgroups per node
    call MPI_Comm_split_type(f4%mpicomm, MPI_COMM_TYPE_SHARED, 0, &
                            MPI_INFO_NULL, host_comm, ierr)
    call MPI_Comm_rank(host_comm, host_rank, ierr)
    call MPI_Comm_size(host_comm, host_size, ierr)

    dev_type = acc_get_device_type()
    num_devices = acc_get_num_devices(dev_type)

    if (num_devices < 1) error stop "No devices available on host"

    if (num_devices < host_size) then
       if (f4%mpirank == 0) then
          write(*, "(A,I0,A,I0,A,I0,A)") "WARNING from ", f4%mpirank, &
               ": more local processes (", host_size, &
               ") than GPUs (", num_devices, ")"
       endif

       ! Round up
       tasks_per_device = (host_size + num_devices - 1)/num_devices

       ! Rank 0 to tasks_per_device-1 on device 0, etc.
       my_device = host_rank/tasks_per_device
    else
       my_device = host_rank
    endif

    call acc_set_device_num(my_device, dev_type)
    call acc_init(dev_type)

    call MPI_Comm_free(host_comm, ierr)

  end subroutine set_openacc_device
#endif

  !> Reset wall clock time measurements
  subroutine f4_reset_wtime(f4)
    type(foap4_t), intent(inout) :: f4

    f4%wtime_t0 = MPI_Wtime()
    f4%wtime_gc_fill_round1 = 0.0_dp
    f4%wtime_gc_fill_round2 = 0.0_dp
    f4%wtime_gc_fill_buff_round1 = 0.0_dp
    f4%wtime_gc_fill_buff_round2 = 0.0_dp
    f4%wtime_adjust_ref_p4est = 0.0_dp
    f4%wtime_adjust_ref_foap4 = 0.0_dp
    f4%wtime_partition = 0.0_dp
    f4%wtime_write_grid = 0.0_dp
    f4%wtime_update_gc_pattern = 0.0_dp
    f4%wtime_exchange_buffers = 0.0_dp
    f4%wtime_flux_fix = 0.0_dp
    f4%wtime_finite_volume = 0.0_dp
  end subroutine f4_reset_wtime

  !> Print wall clock time measurements
  subroutine f4_print_wtime(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: ierr
    integer, parameter           :: n_timers = 12
    real(dp)                     :: local_times(n_timers)
    real(dp)                     :: local_total
    real(dp)                     :: local_fracs(n_timers)
    real(dp)                     :: sum_fracs(n_timers), min_fracs(n_timers)
    real(dp)                     :: max_fracs(n_timers)
    real(dp)                     :: avg_frac
    real(dp)                     :: send_total(1), recv_total(1)
    real(dp)                     :: local_sum_frac(1)
    character(len=25)            :: timer_names(n_timers)
    integer                      :: i

    timer_names(1)  = "gc_fill_round1"
    timer_names(2)  = "gc_fill_round2"
    timer_names(3)  = "gc_fill_buff_round1"
    timer_names(4)  = "gc_fill_buff_round2"
    timer_names(5)  = "adjust_ref_p4est"
    timer_names(6)  = "adjust_ref_foap4"
    timer_names(7)  = "partition"
    timer_names(8)  = "write_grid"
    timer_names(9)  = "update_gc_pattern"
    timer_names(10) = "exchange_buffers"
    timer_names(11) = "flux_fix"
    timer_names(12) = "finite_volume"

    local_total = MPI_Wtime() - f4%wtime_t0

    local_times(1)  = f4%wtime_gc_fill_round1
    local_times(2)  = f4%wtime_gc_fill_round2
    local_times(3)  = f4%wtime_gc_fill_buff_round1
    local_times(4)  = f4%wtime_gc_fill_buff_round2
    local_times(5)  = f4%wtime_adjust_ref_p4est
    local_times(6)  = f4%wtime_adjust_ref_foap4
    local_times(7)  = f4%wtime_partition
    local_times(8)  = f4%wtime_write_grid
    local_times(9)  = f4%wtime_update_gc_pattern
    local_times(10) = f4%wtime_exchange_buffers
    local_times(11) = f4%wtime_flux_fix
    local_times(12) = f4%wtime_finite_volume

    if (local_total > 0.0_dp) then
       local_fracs = local_times * (1e2_dp / local_total)
    else
       local_fracs = 0.0_dp
    end if

    local_sum_frac(1) = sum(local_fracs)
    send_total(1) = local_total

    call MPI_Reduce(local_fracs,   sum_fracs,      n_timers, &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, f4%mpicomm, ierr)
    call MPI_Reduce(local_fracs,   min_fracs,      n_timers, &
         MPI_DOUBLE_PRECISION, MPI_MIN, 0, f4%mpicomm, ierr)
    call MPI_Reduce(local_fracs,   max_fracs,      n_timers, &
         MPI_DOUBLE_PRECISION, MPI_MAX, 0, f4%mpicomm, ierr)
    call MPI_Reduce(send_total,    recv_total,  1, &
         MPI_DOUBLE_PRECISION, MPI_MAX, 0, f4%mpicomm, ierr)

    if (f4%mpirank == 0) then
       write(*, "(A)") ""
       write(*, "(A,E15.6,A)") "total_runtime ", recv_total(1), " s"
       write(*, "(A25,3A11)") "timer                    ", "avg %", "min %", "max %"
       write(*, "(A25,3A11)") repeat("-", 25), ("-----------", i=1,3)
       do i = 1, n_timers
          avg_frac = sum_fracs(i) / f4%mpisize
          write(*, "(A25,3F11.3)") timer_names(i), avg_frac, min_fracs(i), max_fracs(i)
       end do
       write(*, "(A25,3F11.3)") "sum_of_above             ", &
            sum(sum_fracs) / f4%mpisize, sum(min_fracs), sum(max_fracs)
       write(*, "(A)") ""
    end if
  end subroutine f4_print_wtime

  !> Destroy all data for the current mesh
  subroutine f4_destroy(f4)
    type(foap4_t), intent(inout) :: f4

    ! OpenACC - Remove data from device

    ${EXIT_DATA_DELETE('f4%bc_simple_type, f4%bc_simple')}$
    ${EXIT_DATA_DELETE('f4%block_level, f4%block_origin')}$
    ${EXIT_DATA_DELETE('f4%uu, f4%uu_prim, f4%refinement_flags')}$
    ${EXIT_DATA_DELETE('f4%bc_data_ix, f4%bc_data, f4%bc_data_type, f4%bflux_ix, f4%bflux')}$
    ${EXIT_DATA_DELETE('f4%recv_buffer, f4%send_buffer')}$
    ${EXIT_DATA_DELETE('f4%gc_srl_local_iface, f4%gc_srl_from_buf_iface')}$
    ${EXIT_DATA_DELETE('f4%gc_srl_to_buf_iface, f4%gc_f2c_local_iface')}$
    ${EXIT_DATA_DELETE('f4%gc_f2c_from_buf_iface, f4%gc_f2c_to_buf_iface')}$
    ${EXIT_DATA_DELETE('f4%gc_c2f_from_buf_iface, f4%gc_c2f_to_buf_iface')}$
    ${EXIT_DATA_DELETE('f4%gc_phys_iface')}$
    ${EXIT_DATA_DELETE('f4')}$

    ! TODO: delete ghost cell patterns

    call pw_destroy(f4%pw)

    deallocate(f4%var_names)
    deallocate(f4%bc_simple_type)
    deallocate(f4%bc_simple)
    deallocate(f4%block_origin)
    deallocate(f4%block_level)
    deallocate(f4%refinement_flags)
    deallocate(f4%uu)
    deallocate(f4%uu_prim)
    deallocate(f4%bflux_ix)
    deallocate(f4%bflux)
    deallocate(f4%bc_data_ix)
    deallocate(f4%bc_data)
    deallocate(f4%bc_data_type)
    deallocate(f4%recv_buffer)
    deallocate(f4%send_buffer)
    deallocate(f4%recv_offset)
    deallocate(f4%send_offset)

    f4%gc_mesh_revision = -1
  end subroutine f4_destroy

  !> Finalize p4est and MPI
  subroutine f4_finalize(f4)
    type(foap4_t), intent(inout) :: f4
    call pw_finalize(f4%pw)
  end subroutine f4_finalize

  !> Construct a brick of blocks
  subroutine f4_construct_brick(f4, trees_per_dim, tree_length, bx, n_gc, &
       n_vars, var_names, var_temporal, n_temporal_states, periodic, &
       min_level, max_blocks, bc_type, bc_value)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: trees_per_dim(ndim) !< How many trees per dimension
    real(dp), intent(in)         :: tree_length(ndim)   !< Length of each tree
    integer, intent(in)          :: bx(ndim)            !< Number of cells in a block
    integer, intent(in)          :: n_gc !< Number of ghost cells
    integer, intent(in)          :: n_vars !< Number of variables
    !< Variable names. Temporal variables should be at the end.
    character(len=*), intent(in) :: var_names(n_vars)
    !> Whether the variables requires several temporal states
    logical, intent(in)          :: var_temporal(n_vars)
    !> Number of temporal states per variable
    integer, intent(in)          :: n_temporal_states
    logical, intent(in)          :: periodic(ndim) !< Periodic flag per dim
    integer, intent(in)          :: min_level   !< Refine up to this level
    integer, intent(in)          :: max_blocks  !< Maximum number of blocks
    integer, intent(in)          :: bc_type !< Default physical boundary type
    real(dp), intent(in)         :: bc_value !< Default physical boundary value
    integer                      :: i, periodic_as_int(ndim)
    real(dp)                     :: t0, t1

    ! Reset timers
    call f4_reset_wtime(f4)

    t0 = MPI_Wtime()

    if (any(bx /= bx(1))) error stop "TODO: unequal bx(:) not yet supported"
    if (any(bx < n_gc)) error stop "Cannot have any(bx < n_gc)"
    if (any(iand(bx, 1) == 1)) error stop "All bx have to be even"

    f4%bx          = bx
    f4%hbx         = bx/2
    f4%n_gc        = n_gc
    f4%ilo         = 1 - n_gc
    f4%ihi         = bx + n_gc
    f4%max_blocks  = max_blocks
    f4%tree_length = tree_length

    f4%n_vars = n_vars
    f4%n_temporal_states = n_temporal_states

    where (periodic)
       periodic_as_int = 1
    elsewhere
       periodic_as_int = 0
    end where

    allocate(f4%var_names(f4%n_vars))
    do i = 1, f4%n_vars
       f4%var_names(i) = var_names(i)
    end do

    allocate(f4%bc_simple_type(f4%n_vars, 0:2*ndim-1))
    allocate(f4%bc_simple(f4%n_vars, 0:2*ndim-1))
    f4%bc_simple_type(:, :) = bc_type
    f4%bc_simple(:, :) = real(bc_value, fp)

    do i = 0, P4EST_MAXLEVEL-1
       f4%dr_level(:, i) = (tree_length/bx) * 0.5**i
    end do

    call pw_set_connectivity_brick(f4%pw, trees_per_dim, periodic_as_int, &
         min_level, 1, max_blocks)

    ! When doing mesh refinement, the extra storage for temporal states is
    ! used as a copy buffer. The level is also used afterwards.
    allocate(f4%block_level(max_blocks * n_temporal_states))

    allocate(f4%block_origin(ndim, max_blocks))
    allocate(f4%refinement_flags(max_blocks))
    allocate(f4%bflux_ix(0:2*NDIM-1, max_blocks))
    allocate(f4%bc_data_ix(0:2*NDIM-1, max_blocks))

#:if NDIM == 2
    allocate(f4%uu(1-n_gc:bx(1)+n_gc, 1-n_gc:bx(2)+n_gc, &
         f4%n_vars, max_blocks * n_temporal_states))
    allocate(f4%uu_prim(1-n_gc:bx(1)+n_gc, 1-n_gc:bx(2)+n_gc, &
         f4%n_vars, max_blocks))
    allocate(f4%bflux(bx(1), n_vars, 0))
    allocate(f4%bc_data(bx(1), n_vars, 0))
    allocate(f4%bc_data_type(bx(1), n_vars, 0))
    f4%gc_data_size = f4%bx(1) * f4%n_gc
    f4%gc_data_size_c2f = (f4%bx(1)/2) * f4%n_gc
    f4%gc_data_size_fluxfix = f4%bx(1)/2
#:elif NDIM == 3
    allocate(f4%uu(1-n_gc:bx(1)+n_gc, 1-n_gc:bx(2)+n_gc, 1-n_gc:bx(3)+n_gc, &
         f4%n_vars, max_blocks * n_temporal_states))
    allocate(f4%uu_prim(1-n_gc:bx(1)+n_gc, 1-n_gc:bx(2)+n_gc, 1-n_gc:bx(3)+n_gc, &
         f4%n_vars, max_blocks))
    allocate(f4%bflux(bx(1), bx(1), n_vars, 0))
    allocate(f4%bc_data(bx(1), bx(1), n_vars, 0))
    allocate(f4%bc_data_type(bx(1), bx(1), n_vars, 0))
    f4%gc_data_size = f4%bx(1)**2 * f4%n_gc
    f4%gc_data_size_c2f = (f4%bx(1)/2)**2 * f4%n_gc
    f4%gc_data_size_fluxfix = (f4%bx(1)/2)**2
#:endif

    f4%uu = 0.0_dp
    f4%bflux = 0.0_dp

    ! Maximum size of recv/send buffer
    if (f4%mpisize > 1) then
       i = max_blocks * 2 * NDIM * f4%gc_data_size * f4%n_vars
    else
       i = 0
    end if

    allocate(f4%recv_buffer(i))
    allocate(f4%send_buffer(i))
    allocate(f4%recv_offset(0:f4%mpisize))
    allocate(f4%send_offset(0:f4%mpisize))

    ! OpenACC - Copy data structure and create allocatable components
    ${ENTER_DATA_COPYIN('f4')}$
    ${ENTER_DATA_COPYIN('f4%bc_simple_type, f4%bc_simple')}$
    ${ENTER_DATA_CREATE('f4%block_level, f4%block_origin')}$
    ${ENTER_DATA_CREATE('f4%uu, f4%uu_prim, f4%refinement_flags')}$
    ${ENTER_DATA_CREATE('f4%bc_data_ix, f4%bc_data, f4%bc_data_type')}$
    ${ENTER_DATA_CREATE('f4%bflux_ix, f4%bflux')}$
    ${ENTER_DATA_CREATE('f4%recv_buffer, f4%send_buffer')}$
    ${ENTER_DATA_CREATE('f4%gc_srl_local_iface, f4%gc_srl_from_buf_iface')}$
    ${ENTER_DATA_CREATE('f4%gc_srl_to_buf_iface, f4%gc_f2c_local_iface')}$
    ${ENTER_DATA_CREATE('f4%gc_f2c_from_buf_iface, f4%gc_f2c_to_buf_iface')}$
    ${ENTER_DATA_CREATE('f4%gc_c2f_from_buf_iface, f4%gc_c2f_to_buf_iface')}$
    ${ENTER_DATA_CREATE('f4%gc_phys_iface')}$

    call f4_set_quadrants(f4)
    call update_ghostcell_pattern(f4)
    call set_face_data_storage(f4)
    if (associated(f4%bc_callback)) call f4%bc_callback(f4)

    t1 = MPI_Wtime()

  end subroutine f4_construct_brick

  !> Allocate storage for fluxes at refinement boundaries and for boundary
  !> condition data, and set indices into this storage.
  !> Note: this routine requires that update_ghostcell_pattern has been called
  subroutine set_face_data_storage(f4)
    type(foap4_t), intent(inout) :: f4

    integer :: n, face, ix, i_block, i_coarse, i_fine, offset(NDIM-1)
    integer :: n_face_bc, n_face_rb

    ! Number of faces next to physical boundary
    n_face_bc = f4%gc_phys_iface(2*NDIM) - 1

    ! Number of faces next to refinement boundary. For the local ones, also
    ! add storage for the coarse side. Since it is not guaranteed that the
    ! fine sides are all on the same MPI rank, add an equal number to be safe
    n_face_rb = 2 * (f4%gc_f2c_local_iface(2*NDIM) - 1) + &
         f4%gc_f2c_to_buf_iface(2*NDIM) - 1 + &
         f4%gc_c2f_from_buf_iface(2*NDIM) - 1

    if (n_face_bc > size(f4%bc_data, NDIM+1)) then
       ! Resize storage, and reserve extra space
       ${EXIT_DATA_DELETE('f4%bc_data, f4%bc_data_type')}$
       deallocate(f4%bc_data, f4%bc_data_type)
#:if NDIM == 2
       allocate(f4%bc_data(f4%bx(1), f4%n_vars, 2*n_face_bc))
       allocate(f4%bc_data_type(f4%bx(1), f4%n_vars, 2*n_face_bc))
#:elif NDIM == 3
       allocate(f4%bc_data(f4%bx(1), f4%bx(1), f4%n_vars, 2*n_face_bc))
       allocate(f4%bc_data_type(f4%bx(1), f4%bx(1), f4%n_vars, 2*n_face_bc))
#:endif
       ${ENTER_DATA_CREATE('f4%bc_data, f4%bc_data_type')}$
    end if

    if (n_face_rb > size(f4%bflux, NDIM+1)) then
       ! Resize storage, and reserve extra space
       ${EXIT_DATA_DELETE('f4%bflux')}$
       deallocate(f4%bflux)
#:if NDIM == 2
       allocate(f4%bflux(f4%bx(1), f4%n_vars, 2*n_face_rb))
#:elif NDIM == 3
       allocate(f4%bflux(f4%bx(1), f4%bx(1), f4%n_vars, 2*n_face_rb))
#:endif
       ${ENTER_DATA_CREATE('f4%bflux')}$
    end if

    ! Set indices into boundary condition data. No boundary is set to zero,
    ! while boundaries are set to a negative index. This index should be
    ! turned positive where the boundary data should be used.
    ix = 0
    f4%bc_data_ix(:, 1:f4%n_blocks) = 0
    do face = 0, 2*NDIM-1
       do n = f4%gc_phys_iface(face), f4%gc_phys_iface(face+1)-1
          i_block = f4%gc_phys(n) + 1
          ix = ix + 1
          f4%bc_data_ix(face, i_block) = -ix
       end do
    end do

    ! Set indices into bflux array (local case)
    ix = 0
    f4%bflux_ix(:, 1:f4%n_blocks) = 0
    do face = 0, 2*NDIM-1
       do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
          i_fine = f4%gc_f2c_local(1, n) + 1
          i_coarse = f4%gc_f2c_local(2, n) + 1
#:if NDIM == 2
          offset(1) = f4%gc_f2c_local(3, n) ! offset
#:elif NDIM == 3
          offset(1:2) = f4%gc_f2c_local(3:4, n) ! offset
#:endif

          ix = ix + 1
          f4%bflux_ix(face, i_fine) = ix

          ! Add index for coarse block if not set already. Note that the face
          ! is reversed.
          if (f4%bflux_ix(face_swap(face), i_coarse) == 0) then
             ix = ix + 1
             f4%bflux_ix(face_swap(face), i_coarse) = ix
          end if
       end do
    end do

    ! Case where coarse block is on other MPI rank
    do face = 0, 2*NDIM-1
       do n = f4%gc_f2c_to_buf_iface(face), f4%gc_f2c_to_buf_iface(face+1)-1
          i_fine = f4%gc_f2c_to_buf(1, n) + 1
          ix = ix + 1
          f4%bflux_ix(face, i_fine) = ix
       end do
    end do

    ! Case where fine blocks are on other MPI rank
    do face = 0, 2*NDIM-1
       do n = f4%gc_c2f_from_buf_iface(face), f4%gc_c2f_from_buf_iface(face+1)-1
          i_coarse = f4%gc_c2f_from_buf(1, n) + 1
          ix = ix + 1
          f4%bflux_ix(face, i_coarse) = ix
       end do
    end do

    ${UPDATE_DEVICE('f4%bc_data_ix(:, 1:f4%n_blocks)')}$
    ${UPDATE_DEVICE('f4%bflux_ix(:, 1:f4%n_blocks)')}$

  end subroutine set_face_data_storage

  !> Set a default physical boundary conditions for a variable for a given
  !> face direction
  subroutine f4_set_bc_scalar(f4, ivar, iface, bc_type, bc_value)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: ivar
    integer, intent(in)          :: iface !< Range 0 - 2*NDIM-1
    integer, intent(in)          :: bc_type
    real(dp), intent(in)         :: bc_value

    if (iface < 0 .or. iface > 2*NDIM-1) &
         error stop "Must have 0 <= iface <= 2*NDIM-1"
    if (ivar < 1 .or. ivar > f4%n_vars) &
         error stop "Must have 1 <= ivar <= n_vars"
    if (bc_type > f4_bc_fixed_value) &
         error stop "Unsupported boundary condition type"

    f4%bc_simple_type(ivar, iface) = bc_type
    f4%bc_simple(ivar, iface) = real(bc_value, fp)
    ${UPDATE_DEVICE('f4%bc_simple_type(ivar, iface), f4%bc_simple(ivar, iface)')}$

  end subroutine f4_set_bc_scalar

  !> Return the mesh revision number
  pure integer function f4_get_mesh_revision(f4)
    type(foap4_t), intent(in) :: f4
    f4_get_mesh_revision = pw_get_mesh_revision(f4%pw)
  end function f4_get_mesh_revision

  !> Return the number of blocks on this MPI rank
  pure integer function f4_get_num_local_blocks(f4)
    type(foap4_t), intent(in) :: f4
    f4_get_num_local_blocks = pw_get_num_local_quadrants(f4%pw)
  end function f4_get_num_local_blocks

  !> Return the number of blocks on all MPI ranks together
  pure integer function f4_get_num_global_blocks(f4)
    type(foap4_t), intent(in) :: f4
    f4_get_num_global_blocks = pw_get_num_global_quadrants(f4%pw)
  end function f4_get_num_global_blocks

  !> Set the number of blocks, their origins and their refinement levels
  subroutine f4_set_quadrants(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, n_blocks

    if (.not. allocated(f4%block_origin)) error stop "block_origin not allocated"
    if (.not. allocated(f4%block_level)) error stop "block_level not allocated"

    n_blocks = f4_get_num_local_blocks(f4)
    if (n_blocks > f4%max_blocks) then
       write(error_unit, "(A,I0,A,I0)") "ERROR: n_blocks = ", n_blocks, &
            ", max_blocks = ", f4%max_blocks
       error stop "n_blocks > max_blocks"
    end if

    f4%n_blocks = n_blocks

    call pw_get_quadrants(f4%pw, f4%n_blocks, &
         f4%block_origin(:, 1:f4%n_blocks), &
         f4%block_level(1:f4%n_blocks))

    do n = 1, f4%n_blocks
       f4%block_origin(:, n) = f4%block_origin(:, n) * f4%tree_length
    end do

    ! OpenACC - synchronize block information to device
    ${UPDATE_DEVICE('f4%n_blocks, f4%block_origin(:, 1:f4%n_blocks)')}$
    ${UPDATE_DEVICE('f4%block_level(1:f4%n_blocks)')}$
  end subroutine f4_set_quadrants

  !> Get the global highest refinement level (on all MPI ranks)
  subroutine f4_get_global_highest_level(f4, max_level)
    type(foap4_t), intent(in) :: f4
    integer, intent(out)      :: max_level
    integer                   :: ierror

    max_level = pw_get_highest_local_level(f4%pw)

    call MPI_Allreduce(MPI_IN_PLACE, max_level, 1, MPI_INTEGER, &
         MPI_MAX, f4%mpicomm, ierror)
  end subroutine f4_get_global_highest_level


  !> Return the coordinates at the center of a grid cell
  pure function f4_cell_coord(f4, i_block, ${IJK}$) result(rr)
    @{ROUTINE_SEQ()}@
    type(foap4_t), intent(in) :: f4
    integer, intent(in)       :: i_block, ${IJK}$
    real(dp)                  :: rr(NDIM), dr(NDIM)

    dr = f4%dr_level(:, f4%block_level(i_block))
    rr(1) = f4%block_origin(1, i_block) + dr(1) * (i - 0.5_dp)
    rr(2) = f4%block_origin(2, i_block) + dr(2) * (j - 0.5_dp)
#:if NDIM == 3
    rr(3) = f4%block_origin(3, i_block) + dr(3) * (k - 0.5_dp)
#:endif
  end function f4_cell_coord

  !> Return the coordinates at the center of a cell face on a block
#:if NDIM == 2
  pure function f4_block_face_coord(f4, i_block, face, i) result(rr)
    @{ROUTINE_SEQ()}@
    type(foap4_t), intent(in) :: f4
    integer, intent(in)       :: i_block
    integer, intent(in)       :: face
    integer, intent(in)       :: i
    real(dp)                  :: rr(NDIM), dr(NDIM)

    dr = f4%dr_level(:, f4%block_level(i_block))
    rr = f4%block_origin(:, i_block)

    select case (face)
    case (0)
       rr(2) = rr(2) + dr(2) * (i - 0.5_dp)
    case (1)
       rr(1) = rr(1) + dr(1) * f4%bx(1)
       rr(2) = rr(2) + dr(2) * (i - 0.5_dp)
    case (2)
       rr(1) = rr(1) + dr(1) * (i - 0.5_dp)
    case (3)
       rr(1) = rr(1) + dr(1) * (i - 0.5_dp)
       rr(2) = rr(2) + dr(2) * f4%bx(2)
    end select
  end function f4_block_face_coord
#:elif NDIM == 3
  pure function f4_block_face_coord(f4, i_block, face, i, j) result(rr)
    @{ROUTINE_SEQ()}@
    type(foap4_t), intent(in) :: f4
    integer, intent(in)       :: i_block
    integer, intent(in)       :: face
    integer, intent(in)       :: i, j
    real(dp)                  :: rr(NDIM), dr(NDIM)

    dr = f4%dr_level(:, f4%block_level(i_block))
    rr = f4%block_origin(:, i_block)

    select case (face)
    case (0)
       rr(2) = rr(2) + dr(2) * (i - 0.5_dp)
       rr(3) = rr(3) + dr(3) * (j - 0.5_dp)
    case (1)
       rr(1) = rr(1) + dr(1) * f4%bx(1)
       rr(2) = rr(2) + dr(2) * (i - 0.5_dp)
       rr(3) = rr(3) + dr(3) * (j - 0.5_dp)
    case (2)
       rr(1) = rr(1) + dr(1) * (i - 0.5_dp)
       rr(3) = rr(3) + dr(3) * (j - 0.5_dp)
    case (3)
       rr(1) = rr(1) + dr(1) * (i - 0.5_dp)
       rr(2) = rr(2) + dr(2) * f4%bx(2)
       rr(3) = rr(3) + dr(3) * (j - 0.5_dp)
    case (4)
       rr(1) = rr(1) + dr(1) * (i - 0.5_dp)
       rr(2) = rr(2) + dr(2) * (j - 0.5_dp)
    case (5)
       rr(1) = rr(1) + dr(1) * (i - 0.5_dp)
       rr(2) = rr(2) + dr(2) * (j - 0.5_dp)
       rr(3) = rr(3) + dr(3) * f4%bx(3)
    end select
  end function f4_block_face_coord
#:endif

  !> Update the information required to update ghost cells
  subroutine update_ghostcell_pattern(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: mesh_revision
    type(bnd_face_t), pointer    :: bnd_face(:)
    type(c_ptr)                  :: tmp
    integer                      :: n_faces
    real(dp)                     :: t0, t1

    mesh_revision = pw_get_mesh_revision(f4%pw)
    if (mesh_revision == f4%gc_mesh_revision) return

    t0 = MPI_Wtime()
    call pw_get_all_faces(f4%pw, n_faces, tmp)

    call c_f_pointer(tmp, bnd_face, shape=[n_faces])

    call set_ghost_cell_pattern(f4, size(bnd_face), bnd_face, &
         f4%mpirank, f4%mpisize)
    f4%gc_mesh_revision = mesh_revision
    t1 = MPI_Wtime()
    f4%wtime_update_gc_pattern = f4%wtime_update_gc_pattern + t1 - t0
  end subroutine update_ghostcell_pattern

  !> Store the information required to update ghost cells
  subroutine set_ghost_cell_pattern(f4, n_faces, bnd_face, mpirank, mpisize)
    type(foap4_t), intent(inout)    :: f4
    integer, intent(in)             :: n_faces
    type(bnd_face_t), intent(inout) :: bnd_face(n_faces)
    integer, intent(in)             :: mpirank
    integer, intent(in)             :: mpisize

    integer :: n, rank, i, j
    integer :: i_phys
    integer :: i_same(0:mpisize-1)
    integer :: i_c2f(0:mpisize-1)
    integer :: i_f2c(0:mpisize-1)
    integer :: i_buf_recv(0:mpisize-1)
    integer :: i_buf_send(0:mpisize-1)
    integer :: i_srl_nonlocal
    integer :: i_f2c_to_buf, i_f2c_from_buf
    integer :: i_c2f_to_buf, i_c2f_from_buf

    type(int_array_t) :: same_ix(0:mpisize-1)
    type(int_array_t) :: c2f_ix(0:mpisize-1)
    type(int_array_t) :: f2c_ix(0:mpisize-1)
    type(int_array_t) :: phys_ix
    type(int_array_t) :: ix_send

    type(int_array_t) :: all_srl_from_buf
    type(int_array_t) :: all_srl_to_buf
    type(int_array_t) :: all_f2c_from_buf
    type(int_array_t) :: all_f2c_to_buf
    type(int_array_t) :: all_c2f_from_buf
    type(int_array_t) :: all_c2f_to_buf

    logical, parameter :: recv = .true., send = .false.

    i_same = 0
    i_c2f  = 0
    i_f2c  = 0
    i_phys = 0

    ! Count different communication patterns
    do n = 1, n_faces
       rank = bnd_face(n)%other_proc

       if (bnd_face(n)%face_type == FACE_SAME_LEVEL) then
          i_same(rank) = i_same(rank) + 1
       else if (bnd_face(n)%face_type == FACE_BOUNDARY) then
          i_phys = i_phys + 1
       else if (bnd_face(n)%face_type == FACE_COARSE_TO_FINE) then
          i_c2f(rank) = i_c2f(rank) + 1
       else if (bnd_face(n)%face_type == FACE_FINE_TO_COARSE) then
          i_f2c(rank) = i_f2c(rank) + 1
       else
          error stop "Unknown face type"
       end if
    end do

    ! Determine (cumulatively) how much to receive from and send to each rank
    if (.not. allocated(f4%gc_recv_offset)) then
       allocate(f4%gc_recv_offset(0:mpisize))
       allocate(f4%gc_send_offset(0:mpisize))
       allocate(f4%gc_recv_offset_c2f(0:mpisize))
       allocate(f4%gc_send_offset_c2f(0:mpisize))
       allocate(f4%gc_recv_offset_fluxfix(0:mpisize))
       allocate(f4%gc_send_offset_fluxfix(0:mpisize))
    end if

    f4%gc_recv_offset(0) = 0
    f4%gc_send_offset(0) = 0
    f4%gc_recv_offset_c2f(0) = 0
    f4%gc_send_offset_c2f(0) = 0
    f4%gc_recv_offset_fluxfix(0) = 0
    f4%gc_send_offset_fluxfix(0) = 0

    ! Determine offsets for sending and receiving with other ranks
    do rank = 1, mpisize
       if (rank-1 == mpirank) then
          f4%gc_recv_offset(rank) = f4%gc_recv_offset(rank-1)
          f4%gc_send_offset(rank) = f4%gc_send_offset(rank-1)

          f4%gc_recv_offset_c2f(rank) = f4%gc_recv_offset_c2f(rank-1)
          f4%gc_send_offset_c2f(rank) = f4%gc_send_offset_c2f(rank-1)

          f4%gc_recv_offset_fluxfix(rank) = f4%gc_recv_offset_fluxfix(rank-1)
          f4%gc_send_offset_fluxfix(rank) = f4%gc_send_offset_fluxfix(rank-1)
       else
          f4%gc_recv_offset(rank) = f4%gc_recv_offset(rank-1) + &
               f4%gc_data_size * i_same(rank-1) + &
               f4%gc_data_size_c2f * i_c2f(rank-1)
          f4%gc_send_offset(rank) = f4%gc_send_offset(rank-1) + &
               f4%gc_data_size * i_same(rank-1) + &
               f4%gc_data_size_c2f * i_f2c(rank-1)

          ! In a second round of communication, handle the fine side of
          ! refinement boundaries
          f4%gc_recv_offset_c2f(rank) = f4%gc_recv_offset_c2f(rank-1) + &
               f4%gc_data_size * i_f2c(rank-1)
          f4%gc_send_offset_c2f(rank) = f4%gc_send_offset_c2f(rank-1) + &
               f4%gc_data_size * i_c2f(rank-1)

          f4%gc_recv_offset_fluxfix(rank) = f4%gc_recv_offset_fluxfix(rank-1) + &
               f4%gc_data_size_fluxfix * i_c2f(rank-1)
          f4%gc_send_offset_fluxfix(rank) = f4%gc_send_offset_fluxfix(rank-1) + &
               f4%gc_data_size_fluxfix * i_f2c(rank-1)
       end if
    end do

    if (allocated(f4%gc_srl_local)) then
       ! OpenACC - deallocate arrays
       ${EXIT_DATA_DELETE('f4%gc_srl_local, f4%gc_srl_from_buf, f4%gc_srl_to_buf')}$
       ${EXIT_DATA_DELETE('f4%gc_f2c_local, f4%gc_f2c_from_buf, f4%gc_f2c_to_buf')}$
       ${EXIT_DATA_DELETE('f4%gc_c2f_from_buf, f4%gc_c2f_to_buf, f4%gc_phys')}$
       ${EXIT_DATA_DELETE('f4%gc_f2c_to_buf_fluxfix, f4%gc_c2f_from_buf_fluxfix')}$

       deallocate(f4%gc_srl_local, f4%gc_srl_from_buf, f4%gc_srl_to_buf, &
            f4%gc_f2c_local, f4%gc_f2c_from_buf, f4%gc_f2c_to_buf, &
            f4%gc_c2f_from_buf, f4%gc_c2f_to_buf, f4%gc_phys, &
            f4%gc_f2c_to_buf_fluxfix, f4%gc_c2f_from_buf_fluxfix)
    end if

    ! Local ghost cell exchange at the same level
    allocate(f4%gc_srl_local(2, i_same(mpirank)))

    ! Local ghost cell exchange at refinement boundaries
    allocate(f4%gc_f2c_local(2+NDIM-1, i_f2c(mpirank)))

    ! Physical boundaries
    allocate(f4%gc_phys(i_phys))
    allocate(phys_ix%i(i_phys))

    ! To store indices for different types of face boundaries
    do rank = 0, mpisize - 1
       allocate(same_ix(rank)%i(i_same(rank)))
       allocate(c2f_ix(rank)%i(i_c2f(rank)))
       allocate(f2c_ix(rank)%i(i_f2c(rank)))
    end do

    i_same = 0
    i_phys = 0
    i_c2f  = 0
    i_f2c  = 0

    ! Store indices of different face types
    do n = 1, n_faces
       rank = bnd_face(n)%other_proc

       if (bnd_face(n)%face_type == FACE_SAME_LEVEL) then
          i_same(rank) = i_same(rank) + 1
          same_ix(rank)%i(i_same(rank)) = n
       else if (bnd_face(n)%face_type == FACE_BOUNDARY) then
          i_phys = i_phys + 1
          phys_ix%i(i_phys) = n
       else if (bnd_face(n)%face_type == FACE_COARSE_TO_FINE) then
          i_c2f(rank) = i_c2f(rank) + 1
          c2f_ix(rank)%i(i_c2f(rank)) = n
       else if (bnd_face(n)%face_type == FACE_FINE_TO_COARSE) then
          i_f2c(rank) = i_f2c(rank) + 1
          f2c_ix(rank)%i(i_f2c(rank)) = n
       else
          error stop "Unknown face type"
       end if
    end do

    ! Sort local faces by face direction
    call sort_by_face(phys_ix, n_faces, bnd_face, f4%gc_phys_iface)
    call sort_by_face(same_ix(mpirank), n_faces, bnd_face, f4%gc_srl_local_iface)
    call sort_by_face(f2c_ix(mpirank), n_faces, bnd_face, f4%gc_f2c_local_iface)

    do n = 1, i_phys
       i = phys_ix%i(n)
       f4%gc_phys(n) = bnd_face(i)%quadid(1)
    end do

    do n = 1, i_same(mpirank)
       i = same_ix(mpirank)%i(n)
       f4%gc_srl_local(:, n) = bnd_face(i)%quadid(1:2)
    end do

    do n = 1, i_f2c(mpirank)
       i = f2c_ix(mpirank)%i(n)
       f4%gc_f2c_local(:, n) = [bnd_face(i)%quadid(1:2), bnd_face(i)%offset]
    end do

    ! Convert faces 0 to 1 and 2 to 3
    do n = f4%gc_srl_local_iface(0), f4%gc_srl_local_iface(1)-1
       f4%gc_srl_local(:, n) = f4%gc_srl_local([2, 1], n)
    end do
    do n = f4%gc_srl_local_iface(2), f4%gc_srl_local_iface(3)-1
       f4%gc_srl_local(:, n) = f4%gc_srl_local([2, 1], n)
    end do
    f4%gc_srl_local_iface(1) = f4%gc_srl_local_iface(0)
    f4%gc_srl_local_iface(3) = f4%gc_srl_local_iface(2)

#:if NDIM == 3
    ! Convert faces 4 to 5
    do n = f4%gc_srl_local_iface(4), f4%gc_srl_local_iface(5)-1
       f4%gc_srl_local(:, n) = f4%gc_srl_local([2, 1], n)
    end do
    f4%gc_srl_local_iface(5) = f4%gc_srl_local_iface(4)
#:endif

    ! Non-local ghost cell exchange at the same level
    n = sum(i_same) - i_same(mpirank)
    allocate(f4%gc_srl_from_buf(2, n))
    allocate(f4%gc_srl_to_buf(2, n))
    allocate(all_srl_from_buf%i(n))
    allocate(all_srl_to_buf%i(n))

    ! Non-local ghost cell exchange from fine to coarse
    n = sum(i_f2c) - i_f2c(mpirank)
    allocate(f4%gc_f2c_from_buf(2, n))
    allocate(f4%gc_f2c_to_buf(2, n))
    allocate(f4%gc_f2c_to_buf_fluxfix(2, n))
    allocate(all_f2c_from_buf%i(n))
    allocate(all_f2c_to_buf%i(n))

    ! Non-local ghost cell exchange from coarse to fine
    n = sum(i_c2f) - i_c2f(mpirank)
    allocate(f4%gc_c2f_from_buf(2+NDIM-1, n))
    allocate(f4%gc_c2f_from_buf_fluxfix(2+NDIM-1, n))
    allocate(f4%gc_c2f_to_buf(2+NDIM-1, n))
    allocate(all_c2f_from_buf%i(n))
    allocate(all_c2f_to_buf%i(n))

    i_srl_nonlocal = 0
    i_f2c_to_buf = 0
    i_f2c_from_buf = 0
    i_c2f_to_buf = 0
    i_c2f_from_buf = 0

    ! Determine buffer locations for each send and receive with another rank.
    ! For each type of face, sort the data exchanged with another rank by
    ! block, face and offset (in case of hanging faces).
    do rank = 0, mpisize - 1
       if (rank == mpirank) cycle

       ! Incremental offsets for sending and receiving data
       i_buf_recv(rank) = f4%gc_recv_offset(rank)
       i_buf_send(rank) = f4%gc_send_offset(rank)

       ! Boundaries at the same level
       ix_send = same_ix(rank)
       call sort_for_recv_or_send(same_ix(rank), n_faces, bnd_face, recv)
       call sort_for_recv_or_send(ix_send, n_faces, bnd_face, send)

       do n = 1, size(same_ix(rank)%i)
          i = same_ix(rank)%i(n) ! Index in bnd_face array
          j = ix_send%i(n)       ! Index in bnd_face array
          i_srl_nonlocal = i_srl_nonlocal + 1
          all_srl_from_buf%i(i_srl_nonlocal) = i
          all_srl_to_buf%i(i_srl_nonlocal) = j

          bnd_face(i)%ibuf_recv = i_buf_recv(rank)
          bnd_face(j)%ibuf_send = i_buf_send(rank)
          i_buf_recv(rank) = i_buf_recv(rank) + f4%gc_data_size
          i_buf_send(rank) = i_buf_send(rank) + f4%gc_data_size
       end do

       ! Sending from fine to coarse
       call sort_for_recv_or_send(f2c_ix(rank), n_faces, bnd_face, send)

       do n = 1, size(f2c_ix(rank)%i)
          i = f2c_ix(rank)%i(n) ! Index in bnd_face array
          i_f2c_to_buf = i_f2c_to_buf + 1
          all_f2c_to_buf%i(i_f2c_to_buf) = i

          bnd_face(i)%ibuf_send = i_buf_send(rank)
          i_buf_send(rank) = i_buf_send(rank) + f4%gc_data_size_c2f
       end do

       ! Receiving coarse from fine
       call sort_for_recv_or_send(c2f_ix(rank), n_faces, bnd_face, recv)

       do n = 1, size(c2f_ix(rank)%i)
          i = c2f_ix(rank)%i(n) ! Index in bnd_face array
          i_c2f_from_buf = i_c2f_from_buf + 1
          all_c2f_from_buf%i(i_c2f_from_buf) = i

          bnd_face(i)%ibuf_recv = i_buf_recv(rank)
          i_buf_recv(rank) = i_buf_recv(rank) + f4%gc_data_size_c2f
       end do

       ! After the above ghost cells have been updated, we can handle the fine
       ! side of refinement boundaries. This involves a new round of
       ! communication, so reset buffer offsets.
       i_buf_recv(rank) = f4%gc_recv_offset_c2f(rank)
       i_buf_send(rank) = f4%gc_send_offset_c2f(rank)

       ! Sending from coarse to fine
       call sort_for_recv_or_send(c2f_ix(rank), n_faces, bnd_face, send)

       do n = 1, size(c2f_ix(rank)%i)
          i = c2f_ix(rank)%i(n) ! Index in bnd_face array
          i_c2f_to_buf = i_c2f_to_buf + 1
          all_c2f_to_buf%i(i_c2f_to_buf) = i

          bnd_face(i)%ibuf_send = i_buf_send(rank)
          i_buf_send(rank) = i_buf_send(rank) + f4%gc_data_size
       end do

       ! Receiving fine from coarse
       call sort_for_recv_or_send(f2c_ix(rank), n_faces, bnd_face, recv)

       do n = 1, size(f2c_ix(rank)%i)
          i = f2c_ix(rank)%i(n) ! Index in bnd_face array
          i_f2c_from_buf = i_f2c_from_buf + 1
          all_f2c_from_buf%i(i_f2c_from_buf) = i

          bnd_face(i)%ibuf_recv = i_buf_recv(rank)
          i_buf_recv(rank) = i_buf_recv(rank) + f4%gc_data_size
       end do

    end do

    ! The buffer locations have been determined, sort over *all* ranks by face
    call sort_by_face(all_srl_from_buf, n_faces, bnd_face, f4%gc_srl_from_buf_iface)
    call sort_by_face(all_srl_to_buf, n_faces, bnd_face, f4%gc_srl_to_buf_iface)
    call sort_by_face(all_f2c_from_buf, n_faces, bnd_face, f4%gc_f2c_from_buf_iface)
    call sort_by_face(all_f2c_to_buf, n_faces, bnd_face, f4%gc_f2c_to_buf_iface)
    call sort_by_face(all_c2f_from_buf, n_faces, bnd_face, f4%gc_c2f_from_buf_iface)
    call sort_by_face(all_c2f_to_buf, n_faces, bnd_face, f4%gc_c2f_to_buf_iface)

    ! Extract ghost cell patterns
    do n = 1, i_srl_nonlocal
       i = all_srl_from_buf%i(n)
       j = all_srl_to_buf%i(n)
       f4%gc_srl_from_buf(:, n) = [bnd_face(i)%quadid(1), bnd_face(i)%ibuf_recv]
       f4%gc_srl_to_buf(:, n) = [bnd_face(j)%quadid(1), bnd_face(j)%ibuf_send]
    end do

    do n = 1, i_f2c_to_buf
       i = all_f2c_to_buf%i(n)
       f4%gc_f2c_to_buf(:, n) = [bnd_face(i)%quadid(1), bnd_face(i)%ibuf_send]
    end do

    do n = 1, i_f2c_from_buf
       i = all_f2c_from_buf%i(n)
       f4%gc_f2c_from_buf(:, n) = [bnd_face(i)%quadid(1), bnd_face(i)%ibuf_recv]
    end do

    do n = 1, i_c2f_from_buf
       i = all_c2f_from_buf%i(n)
       f4%gc_c2f_from_buf(:, n) = [bnd_face(i)%quadid(1), &
            bnd_face(i)%offset, bnd_face(i)%ibuf_recv]
    end do

    do n = 1, i_c2f_to_buf
       i = all_c2f_to_buf%i(n)
       f4%gc_c2f_to_buf(:, n) = [bnd_face(i)%quadid(1), &
            bnd_face(i)%offset, bnd_face(i)%ibuf_send]
    end do

    ! Repeat the above once more, but now for flux fixing. The code below
    ! could be merged with the code above if we could store two buffer indices
    ! in the bnd_face array.
    do rank = 0, mpisize - 1
       if (rank == mpirank) cycle

       ! Incremental offsets for sending and receiving data
       i_buf_recv(rank) = f4%gc_recv_offset_fluxfix(rank)
       i_buf_send(rank) = f4%gc_send_offset_fluxfix(rank)

       ! Sending from fine to coarse
       call sort_for_recv_or_send(f2c_ix(rank), n_faces, bnd_face, send)

       do n = 1, size(f2c_ix(rank)%i)
          i = f2c_ix(rank)%i(n) ! Index in bnd_face array
          bnd_face(i)%ibuf_send = i_buf_send(rank)
          i_buf_send(rank) = i_buf_send(rank) + f4%gc_data_size_fluxfix
       end do

       ! Receiving coarse from fine
       call sort_for_recv_or_send(c2f_ix(rank), n_faces, bnd_face, recv)

       do n = 1, size(c2f_ix(rank)%i)
          i = c2f_ix(rank)%i(n) ! Index in bnd_face array
          bnd_face(i)%ibuf_recv = i_buf_recv(rank)
          i_buf_recv(rank) = i_buf_recv(rank) + f4%gc_data_size_fluxfix
       end do
    end do

    do n = 1, i_f2c_to_buf
       i = all_f2c_to_buf%i(n)
       f4%gc_f2c_to_buf_fluxfix(:, n) = [bnd_face(i)%quadid(1), bnd_face(i)%ibuf_send]
    end do

    do n = 1, i_c2f_from_buf
       i = all_c2f_from_buf%i(n)
       f4%gc_c2f_from_buf_fluxfix(:, n) = [bnd_face(i)%quadid(1), &
            bnd_face(i)%offset, bnd_face(i)%ibuf_recv]
    end do

    ! OpenACC - copy/sync data to device

    ${ENTER_DATA_COPYIN('f4%gc_srl_local, f4%gc_srl_from_buf, f4%gc_srl_to_buf')}$
    ${ENTER_DATA_COPYIN('f4%gc_f2c_local, f4%gc_f2c_from_buf, f4%gc_f2c_to_buf')}$
    ${ENTER_DATA_COPYIN('f4%gc_c2f_from_buf, f4%gc_c2f_to_buf, f4%gc_phys')}$
    ${ENTER_DATA_COPYIN('f4%gc_f2c_to_buf_fluxfix, f4%gc_c2f_from_buf_fluxfix')}$

    ${UPDATE_DEVICE('f4%gc_srl_local_iface, f4%gc_srl_from_buf_iface')}$
    ${UPDATE_DEVICE('f4%gc_srl_to_buf_iface, f4%gc_f2c_local_iface')}$
    ${UPDATE_DEVICE('f4%gc_f2c_from_buf_iface, f4%gc_f2c_to_buf_iface')}$
    ${UPDATE_DEVICE('f4%gc_c2f_from_buf_iface, f4%gc_c2f_to_buf_iface')}$
    ${UPDATE_DEVICE('f4%gc_phys_iface')}$

  end subroutine set_ghost_cell_pattern

  !> Sort index array by face direction
  subroutine sort_by_face(ix, n_bnd_face, bnd_face, iface)
    type(int_array_t), intent(inout) :: ix !< Index array
    integer, intent(in)              :: n_bnd_face
    type(bnd_face_t), intent(in)     :: bnd_face(n_bnd_face)
    !> Start index for each face direction
    integer, intent(out)             :: iface(0:2*NDIM)
    integer                          :: face_count(0:2*NDIM-1)
    integer                          :: face_offset(0:2*NDIM-1)
    type(int_array_t)                :: ix_sorted
    integer                          :: n, face

    allocate(ix_sorted%i(size(ix%i)))

    ! Count number of each face type
    face_count(:) = 0
    do n = 1, size(ix%i)
       face = bnd_face(ix%i(n))%face
       face_count(face) = face_count(face) + 1
    end do

    ! Determine initial offset in index array
    face_offset(0) = 0
    do n = 1, 2*NDIM-1
       face_offset(n) = face_offset(n-1) + face_count(n-1)
    end do

    do n = 1, size(ix%i)
       face = bnd_face(ix%i(n))%face
       face_offset(face) = face_offset(face) + 1
       ix_sorted%i(face_offset(face)) = ix%i(n)
    end do

    ix%i(:) = ix_sorted%i
    call face_count_to_iface(face_count, iface)

  end subroutine sort_by_face

  !> Determine start index for each face direction
  subroutine face_count_to_iface(face_count, iface)
    integer, intent(in)  :: face_count(0:2*NDIM-1)
    integer, intent(out) :: iface(0:2*NDIM)
    integer              :: n

    iface(0) = 1
    do n = 1, 2*NDIM
       iface(n) = iface(n-1) + face_count(n-1)
    end do
  end subroutine face_count_to_iface

  !> Sort face boundaries for receiving or sending
  subroutine sort_for_recv_or_send(ix, n_bnd_face, bnd_face, recv)
    type(int_array_t), intent(inout) :: ix
    integer, intent(in)              :: n_bnd_face
    type(bnd_face_t), intent(in)     :: bnd_face(n_bnd_face)
    logical, intent(in)              :: recv
    integer, parameter               :: QSORT_THRESHOLD = 32
    integer                          :: array_size

    ! Fast in-line QSORT+INSERTION SORT for Fortran.
    ! Author: Joseph M. Krahn

    ! Generate a custom array sort procedure for a specific type,
    ! without the comparison-callback overhead of a generic sort procedure.
    ! This is essentially the same as an in-line optimization, which generally
    ! is not feasible for a library-based generic sort procedure.
    !
    ! NOTES:
    ! The procedure uses a optimized combination of QSORT and INSERTION
    ! sorting. The algorithm is based on code used in GLIBC. 
    ! A stack is used in place of recursive calls. The stack size must
    ! be at least as big as the number of bits in the largest array index.
    !
    ! Sorting vectors of a multidimensional allocatable array can be
    ! VERY slow. In this case, or with large derived types, it is better
    ! to sort a simple derived type of key/index pairs, then reorder
    ! tha actual data using the sorted indices.
    !
    !---------------------------------------------------------------------
    integer :: stack_top, right_size, left_size
    integer :: mid, left, right, low, high

    ! A stack of 32 can handle the entire extent of a 32-bit
    ! index, so this value is fixed. If you have 64-bit indexed
    ! arrays, which might contain more thant 2^32 elements, this
    ! should be set to 64.
    integer, parameter :: QSORT_STACK_SIZE = 32
    type qsort_stack
       integer :: low, high
    end type qsort_stack

    type(qsort_stack) :: stack(QSORT_STACK_SIZE)

    call init()

    if (array_size > QSORT_THRESHOLD) then
       low = 1
       high = array_size
       stack_top = 0

       QSORT_LOOP: do
          mid = (low + high)/2
          if (LESS_THAN (mid, low)) then
             call SWAP(mid,low)
          end if
          if (LESS_THAN (high, mid)) then
             call SWAP(high,mid)
             if (LESS_THAN (mid, low)) then
                call SWAP(mid,low)
             end if
          end if
          left  = low + 1
          right = high - 1

          COLLAPSE_WALLS: do
             do while (LESS_THAN (left, mid))
                left=left+1
             end do
             do while (LESS_THAN (mid, right))
                right=right-1
             end do
             if (left < right) then
                call SWAP(left,right)
                if (mid == left) then
                   mid = right
                else if (mid == right) then
                   mid = left
                end if
                left=left+1
                right=right-1
             else
                if (left == right) then
                   left=left+1
                   right=right-1
                end if
                exit COLLAPSE_WALLS
             end if
          end do COLLAPSE_WALLS

          ! Set up indices for the next iteration.
          ! Determine left and right partition sizes.
          ! Defer partitions smaller than the QSORT_THRESHOLD.
          ! If both partitions are significant,
          ! push the larger one onto the stack.
          right_size = right - low
          left_size = high - left
          if (right_size <= QSORT_THRESHOLD) then
             if (left_size <= QSORT_THRESHOLD) then
                ! Ignore both small partitions: Pop a partition or exit.
                if (stack_top<1) exit QSORT_LOOP
                low=stack(stack_top)%low; high=stack(stack_top)%high
                stack_top=stack_top-1
             else
                ! Ignore small left partition.
                low = left
             end if
          else if (left_size <= QSORT_THRESHOLD) then
             ! Ignore small right partition.
             high = right
          else if (right_size > left_size) then
             ! Push larger left partition indices.
             stack_top=stack_top+1
             stack(stack_top)=qsort_stack(low,right)
             low = left
          else
             ! Push larger right partition indices.
             stack_top=stack_top+1
             stack(stack_top)=qsort_stack(left,high)
             high = right
          end if
       end do QSORT_LOOP
    end if ! (array_size > QSORT_THRESHOLD)

    ! Sort the remaining small partitions using insertion sort,
    ! which should be faster for partitions smaller than the
    ! appropriate QSORT_THRESHOLD.

    ! First, find smallest element in first QSORT_THRESHOLD and
    ! place it at the array's beginning. This places a lower
    ! bound 'guard' position, and speeds up the inner loop
    ! below, because it will not need a lower-bound test.
    low = 1
    high = array_size

    ! left is the MIN_LOC index here:
    left=low
    do right = low+1, min(low+QSORT_THRESHOLD,high)
       if (LESS_THAN(right,left)) left=right
    end do
    if (left/=low) call SWAP(left,low)

    ! Insertion sort, from left to right.
    ! (assuming that the left is the lowest numbered index)
    INSERTION_SORT: do right = low+2,high
       left=right-1
       if (LESS_THAN(right,left)) then
          do while (LESS_THAN(right,left-1))
             left=left-1
          end do
          call RSHIFT(left,right)
       end if
    end do INSERTION_SORT

  contains

    pure logical function less_than(ia, ib)
      integer, intent(in) :: ia, ib
      integer             :: a, b

      a = ix%i(ia)
      b = ix%i(ib)

      if (recv) then
         ! Order by face and quadid of 'this' side
         if (bnd_face(a)%face /= bnd_face(b)%face) then
            less_than = (bnd_face(a)%face < bnd_face(b)%face)
         else if (bnd_face(a)%quadid(1) /= bnd_face(b)%quadid(1)) then
            less_than = (bnd_face(a)%quadid(1) < bnd_face(b)%quadid(1))
#:if NDIM == 2
         else
            ! For hanging faces, sort by offset
            less_than = (bnd_face(a)%offset(1) < bnd_face(b)%offset(1))
#:elif NDIM == 3
         else if (bnd_face(a)%offset(1) /= bnd_face(b)%offset(1)) then
            ! For hanging faces, sort by offset
            less_than = (bnd_face(a)%offset(1) < bnd_face(b)%offset(1))
         else
            less_than = (bnd_face(a)%offset(2) < bnd_face(b)%offset(2))
#:endif
         end if
      else
         ! Order by face and quadid of 'other' side
         if (bnd_face(a)%face /= bnd_face(b)%face) then
            ! Note that the face is swapped for the receiving side
            less_than = (face_swap(bnd_face(a)%face) < &
                 face_swap(bnd_face(b)%face))
         else if (bnd_face(a)%quadid(2) /= bnd_face(b)%quadid(2)) then
            less_than = (bnd_face(a)%quadid(2) < bnd_face(b)%quadid(2))
#:if NDIM == 2
         else
            less_than = (bnd_face(a)%offset(1) < bnd_face(b)%offset(1))
#:elif NDIM == 3
         else if (bnd_face(a)%offset(1) /= bnd_face(b)%offset(1)) then
            less_than = (bnd_face(a)%offset(1) < bnd_face(b)%offset(1))
         else
            less_than = (bnd_face(a)%offset(2) < bnd_face(b)%offset(2))
#:endif
         end if
      end if
    end function less_than

    subroutine init()
      array_size = size(ix%i)
    end subroutine init

    subroutine swap(a, b)
      integer, intent(in) :: a, b
      ix%i([a, b]) = ix%i([b, a])
    end subroutine swap

    subroutine rshift(left, right)
      integer, intent(in) :: left, right
      integer             :: tmp

      tmp                = ix%i(right)
      ix%i(left+1:right) = ix%i(left:right-1)
      ix%i(left)         = tmp
    end subroutine rshift

  end subroutine sort_for_recv_or_send

  !> Fill buffers for round one of the ghost cell exchange, which excludes the
  !> fine side of coarse-to-fine boundaries
  subroutine fill_ghostcell_buffers_round_one(f4, n_vars, i_vars, block_offset)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer, intent(in)          :: block_offset
    integer                      :: n, ivar, iv
    integer                      :: iq, i_buf, i_buf0
#:if NDIM   == 2
    integer                      :: i, j, i_f, j_f
#:elif NDIM == 3
    integer                      :: i, j, k, i_f, j_f, k_f
#:endif

    if (maxval(f4%gc_send_offset) * n_vars > size(f4%send_buffer)) then
       write(error_unit, *) "maxval(f4%gc_send_offset): ", maxval(f4%gc_send_offset)
       write(error_unit, *) "n_vars: ", n_vars
       write(error_unit, *) "size(f4%send_buffer): ", size(f4%send_buffer)
       error stop "send buffer too small"
    end if

#:def fyp_srl_to_buf(face, ilim, jlim, klim=None, i0=0, j0=0, k0=0)
    ${LOOP('collapse(NDIM+2) private(iq, i_buf0, ivar, i_buf)')}$
    do n = f4%gc_srl_to_buf_iface(${face}$), f4%gc_srl_to_buf_iface(${face}$+1)-1
#:if NDIM == 2
       do iv = 1, n_vars
          do j = 1, ${jlim}$
             do i = 1, ${ilim}$
                iq = f4%gc_srl_to_buf(1, n) + 1 + block_offset
                i_buf0 = f4%gc_srl_to_buf(2, n) * n_vars
                ivar = i_vars(iv)
                i_buf = i_buf0 + ix_offset3(iv, j, i, ${jlim}$, ${ilim}$) + 1
                f4%send_buffer(i_buf) = f4%uu(${i0}$+i, ${j0}$+j, ivar, iq)
             end do
          end do
       end do
#:elif NDIM == 3
       do iv = 1, n_vars
          do k = 1, ${klim}$
             do j = 1, ${jlim}$
                do i = 1, ${ilim}$
                   iq = f4%gc_srl_to_buf(1, n) + 1 + block_offset
                   i_buf0 = f4%gc_srl_to_buf(2, n) * n_vars
                   ivar = i_vars(iv)
                   i_buf = i_buf0 + ix_offset4(iv, k, j, i, &
                        ${klim}$, ${jlim}$, ${ilim}$) + 1
                   f4%send_buffer(i_buf) = f4%uu(${i0}$+i, ${j0}$+j, ${k0}$+k, ivar, iq)
                end do
             end do
          end do
       end do
#:endif
    end do
#:enddef

    ! Nonlocal fine-to-coarse boundaries, fill buffer for coarse side

#:def fyp_f2c_to_buf(face, ilim, jlim, klim=None, i0=0, j0=0, k0=0)
#:if NDIM == 2
    ${LOOP('collapse(NDIM+1) private(ivar, j_f, i_f, i_buf, iq, i_buf0)')}$
    do n = f4%gc_f2c_to_buf_iface(${face}$), f4%gc_f2c_to_buf_iface(${face}$+1)-1
       do iv = 1, n_vars
          do j = 1, ${jlim}$
             do i = 1, ${ilim}$
                ivar = i_vars(iv)
                iq = f4%gc_f2c_to_buf(1, n) + 1 + block_offset ! fine block
                i_buf0 = f4%gc_f2c_to_buf(2, n) * n_vars

                j_f = ${j0}$ + 2 * j - 1
                i_f = ${i0}$ + 2 * i - 1
                i_buf = i_buf0 + ix_offset3(iv, j, i, ${jlim}$, ${ilim}$) + 1
                f4%send_buffer(i_buf) = 0.25_fp * ( &
                     f4%uu(i_f,   j_f,   ivar, iq) + &
                     f4%uu(i_f+1, j_f,   ivar, iq) + &
                     f4%uu(i_f,   j_f+1, ivar, iq) + &
                     f4%uu(i_f+1, j_f+1, ivar, iq) )
             end do
          end do
       end do
    end do
#:elif NDIM == 3
    ${LOOP('collapse(NDIM+1) private(ivar, k_f, j_f, i_f, i_buf, iq, i_buf0)')}$
    do n = f4%gc_f2c_to_buf_iface(${face}$), f4%gc_f2c_to_buf_iface(${face}$+1)-1
       do iv = 1, n_vars
          do k = 1, ${klim}$
             do j = 1, ${jlim}$
                do i = 1, ${ilim}$
                   ivar = i_vars(iv)
                   iq = f4%gc_f2c_to_buf(1, n) + 1 + block_offset ! fine block
                   i_buf0 = f4%gc_f2c_to_buf(2, n) * n_vars

                   k_f = ${k0}$ + 2 * k - 1
                   j_f = ${j0}$ + 2 * j - 1
                   i_f = ${i0}$ + 2 * i - 1
                   i_buf = i_buf0 + ix_offset4(iv, k, j, i, &
                        ${klim}$, ${jlim}$, ${ilim}$) + 1
                   f4%send_buffer(i_buf) = 0.125_fp * ( &
                        f4%uu(i_f,   j_f,   k_f,   ivar, iq) + &
                        f4%uu(i_f+1, j_f,   k_f,   ivar, iq) + &
                        f4%uu(i_f,   j_f+1, k_f,   ivar, iq) + &
                        f4%uu(i_f+1, j_f+1, k_f,   ivar, iq) + &
                        f4%uu(i_f,   j_f,   k_f+1, ivar, iq) + &
                        f4%uu(i_f+1, j_f,   k_f+1, ivar, iq) + &
                        f4%uu(i_f,   j_f+1, k_f+1, ivar, iq) + &
                        f4%uu(i_f+1, j_f+1, k_f+1, ivar, iq) )
                end do
             end do
          end do
       end do
    end do
#:endif
#:enddef

    ${PARALLEL(COPYIN('i_vars'))}$ ${DEFAULT_PRESENT()}$

#:if NDIM == 2
    @:fyp_srl_to_buf(0, f4%n_gc, f4%bx(2))
    @:fyp_srl_to_buf(1, f4%n_gc, f4%bx(2), i0=f4%bx(1)-f4%n_gc)
    @:fyp_srl_to_buf(2, f4%bx(1), f4%n_gc)
    @:fyp_srl_to_buf(3, f4%bx(1), f4%n_gc, j0=f4%bx(2)-f4%n_gc)
#:elif NDIM == 3
    @:fyp_srl_to_buf(0, f4%n_gc, f4%bx(2), f4%bx(3))
    @:fyp_srl_to_buf(1, f4%n_gc, f4%bx(2), f4%bx(3), i0=f4%bx(1)-f4%n_gc)
    @:fyp_srl_to_buf(2, f4%bx(1), f4%n_gc, f4%bx(3))
    @:fyp_srl_to_buf(3, f4%bx(1), f4%n_gc, f4%bx(3), j0=f4%bx(2)-f4%n_gc)
    @:fyp_srl_to_buf(4, f4%bx(1), f4%bx(2), f4%n_gc)
    @:fyp_srl_to_buf(5, f4%bx(1), f4%bx(2), f4%n_gc, k0=f4%bx(3)-f4%n_gc)
#:endif

#:if NDIM == 2
    @:fyp_f2c_to_buf(0, f4%n_gc, f4%hbx(2))
    @:fyp_f2c_to_buf(1, f4%n_gc, f4%hbx(2), i0=f4%bx(1) - 2*f4%n_gc)
    @:fyp_f2c_to_buf(2, f4%hbx(1), f4%n_gc)
    @:fyp_f2c_to_buf(3, f4%hbx(1), f4%n_gc, j0=f4%bx(2) - 2*f4%n_gc)
#:elif NDIM == 3
    @:fyp_f2c_to_buf(0, f4%n_gc, f4%hbx(2), f4%hbx(3))
    @:fyp_f2c_to_buf(1, f4%n_gc, f4%hbx(2), f4%hbx(3), i0=f4%bx(1) - 2*f4%n_gc)
    @:fyp_f2c_to_buf(2, f4%hbx(1), f4%n_gc, f4%hbx(3))
    @:fyp_f2c_to_buf(3, f4%hbx(1), f4%n_gc, f4%hbx(3), j0=f4%bx(2) - 2*f4%n_gc)
    @:fyp_f2c_to_buf(4, f4%hbx(1), f4%hbx(2), f4%n_gc)
    @:fyp_f2c_to_buf(5, f4%hbx(1), f4%hbx(2), f4%n_gc, k0=f4%bx(3) - 2*f4%n_gc)
#:endif
    ${END_PARALLEL()}$

    f4%recv_offset(:) = f4%gc_recv_offset * n_vars
    f4%send_offset(:) = f4%gc_send_offset * n_vars

  end subroutine fill_ghostcell_buffers_round_one

  !> Fill buffers for the fine side of coarse-to-fine boundaries
  subroutine fill_ghostcell_buffers_round_two(f4, n_vars, i_vars, block_offset)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer, intent(in)          :: block_offset
    integer                      :: i, j, n, ivar, iv
    integer                      :: i_c, j_c
    integer                      :: iq, i_buf, i_buf0, offset(NDIM-1)
    integer                      :: half_n_gc
    logical                      :: odd_n_gc
    real(fp)                     :: fine(2**NDIM)
#:if NDIM == 3
    integer                      :: k, k_c
#:endif

    ! Update send/recv offsets
    f4%recv_offset(:) = f4%gc_recv_offset_c2f * n_vars
    f4%send_offset(:) = f4%gc_send_offset_c2f * n_vars

    ! If nothing to do, save time by not starting parallel region
    if (f4%gc_c2f_to_buf_iface(2*NDIM) == 1) return

    if (maxval(f4%gc_send_offset_c2f) * n_vars > size(f4%send_buffer)) then
       write(error_unit, *) "maxval(f4%gc_send_offset_c2f): ", &
            maxval(f4%gc_send_offset_c2f)
       write(error_unit, *) "n_vars: ", n_vars
       write(error_unit, *) "size(f4%send_buffer): ", size(f4%send_buffer)
       error stop "send buffer too small"
    end if

    half_n_gc = f4%n_gc/2 ! Round down
    odd_n_gc  = (iand(f4%n_gc, 1) == 1)

#:if NDIM == 2
#:def fyp_c2f_to_buf(face, ilim='f4%hbx(1)', jlim='f4%hbx(2)', &
    &ic0=0, jc0=0)
    ${LOOP('private(iq, offset, i_buf0)')}$
    do n = f4%gc_c2f_to_buf_iface(${face}$), f4%gc_c2f_to_buf_iface(${face}$+1)-1
       iq = f4%gc_c2f_to_buf(1, n) + 1 + block_offset ! coarse block
       offset(1) = f4%gc_c2f_to_buf(2, n)
       i_buf0 = f4%gc_c2f_to_buf(3, n) * n_vars

       ${LOOP('collapse(3) private(ivar, j_c, i_c, i_buf, fine)')}$
       do iv = 1, n_vars
          do j = 1, ${jlim}$
             do i = 1,  ${ilim}$
                ivar = i_vars(iv)
                i_c = ${ic0}$ + i
                j_c = ${jc0}$ + j
                call prolong(f4%gminmod_theta_prolong, &
                     f4%uu(i_c, j_c, ivar, iq), &
                     f4%uu(i_c-1, j_c, ivar, iq), &
                     f4%uu(i_c+1, j_c, ivar, iq), &
                     f4%uu(i_c, j_c-1, ivar, iq), &
                     f4%uu(i_c, j_c+1, ivar, iq), fine)

                i_buf = i_buf0 + 4 * ix_offset3(iv, j, i, ${jlim}$, ${ilim}$)
                f4%send_buffer(i_buf+1:i_buf+4) = fine
             end do
          end do
       end do

       i_buf0 = i_buf0 + 4 * n_vars * ${jlim}$ * ${ilim}$

       if (odd_n_gc) then
          ${LOOP('collapse(2) private(ivar, i_c, j_c, i_buf, fine)')}$
          do iv = 1, n_vars
#:if face == '0'
             do j = 1, ${jlim}$
                ! i = half_n_gc + 1
                ivar = i_vars(iv)
                i_c = ${ic0}$ + half_n_gc + 1
                j_c = ${jc0}$ + j

                call prolong(f4%gminmod_theta_prolong, &
                     f4%uu(i_c, j_c, ivar, iq), &
                     f4%uu(i_c-1, j_c, ivar, iq), &
                     f4%uu(i_c+1, j_c, ivar, iq), &
                     f4%uu(i_c, j_c-1, ivar, iq), &
                     f4%uu(i_c, j_c+1, ivar, iq), fine)

                i_buf = i_buf0 + 2 * ix_offset2(iv, j, ${jlim}$)
                f4%send_buffer(i_buf+1) = fine(1)
                f4%send_buffer(i_buf+2) = fine(3)
             end do
#:elif face == '1'
             do j = 1, ${jlim}$
                ! i = 0
                ivar = i_vars(iv)
                i_c = ${ic0}$ + 0
                j_c = ${jc0}$ + j

                call prolong(f4%gminmod_theta_prolong, &
                     f4%uu(i_c, j_c, ivar, iq), &
                     f4%uu(i_c-1, j_c, ivar, iq), &
                     f4%uu(i_c+1, j_c, ivar, iq), &
                     f4%uu(i_c, j_c-1, ivar, iq), &
                     f4%uu(i_c, j_c+1, ivar, iq), fine)

                i_buf = i_buf0 + 2 * ix_offset2(iv, j, ${jlim}$)
                f4%send_buffer(i_buf+1) = fine(2)
                f4%send_buffer(i_buf+2) = fine(4)
             end do
#:elif face == '2'
             do i = 1, ${ilim}$
                ! j = half_n_gc + 1
                ivar = i_vars(iv)
                i_c = ${ic0}$ + i
                j_c = ${jc0}$ + half_n_gc + 1

                call prolong(f4%gminmod_theta_prolong, &
                     f4%uu(i_c, j_c, ivar, iq), &
                     f4%uu(i_c-1, j_c, ivar, iq), &
                     f4%uu(i_c+1, j_c, ivar, iq), &
                     f4%uu(i_c, j_c-1, ivar, iq), &
                     f4%uu(i_c, j_c+1, ivar, iq), fine)

                i_buf = i_buf0 + 2 * ix_offset2(iv, i, ${ilim}$)
                f4%send_buffer(i_buf+1) = fine(1)
                f4%send_buffer(i_buf+2) = fine(2)
             end do
#:elif face == '3'
             do i = 1, ${ilim}$
                ! j = 0
                ivar = i_vars(iv)
                i_c = ${ic0}$ + i
                j_c = ${jc0}$ + 0

                call prolong(f4%gminmod_theta_prolong, &
                     f4%uu(i_c, j_c, ivar, iq), &
                     f4%uu(i_c-1, j_c, ivar, iq), &
                     f4%uu(i_c+1, j_c, ivar, iq), &
                     f4%uu(i_c, j_c-1, ivar, iq), &
                     f4%uu(i_c, j_c+1, ivar, iq), fine)

                i_buf = i_buf0 + 2 * ix_offset2(iv, i, ${ilim}$)
                f4%send_buffer(i_buf+1) = fine(3)
                f4%send_buffer(i_buf+2) = fine(4)
             end do
#:endif
          end do
       end if
    end do
#:enddef
#:elif NDIM == 3
#:def fyp_c2f_to_buf(face, ilim='f4%hbx(1)', jlim='f4%hbx(2)', &
    klim='f4%hbx(3)', ic0=0, jc0=0, kc0=0)
    ${LOOP('private(iq, offset, i_buf0)')}$
    do n = f4%gc_c2f_to_buf_iface(${face}$), f4%gc_c2f_to_buf_iface(${face}$+1)-1
       iq = f4%gc_c2f_to_buf(1, n) + 1 + block_offset ! coarse block
       offset(1:2) = f4%gc_c2f_to_buf(2:3, n)
       i_buf0 = f4%gc_c2f_to_buf(4, n) * n_vars

       ${LOOP('collapse(4) private(ivar, k_c, j_c, i_c, i_buf, fine)')}$
       do iv = 1, n_vars
          do k = 1, ${klim}$
             do j = 1, ${jlim}$
                do i = 1,  ${ilim}$
                   ivar = i_vars(iv)
                   i_c = ${ic0}$ + i
                   j_c = ${jc0}$ + j
                   k_c = ${kc0}$ + k

                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, k_c, ivar, iq), &
                        f4%uu(i_c-1, j_c, k_c, ivar, iq), &
                        f4%uu(i_c+1, j_c, k_c, ivar, iq), &
                        f4%uu(i_c, j_c-1, k_c, ivar, iq), &
                        f4%uu(i_c, j_c+1, k_c, ivar, iq), &
                        f4%uu(i_c, j_c, k_c-1, ivar, iq), &
                        f4%uu(i_c, j_c, k_c+1, ivar, iq), &
                        fine)

                   i_buf = i_buf0 + 8 * &
                        ix_offset4(iv, k, j, i, ${klim}$, ${jlim}$, ${ilim}$)
                   f4%send_buffer(i_buf+1:i_buf+8) = fine
                end do
             end do
          end do
       end do

       i_buf0 = i_buf0 + 8 * n_vars * ${klim}$ * ${jlim}$ * ${ilim}$

       if (odd_n_gc) then
          ${LOOP('collapse(3) private(ivar, i_c, j_c, k_c, i_buf, fine)')}$
          do iv = 1, n_vars
#:if face == '0'
             do k = 1, ${klim}$
                do j = 1, ${jlim}$
                   ! i = half_n_gc + 1
                   ivar = i_vars(iv)
                   i_c = ${ic0}$ + half_n_gc + 1
                   j_c = ${jc0}$ + j
                   k_c = ${kc0}$ + k

                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, k_c, ivar, iq), &
                        f4%uu(i_c-1, j_c, k_c, ivar, iq), &
                        f4%uu(i_c+1, j_c, k_c, ivar, iq), &
                        f4%uu(i_c, j_c-1, k_c, ivar, iq), &
                        f4%uu(i_c, j_c+1, k_c, ivar, iq), &
                        f4%uu(i_c, j_c, k_c-1, ivar, iq), &
                        f4%uu(i_c, j_c, k_c+1, ivar, iq), &
                        fine)
                   i_buf = i_buf0 + 4 * ix_offset3(iv, k, j, n_vars, ${klim}$)
                   f4%send_buffer(i_buf+1) = fine(1)
                   f4%send_buffer(i_buf+2) = fine(3)
                   f4%send_buffer(i_buf+3) = fine(5)
                   f4%send_buffer(i_buf+4) = fine(7)
                end do
             end do
#:elif face == '1'
             do k = 1, ${klim}$
                do j = 1, ${jlim}$
                   ! i = 0
                   ivar = i_vars(iv)
                   i_c = ${ic0}$ + 0
                   j_c = ${jc0}$ + j
                   k_c = ${kc0}$ + k

                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, k_c, ivar, iq), &
                        f4%uu(i_c-1, j_c, k_c, ivar, iq), &
                        f4%uu(i_c+1, j_c, k_c, ivar, iq), &
                        f4%uu(i_c, j_c-1, k_c, ivar, iq), &
                        f4%uu(i_c, j_c+1, k_c, ivar, iq), &
                        f4%uu(i_c, j_c, k_c-1, ivar, iq), &
                        f4%uu(i_c, j_c, k_c+1, ivar, iq), &
                        fine)

                   i_buf = i_buf0 + 4 * ix_offset3(iv, k, j, n_vars, ${klim}$)
                   f4%send_buffer(i_buf+1) = fine(2)
                   f4%send_buffer(i_buf+2) = fine(4)
                   f4%send_buffer(i_buf+3) = fine(6)
                   f4%send_buffer(i_buf+4) = fine(8)
                end do
             end do
#:elif face == '2'
             do k = 1, ${klim}$
                do i = 1, ${ilim}$
                   ! j = half_n_gc + 1
                   ivar = i_vars(iv)
                   i_c = ${ic0}$ + i
                   j_c = ${jc0}$ + half_n_gc + 1
                   k_c = ${kc0}$ + k

                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, k_c, ivar, iq), &
                        f4%uu(i_c-1, j_c, k_c, ivar, iq), &
                        f4%uu(i_c+1, j_c, k_c, ivar, iq), &
                        f4%uu(i_c, j_c-1, k_c, ivar, iq), &
                        f4%uu(i_c, j_c+1, k_c, ivar, iq), &
                        f4%uu(i_c, j_c, k_c-1, ivar, iq), &
                        f4%uu(i_c, j_c, k_c+1, ivar, iq), &
                        fine)

                   i_buf = i_buf0 + 4 * ix_offset3(iv, k, i, n_vars, ${klim}$)
                   f4%send_buffer(i_buf+1) = fine(1)
                   f4%send_buffer(i_buf+2) = fine(2)
                   f4%send_buffer(i_buf+3) = fine(5)
                   f4%send_buffer(i_buf+4) = fine(6)
                end do
             end do
#:elif face == '3'
             do k = 1, ${klim}$
                do i = 1, ${ilim}$
                   ! j = 0
                   ivar = i_vars(iv)
                   i_c = ${ic0}$ + i
                   j_c = ${jc0}$ + 0
                   k_c = ${kc0}$ + k

                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, k_c, ivar, iq), &
                        f4%uu(i_c-1, j_c, k_c, ivar, iq), &
                        f4%uu(i_c+1, j_c, k_c, ivar, iq), &
                        f4%uu(i_c, j_c-1, k_c, ivar, iq), &
                        f4%uu(i_c, j_c+1, k_c, ivar, iq), &
                        f4%uu(i_c, j_c, k_c-1, ivar, iq), &
                        f4%uu(i_c, j_c, k_c+1, ivar, iq), &
                        fine)

                   i_buf = i_buf0 + 4 * ix_offset3(iv, k, i, n_vars, ${klim}$)
                   f4%send_buffer(i_buf+1) = fine(3)
                   f4%send_buffer(i_buf+2) = fine(4)
                   f4%send_buffer(i_buf+3) = fine(7)
                   f4%send_buffer(i_buf+4) = fine(8)
                end do
             end do
#:elif face == '4'
             do j = 1, ${jlim}$
                do i = 1, ${ilim}$
                   ! k = half_n_gc + 1
                   ivar = i_vars(iv)
                   i_c = ${ic0}$ + i
                   j_c = ${jc0}$ + j
                   k_c = ${kc0}$ + half_n_gc + 1

                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, k_c, ivar, iq), &
                        f4%uu(i_c-1, j_c, k_c, ivar, iq), &
                        f4%uu(i_c+1, j_c, k_c, ivar, iq), &
                        f4%uu(i_c, j_c-1, k_c, ivar, iq), &
                        f4%uu(i_c, j_c+1, k_c, ivar, iq), &
                        f4%uu(i_c, j_c, k_c-1, ivar, iq), &
                        f4%uu(i_c, j_c, k_c+1, ivar, iq), &
                        fine)

                   i_buf = i_buf0 + 4 * ix_offset3(iv, j, i, n_vars, ${jlim}$)
                   f4%send_buffer(i_buf+1) = fine(1)
                   f4%send_buffer(i_buf+2) = fine(2)
                   f4%send_buffer(i_buf+3) = fine(3)
                   f4%send_buffer(i_buf+4) = fine(4)
                end do
             end do

#:elif face == '5'
             do j = 1, ${jlim}$
                do i = 1, ${ilim}$
                   ! k = 0
                   ivar = i_vars(iv)
                   i_c = ${ic0}$ + i
                   j_c = ${jc0}$ + j
                   k_c = ${kc0}$ + 0

                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, k_c, ivar, iq), &
                        f4%uu(i_c-1, j_c, k_c, ivar, iq), &
                        f4%uu(i_c+1, j_c, k_c, ivar, iq), &
                        f4%uu(i_c, j_c-1, k_c, ivar, iq), &
                        f4%uu(i_c, j_c+1, k_c, ivar, iq), &
                        f4%uu(i_c, j_c, k_c-1, ivar, iq), &
                        f4%uu(i_c, j_c, k_c+1, ivar, iq), &
                        fine)

                   i_buf = i_buf0 + 4 * ix_offset3(iv, j, i, n_vars, ${jlim}$)
                   f4%send_buffer(i_buf+1) = fine(5)
                   f4%send_buffer(i_buf+2) = fine(6)
                   f4%send_buffer(i_buf+3) = fine(7)
                   f4%send_buffer(i_buf+4) = fine(8)
                end do
             end do
#:endif
          end do
       end if
    end do
#:enddef
#:endif

    ${PARALLEL(COPYIN('i_vars'))}$ ${DEFAULT_PRESENT()}$

#:if NDIM == 2
    @:fyp_c2f_to_buf(0, ilim=half_n_gc, jc0=offset(1)*f4%hbx(2))
    @:fyp_c2f_to_buf(1, ilim=half_n_gc, ic0=f4%bx(1)-half_n_gc, &
         &jc0=offset(1)*f4%hbx(2))
    @:fyp_c2f_to_buf(2, jlim=half_n_gc, ic0=offset(1)*f4%hbx(1))
    @:fyp_c2f_to_buf(3, jlim=half_n_gc, ic0=offset(1)*f4%hbx(1), &
         &jc0=f4%bx(2)-half_n_gc)
#:elif NDIM == 3
    @:fyp_c2f_to_buf(0, ilim=half_n_gc, jc0=offset(1)*f4%hbx(2), &
         &kc0=offset(2)*f4%hbx(3))
    @:fyp_c2f_to_buf(1, ilim=half_n_gc, ic0=f4%bx(1)-half_n_gc, &
         &jc0=offset(1)*f4%hbx(2), kc0=offset(2)*f4%hbx(3))
    @:fyp_c2f_to_buf(2, jlim=half_n_gc, ic0=offset(1)*f4%hbx(1), &
         &kc0=offset(2)*f4%hbx(3))
    @:fyp_c2f_to_buf(3, jlim=half_n_gc, ic0=offset(1)*f4%hbx(1), &
         &jc0=f4%bx(2)-half_n_gc, kc0=offset(2)*f4%hbx(3))
    @:fyp_c2f_to_buf(4, klim=half_n_gc, ic0=offset(1)*f4%hbx(1), &
         &jc0=offset(2)*f4%hbx(2))
    @:fyp_c2f_to_buf(5, klim=half_n_gc, ic0=offset(1)*f4%hbx(1), &
         &jc0=offset(2)*f4%hbx(2), kc0=f4%bx(3)-half_n_gc)
#:endif

    ${END_PARALLEL()}$

  end subroutine fill_ghostcell_buffers_round_two

  !> Exchange the receive and send buffers according to the specified offsets
  !> per MPI rank
  subroutine f4_exchange_buffers(f4)
    type(foap4_t), intent(inout) :: f4
    type(MPI_Request)            :: send_req(f4%max_requests)
    type(MPI_Request)            :: recv_req(f4%max_requests)
    integer                      :: n_send, n_recv, ierr, rank
    integer(MPI_COUNT_KIND)      :: ilo, ihi
    integer, parameter           :: tag = 0

    n_send = 0
    n_recv = 0

    ! Use device pointers for device data
    do rank = 0, f4%mpisize - 1
       ilo = f4%send_offset(rank) + 1
       ihi = f4%send_offset(rank+1)

       if (ihi >= ilo) then
          call mpi_isend_wrapper(f4%send_buffer(ilo:ihi), ihi-ilo+1, &
               rank, tag, f4%mpicomm, send_req, n_send, ierr)
       end if

       ilo = f4%recv_offset(rank) + 1
       ihi = f4%recv_offset(rank+1)

       if (ihi >= ilo) then
          call mpi_irecv_wrapper(f4%recv_buffer(ilo:ihi), ihi-ilo+1, &
               rank, tag, f4%mpicomm, recv_req, n_recv, ierr)
       end if
    end do

    if (n_send > 0) then
       call MPI_Waitall(n_send, send_req(1:n_send), MPI_STATUSES_IGNORE, ierr)
    end if
    if (n_recv > 0) then
       call MPI_Waitall(n_recv, recv_req(1:n_recv), MPI_STATUSES_IGNORE, ierr)
    end if
  end subroutine f4_exchange_buffers

  !> Wrapper for MPI_Isend. This was done because of problems with the
  !> use_device clause combined with an index offset and NVHPC-v25. With this
  !> wrapper, there is no offset in "buf". The large-count interface
  !> (mpi_isend_c) is not widely supported, but can replace the use of
  !> multiple messages in the future.
  subroutine mpi_isend_wrapper(buf, count, dest, tag, comm, requests, nreqs, ierror)
    integer(MPI_COUNT_KIND), intent(in) :: count
    integer, intent(in)                 :: dest, tag
    real(fp), intent(in)                :: buf(count)
    type(mpi_comm), intent(in)          :: comm
    type(mpi_request), intent(inout)    :: requests(:)
    integer, intent(inout)              :: nreqs
    integer, optional, intent(out)      :: ierror

    integer(MPI_COUNT_KIND) :: offset, remaining, chunk_size
    integer :: chunk_count, chunk_tag, ierr, chunk_index

    offset = 0
    remaining = count
    chunk_index = 0

    do while (remaining > 0)
      chunk_size = min(remaining, f4_mpi_max_count)
      chunk_count = int(chunk_size)
      chunk_tag = tag + chunk_index
      chunk_index = chunk_index + 1
      nreqs = nreqs + 1

      if (nreqs > size(requests)) then
         write(error_unit, *) "Too many requests"
         write(error_unit, *) "nreqs: ", nreqs, " size(requests): ", size(requests)
         write(error_unit, *) "Increase f4%max_requests or f4_mpi_max_count"
         error stop "increase f4%max_requests"
      end if

      ${HOST_DATA_USE_DEVICE('buf')}$
      call mpi_isend(buf(offset + 1), chunk_count, f4_mpi_float, &
                     dest, chunk_tag, comm, requests(nreqs), ierr)
      ${END_HOST_DATA()}$

      if (ierr /= MPI_SUCCESS) then
        if (present(ierror)) ierror = ierr
        return
      end if

      offset = offset + chunk_size
      remaining = remaining - chunk_size
    end do

    if (present(ierror)) ierror = MPI_SUCCESS
  end subroutine mpi_isend_wrapper

  subroutine mpi_irecv_wrapper(buf, count, source, tag, comm, requests, nreqs, ierror)
    integer(MPI_COUNT_KIND), intent(in) :: count
    integer, intent(in)                 :: source, tag
    real(fp), intent(inout)             :: buf(count)
    type(mpi_comm), intent(in)          :: comm
    type(mpi_request), intent(inout)    :: requests(:)
    integer, intent(inout)              :: nreqs
    integer, optional, intent(out)      :: ierror

    integer(MPI_COUNT_KIND) :: offset, remaining, chunk_size
    integer :: chunk_count, chunk_tag, ierr, chunk_index

    offset = 0
    remaining = count
    chunk_index = 0

    do while (remaining > 0)
      chunk_size = min(remaining, f4_mpi_max_count)
      chunk_count = int(chunk_size)
      chunk_tag = tag + chunk_index
      chunk_index = chunk_index + 1
      nreqs = nreqs + 1

      if (nreqs > size(requests)) then
         write(error_unit, *) "Too many requests"
         write(error_unit, *) "nreqs: ", nreqs, " size(requests): ", size(requests)
         write(error_unit, *) "Increase f4%max_requests or f4_mpi_max_count"
         error stop "increase f4%max_requests"
      end if

      ${HOST_DATA_USE_DEVICE('buf')}$
      call mpi_irecv(buf(offset + 1), chunk_count, f4_mpi_float, &
                     source, chunk_tag, comm, requests(nreqs), ierr)
      ${END_HOST_DATA()}$

      if (ierr /= MPI_SUCCESS) then
        if (present(ierror)) ierror = ierr
        return
      end if

      offset = offset + chunk_size
      remaining = remaining - chunk_size
    end do

    if (present(ierror)) ierror = MPI_SUCCESS
  end subroutine mpi_irecv_wrapper

  !> After buffers have been communicated, handle all ghost cells for "round
  !> one", which excludes coarse-to-fine interpolation
  subroutine fill_ghostcells_round_one(f4, n_vars, i_vars, block_offset)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer, intent(in)          :: block_offset
    integer                      :: n, i, j, iq, jq, i_f, j_f
    integer                      :: i_buf, i_buf0, iv, ivar, i_bc_data
    integer                      :: offset(NDIM-1), bc_type, level
    real(fp)                     :: bc_value, dr(NDIM), slope
#:if NDIM == 3
    integer                      :: k, k_f
#:endif

#:def fyp_srl_local(face, ilim, jlim, klim=None, i0=0, j0=0, k0=0, i1=0, j1=0, k1=0, &
    &i2=0, j2=0, k2=0)
    ${LOOP('collapse(NDIM+2) private(iq, jq, ivar)')}$
    do n = f4%gc_srl_local_iface(${face}$), f4%gc_srl_local_iface(${face}$+1)-1
#:if NDIM == 2
       do iv = 1, n_vars
          do j = 1, ${jlim}$
             do i = 1, ${ilim}$
                iq = f4%gc_srl_local(1, n) + 1 + block_offset
                jq = f4%gc_srl_local(2, n) + 1 + block_offset
                ivar = i_vars(iv)
                f4%uu(${i0}$+i, ${j0}$+j, ivar, iq) = &
                     f4%uu(i, j, ivar, jq)
                f4%uu(${i1}$+i, ${j1}$+j, ivar, jq) = &
                     f4%uu(${i2}$+i, ${j2}$+j, ivar, iq)
             end do
          end do
       end do
#:elif NDIM == 3
       do iv = 1, n_vars
          do k = 1, ${klim}$
             do j = 1, ${jlim}$
                do i = 1, ${ilim}$
                   iq = f4%gc_srl_local(1, n) + 1 + block_offset
                   jq = f4%gc_srl_local(2, n) + 1 + block_offset
                   ivar = i_vars(iv)
                   f4%uu(${i0}$+i, ${j0}$+j, ${k0}$+k, ivar, iq) = &
                        f4%uu(i, j, k, ivar, jq)
                   f4%uu(${i1}$+i, ${j1}$+j, ${k1}$+k, ivar, jq) = &
                        f4%uu(${i2}$+i, ${j2}$+j, ${k2}$+k, ivar, iq)
                end do
             end do
          end do
       end do
#:endif
    end do
#:enddef

    ! Fill physical boundaries

#:def fyp_phys(face, ilim, jlim, klim=None)
    ${LOOP('collapse(NDIM+2) private(ivar, bc_type, bc_value, slope, iq, level, dr, i_bc_data)')}$
    do n = f4%gc_phys_iface(${face}$), f4%gc_phys_iface(${face}$+1)-1
#:if NDIM == 2
       do iv = 1, n_vars
          do j = 1, ${jlim}$
             do i = 1, ${ilim}$
                ivar = i_vars(iv)
                iq    = f4%gc_phys(n) + 1 + block_offset
                level = f4%block_level(iq - block_offset)
                dr    = real(f4%dr_level(:, level), fp)
                ! Do not use block_offset when accessing b.c. data
                i_bc_data = f4%bc_data_ix(${face}$, iq - block_offset)

                if (i_bc_data > 0) then
                   ! Use array value
#:if face in ['0', '1']
                   bc_value = f4%bc_data(j, i_vars(iv), i_bc_data)
                   bc_type = f4%bc_data_type(j, i_vars(iv), i_bc_data)
#:else
                   bc_value = f4%bc_data(i, i_vars(iv), i_bc_data)
                   bc_type = f4%bc_data_type(i, i_vars(iv), i_bc_data)
#:endif
                else
                   ! Use scalar values
                   bc_value = f4%bc_simple(ivar, ${face}$)
                   bc_type = f4%bc_simple_type(ivar, ${face}$)
                end if

                select case (bc_type)
                case (f4_bc_dirichlet)
#:if face == '0'
                   f4%uu(1-i, j, ivar, iq) = &
                        2 * bc_value - f4%uu(i, j, ivar, iq)
#:elif face == '1'
                   f4%uu(f4%bx(1)+i, j, ivar, iq) = &
                        2 * bc_value - f4%uu(f4%bx(1)+1-i, j, ivar, iq)
#:elif face == '2'
                   f4%uu(i, 1-j, ivar, iq) = &
                        2 * bc_value - f4%uu(i, j, ivar, iq)
#:elif face == '3'
                   f4%uu(i, f4%bx(2)+j, ivar, iq) = &
                        2 * bc_value - f4%uu(i, f4%bx(2)+1-j, ivar, iq)
#:endif
                case (f4_bc_neumann)
#:if face == '0'
                   f4%uu(1-i, j, ivar, iq) = f4%uu(i, j, ivar, iq) - &
                        (2*i-1) * dr(1) * bc_value
#:elif face == '1'
                   f4%uu(f4%bx(1)+i, j, ivar, iq) = f4%uu(f4%bx(1)+1-i, j, ivar, iq) + &
                        (2*i-1) * dr(1) * bc_value
#:elif face == '2'
                   f4%uu(i, 1-j, ivar, iq) = f4%uu(i, j, ivar, iq) - &
                        (2*j-1) * dr(2) * bc_value
#:elif face == '3'
                   f4%uu(i, f4%bx(2)+j, ivar, iq) = f4%uu(i, f4%bx(2)+1-j, ivar, iq) + &
                        (2*j-1) * dr(2) * bc_value
#:endif

                case (f4_bc_linear_extrap)
#:if face == '0'
                   slope = (f4%uu(2, j, ivar, iq) - f4%uu(1, j, ivar, iq))
                   f4%uu(1-i, j, ivar, iq) = f4%uu(1, j, ivar, iq) - i * slope
#:elif face == '1'
                   slope = (f4%uu(f4%bx(1), j, ivar, iq) - f4%uu(f4%bx(1)-1, j, ivar, iq))
                   f4%uu(f4%bx(1)+i, j, ivar, iq) = f4%uu(f4%bx(1), j, ivar, iq) + i * slope

#:elif face == '2'
                   slope = (f4%uu(i, 2, ivar, iq) - f4%uu(i, 1, ivar, iq))
                   f4%uu(i, 1-j, ivar, iq) = f4%uu(i, 1, ivar, iq) - j * slope

#:elif face == '3'
                   slope = (f4%uu(i, f4%bx(2), ivar, iq) - f4%uu(i, f4%bx(2)-1, ivar, iq))
                   f4%uu(i, f4%bx(2)+j, ivar, iq) = f4%uu(i, f4%bx(2), ivar, iq) + j * slope
#:endif
                case (f4_bc_fixed_value)
#:if face == '0'
                   f4%uu(1-i, j, ivar, iq) = bc_value
#:elif face == '1'
                   f4%uu(f4%bx(1)+i, j, ivar, iq) = bc_value
#:elif face == '2'
                   f4%uu(i, 1-j, ivar, iq) = bc_value
#:elif face == '3'
                   f4%uu(i, f4%bx(2)+j, ivar, iq) = bc_value
#:endif
                end select
             end do
          end do
       end do
#:elif NDIM == 3
       do iv = 1, n_vars
          do k = 1, ${klim}$
             do j = 1, ${jlim}$
                do i = 1, ${ilim}$
                   ivar = i_vars(iv)
                   iq    = f4%gc_phys(n) + 1 + block_offset
                   level = f4%block_level(iq - block_offset)
                   dr    = real(f4%dr_level(:, level), fp)
                   ! Do not use block_offset when accessing b.c. data
                   i_bc_data = f4%bc_data_ix(${face}$, iq - block_offset)

                   if (i_bc_data > 0) then
                   ! Use array value
#:if face in ['0', '1']
                   bc_value = f4%bc_data(j, k, i_vars(iv), i_bc_data)
                   bc_type = f4%bc_data_type(j, k, i_vars(iv), i_bc_data)
#:elif face in ['2', '3']
                   bc_value = f4%bc_data(i, k, i_vars(iv), i_bc_data)
                   bc_type = f4%bc_data_type(i, k, i_vars(iv), i_bc_data)
#:else
                   bc_value = f4%bc_data(i, j, i_vars(iv), i_bc_data)
                   bc_type = f4%bc_data_type(i, j, i_vars(iv), i_bc_data)
#:endif
                else
                   ! Use stored scalar
                   bc_value = f4%bc_simple(ivar, ${face}$)
                   bc_type = f4%bc_simple_type(ivar, ${face}$)
                end if

                   select case (bc_type)
                   case (f4_bc_dirichlet)
#:if face == '0'
                      f4%uu(1-i, j, k, ivar, iq) = &
                           2 * bc_value - f4%uu(i, j, k, ivar, iq)
#:elif face == '1'
                      f4%uu(f4%bx(1)+i, j, k, ivar, iq) = &
                           2 * bc_value - f4%uu(f4%bx(1)+1-i, j, k, ivar, iq)
#:elif face == '2'
                      f4%uu(i, 1-j, k, ivar, iq) = &
                           2 * bc_value - f4%uu(i, j, k, ivar, iq)
#:elif face == '3'
                      f4%uu(i, f4%bx(2)+j, k, ivar, iq) = &
                           2 * bc_value - f4%uu(i, f4%bx(2)+1-j, k, ivar, iq)
#:elif face == '4'
                      f4%uu(i, j, 1-k, ivar, iq) = &
                           2 * bc_value - f4%uu(i, j, k, ivar, iq)
#:elif face == '5'
                      f4%uu(i, j, f4%bx(3)+k, ivar, iq) = &
                           2 * bc_value - f4%uu(i, j, f4%bx(3)+1-k, ivar, iq)
#:endif
                   case (f4_bc_neumann)
#:if face == '0'
                      f4%uu(1-i, j, k, ivar, iq) = f4%uu(i, j, k, ivar, iq) - &
                           (2*i-1) * dr(1) * bc_value
#:elif face == '1'
                      f4%uu(f4%bx(1)+i, j, k, ivar, iq) = f4%uu(f4%bx(1)+1-i, j, k, ivar, iq) + &
                           (2*i-1) * dr(1) * bc_value
#:elif face == '2'
                      f4%uu(i, 1-j, k, ivar, iq) = f4%uu(i, j, k, ivar, iq) - &
                           (2*j-1) * dr(2) * bc_value
#:elif face == '3'
                      f4%uu(i, f4%bx(2)+j, k, ivar, iq) = f4%uu(i, f4%bx(2)+1-j, k, ivar, iq) + &
                           (2*j-1) * dr(2) * bc_value
#:elif face == '4'
                      f4%uu(i, j, 1-k, ivar, iq) = f4%uu(i, j, k, ivar, iq) - &
                           (2*k-1) * dr(3) * bc_value
#:elif face == '5'
                      f4%uu(i, j, f4%bx(3)+k, ivar, iq) = f4%uu(i, j, f4%bx(3)+1-k, ivar, iq) + &
                           (2*k-1) * dr(3) * bc_value
#:endif

                   case (f4_bc_linear_extrap)
#:if face == '0'
                      slope = (f4%uu(2, j, k, ivar, iq) - f4%uu(1, j, k, ivar, iq))
                      f4%uu(1-i, j, k, ivar, iq) = f4%uu(1, j, k, ivar, iq) - i * slope
#:elif face == '1'
                      slope = (f4%uu(f4%bx(1), j, k, ivar, iq) - f4%uu(f4%bx(1)-1, j, k, ivar, iq))
                      f4%uu(f4%bx(1)+i, j, k, ivar, iq) = f4%uu(f4%bx(1), j, k, ivar, iq) + i * slope

#:elif face == '2'
                      slope = (f4%uu(i, 2, k, ivar, iq) - f4%uu(i, 1, k, ivar, iq))
                      f4%uu(i, 1-j, k, ivar, iq) = f4%uu(i, 1, k, ivar, iq) - j * slope

#:elif face == '3'
                      slope = (f4%uu(i, f4%bx(2), k, ivar, iq) - f4%uu(i, f4%bx(2)-1, k, ivar, iq))
                      f4%uu(i, f4%bx(2)+j, k, ivar, iq) = f4%uu(i, f4%bx(2), k, ivar, iq) + j * slope
#:elif face == '4'
                      slope = (f4%uu(i, j, 2, ivar, iq) - f4%uu(i, j, 1, ivar, iq))
                      f4%uu(i, j, 1-k, ivar, iq) = f4%uu(i, j, 1, ivar, iq) - k * slope
#:elif face == '5'
                      slope = (f4%uu(i, j, f4%bx(3), ivar, iq) - f4%uu(i, j, f4%bx(3)-1, ivar, iq))
                      f4%uu(i, j, f4%bx(3)+k, ivar, iq) = f4%uu(i, j, f4%bx(3), ivar, iq) + k * slope
#:endif
                   case (f4_bc_fixed_value)
#:if face == '0'
                      f4%uu(1-i, j, k, ivar, iq) = bc_value
#:elif face == '1'
                      f4%uu(f4%bx(1)+i, j, k, ivar, iq) = bc_value
#:elif face == '2'
                      f4%uu(i, 1-j, k, ivar, iq) = bc_value
#:elif face == '3'
                      f4%uu(i, f4%bx(2)+j, k, ivar, iq) = bc_value
#:elif face == '4'
                      f4%uu(i, j, 1-k, ivar, iq) = bc_value
#:elif face == '5'
                      f4%uu(i, j, f4%bx(3)+k, ivar, iq) = bc_value
#:endif
                   end select
                end do
             end do
          end do
       end do
#:endif
    end do
#:enddef

#:def fyp_f2c_local(face, ilim, jlim, klim=None, i0=0, j0=0, k0=0, if0=0, jf0=0, kf0=0)
#:if NDIM == 2
    ${LOOP('collapse(NDIM+2) private(i_f, j_f, ivar, iq, jq, offset)')}$
    do n = f4%gc_f2c_local_iface(${face}$), f4%gc_f2c_local_iface(${face}$+1)-1
       do iv = 1, n_vars
          do j = 1, ${jlim}$
             do i = 1, ${ilim}$
                ivar = i_vars(iv)
                iq     = f4%gc_f2c_local(1, n) + 1 + block_offset ! Fine block
                jq     = f4%gc_f2c_local(2, n) + 1 + block_offset ! coarse block
                offset(1) = f4%gc_f2c_local(3, n)     ! offset

                j_f = ${jf0}$ + 2 * j - 1
                i_f = ${if0}$ + 2 * i - 1
                f4%uu(${i0}$+i, ${j0}$+j, ivar, jq) = 0.25_fp * ( &
                     f4%uu(i_f,   j_f,   ivar, iq) + &
                     f4%uu(i_f+1, j_f,   ivar, iq) + &
                     f4%uu(i_f,   j_f+1, ivar, iq) + &
                     f4%uu(i_f+1, j_f+1, ivar, iq) )
             end do
          end do
       end do
    end do
#:elif NDIM == 3
    ${LOOP('collapse(NDIM+2) private(i_f, j_f, k_f, ivar, iq, jq, offset)')}$
    do n = f4%gc_f2c_local_iface(${face}$), f4%gc_f2c_local_iface(${face}$+1)-1
       do iv = 1, n_vars
          do k = 1, ${klim}$
             do j = 1, ${jlim}$
                do i = 1, ${ilim}$
                   ivar = i_vars(iv)
                   iq     = f4%gc_f2c_local(1, n) + 1 + block_offset ! Fine block
                   jq     = f4%gc_f2c_local(2, n) + 1 + block_offset ! coarse block
                   offset(1:2) = f4%gc_f2c_local(3:4, n)     ! offset

                   k_f = ${kf0}$ + 2 * k - 1
                   j_f = ${jf0}$ + 2 * j - 1
                   i_f = ${if0}$ + 2 * i - 1
                   f4%uu(${i0}$+i, ${j0}$+j, ${k0}$+k, ivar, jq) = 0.125_fp * ( &
                        f4%uu(i_f,   j_f,   k_f,   ivar, iq) + &
                        f4%uu(i_f+1, j_f,   k_f,   ivar, iq) + &
                        f4%uu(i_f,   j_f+1, k_f,   ivar, iq) + &
                        f4%uu(i_f+1, j_f+1, k_f,   ivar, iq) + &
                        f4%uu(i_f,   j_f,   k_f+1, ivar, iq) + &
                        f4%uu(i_f+1, j_f,   k_f+1, ivar, iq) + &
                        f4%uu(i_f,   j_f+1, k_f+1, ivar, iq) + &
                        f4%uu(i_f+1, j_f+1, k_f+1, ivar, iq) )
                end do
             end do
          end do
       end do
    end do
#:endif
#:enddef

#:def fyp_srl_from_buf(face, ilim, jlim, klim=None, i0=0, j0=0, k0=0)
    ${LOOP('collapse(NDIM+2) private(ivar, i_buf, iq, i_buf0)')}$
    do n = f4%gc_srl_from_buf_iface(${face}$), f4%gc_srl_from_buf_iface(${face}$+1)-1
#:if NDIM == 2
       do iv = 1, n_vars
          do j = 1, ${jlim}$
             do i = 1, ${ilim}$
                ivar = i_vars(iv)
                iq = f4%gc_srl_from_buf(1, n) + 1 + block_offset
                i_buf0 = f4%gc_srl_from_buf(2, n) * n_vars
                i_buf = i_buf0 + ix_offset3(iv, j, i, ${jlim}$, ${ilim}$)
                f4%uu(${i0}$+i, ${j0}$+j, ivar, iq) = f4%recv_buffer(i_buf+1)
             end do
          end do
       end do
#:elif NDIM == 3
       do iv = 1, n_vars
          do k = 1, ${klim}$
             do j = 1, ${jlim}$
                do i = 1, ${ilim}$
                   ivar = i_vars(iv)
                   iq = f4%gc_srl_from_buf(1, n) + 1 + block_offset
                   i_buf0 = f4%gc_srl_from_buf(2, n) * n_vars
                   i_buf = i_buf0 + ix_offset4(iv, k, j, i, &
                        ${klim}$, ${jlim}$, ${ilim}$)
                   f4%uu(${i0}$+i, ${j0}$+j, ${k0}$+k, ivar, iq) = f4%recv_buffer(i_buf+1)
                end do
             end do
          end do
       end do
#:endif
    end do
#:enddef

#:def fyp_c2f_from_buf(face, ilim, jlim, klim=None, i0=0, j0=0, k0=0)
    ${LOOP('collapse(NDIM+2) private(ivar, i_buf, iq, offset, i_buf0)')}$
    do n = f4%gc_c2f_from_buf_iface(${face}$), f4%gc_c2f_from_buf_iface(${face}$+1)-1
#:if NDIM == 2
       do iv = 1, n_vars
          do j = 1, ${jlim}$
             do i = 1, ${ilim}$
                ivar = i_vars(iv)
                iq     = f4%gc_c2f_from_buf(1, n) + 1 + block_offset ! Coarse block
                offset(1) = f4%gc_c2f_from_buf(2, n)     ! Offset
                i_buf0  = f4%gc_c2f_from_buf(3, n) * n_vars
                i_buf = i_buf0 + ix_offset3(iv, j, i, ${jlim}$, ${ilim}$) + 1
                f4%uu(${i0}$+i, ${j0}$+j, ivar, iq) = &
                     f4%recv_buffer(i_buf)
             end do
          end do
       end do
#:elif NDIM == 3
       do iv = 1, n_vars
          do k = 1, ${klim}$
             do j = 1, ${jlim}$
                do i = 1, ${ilim}$
                   ivar = i_vars(iv)
                   iq     = f4%gc_c2f_from_buf(1, n) + 1 + block_offset ! Coarse block
                   offset(1:2) = f4%gc_c2f_from_buf(2:3, n)     ! Offset
                   i_buf0  = f4%gc_c2f_from_buf(4, n) * n_vars
                   i_buf = i_buf0 + ix_offset4(iv, k, j, i, &
                        ${klim}$, ${jlim}$, ${ilim}$) + 1
                   f4%uu(${i0}$+i, ${j0}$+j, ${k0}$+k, ivar, iq) = &
                        f4%recv_buffer(i_buf)
                end do
             end do
          end do
       end do
#:endif
    end do
#:enddef

    ${PARALLEL(COPYIN('i_vars'))}$ ${DEFAULT_PRESENT()}$

    ! Fill local boundaries at the same refinement level
#:if NDIM == 2
    @:fyp_srl_local(1, f4%n_gc,  f4%bx(2), i0=f4%bx(1), i1=-f4%n_gc, i2=f4%bx(1)-f4%n_gc)
    @:fyp_srl_local(3, f4%bx(1), f4%n_gc,  j0=f4%bx(2), j1=-f4%n_gc, j2=f4%bx(2)-f4%n_gc)
#:elif NDIM == 3
    @:fyp_srl_local(1, f4%n_gc,  f4%bx(2), f4%bx(3), i0=f4%bx(1), i1=-f4%n_gc, i2=f4%bx(1)-f4%n_gc)
    @:fyp_srl_local(3, f4%bx(1), f4%n_gc,  f4%bx(3), j0=f4%bx(2), j1=-f4%n_gc, j2=f4%bx(2)-f4%n_gc)
    @:fyp_srl_local(5, f4%bx(1), f4%bx(2), f4%n_gc,  k0=f4%bx(3), k1=-f4%n_gc, k2=f4%bx(3)-f4%n_gc)
#:endif

    ! Fill ghost cells at physical boundaries
#:if NDIM == 2
    @:fyp_phys(0, f4%n_gc, f4%bx(2))
    @:fyp_phys(1, f4%n_gc, f4%bx(2))
    @:fyp_phys(2, f4%bx(1), f4%n_gc)
    @:fyp_phys(3, f4%bx(1), f4%n_gc)
#:elif NDIM == 3
    @:fyp_phys(0, f4%n_gc, f4%bx(2), f4%bx(3))
    @:fyp_phys(1, f4%n_gc, f4%bx(2), f4%bx(3))
    @:fyp_phys(2, f4%bx(1), f4%n_gc, f4%bx(3))
    @:fyp_phys(3, f4%bx(1), f4%n_gc, f4%bx(3))
    @:fyp_phys(4, f4%bx(1), f4%bx(2), f4%n_gc)
    @:fyp_phys(5, f4%bx(1), f4%bx(2), f4%n_gc)
#:endif

    ! Fill coarse side of local fine-to-coarse refinement boundaries
#:if NDIM == 2
    @:fyp_f2c_local(0, f4%n_gc, f4%hbx(2), i0=f4%bx(1), j0=offset(1)*f4%hbx(2))
    @:fyp_f2c_local(1, f4%n_gc, f4%hbx(2), i0=-f4%n_gc, j0=offset(1)*f4%hbx(2), if0=f4%bx(1) - 2*f4%n_gc)
    @:fyp_f2c_local(2, f4%hbx(1), f4%n_gc, i0=offset(1)*f4%hbx(1), j0=f4%bx(2))
    @:fyp_f2c_local(3, f4%hbx(1), f4%n_gc, i0=offset(1)*f4%hbx(1), j0=-f4%n_gc, jf0=f4%bx(2) -2*f4%n_gc)
#:elif NDIM == 3
    @:fyp_f2c_local(0, f4%n_gc, f4%hbx(2), f4%hbx(3), i0=f4%bx(1), &
         &j0=offset(1)*f4%hbx(2), k0=offset(2)*f4%hbx(3))
    @:fyp_f2c_local(1, f4%n_gc, f4%hbx(2), f4%hbx(3), i0=-f4%n_gc, &
         &j0=offset(1)*f4%hbx(2), k0=offset(2)*f4%hbx(3), if0=f4%bx(1) - 2*f4%n_gc)
    @:fyp_f2c_local(2, f4%hbx(1), f4%n_gc, f4%hbx(3), &
         &i0=offset(1)*f4%hbx(1), j0=f4%bx(2), k0=offset(2)*f4%hbx(3))
    @:fyp_f2c_local(3, f4%hbx(1), f4%n_gc, f4%hbx(3), &
         i0=offset(1)*f4%hbx(1), j0=-f4%n_gc, k0=offset(2)*f4%hbx(3), jf0=f4%bx(2) -2*f4%n_gc)
    @:fyp_f2c_local(4, f4%hbx(1), f4%hbx(2), f4%n_gc, &
         &i0=offset(1)*f4%hbx(1), j0=offset(2)*f4%hbx(2), k0=f4%bx(3))
    @:fyp_f2c_local(5, f4%hbx(1), f4%hbx(2), f4%n_gc, &
         i0=offset(1)*f4%hbx(1), j0=offset(2)*f4%hbx(2), k0=-f4%n_gc, kf0=f4%bx(3) -2*f4%n_gc)
#:endif

    ! Fill ghost cells at the same refinement level from the buffer
#:if NDIM == 2
    @:fyp_srl_from_buf(0, f4%n_gc, f4%bx(2), i0=-f4%n_gc)
    @:fyp_srl_from_buf(1, f4%n_gc, f4%bx(2), i0=f4%bx(1))
    @:fyp_srl_from_buf(2, f4%bx(1), f4%n_gc, j0=-f4%n_gc)
    @:fyp_srl_from_buf(3, f4%bx(1), f4%n_gc, j0=f4%bx(2))
#:elif NDIM == 3
    @:fyp_srl_from_buf(0, f4%n_gc, f4%bx(2), f4%bx(3), i0=-f4%n_gc)
    @:fyp_srl_from_buf(1, f4%n_gc, f4%bx(2), f4%bx(3), i0=f4%bx(1))
    @:fyp_srl_from_buf(2, f4%bx(1), f4%n_gc, f4%bx(3), j0=-f4%n_gc)
    @:fyp_srl_from_buf(3, f4%bx(1), f4%n_gc, f4%bx(3), j0=f4%bx(2))
    @:fyp_srl_from_buf(4, f4%bx(1), f4%bx(2), f4%n_gc, k0=-f4%n_gc)
    @:fyp_srl_from_buf(5, f4%bx(1), f4%bx(2), f4%n_gc, k0=f4%bx(3))
#:endif

    ! Update coarse side from buffers at coarse-to-fine buffer
#:if NDIM == 2
    @:fyp_c2f_from_buf(0, f4%n_gc, f4%hbx(2), i0=-f4%n_gc, j0=offset(1)*f4%hbx(2))
    @:fyp_c2f_from_buf(1, f4%n_gc, f4%hbx(2), i0=f4%bx(1), j0=offset(1)*f4%hbx(2))
    @:fyp_c2f_from_buf(2, f4%hbx(1), f4%n_gc, i0=offset(1)*f4%hbx(1), j0=-f4%n_gc)
    @:fyp_c2f_from_buf(3, f4%hbx(1), f4%n_gc, i0=offset(1)*f4%hbx(1), j0=f4%bx(2))
#:elif NDIM == 3
    @:fyp_c2f_from_buf(0, f4%n_gc, f4%hbx(2), f4%hbx(3), &
         &i0=-f4%n_gc, j0=offset(1)*f4%hbx(2), k0=offset(2)*f4%hbx(3))
    @:fyp_c2f_from_buf(1, f4%n_gc, f4%hbx(2), f4%hbx(3), &
         &i0=f4%bx(1), j0=offset(1)*f4%hbx(2), k0=offset(2)*f4%hbx(3))
    @:fyp_c2f_from_buf(2, f4%hbx(1), f4%n_gc, f4%hbx(3), &
         &i0=offset(1)*f4%hbx(1), j0=-f4%n_gc, k0=offset(2)*f4%hbx(3))
    @:fyp_c2f_from_buf(3, f4%hbx(1), f4%n_gc, f4%hbx(3), &
         &i0=offset(1)*f4%hbx(1), j0=f4%bx(2), k0=offset(2)*f4%hbx(3))
    @:fyp_c2f_from_buf(4, f4%hbx(1), f4%hbx(2), f4%n_gc, &
         &i0=offset(1)*f4%hbx(1), j0=offset(2)*f4%hbx(2), k0=-f4%n_gc)
    @:fyp_c2f_from_buf(5, f4%hbx(1), f4%hbx(2), f4%n_gc, &
         &i0=offset(1)*f4%hbx(1), j0=offset(2)*f4%hbx(2), k0=f4%bx(3))
#:endif

    ${END_PARALLEL()}$

  end subroutine fill_ghostcells_round_one

  !> Handle coarse-to-fine ghost cells on the fine side
  subroutine fill_ghostcells_round_two(f4, n_vars, i_vars, block_offset)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer, intent(in)          :: block_offset
#:if NDIM == 3
    integer                      :: k, k_c, k_f
#:endif
    integer                      :: n, i, j, iq, jq, i_c, j_c, i_f, j_f
    integer                      :: i_buf, i_buf0, iv, ivar
    integer                      :: half_n_gc, offset(NDIM-1)
    logical                      :: odd_n_gc
    real(fp)                     :: fine(2**NDIM)

    ! If nothing to do, save time by not starting parallel region
    if (f4%gc_f2c_local_iface(2*NDIM) == 1 .and. &
         f4%gc_f2c_from_buf_iface(2*NDIM) == 1) return

    half_n_gc = f4%n_gc/2 ! Round down
    odd_n_gc  = (iand(f4%n_gc, 1) == 1)

#:if NDIM == 2
#:def fyp_f2c_local_fine(face, ilim='f4%hbx(1)', jlim='f4%hbx(2)', &
    &ic0=0, jc0=0, if0=0, jf0=0)
    ${LOOP('private(iq, jq, offset)')}$
    do n = f4%gc_f2c_local_iface(${face}$), f4%gc_f2c_local_iface(${face}$+1)-1
       iq     = f4%gc_f2c_local(1, n) + 1 + block_offset ! Fine block
       jq     = f4%gc_f2c_local(2, n) + 1 + block_offset ! coarse block
       offset(1) = f4%gc_f2c_local(3, n)     ! Offset

       ${LOOP('collapse(3) private(ivar, j_f, j_c, i_f, i_c, fine)')}$
       do iv = 1, n_vars
          do j = 1, ${jlim}$
             do i = 1,  ${ilim}$
                ivar = i_vars(iv)
                i_f = ${if0}$ + 2 * i - 1
                j_f = ${jf0}$ + 2 * j - 1
                i_c = ${ic0}$ + i
                j_c = ${jc0}$ + j
                call prolong(f4%gminmod_theta_prolong, &
                     f4%uu(i_c, j_c, ivar, jq), &
                     f4%uu(i_c-1, j_c, ivar, jq), &
                     f4%uu(i_c+1, j_c, ivar, jq), &
                     f4%uu(i_c, j_c-1, ivar, jq), &
                     f4%uu(i_c, j_c+1, ivar, jq), fine)
                f4%uu(i_f  , j_f  , ivar, iq) = fine(1)
                f4%uu(i_f+1, j_f  , ivar, iq) = fine(2)
                f4%uu(i_f  , j_f+1, ivar, iq) = fine(3)
                f4%uu(i_f+1, j_f+1, ivar, iq) = fine(4)
             end do
          end do
       end do

       if (odd_n_gc) then
          ${LOOP('collapse(2) private(ivar, j_f, j_c, i_f, i_c, fine)')}$
          do iv = 1, n_vars
#:if face == '0'
             do j = 1, ${jlim}$
                ! i = 0
                ivar = i_vars(iv)
                i_f = ${if0}$ + 2 * 0 - 1
                j_f = ${jf0}$ + 2 * j - 1
                i_c = ${ic0}$ + 0
                j_c = ${jc0}$ + j

                call prolong(f4%gminmod_theta_prolong, &
                     f4%uu(i_c, j_c, ivar, jq), &
                     f4%uu(i_c-1, j_c, ivar, jq), &
                     f4%uu(i_c+1, j_c, ivar, jq), &
                     f4%uu(i_c, j_c-1, ivar, jq), &
                     f4%uu(i_c, j_c+1, ivar, jq), fine)

                f4%uu(i_f+1, j_f  , ivar, iq) = fine(2)
                f4%uu(i_f+1, j_f+1, ivar, iq) = fine(4)
             end do
#:elif face == '1'
             do j = 1, ${jlim}$
                ! i = half_n_gc + 1
                ivar = i_vars(iv)
                i_f = ${if0}$ + 2 * (half_n_gc + 1) - 1
                j_f = ${jf0}$ + 2 * j - 1
                i_c = ${ic0}$ + (half_n_gc + 1)
                j_c = ${jc0}$ + j

                call prolong(f4%gminmod_theta_prolong, &
                     f4%uu(i_c, j_c, ivar, jq), &
                     f4%uu(i_c-1, j_c, ivar, jq), &
                     f4%uu(i_c+1, j_c, ivar, jq), &
                     f4%uu(i_c, j_c-1, ivar, jq), &
                     f4%uu(i_c, j_c+1, ivar, jq), fine)

                f4%uu(i_f  , j_f  , ivar, iq) = fine(1)
                f4%uu(i_f  , j_f+1, ivar, iq) = fine(3)
             end do
#:elif face == '2'
             do i = 1, ${ilim}$
                ! j = 0
                ivar = i_vars(iv)
                i_f = ${if0}$ + 2 * i - 1
                j_f = ${jf0}$ + 2 * 0 - 1
                i_c = ${ic0}$ + i
                j_c = ${jc0}$ + 0

                call prolong(f4%gminmod_theta_prolong, &
                     f4%uu(i_c, j_c, ivar, jq), &
                     f4%uu(i_c-1, j_c, ivar, jq), &
                     f4%uu(i_c+1, j_c, ivar, jq), &
                     f4%uu(i_c, j_c-1, ivar, jq), &
                     f4%uu(i_c, j_c+1, ivar, jq), fine)

                f4%uu(i_f  , j_f+1, ivar, iq) = fine(3)
                f4%uu(i_f+1, j_f+1, ivar, iq) = fine(4)
             end do
#:elif face == '3'
             do i = 1, ${ilim}$
                ! j = half_n_gc + 1
                ivar = i_vars(iv)
                i_f = ${if0}$ + 2 * i - 1
                j_f = ${jf0}$ + 2 * (half_n_gc + 1) - 1
                i_c = ${ic0}$ + i
                j_c = ${jc0}$ + (half_n_gc + 1)

                call prolong(f4%gminmod_theta_prolong, &
                     f4%uu(i_c, j_c, ivar, jq), &
                     f4%uu(i_c-1, j_c, ivar, jq), &
                     f4%uu(i_c+1, j_c, ivar, jq), &
                     f4%uu(i_c, j_c-1, ivar, jq), &
                     f4%uu(i_c, j_c+1, ivar, jq), fine)

                f4%uu(i_f  , j_f  , ivar, iq) = fine(1)
                f4%uu(i_f+1, j_f  , ivar, iq) = fine(2)
             end do
#:endif
          end do
       end if
    end do
#:enddef
#:elif NDIM == 3
#:def fyp_f2c_local_fine(face, ilim='f4%hbx(1)', jlim='f4%hbx(2)', &
    &klim='f4%hbx(3)', ic0=0, jc0=0, kc0=0, if0=0, jf0=0, kf0=0)
    ${LOOP('private(iq, jq, offset)')}$
    do n = f4%gc_f2c_local_iface(${face}$), f4%gc_f2c_local_iface(${face}$+1)-1
       iq     = f4%gc_f2c_local(1, n) + 1 + block_offset ! Fine block
       jq     = f4%gc_f2c_local(2, n) + 1 + block_offset ! coarse block
       offset(1:2) = f4%gc_f2c_local(3:4, n)     ! Offset

       ${LOOP('collapse(4) private(ivar, k_f, k_c, j_f, j_c, i_f, i_c, fine)')}$
       do iv = 1, n_vars
          do k = 1, ${klim}$
             do j = 1, ${jlim}$
                do i = 1,  ${ilim}$
                   ivar = i_vars(iv)
                   i_f = ${if0}$ + 2 * i - 1
                   j_f = ${jf0}$ + 2 * j - 1
                   k_f = ${kf0}$ + 2 * k - 1
                   i_c = ${ic0}$ + i
                   j_c = ${jc0}$ + j
                   k_c = ${kc0}$ + k
                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, k_c, ivar, jq), &
                        f4%uu(i_c-1, j_c, k_c, ivar, jq), &
                        f4%uu(i_c+1, j_c, k_c, ivar, jq), &
                        f4%uu(i_c, j_c-1, k_c, ivar, jq), &
                        f4%uu(i_c, j_c+1, k_c, ivar, jq), &
                        f4%uu(i_c, j_c, k_c-1, ivar, jq), &
                        f4%uu(i_c, j_c, k_c+1, ivar, jq), &
                        fine)
                   f4%uu(i_f  , j_f  , k_f,   ivar, iq) = fine(1)
                   f4%uu(i_f+1, j_f  , k_f,   ivar, iq) = fine(2)
                   f4%uu(i_f  , j_f+1, k_f,   ivar, iq) = fine(3)
                   f4%uu(i_f+1, j_f+1, k_f,   ivar, iq) = fine(4)
                   f4%uu(i_f  , j_f  , k_f+1, ivar, iq) = fine(5)
                   f4%uu(i_f+1, j_f  , k_f+1, ivar, iq) = fine(6)
                   f4%uu(i_f  , j_f+1, k_f+1, ivar, iq) = fine(7)
                   f4%uu(i_f+1, j_f+1, k_f+1, ivar, iq) = fine(8)
                end do
             end do
          end do
       end do

       if (odd_n_gc) then
          ${LOOP('collapse(3) private(ivar, k_f, k_c, j_f, j_c, i_f, i_c, fine)')}$
          do iv = 1, n_vars
#:if face == '0'
             do k = 1, ${klim}$
                do j = 1, ${jlim}$
                   ! i = 0
                   ivar = i_vars(iv)
                   i_f = ${if0}$ + 2 * 0 - 1
                   j_f = ${jf0}$ + 2 * j - 1
                   k_f = ${kf0}$ + 2 * k - 1
                   i_c = ${ic0}$ + 0
                   j_c = ${jc0}$ + j
                   k_c = ${kc0}$ + k

                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, k_c, ivar, jq), &
                        f4%uu(i_c-1, j_c, k_c, ivar, jq), &
                        f4%uu(i_c+1, j_c, k_c, ivar, jq), &
                        f4%uu(i_c, j_c-1, k_c, ivar, jq), &
                        f4%uu(i_c, j_c+1, k_c, ivar, jq), &
                        f4%uu(i_c, j_c, k_c-1, ivar, jq), &
                        f4%uu(i_c, j_c, k_c+1, ivar, jq), &
                        fine)

                   f4%uu(i_f+1, j_f  , k_f, ivar, iq) = fine(2)
                   f4%uu(i_f+1, j_f+1, k_f, ivar, iq) = fine(4)
                   f4%uu(i_f+1, j_f  , k_f+1, ivar, iq) = fine(6)
                   f4%uu(i_f+1, j_f+1, k_f+1, ivar, iq) = fine(8)
                end do
             end do
#:elif face == '1'
             do k = 1, ${klim}$
                do j = 1, ${jlim}$
                   ! i = half_n_gc + 1
                   ivar = i_vars(iv)
                   i_f = ${if0}$ + 2 * (half_n_gc + 1) - 1
                   j_f = ${jf0}$ + 2 * j - 1
                   k_f = ${kf0}$ + 2 * k - 1
                   i_c = ${ic0}$ + (half_n_gc + 1)
                   j_c = ${jc0}$ + j
                   k_c = ${kc0}$ + k

                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, k_c, ivar, jq), &
                        f4%uu(i_c-1, j_c, k_c, ivar, jq), &
                        f4%uu(i_c+1, j_c, k_c, ivar, jq), &
                        f4%uu(i_c, j_c-1, k_c, ivar, jq), &
                        f4%uu(i_c, j_c+1, k_c, ivar, jq), &
                        f4%uu(i_c, j_c, k_c-1, ivar, jq), &
                        f4%uu(i_c, j_c, k_c+1, ivar, jq), &
                        fine)

                   f4%uu(i_f  , j_f  , k_f, ivar, iq) = fine(1)
                   f4%uu(i_f  , j_f+1, k_f, ivar, iq) = fine(3)
                   f4%uu(i_f  , j_f  , k_f+1, ivar, iq) = fine(5)
                   f4%uu(i_f  , j_f+1, k_f+1, ivar, iq) = fine(7)
                end do
             end do
#:elif face == '2'
             do k = 1, ${klim}$
                do i = 1, ${ilim}$
                   ! j = 0
                   ivar = i_vars(iv)
                   i_f = ${if0}$ + 2 * i - 1
                   j_f = ${jf0}$ + 2 * 0 - 1
                   k_f = ${kf0}$ + 2 * k - 1
                   i_c = ${ic0}$ + i
                   j_c = ${jc0}$ + 0
                   k_c = ${kc0}$ + k

                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, k_c, ivar, jq), &
                        f4%uu(i_c-1, j_c, k_c, ivar, jq), &
                        f4%uu(i_c+1, j_c, k_c, ivar, jq), &
                        f4%uu(i_c, j_c-1, k_c, ivar, jq), &
                        f4%uu(i_c, j_c+1, k_c, ivar, jq), &
                        f4%uu(i_c, j_c, k_c-1, ivar, jq), &
                        f4%uu(i_c, j_c, k_c+1, ivar, jq), &
                        fine)

                   f4%uu(i_f  , j_f+1, k_f, ivar, iq) = fine(3)
                   f4%uu(i_f+1, j_f+1, k_f, ivar, iq) = fine(4)
                   f4%uu(i_f  , j_f+1, k_f+1, ivar, iq) = fine(7)
                   f4%uu(i_f+1, j_f+1, k_f+1, ivar, iq) = fine(8)
                end do
             end do
#:elif face == '3'
             do k = 1, ${klim}$
                do i = 1, ${ilim}$
                   ! j = half_n_gc + 1
                   ivar = i_vars(iv)
                   i_f = ${if0}$ + 2 * i - 1
                   j_f = ${jf0}$ + 2 * (half_n_gc + 1) - 1
                   k_f = ${kf0}$ + 2 * k - 1
                   i_c = ${ic0}$ + i
                   j_c = ${jc0}$ + (half_n_gc + 1)
                   k_c = ${kc0}$ + k

                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, k_c, ivar, jq), &
                        f4%uu(i_c-1, j_c, k_c, ivar, jq), &
                        f4%uu(i_c+1, j_c, k_c, ivar, jq), &
                        f4%uu(i_c, j_c-1, k_c, ivar, jq), &
                        f4%uu(i_c, j_c+1, k_c, ivar, jq), &
                        f4%uu(i_c, j_c, k_c-1, ivar, jq), &
                        f4%uu(i_c, j_c, k_c+1, ivar, jq), &
                        fine)

                   f4%uu(i_f  , j_f  , k_f, ivar, iq) = fine(1)
                   f4%uu(i_f+1, j_f  , k_f, ivar, iq) = fine(2)
                   f4%uu(i_f  , j_f  , k_f+1, ivar, iq) = fine(5)
                   f4%uu(i_f+1, j_f  , k_f+1, ivar, iq) = fine(6)
                end do
             end do
#:elif face == '4'
             do j = 1, ${jlim}$
                do i = 1, ${ilim}$
                   ! k = 0
                   ivar = i_vars(iv)
                   i_f = ${if0}$ + 2 * i - 1
                   j_f = ${jf0}$ + 2 * j - 1
                   k_f = ${kf0}$ + 2 * 0 - 1
                   i_c = ${ic0}$ + i
                   j_c = ${jc0}$ + j
                   k_c = ${kc0}$ + 0

                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, k_c, ivar, jq), &
                        f4%uu(i_c-1, j_c, k_c, ivar, jq), &
                        f4%uu(i_c+1, j_c, k_c, ivar, jq), &
                        f4%uu(i_c, j_c-1, k_c, ivar, jq), &
                        f4%uu(i_c, j_c+1, k_c, ivar, jq), &
                        f4%uu(i_c, j_c, k_c-1, ivar, jq), &
                        f4%uu(i_c, j_c, k_c+1, ivar, jq), &
                        fine)

                   f4%uu(i_f  , j_f  , k_f+1, ivar, iq) = fine(5)
                   f4%uu(i_f+1, j_f  , k_f+1, ivar, iq) = fine(6)
                   f4%uu(i_f  , j_f+1, k_f+1, ivar, iq) = fine(7)
                   f4%uu(i_f+1, j_f+1, k_f+1, ivar, iq) = fine(8)
                end do
             end do
#:elif face == '5'
             do j = 1, ${jlim}$
                do i = 1, ${ilim}$
                   ! k = half_n_gc + 1
                   ivar = i_vars(iv)
                   i_f = ${if0}$ + 2 * i - 1
                   j_f = ${jf0}$ + 2 * j - 1
                   k_f = ${kf0}$ + 2 * (half_n_gc + 1) - 1
                   i_c = ${ic0}$ + i
                   j_c = ${jc0}$ + j
                   k_c = ${kc0}$ + (half_n_gc + 1)

                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, k_c, ivar, jq), &
                        f4%uu(i_c-1, j_c, k_c, ivar, jq), &
                        f4%uu(i_c+1, j_c, k_c, ivar, jq), &
                        f4%uu(i_c, j_c-1, k_c, ivar, jq), &
                        f4%uu(i_c, j_c+1, k_c, ivar, jq), &
                        f4%uu(i_c, j_c, k_c-1, ivar, jq), &
                        f4%uu(i_c, j_c, k_c+1, ivar, jq), &
                        fine)

                   f4%uu(i_f  , j_f  , k_f, ivar, iq) = fine(1)
                   f4%uu(i_f+1, j_f  , k_f, ivar, iq) = fine(2)
                   f4%uu(i_f  , j_f+1, k_f, ivar, iq) = fine(3)
                   f4%uu(i_f+1, j_f+1, k_f, ivar, iq) = fine(4)
                end do
             end do
#:endif
          end do
       end if
    end do
#:enddef
#:endif

#:if NDIM == 2
#:def fyp_f2c_from_buf(face, ilim='f4%hbx(1)', jlim='f4%hbx(2)', if0=0, jf0=0)
    ${LOOP('private(iq, i_buf0)')}$
    do n = f4%gc_f2c_from_buf_iface(${face}$), f4%gc_f2c_from_buf_iface(${face}$+1)-1
       iq    = f4%gc_f2c_from_buf(1, n) + 1 + block_offset ! Fine block
       i_buf0 = f4%gc_f2c_from_buf(2, n) * n_vars

       ${LOOP('collapse(3) private(ivar, j_f, i_f, i_buf)')}$
       do iv = 1, n_vars
          do j = 1, ${jlim}$
             do i = 1, ${ilim}$
                ivar = i_vars(iv)
                i_f = ${if0}$ + 2 * i - 1
                j_f = ${jf0}$ + 2 * j - 1

                i_buf = i_buf0 + 4 * ix_offset3(iv, j, i, ${jlim}$, ${ilim}$)
                f4%uu(i_f  , j_f  , ivar, iq) = f4%recv_buffer(i_buf+1)
                f4%uu(i_f+1, j_f  , ivar, iq) = f4%recv_buffer(i_buf+2)
                f4%uu(i_f  , j_f+1, ivar, iq) = f4%recv_buffer(i_buf+3)
                f4%uu(i_f+1, j_f+1, ivar, iq) = f4%recv_buffer(i_buf+4)
             end do
          end do
       end do

       if (odd_n_gc) then
          i_buf0 = i_buf0 + 4 * n_vars * ${jlim}$ * ${ilim}$
#:if face == '0'
          ${LOOP('collapse(2) private(ivar, i_f, j_f, i_buf)')}$
          do iv = 1, n_vars
             do j = 1, ${jlim}$
                ivar = i_vars(iv)
                i_f = -f4%n_gc + 1
                j_f = 2 * j - 1

                i_buf = i_buf0 + 2 * ix_offset2(iv, j, ${jlim}$)
                f4%uu(i_f, j_f  , ivar, iq) = f4%recv_buffer(i_buf+1)
                f4%uu(i_f, j_f+1, ivar, iq) = f4%recv_buffer(i_buf+2)
             end do
          end do
#:elif face == '1'
          ${LOOP('collapse(2) private(ivar, i_f, j_f, i_buf)')}$
          do iv = 1, n_vars
             do j = 1, ${jlim}$
                ivar = i_vars(iv)
                i_f = f4%bx(1) + f4%n_gc
                j_f = 2 * j - 1

                i_buf = i_buf0 + 2 * ix_offset2(iv, j, ${jlim}$)
                f4%uu(i_f, j_f  , ivar, iq) = f4%recv_buffer(i_buf+1)
                f4%uu(i_f, j_f+1, ivar, iq) = f4%recv_buffer(i_buf+2)
             end do
          end do
#:elif face == '2'
          ${LOOP('collapse(2) private(ivar, i_f, j_f, i_buf)')}$
          do iv = 1, n_vars
             do i = 1, ${ilim}$
                ivar = i_vars(iv)
                i_f = 2 * i - 1
                j_f = -f4%n_gc + 1

                i_buf = i_buf0 + 2 * ix_offset2(iv, i, ${ilim}$)
                f4%uu(i_f  , j_f  , ivar, iq) = f4%recv_buffer(i_buf+1)
                f4%uu(i_f+1, j_f  , ivar, iq) = f4%recv_buffer(i_buf+2)
             end do
          end do
#:elif face == '3'
          ${LOOP('collapse(2) private(ivar, i_f, j_f, i_buf)')}$
          do iv = 1, n_vars
             do i = 1, ${ilim}$
                ivar = i_vars(iv)
                i_f = 2 * i - 1
                j_f = f4%bx(2) + f4%n_gc

                i_buf = i_buf0 + 2 * ix_offset2(iv, i, ${ilim}$)
                f4%uu(i_f  , j_f  , ivar, iq) = f4%recv_buffer(i_buf+1)
                f4%uu(i_f+1, j_f  , ivar, iq) = f4%recv_buffer(i_buf+2)
             end do
          end do
#:endif
       end if
    end do
#:enddef
#:elif NDIM == 3
#:def fyp_f2c_from_buf(face, ilim='f4%hbx(1)', jlim='f4%hbx(2)', &
    & klim='f4%hbx(3)', if0=0, jf0=0, kf0=0)
    ${LOOP('private(iq, i_buf0)')}$
    do n = f4%gc_f2c_from_buf_iface(${face}$), f4%gc_f2c_from_buf_iface(${face}$+1)-1
       iq    = f4%gc_f2c_from_buf(1, n) + 1 + block_offset ! Fine block
       i_buf0 = f4%gc_f2c_from_buf(2, n) * n_vars

       ${LOOP('collapse(4) private(ivar, k_f, j_f, i_f, i_buf)')}$
       do iv = 1, n_vars
          do k = 1, ${klim}$
             do j = 1, ${jlim}$
                do i = 1, ${ilim}$
                   ivar = i_vars(iv)
                   i_f = ${if0}$ + 2 * i - 1
                   j_f = ${jf0}$ + 2 * j - 1
                   k_f = ${kf0}$ + 2 * k - 1

                   i_buf = i_buf0 + 8 * ix_offset4(iv, k, j, i, &
                        ${klim}$, ${jlim}$, ${ilim}$)
                   f4%uu(i_f  , j_f  , k_f, ivar, iq) = f4%recv_buffer(i_buf+1)
                   f4%uu(i_f+1, j_f  , k_f, ivar, iq) = f4%recv_buffer(i_buf+2)
                   f4%uu(i_f  , j_f+1, k_f, ivar, iq) = f4%recv_buffer(i_buf+3)
                   f4%uu(i_f+1, j_f+1, k_f, ivar, iq) = f4%recv_buffer(i_buf+4)
                   f4%uu(i_f  , j_f  , k_f+1, ivar, iq) = f4%recv_buffer(i_buf+5)
                   f4%uu(i_f+1, j_f  , k_f+1, ivar, iq) = f4%recv_buffer(i_buf+6)
                   f4%uu(i_f  , j_f+1, k_f+1, ivar, iq) = f4%recv_buffer(i_buf+7)
                   f4%uu(i_f+1, j_f+1, k_f+1, ivar, iq) = f4%recv_buffer(i_buf+8)
                end do
             end do
          end do
       end do

       if (odd_n_gc) then
          i_buf0 = i_buf0 + 8 * n_vars * ${klim}$ * ${jlim}$ * ${ilim}$
#:if face in ['0', '1']
          ${LOOP('collapse(3) private(ivar, k_f, j_f, i_f, i_buf)')}$
          do iv = 1, n_vars
             do k = 1, ${klim}$
                do j = 1, ${jlim}$
                   ivar = i_vars(iv)
#:if face == '0'
                   i_f = -f4%n_gc + 1
#:else
                   i_f = f4%bx(1) + f4%n_gc
#:endif
                   j_f = 2 * j - 1
                   k_f = 2 * k - 1

                   i_buf = i_buf0 + 4 * ix_offset3(iv, k, j, ${klim}$, ${jlim}$)
                   f4%uu(i_f, j_f  , k_f, ivar, iq) = f4%recv_buffer(i_buf+1)
                   f4%uu(i_f, j_f+1, k_f, ivar, iq) = f4%recv_buffer(i_buf+2)
                   f4%uu(i_f, j_f  , k_f+1, ivar, iq) = f4%recv_buffer(i_buf+3)
                   f4%uu(i_f, j_f+1, k_f+1, ivar, iq) = f4%recv_buffer(i_buf+4)
                end do
             end do
          end do
#:elif face in ['2', '3']
          ${LOOP('collapse(3) private(ivar, k_f, j_f, i_f, i_buf)')}$
          do iv = 1, n_vars
             do k = 1, ${klim}$
                do i = 1, ${ilim}$
                   ivar = i_vars(iv)
                   i_f = 2 * i - 1
#:if face == '2'
                   j_f = -f4%n_gc + 1
#:else
                   j_f = f4%bx(2) + f4%n_gc
#:endif
                   k_f = 2 * k - 1

                   i_buf = i_buf0 + 4 * ix_offset3(iv, k, i, ${klim}$, ${ilim}$)
                   f4%uu(i_f  , j_f  , k_f, ivar, iq) = f4%recv_buffer(i_buf+1)
                   f4%uu(i_f+1, j_f  , k_f, ivar, iq) = f4%recv_buffer(i_buf+2)
                   f4%uu(i_f  , j_f  , k_f+1, ivar, iq) = f4%recv_buffer(i_buf+3)
                   f4%uu(i_f+1, j_f  , k_f+1, ivar, iq) = f4%recv_buffer(i_buf+4)
                end do
             end do
          end do
#:elif face in ['4', '5']
          ${LOOP('collapse(3) private(ivar, k_f, j_f, i_f, i_buf)')}$
          do iv = 1, n_vars
             do j = 1, ${jlim}$
                do i = 1, ${ilim}$
                   ivar = i_vars(iv)
                   i_f = 2 * i - 1
                   j_f = 2 * j - 1
#:if face == '4'
                   k_f = -f4%n_gc + 1
#:else
                   k_f = f4%bx(3) + f4%n_gc
#:endif

                   i_buf = i_buf0 + 4 * ix_offset3(iv, j, i, ${jlim}$, ${ilim}$)
                   f4%uu(i_f  , j_f  , k_f, ivar, iq) = f4%recv_buffer(i_buf+1)
                   f4%uu(i_f+1, j_f  , k_f, ivar, iq) = f4%recv_buffer(i_buf+2)
                   f4%uu(i_f  , j_f+1, k_f, ivar, iq) = f4%recv_buffer(i_buf+3)
                   f4%uu(i_f+1, j_f+1, k_f, ivar, iq) = f4%recv_buffer(i_buf+4)
                end do
             end do
          end do
#:endif
       end if
    end do
#:enddef
#:endif

    ${PARALLEL(COPYIN('i_vars'))}$ ${DEFAULT_PRESENT()}$

    ! ----------------------------------------
    ! Fill fine side of local coarse-to-fine boundaries
    ! ----------------------------------------

#:if NDIM == 2
    @:fyp_f2c_local_fine(0, ilim=half_n_gc, ic0=f4%bx(1)-half_n_gc, &
         &jc0=offset(1)*f4%hbx(2), if0=-2*half_n_gc)
    @:fyp_f2c_local_fine(1, ilim=half_n_gc, &
         &jc0=offset(1)*f4%hbx(2), if0=f4%bx(1))
    @:fyp_f2c_local_fine(2, jlim=half_n_gc, ic0=offset(1)*f4%hbx(1), &
         &jc0=f4%bx(2)-half_n_gc, jf0=-2*half_n_gc)
    @:fyp_f2c_local_fine(3, jlim=half_n_gc, &
         &ic0=offset(1)*f4%hbx(1), jf0=f4%bx(2))
#:elif NDIM == 3
    @:fyp_f2c_local_fine(0, ilim=half_n_gc, ic0=f4%bx(1)-half_n_gc, &
         &jc0=offset(1)*f4%hbx(2), kc0=offset(2)*f4%hbx(3), if0=-2*half_n_gc)
    @:fyp_f2c_local_fine(1, ilim=half_n_gc, &
         &jc0=offset(1)*f4%hbx(2), kc0=offset(2)*f4%hbx(3), if0=f4%bx(1))
    @:fyp_f2c_local_fine(2, jlim=half_n_gc, ic0=offset(1)*f4%hbx(1), &
         &jc0=f4%bx(2)-half_n_gc, kc0=offset(2)*f4%hbx(3), jf0=-2*half_n_gc)
    @:fyp_f2c_local_fine(3, jlim=half_n_gc, &
         &ic0=offset(1)*f4%hbx(1), kc0=offset(2)*f4%hbx(3), jf0=f4%bx(2))
    @:fyp_f2c_local_fine(4, klim=half_n_gc, ic0=offset(1)*f4%hbx(1), &
         &jc0=offset(2)*f4%hbx(2), kc0=f4%bx(3)-half_n_gc, kf0=-2*half_n_gc)
    @:fyp_f2c_local_fine(5, klim=half_n_gc, &
         &ic0=offset(1)*f4%hbx(1), jc0=offset(2)*f4%hbx(2), kf0=f4%bx(3))
#:endif

    ! ----------------------------------------
    ! Fill fine side of nonlocal coarse-to-fine boundaries
    ! ----------------------------------------

#:if NDIM == 2
    @:fyp_f2c_from_buf(0, ilim=half_n_gc, if0=-2*half_n_gc)
    @:fyp_f2c_from_buf(1, ilim=half_n_gc, if0=f4%bx(1))
    @:fyp_f2c_from_buf(2, jlim=half_n_gc, jf0=-2*half_n_gc)
    @:fyp_f2c_from_buf(3, jlim=half_n_gc, jf0=f4%bx(2))
#:elif NDIM == 3
    @:fyp_f2c_from_buf(0, ilim=half_n_gc, if0=-2*half_n_gc)
    @:fyp_f2c_from_buf(1, ilim=half_n_gc, if0=f4%bx(1))
    @:fyp_f2c_from_buf(2, jlim=half_n_gc, jf0=-2*half_n_gc)
    @:fyp_f2c_from_buf(3, jlim=half_n_gc, jf0=f4%bx(2))
    @:fyp_f2c_from_buf(4, klim=half_n_gc, kf0=-2*half_n_gc)
    @:fyp_f2c_from_buf(5, klim=half_n_gc, kf0=f4%bx(2))
#:endif
    ${END_PARALLEL()}$

  end subroutine fill_ghostcells_round_two

  !> Update ghost cells for selected variables
  subroutine f4_update_ghostcells(f4, n_vars, i_vars, block_offset)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer, intent(in)          :: block_offset
    real(dp)                     :: t0, t1

    t1 = MPI_Wtime()
    call fill_ghostcell_buffers_round_one(f4, n_vars, i_vars, block_offset)
    t0 = MPI_Wtime()
    f4%wtime_gc_fill_buff_round1 = f4%wtime_gc_fill_buff_round1 + t0 - t1

    call f4_exchange_buffers(f4)
    t1 = MPI_Wtime()
    f4%wtime_exchange_buffers = f4%wtime_exchange_buffers + t1 - t0

    call fill_ghostcells_round_one(f4, n_vars, i_vars, block_offset)

    t0 = MPI_Wtime()
    f4%wtime_gc_fill_round1 = f4%wtime_gc_fill_round1 + t0 - t1

    ! Do coarse-to-fine refinement boundaries last, so that ghost cells
    ! required for interpolation have been filled
    call fill_ghostcell_buffers_round_two(f4, n_vars, i_vars, block_offset)

    t1 = MPI_Wtime()
    f4%wtime_gc_fill_buff_round2 = f4%wtime_gc_fill_buff_round2 + t1 - t0

    call f4_exchange_buffers(f4)

    t0 = MPI_Wtime()
    f4%wtime_exchange_buffers = f4%wtime_exchange_buffers + t0 - t1

    call fill_ghostcells_round_two(f4, n_vars, i_vars, block_offset)

    t1 = MPI_Wtime()
    f4%wtime_gc_fill_round2 = f4%wtime_gc_fill_round2 + t1 - t0

  end subroutine f4_update_ghostcells

  !> Refine the mesh according to f4%refinement_flags
  subroutine f4_adjust_refinement(f4, load_imbalance_threshold)
    type(foap4_t), intent(inout) :: f4
    !> Threshold for load imbalance, for example 1.1 could be reasonable
    real(dp), intent(in)         :: load_imbalance_threshold
    integer                      :: n_blocks_new, n_blocks_old, iv
    integer                      :: has_changed, i_srl, i_refine, i_coarsen
    integer, allocatable         :: srl(:, :), refine(:, :), coarsen(:, :)
    integer                      :: n, i, j, i_c, j_c, i_f, j_f, n_old
    integer                      :: offset_copy
    integer                      :: i_from, i_to, i_ch
    real(dp)                     :: t0, t1, load_imbalance
    real(fp)                     :: fine(2**NDIM)
#:if NDIM == 3
    integer                      :: k, k_c, k_f
#:endif

    t0 = MPI_Wtime()

    n_blocks_old = f4%n_blocks
    call pw_adjust_refinement(f4%pw, f4%n_blocks, &
         f4%refinement_flags(1:f4%n_blocks), has_changed)

    t1 = MPI_Wtime()
    f4%wtime_adjust_ref_p4est = f4%wtime_adjust_ref_p4est + t1 - t0

    if (has_changed == 0) return

    n_blocks_new = pw_get_num_local_quadrants(f4%pw)

    call copy_blocks_to_end(f4, n_blocks_old, n_blocks_new, offset_copy)
    call f4_set_quadrants(f4)

    allocate(srl(2, n_blocks_new), refine(2, n_blocks_new), coarsen(2, n_blocks_new))

    i_srl     = 0
    i_refine  = 0
    i_coarsen = 0
    n_old     = offset_copy + 1
    n         = 1
    do while (n <= n_blocks_new)
       select case (f4%block_level(n) - f4%block_level(n_old))
       case (0)
          ! Same refinement level
          i_srl = i_srl + 1
          srl(:, i_srl) = [n_old, n]
          n = n + 1
          n_old = n_old + 1
       case (1)
          ! Block has been refined
          i_refine = i_refine + 1
          refine(:, i_refine) = [n_old, n]
          n = n + 2**NDIM
          n_old = n_old + 1
       case (-1)
          ! Block has been coarsened
          i_coarsen = i_coarsen + 1
          coarsen(:, i_coarsen) = [n_old, n]
          n = n + 1
          n_old = n_old + 2**NDIM
       case default
          error stop "Refinement: difference in levels > 1"
       end select
    end do

    if (n_old /= offset_copy + n_blocks_old + 1) &
         error stop "Refinement: loops did not end simultaneously"

    t1 = MPI_Wtime()

    ${ENTER_DATA_COPYIN('srl, refine, coarsen')}$

    ! Copy on device
    ${PARALLEL_LOOP('collapse(NDIM+2) private(i_from, i_to)')}$ ${DEFAULT_PRESENT()}$ ${ASYNC()}$
    do n = 1, i_srl
       do iv = 1, f4%n_vars
          do @{KJI_LOOP_array_to_array(f4%ilo, f4%ihi)}@
             i_from = srl(1, n)
             i_to = srl(2, n)
             f4%uu(${IJK}$, iv, i_to) = f4%uu(${IJK}$, iv, i_from)
          end do; ${KJI_CLOSE_LOOP}$
       end do
    end do

    ! Refine on device
#:if NDIM == 2
    ${PARALLEL_LOOP('collapse(NDIM+3) private(i_from, i_to, j_c, j_f, i_c, i_f, fine)')}$ ${DEFAULT_PRESENT()}$ ${ASYNC()}$
    do n = 1, i_refine
       do i_ch = 1, 2**NDIM
          do iv = 1, f4%n_vars
             do j = 1, f4%hbx(2)
                do i = 1, f4%hbx(1)
                   i_from = refine(1, n)
                   i_to = refine(2, n)
                   j_c = j + f4_child_offset(2, i_ch) * f4%hbx(2)
                   i_c = i + f4_child_offset(1, i_ch) * f4%hbx(1)
                   j_f = 2 * j - 1
                   i_f = 2 * i - 1

                   call prolong(f4%gminmod_theta_prolong, &
                        f4%uu(i_c, j_c, iv, i_from), &
                        f4%uu(i_c-1, j_c, iv, i_from), &
                        f4%uu(i_c+1, j_c, iv, i_from), &
                        f4%uu(i_c, j_c-1, iv, i_from), &
                        f4%uu(i_c, j_c+1, iv, i_from), &
                        fine)

                   f4%uu(i_f,   j_f, iv, i_to+i_ch-1)   = fine(1)
                   f4%uu(i_f+1, j_f, iv, i_to+i_ch-1)   = fine(2)
                   f4%uu(i_f,   j_f+1, iv, i_to+i_ch-1) = fine(3)
                   f4%uu(i_f+1, j_f+1, iv, i_to+i_ch-1) = fine(4)
                end do
             end do
          end do
       end do
    end do
#:elif NDIM == 3
    ${PARALLEL_LOOP('collapse(NDIM+3) private(i_from, i_to, k_c, k_f, j_c, j_f, i_c, i_f, fine)')}$ ${DEFAULT_PRESENT()}$ ${ASYNC()}$
    do n = 1, i_refine
       do i_ch = 1, 2**NDIM
          do iv = 1, f4%n_vars
             do k = 1, f4%hbx(3)
                do j = 1, f4%hbx(2)
                   do i = 1, f4%hbx(1)
                      i_from = refine(1, n)
                      i_to = refine(2, n)
                      k_c = k + f4_child_offset(3, i_ch) * f4%hbx(3)
                      j_c = j + f4_child_offset(2, i_ch) * f4%hbx(2)
                      i_c = i + f4_child_offset(1, i_ch) * f4%hbx(1)
                      k_f = 2 * k - 1
                      j_f = 2 * j - 1
                      i_f = 2 * i - 1

                      call prolong(f4%gminmod_theta_prolong, &
                           f4%uu(i_c, j_c, k_c, iv, i_from), &
                           f4%uu(i_c-1, j_c, k_c, iv, i_from), &
                           f4%uu(i_c+1, j_c, k_c, iv, i_from), &
                           f4%uu(i_c, j_c-1, k_c, iv, i_from), &
                           f4%uu(i_c, j_c+1, k_c, iv, i_from), &
                           f4%uu(i_c, j_c, k_c-1, iv, i_from), &
                           f4%uu(i_c, j_c, k_c+1, iv, i_from), &
                           fine)

                      f4%uu(i_f  , j_f  , k_f,   iv, i_to+i_ch-1) = fine(1)
                      f4%uu(i_f+1, j_f  , k_f,   iv, i_to+i_ch-1) = fine(2)
                      f4%uu(i_f  , j_f+1, k_f,   iv, i_to+i_ch-1) = fine(3)
                      f4%uu(i_f+1, j_f+1, k_f,   iv, i_to+i_ch-1) = fine(4)
                      f4%uu(i_f  , j_f  , k_f+1, iv, i_to+i_ch-1) = fine(5)
                      f4%uu(i_f+1, j_f  , k_f+1, iv, i_to+i_ch-1) = fine(6)
                      f4%uu(i_f  , j_f+1, k_f+1, iv, i_to+i_ch-1) = fine(7)
                      f4%uu(i_f+1, j_f+1, k_f+1, iv, i_to+i_ch-1) = fine(8)
                   end do
                end do
             end do
          end do
       end do
    end do
#:endif

    ! Coarsen on device
    #:if NDIM == 2
    ${PARALLEL_LOOP('collapse(NDIM+3) private(i_from, i_to, j_c, j_f, i_c, i_f)')}$ ${DEFAULT_PRESENT()}$ ${ASYNC()}$
    do n = 1, i_coarsen
       do i_ch = 1, 2**NDIM
          do iv = 1, f4%n_vars
             do j = 1, f4%hbx(2)
                do i = 1, f4%hbx(1)
                   i_from = coarsen(1, n)
                   i_to = coarsen(2, n)
                   j_c = j + f4_child_offset(2, i_ch) * f4%hbx(2)
                   j_f = 2 * j - 1
                   i_c = i + f4_child_offset(1, i_ch) * f4%hbx(1)
                   i_f = 2 * i - 1

                   f4%uu(i_c, j_c, iv, i_to) = 0.25_fp * (&
                        f4%uu(i_f,   j_f, iv, i_from+i_ch-1) + &
                        f4%uu(i_f+1, j_f, iv, i_from+i_ch-1) + &
                        f4%uu(i_f,   j_f+1, iv, i_from+i_ch-1) + &
                        f4%uu(i_f+1, j_f+1, iv, i_from+i_ch-1))
                end do
             end do
          end do
       end do
    end do
#:elif NDIM == 3
    ${PARALLEL_LOOP('collapse(NDIM+3) private(i_from, i_to, k_c, k_f, j_c, j_f, i_c, i_f)')}$ ${DEFAULT_PRESENT()}$ ${ASYNC()}$
    do n = 1, i_coarsen
       do i_ch = 1, 2**NDIM
          do iv = 1, f4%n_vars
             do k = 1, f4%hbx(3)
                do j = 1, f4%hbx(2)
                   do i = 1, f4%hbx(1)
                      i_from = coarsen(1, n)
                      i_to = coarsen(2, n)
                      k_c = k + f4_child_offset(3, i_ch) * f4%hbx(3)
                      k_f = 2 * k - 1
                      j_c = j + f4_child_offset(2, i_ch) * f4%hbx(2)
                      j_f = 2 * j - 1
                      i_c = i + f4_child_offset(1, i_ch) * f4%hbx(1)
                      i_f = 2 * i - 1

                      f4%uu(i_c, j_c, k_c, iv, i_to) = 0.125_fp * (&
                           f4%uu(i_f  , j_f  , k_f,   iv, i_from+i_ch-1) + &
                           f4%uu(i_f+1, j_f  , k_f,   iv, i_from+i_ch-1) + &
                           f4%uu(i_f  , j_f+1, k_f,   iv, i_from+i_ch-1) + &
                           f4%uu(i_f+1, j_f+1, k_f,   iv, i_from+i_ch-1) + &
                           f4%uu(i_f  , j_f  , k_f+1, iv, i_from+i_ch-1) + &
                           f4%uu(i_f+1, j_f  , k_f+1, iv, i_from+i_ch-1) + &
                           f4%uu(i_f  , j_f+1, k_f+1, iv, i_from+i_ch-1) + &
                           f4%uu(i_f+1, j_f+1, k_f+1, iv, i_from+i_ch-1))
                   end do
                end do
             end do
          end do
       end do
    end do
#:endif

    ${WAIT()}$
    ${EXIT_DATA_DELETE('srl, refine, coarsen')}$

    t0 = MPI_Wtime()
    f4%wtime_adjust_ref_foap4 = f4%wtime_adjust_ref_foap4 + t0 - t1

    call f4_get_load_imbalance(f4, load_imbalance)

    if (load_imbalance > load_imbalance_threshold) call f4_partition(f4)

    t0 = MPI_Wtime()
    call update_ghostcell_pattern(f4)
    call set_face_data_storage(f4)
    if (associated(f4%bc_callback)) call f4%bc_callback(f4)
    t1 = MPI_Wtime()
    f4%wtime_adjust_ref_foap4 = f4%wtime_adjust_ref_foap4 + t1 - t0

  end subroutine f4_adjust_refinement

  !> Temporarily copy blocks to the end of the block array
  subroutine copy_blocks_to_end(f4, n_blocks_old, n_blocks_new, offset_copy)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_blocks_old
    integer, intent(in)          :: n_blocks_new
    integer, intent(out)         :: offset_copy
    integer                      :: n, ${IJK}$, iv

    offset_copy = max(n_blocks_old, n_blocks_new)

    if (offset_copy + n_blocks_old > size(f4%uu, NDIM+2)) then
       write(error_unit, "(A,I0,A,I0)") "ERROR: allocated_blocks = ", &
            size(f4%uu, NDIM+2), ", copy requires ", offset_copy + n_blocks_old
       error stop "Not enough block memory for copying"
    end if

    ! Copy block level on host
    do n = 1, n_blocks_old
       f4%block_level(offset_copy+n) = f4%block_level(n)
    end do

    ! Copy block solution data on device
    ${PARALLEL_LOOP('collapse(NDIM+2)')}$ ${DEFAULT_PRESENT()}$
    do n = 1, n_blocks_old
       do iv = 1, f4%n_vars
          do @{KJI_LOOP_array_to_array(f4%ilo, f4%ihi)}@
             f4%uu(${IJK}$, iv, offset_copy+n) = f4%uu(${IJK}$, iv, n)
          end do; ${KJI_CLOSE_LOOP}$
       end do
    end do
  end subroutine copy_blocks_to_end

  !> Method for prolongation (interpolation) of a coarse block to its children
#:if NDIM == 2
  subroutine prolong(theta, center, xlo, xhi, ylo, yhi, fine)
    @{ROUTINE_SEQ()}@
    real(fp), intent(in)  :: theta
    real(fp), intent(in)  :: center ! Center value
    real(fp), intent(in)  :: xlo, xhi ! x-neighbors (-1, +1)
    real(fp), intent(in)  :: ylo, yhi ! y-neighbors (-1, +1)
    real(fp), intent(out) :: fine(2**NDIM)
    real(fp)              :: f(0:NDIM), slopes_a(NDIM), slopes_b(NDIM)

    f(0) = center
    slopes_a(1) = center - xlo
    slopes_a(2) = center - ylo
    slopes_b(1) = xhi - center
    slopes_b(2) = yhi - center

    f(1:) = 0.25_fp * limiter_gminmod(slopes_a, slopes_b, theta)

    fine(1) = f(0) - f(1) - f(2)
    fine(2) = f(0) + f(1) - f(2)
    fine(3) = f(0) - f(1) + f(2)
    fine(4) = f(0) + f(1) + f(2)
  end subroutine prolong
#:elif NDIM == 3
  subroutine prolong(theta, center, xlo, xhi, ylo, yhi, zlo, zhi, fine)
    @{ROUTINE_SEQ()}@
    real(fp), intent(in)  :: theta
    real(fp), intent(in)  :: center   ! Center value
    real(fp), intent(in)  :: xlo, xhi ! x-neighbors (-1, +1)
    real(fp), intent(in)  :: ylo, yhi ! y-neighbors (-1, +1)
    real(fp), intent(in)  :: zlo, zhi ! z-neighbors (-1, +1)
    real(fp), intent(out) :: fine(2**NDIM)
    real(fp)              :: f(0:NDIM), slopes_a(NDIM), slopes_b(NDIM)

    f(0) = center
    slopes_a(1) = center - xlo
    slopes_a(2) = center - ylo
    slopes_a(3) = center - zlo
    slopes_b(1) = xhi - center
    slopes_b(2) = yhi - center
    slopes_b(3) = zhi - center

    f(1:) = 0.25_fp * limiter_gminmod(slopes_a, slopes_b, theta)

    fine(1) = f(0) - f(1) - f(2) - f(3)
    fine(2) = f(0) + f(1) - f(2) - f(3)
    fine(3) = f(0) - f(1) + f(2) - f(3)
    fine(4) = f(0) + f(1) + f(2) - f(3)
    fine(5) = f(0) - f(1) - f(2) + f(3)
    fine(6) = f(0) + f(1) - f(2) + f(3)
    fine(7) = f(0) - f(1) + f(2) + f(3)
    fine(8) = f(0) + f(1) + f(2) + f(3)
  end subroutine prolong
#:endif

  !> Generalized minmod limiter. The parameter theta controls how dissipative
  !> the limiter is, with 1 corresponding to the minmod limiter and 2 to the MC
  !> limiter.
  elemental function limiter_gminmod(a, b, theta) result(phi)
    @{ROUTINE_SEQ()}@
    real(fp), intent(in) :: a, b, theta
    real(fp)             :: phi

    if (a * b > 0) then
       phi = sign(min(theta*abs(a), theta*abs(b), 0.5_fp*abs(a+b)), a)
    else
       phi = 0.0_fp
    end if
  end function limiter_gminmod

  !> Get load imbalance, defined as the ratio of max_blocks/avg_blocks,
  !> normalized by the minimum achievable imbalance for the given block count
  !> and number of ranks. Returns 1.0 for a perfect distribution.
  subroutine f4_get_load_imbalance(f4, load_imbalance)
    type(foap4_t), intent(in) :: f4
    real(dp), intent(out)     :: load_imbalance
    integer(c_int64_t)        :: gfq(0:f4%mpisize)
    integer(c_int64_t)        :: blocks_per_rank(f4%mpisize)
    integer(c_int64_t)        :: max_blocks, total_blocks
    integer(c_int64_t)        :: optimal_max_blocks

    call pw_get_global_filling_curve(f4%pw, gfq)

    total_blocks = gfq(f4%mpisize) - gfq(0)
    blocks_per_rank(:) = gfq(1:f4%mpisize) - gfq(0:f4%mpisize-1)
    max_blocks = maxval(blocks_per_rank)

    if (total_blocks > 0) then
       ! The best possible distribution assigns ceil(total/nprocs) to some
       ! ranks and floor(total/nprocs) to the rest
       optimal_max_blocks = (total_blocks + f4%mpisize - 1) / f4%mpisize
       load_imbalance = real(max_blocks, dp) / optimal_max_blocks
    else
       load_imbalance = 1.0_dp
    end if

  end subroutine f4_get_load_imbalance

  !> Partition the blocks over the MPI ranks
  subroutine f4_partition(f4)
    type(foap4_t), intent(inout)   :: f4
    integer                        :: n_changed_global
    integer                        :: n_blocks_old, n_blocks_new
    integer(c_int64_t)             :: gfq_old(0:f4%mpisize)
    integer(c_int64_t)             :: gfq_new(0:f4%mpisize)
    integer                        :: offset_copy, ierr
    integer(MPI_COUNT_KIND)        :: count, dsize
    integer, parameter             :: tag = 0
    integer                        :: gfq_src(0:f4%mpisize)
    integer                        :: gfq_dest(0:f4%mpisize)
    integer                        :: dest_begin, dest_end
    integer                        :: src_begin, src_end
    integer                        :: gend, gbegin, n_blocks_transfer
    integer                        :: first_sender, last_sender
    integer                        :: first_receiver, last_receiver
    integer                        :: n_recv, n_send
    integer                        :: rank, block_ix
    type(MPI_Request)              :: send_req(f4%max_requests)
    type(MPI_Request)              :: recv_req(f4%max_requests)
    real(dp)                       :: t0, t1

    t0 = MPI_Wtime()
    n_blocks_old = f4%n_blocks
    call pw_partition(f4%pw, n_changed_global, gfq_old, gfq_new)

    ! No need to do anything if the global number of blocks shipped is zero
    if (n_changed_global == 0) return

    ! Convert to one-based indexing and to default integer type
    gfq_src(:) = int(gfq_old(:) + 1)
    gfq_dest(:) = int(gfq_new(:) + 1)

    ! Copy blocks
    n_blocks_new = pw_get_num_local_quadrants(f4%pw)
    call copy_blocks_to_end(f4, n_blocks_old, n_blocks_new, offset_copy)

    ! Size of a block
    dsize = product(f4%bx + 2*f4%n_gc) * f4%n_vars

    ! First and last block this rank owns after partitioning
    dest_begin = gfq_dest(f4%mpirank)
    dest_end   = gfq_dest(f4%mpirank + 1)

    ! First and last block this rank owned before partitioning
    src_begin = gfq_src(f4%mpirank)
    src_end   = gfq_src(f4%mpirank + 1)

    n_recv = 0
    n_send = 0

    if (dest_end > dest_begin) then
       ! Find index of first/last sender, -1 to get zero-based index
       first_sender = find_bracket(f4%mpisize+1, gfq_src, dest_begin) - 1
       last_sender = find_bracket(f4%mpisize+1, gfq_src, dest_end-1) - 1
       gend = dest_begin
       block_ix = 1

       do rank = first_sender, last_sender
          gbegin = gend
          gend = min(dest_end, gfq_src(rank + 1))
          n_blocks_transfer = gend - gbegin

          if (n_blocks_transfer > 0) then
             count = dsize * n_blocks_transfer
             call mpi_irecv_wrapper(f4%uu(@{DTIMES(:)}@, :, block_ix), &
                  count, rank, tag, f4%mpicomm, recv_req, n_recv, ierr)
             block_ix = block_ix + n_blocks_transfer
          end if
       end do
    end if

    if (src_end > src_begin) then
       ! Find index of first/last receiver, -1 to get zero-based index
       first_receiver = find_bracket(f4%mpisize+1, gfq_dest, src_begin) - 1
       last_receiver = find_bracket(f4%mpisize+1, gfq_dest, src_end-1) - 1
       gend = src_begin
       block_ix = offset_copy + 1

       do rank = first_receiver, last_receiver
          gbegin = gend
          gend = min(src_end, gfq_dest(rank + 1))
          n_blocks_transfer = gend - gbegin

          if (n_blocks_transfer > 0) then
             count = dsize * n_blocks_transfer
             call mpi_isend_wrapper(f4%uu(@{DTIMES(:)}@, :, block_ix), &
                  count, rank, tag, f4%mpicomm, send_req, n_send, ierr)
             block_ix = block_ix + n_blocks_transfer
          end if
       end do
    end if

    if (n_recv > 0) then
       call MPI_Waitall(n_recv, recv_req(1:n_recv), MPI_STATUSES_IGNORE, ierr)
    end if
    if (n_send > 0) then
       call MPI_Waitall(n_send, send_req(1:n_send), MPI_STATUSES_IGNORE, ierr)
    end if

    call f4_set_quadrants(f4)

    t1 = MPI_Wtime()
    f4%wtime_partition = f4%wtime_partition + t1 - t0
  end subroutine f4_partition

  !> Performs a binary search for the index ix such that:
  !> array(ix) <= key < array(ix+1)
  !>
  !> Returns:
  !>   ix    : the bracketing index (1 <= ix < n-1)
  !>   ix = 0  if key <  array(1)
  !>   ix = -1 if key >  array(n)
  !>   ix = -2 if no valid bracketing index was found (shouldn't happen if key in bounds)
  pure function find_bracket(n, array, key) result(ix)
    integer, intent(in) :: n
    integer, intent(in) :: array(n)
    integer, intent(in) :: key
    integer             :: ix
    integer             :: i_min, i_max, i_middle

    if (n < 2) error stop "n < 2"

    if (key < array(1)) then
       error stop "key < array(1)"
    else if (key > array(n)) then
       error stop "key > array(n)"
    end if

    ! Binary search
    i_min = 1
    i_max = n - 1

    do while (i_min <= i_max)
       i_middle = i_min + ishft(i_max - i_min, -1)  ! midpoint = (i_min + i_max) / 2

       if (key < array(i_middle)) then
          i_max = i_middle - 1
       else if (key >= array(i_middle + 1)) then
          i_min = i_middle + 1
       else
          ix = i_middle
          return
       end if
    end do

    error stop "No index found, is the array sorted?"
  end function find_bracket

  !> Fill buffers for flux fixing with average flux from fine side
  subroutine fixflux_to_buf(f4, n_vars, i_vars)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)

    integer  :: i_fine, i_bflux, n, i, i_f
    integer  :: iv, ivar, i_buf0, i_buf
#:if NDIM == 3
    integer :: j, j_f
#:endif

    ! Early exit if nothing to do
    if (f4%gc_f2c_to_buf_iface(2*NDIM) == 1) return

#:def fyp_fixflux_to_buf(face, ilim, jlim=None)
#:if NDIM == 2
    ${LOOP('collapse(NDIM+1) private(i_f, ivar, i_buf, i_fine, i_buf0, i_bflux)')}$
    do n = f4%gc_f2c_to_buf_iface(${face}$), f4%gc_f2c_to_buf_iface(${face}$+1)-1
       do iv = 1, n_vars
          do i = 1, ${ilim}$
             ivar = i_vars(iv)
             i_fine = f4%gc_f2c_to_buf_fluxfix(1, n) + 1 ! Fine block
             i_buf0 = f4%gc_f2c_to_buf_fluxfix(2, n) * n_vars
             i_bflux = f4%bflux_ix(${face}$, i_fine)
             i_f = 2 * i - 1
             i_buf = i_buf0 + ix_offset2(iv, i, ${ilim}$) + 1

             f4%send_buffer(i_buf) = 0.5_fp * ( &
                  f4%bflux(i_f, ivar, i_bflux) + &
                  f4%bflux(i_f+1, ivar, i_bflux))
          end do
       end do
    end do
#:elif NDIM == 3
    ${LOOP('collapse(NDIM+1) private(i_f, j_f, ivar, i_buf, i_fine, i_buf0, i_bflux)')}$
    do n = f4%gc_f2c_to_buf_iface(${face}$), f4%gc_f2c_to_buf_iface(${face}$+1)-1
       do iv = 1, n_vars
          do j = 1, ${jlim}$
             do i = 1, ${ilim}$
                ivar = i_vars(iv)
                i_fine = f4%gc_f2c_to_buf_fluxfix(1, n) + 1 ! Fine block
                i_buf0 = f4%gc_f2c_to_buf_fluxfix(2, n) * n_vars
                i_bflux = f4%bflux_ix(${face}$, i_fine)
                i_f = 2 * i - 1
                j_f = 2 * j - 1

                i_buf = i_buf0 + ix_offset3(iv, j, i, ${jlim}$, ${ilim}$) + 1

                f4%send_buffer(i_buf) = 0.25_fp * ( &
                     f4%bflux(i_f, j_f, ivar, i_bflux) + &
                     f4%bflux(i_f+1, j_f, ivar, i_bflux) + &
                     f4%bflux(i_f, j_f+1, ivar, i_bflux) + &
                     f4%bflux(i_f+1, j_f+1, ivar, i_bflux))
             end do
          end do
       end do
    end do
#:endif
#:enddef

    ${PARALLEL(COPYIN('i_vars'))}$ ${DEFAULT_PRESENT()}$

#:if NDIM == 2
    @:fyp_fixflux_to_buf(0, f4%hbx(2))
    @:fyp_fixflux_to_buf(1, f4%hbx(2))
    @:fyp_fixflux_to_buf(2, f4%hbx(1))
    @:fyp_fixflux_to_buf(3, f4%hbx(1))
#:elif NDIM == 3
    @:fyp_fixflux_to_buf(0, f4%hbx(2), f4%hbx(3))
    @:fyp_fixflux_to_buf(1, f4%hbx(2), f4%hbx(3))
    @:fyp_fixflux_to_buf(2, f4%hbx(1), f4%hbx(3))
    @:fyp_fixflux_to_buf(3, f4%hbx(1), f4%hbx(3))
    @:fyp_fixflux_to_buf(4, f4%hbx(1), f4%hbx(2))
    @:fyp_fixflux_to_buf(5, f4%hbx(1), f4%hbx(2))
#:endif

    ${END_PARALLEL()}$

  end subroutine fixflux_to_buf

  !> Correct solution on coarse side of refinement boundaries, taking into
  !> account the difference in the fine and coarse flux at this boundary
  subroutine fixflux_correct(f4, n_vars, i_vars, s_out)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    !> Block offset for temporal state to modify
    integer, intent(in)          :: s_out

    integer  :: i_coarse, i_fine, i_bflux, i_bflux_fine, i_bflux_coarse
    integer  :: n, i, i_f, i_c
    integer  :: iv, ivar, i_buf0, i_buf
    integer  :: offset(NDIM-1)
    real(fp) :: flux_diff, fac
#:if NDIM == 3
    integer :: j, j_c, j_f
#:endif

    ! Early exit if nothing to do
    if (f4%gc_c2f_from_buf_iface(2*NDIM) == 1 .and. &
         f4%gc_f2c_local_iface(2*NDIM) == 1) return

#:if NDIM == 2
#:def fyp_fixflux_from_buf(face, ilim, ix, sign)
    ${LOOP('collapse(3) private(i_c, ivar, i_buf, flux_diff, i_coarse, offset, i_buf0, fac, i_bflux)')}$
    do n = f4%gc_c2f_from_buf_iface(${face}$), f4%gc_c2f_from_buf_iface(${face}$+1)-1
       do iv = 1, n_vars
          do i = 1, ${ilim}$
             ivar = i_vars(iv)
             i_coarse = f4%gc_c2f_from_buf_fluxfix(1, n) + 1 + s_out ! Coarse block
             offset(1)   = f4%gc_c2f_from_buf_fluxfix(2, n)
             i_buf0   = f4%gc_c2f_from_buf_fluxfix(3, n) * n_vars
             i_bflux = f4%bflux_ix(${face}$, i_coarse-s_out)
             fac      = real(${sign}$ / &
                  f4%dr_level(f4_face_dim(${face}$), f4%block_level(i_coarse-s_out)), fp)
             i_c = i + offset(1) * ${ilim}$

             i_buf = i_buf0 + ix_offset2(iv, i, ${ilim}$) + 1
             flux_diff = fac * ( &
                  f4%bflux(i_c, ivar, i_bflux) - &
                  f4%recv_buffer(i_buf))

             ! Correct solution on coarse side. Prevent a race condition with
             ! the atomic statement.
             ${ATOMIC()}$
             f4%uu(${ix}$, ivar, i_coarse) = &
                  f4%uu(${ix}$, ivar, i_coarse) + flux_diff
          end do
       end do
    end do
#:enddef

#:def fyp_fixflux_local(face, oface, ilim, ix, sign)
    ${LOOP('collapse(3) private(i_c, i_f, ivar, flux_diff, i_coarse, i_fine, offset, fac, i_bflux_fine, i_bflux_coarse)')}$
    do n = f4%gc_f2c_local_iface(${face}$), f4%gc_f2c_local_iface(${face}$+1)-1
       do iv = 1, n_vars
          do i = 1, ${ilim}$
             ivar = i_vars(iv)
             i_fine   = f4%gc_f2c_local(1, n) + 1 + s_out ! Fine block
             i_coarse = f4%gc_f2c_local(2, n) + 1 + s_out ! coarse block
             offset(1) = f4%gc_f2c_local(3, n)  ! offset
             i_bflux_fine = f4%bflux_ix(${face}$, i_fine-s_out)
             i_bflux_coarse = f4%bflux_ix(${oface}$, i_coarse-s_out)
             fac = real(${sign}$ / &
                  f4%dr_level(f4_face_dim(${face}$), f4%block_level(i_coarse-s_out)), fp)
             i_f = 2 * i - 1
             i_c = i + offset(1) * ${ilim}$

             ! Difference between coarse and fine side
             flux_diff = fac * ( &
                  f4%bflux(i_c, ivar, i_bflux_coarse) - 0.5_fp * ( &
                  f4%bflux(i_f, ivar, i_bflux_fine) + &
                  f4%bflux(i_f+1, ivar, i_bflux_fine)))

             ! Correct solution on coarse side. Prevent a race condition with
             ! the atomic statement.
             ${ATOMIC()}$
             f4%uu(${ix}$, ivar, i_coarse) = &
                  f4%uu(${ix}$, ivar, i_coarse) + flux_diff
          end do
       end do
    end do
#:enddef

#:elif NDIM == 3
#:def fyp_fixflux_from_buf(face, ilim, jlim, ix, sign)
    ${LOOP('collapse(4) private(i_c, j_c, ivar, i_buf, flux_diff, i_coarse, offset, i_buf0, fac, i_bflux)')}$
    do n = f4%gc_c2f_from_buf_iface(${face}$), f4%gc_c2f_from_buf_iface(${face}$+1)-1
       do iv = 1, n_vars
          do j = 1, ${jlim}$
             do i = 1, ${ilim}$
                ivar = i_vars(iv)
                i_coarse = f4%gc_c2f_from_buf_fluxfix(1, n) + 1 + s_out ! Coarse block
                offset(1:2) = f4%gc_c2f_from_buf_fluxfix(2:3, n)
                i_buf0   = f4%gc_c2f_from_buf_fluxfix(4, n) * n_vars
                i_bflux = f4%bflux_ix(${face}$, i_coarse-s_out)
                fac      = real(${sign}$ / &
                     f4%dr_level(f4_face_dim(${face}$), f4%block_level(i_coarse-s_out)), fp)
                i_c = i + offset(1) * ${ilim}$
                j_c = j + offset(2) * ${jlim}$

                i_buf = i_buf0 + ix_offset3(iv, j, i, ${jlim}$, ${ilim}$) + 1
                flux_diff = fac * ( &
                     f4%bflux(i_c, j_c, ivar, i_bflux) - &
                     f4%recv_buffer(i_buf))

                ! Correct solution on coarse side. Prevent a race condition with
                ! the atomic statement.
                ${ATOMIC()}$
                f4%uu(${ix}$, ivar, i_coarse) = &
                     f4%uu(${ix}$, ivar, i_coarse) + flux_diff
             end do
          end do
       end do
    end do
#:enddef

#:def fyp_fixflux_local(face, oface, ilim, jlim, ix, sign)
    ${LOOP('collapse(4) private(i_c, i_f, j_c, j_f, ivar, flux_diff, i_coarse, i_fine, offset, fac, i_bflux_fine, i_bflux_coarse)')}$
    do n = f4%gc_f2c_local_iface(${face}$), f4%gc_f2c_local_iface(${face}$+1)-1
       do iv = 1, n_vars
          do j = 1, ${jlim}$
             do i = 1, ${ilim}$
                ivar = i_vars(iv)
                i_fine   = f4%gc_f2c_local(1, n) + 1 + s_out ! Fine block
                i_coarse = f4%gc_f2c_local(2, n) + 1 + s_out ! coarse block
                offset(1:2) = f4%gc_f2c_local(3:4, n)   ! offset
                i_bflux_fine = f4%bflux_ix(${face}$, i_fine-s_out)
                i_bflux_coarse = f4%bflux_ix(${oface}$, i_coarse-s_out)
                fac      = real(${sign}$ / &
                     f4%dr_level(f4_face_dim(${face}$), f4%block_level(i_coarse-s_out)), fp)

                i_f = 2 * i - 1
                j_f = 2 * j - 1
                i_c = i + offset(1) * ${ilim}$
                j_c = j + offset(2) * ${jlim}$

                ! Difference between coarse and fine side
                flux_diff = fac * ( &
                     f4%bflux(i_c, j_c, ivar, i_bflux_coarse) - 0.25_fp * ( &
                     f4%bflux(i_f, j_f, ivar, i_bflux_fine) + &
                     f4%bflux(i_f+1, j_f, ivar, i_bflux_fine) + &
                     f4%bflux(i_f, j_f+1, ivar, i_bflux_fine) + &
                     f4%bflux(i_f+1, j_f+1, ivar, i_bflux_fine)))

                ! Correct solution on coarse side. Prevent a race condition with
                ! the atomic statement.
                ${ATOMIC()}$
                f4%uu(${ix}$, ivar, i_coarse) = &
                     f4%uu(${ix}$, ivar, i_coarse) + flux_diff
             end do
          end do
       end do
    end do
#:enddef
#:endif

    ! Correct solution on coarse side of non-local refinement boundaries

    ${PARALLEL(COPYIN('i_vars'))}$ ${DEFAULT_PRESENT()}$
#:if NDIM == 2
    @:fyp_fixflux_from_buf(0, f4%hbx(2), {1, i_c}, -1)
    @:fyp_fixflux_from_buf(1, f4%hbx(2), {f4%bx(1), i_c}, 1)
    @:fyp_fixflux_from_buf(2, f4%hbx(1), {i_c, 1}, -1)
    @:fyp_fixflux_from_buf(3, f4%hbx(1), {i_c, f4%bx(2)}, 1)
#:elif NDIM == 3
    @:fyp_fixflux_from_buf(0, f4%hbx(2), f4%hbx(3), {1, i_c, j_c}, -1)
    @:fyp_fixflux_from_buf(1, f4%hbx(2), f4%hbx(3), {f4%bx(1), i_c, j_c}, 1)
    @:fyp_fixflux_from_buf(2, f4%hbx(1), f4%hbx(3), {i_c, 1, j_c}, -1)
    @:fyp_fixflux_from_buf(3, f4%hbx(1), f4%hbx(3), {i_c, f4%bx(2), j_c}, 1)
    @:fyp_fixflux_from_buf(4, f4%hbx(2), f4%hbx(3), {i_c, j_c, 1}, -1)
    @:fyp_fixflux_from_buf(5, f4%hbx(2), f4%hbx(3), {i_c, j_c, f4%bx(3)}, 1)
#:endif

    ! Local refinement boundaries
#:if NDIM == 2
    @:fyp_fixflux_local(0, 1, f4%hbx(2), {f4%bx(1), i_c}, 1)
    @:fyp_fixflux_local(1, 0, f4%hbx(2), {1, i_c}, -1)
    @:fyp_fixflux_local(2, 3, f4%hbx(1), {i_c, f4%bx(2)}, 1)
    @:fyp_fixflux_local(3, 2, f4%hbx(1), {i_c, 1}, -1)
#:elif NDIM == 3
    @:fyp_fixflux_local(0, 1, f4%hbx(2), f4%hbx(3), {f4%bx(1), i_c, j_c}, 1)
    @:fyp_fixflux_local(1, 0, f4%hbx(2), f4%hbx(3), {1, i_c, j_c}, -1)
    @:fyp_fixflux_local(2, 3, f4%hbx(1), f4%hbx(3), {i_c, f4%bx(2), j_c}, 1)
    @:fyp_fixflux_local(3, 2, f4%hbx(1), f4%hbx(3), {i_c, 1, j_c}, -1)
    @:fyp_fixflux_local(4, 5, f4%hbx(1), f4%hbx(2), {i_c, j_c, f4%bx(3)}, 1)
    @:fyp_fixflux_local(5, 4, f4%hbx(1), f4%hbx(2), {i_c, j_c, 1}, -1)
#:endif
    ${END_PARALLEL()}$

  end subroutine fixflux_correct

  !> Correct fluxes at refinement boundaries. Use the average flux from the
  !> fine side on the coarse side.
  subroutine f4_fix_c2f_flux(f4, n_vars, i_vars, s_out)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer, intent(in)          :: s_out
    real(dp)                     :: t0, t1

    t0 = MPI_Wtime()
    call fixflux_to_buf(f4, n_vars, i_vars)

    t1 = MPI_Wtime()
    f4%wtime_flux_fix = f4%wtime_flux_fix + t1 - t0

    ! Update send/recv offsets
    f4%recv_offset(:) = f4%gc_recv_offset_fluxfix * n_vars
    f4%send_offset(:) = f4%gc_send_offset_fluxfix * n_vars
    call f4_exchange_buffers(f4)

    t0 = MPI_Wtime()
    f4%wtime_exchange_buffers = f4%wtime_exchange_buffers + t0 - t1

    call fixflux_correct(f4, n_vars, i_vars, s_out)

    t1 = MPI_Wtime()
    f4%wtime_flux_fix = f4%wtime_flux_fix + t1 - t0

  end subroutine f4_fix_c2f_flux

  !> Compute sum of a variable
  subroutine f4_compute_sum(f4, i_var, var_sum)
    type(foap4_t), intent(in) :: f4
    integer, intent(in)       :: i_var
    real(dp), intent(out)     :: var_sum
    integer                   :: level, ${IJK}$, n, ierror
    real(dp)                  :: dvol

    var_sum = 0.0_dp

    ${PARALLEL_LOOP('collapse(NDIM+1) private(level, dvol) reduction(+:var_sum)')}$ ${DEFAULT_PRESENT()}$
    do n = 1, f4%n_blocks
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          level = f4%block_level(n)
#:if NDIM == 2
          dvol = f4%dr_level(1, level) * f4%dr_level(2, level)
#:elif NDIM == 3
          dvol = f4%dr_level(1, level) * f4%dr_level(2, level) * &
               f4%dr_level(3, level)
#:endif
          var_sum = var_sum + f4%uu(${IJK}$, i_var, n) * dvol
       end do; ${KJI_CLOSE_LOOP}$
    end do

    call MPI_Allreduce(MPI_IN_PLACE, var_sum, 1, MPI_DOUBLE_PRECISION, &
         MPI_SUM, f4%mpicomm, ierror)

  end subroutine f4_compute_sum

  !> Compute max of a variable
  subroutine f4_compute_max(f4, i_var, var_max)
    type(foap4_t), intent(in) :: f4
    integer, intent(in)       :: i_var
    real(dp), intent(out)     :: var_max
    integer                   :: ${IJK}$, n, ierror

    var_max = -huge(1.0_dp)

    ${PARALLEL_LOOP('collapse(NDIM+1) reduction(max:var_max)')}$ ${DEFAULT_PRESENT()}$
    do n = 1, f4%n_blocks
       do @{KJI_LOOP_1_to_array(f4%bx)}@
            var_max = max(var_max, f4%uu(${IJK}$, i_var, n))
       end do; ${KJI_CLOSE_LOOP}$
    end do

    call MPI_Allreduce(MPI_IN_PLACE, var_max, 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, f4%mpicomm, ierror)

  end subroutine f4_compute_max

  !> Compute index offset for indexing in 4D array shaped (*, n2, n3, n4)
  pure integer function ix_offset4(i1, i2, i3, i4, n2, n3, n4)
    @{ROUTINE_SEQ()}@
    integer, intent(in) :: i1, i2, i3, i4, n2, n3, n4
    ix_offset4 = &
         (i1 - 1) * n2 * n3 * n4 + &
         (i2 - 1) * n3 * n4 + &
         (i3 - 1) * n4 + &
         (i4 - 1)
  end function ix_offset4

  !> Compute index offset for indexing in 3D array shaped (*, n2, n3)
  pure integer function ix_offset3(i1, i2, i3, n2, n3)
    @{ROUTINE_SEQ()}@
    integer, intent(in) :: i1, i2, i3, n2, n3
    ix_offset3 = (i1 - 1) * n2 * n3 + (i2 - 1) * n3 + (i3 - 1)
  end function ix_offset3

  !> Compute index offset for indexing in 2D array shaped (*, n2)
  pure integer function ix_offset2(i1, i2, n2)
    @{ROUTINE_SEQ()}@
    integer, intent(in) :: i1, i2, n2
    ix_offset2 = (i1 - 1) * n2 + (i2 - 1)
  end function ix_offset2

end module m_foap4_${NDIM}$d
