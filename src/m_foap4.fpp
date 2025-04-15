!> Foap4 stands for "Fortran OpenACC p4est", and implements MPI-parallel
!> quadtree/octree adaptive mesh refinement with support for OpenACC. Together
!> with p4est_wrapper.c, this module implements the required data structures and
!> methods.
!>
!> Author(s): Jannis Teunissen
!> 
module m_foap4
  use, intrinsic :: iso_c_binding
  use mpi_f08

  implicit none
  private

  integer, parameter, private :: dp = kind(0.0d0)

  !> Maximum refinement level in p4est
  integer, parameter :: P4EST_MAXLEVEL = 30

  !> The opposite of the faces 0-3
  integer, parameter :: face_swap(0:3) = [1, 0, 3, 2]

  !> The offset of the children
  integer, parameter :: child_offset(2, 4) = reshape([0,0,1,0,0,1,1,1], [2,4])
  !$acc declare create(child_offset)

  !> Value indicating a physical boundary at a block face
  integer, parameter :: FACE_BOUNDARY = 0

  !> Value indicating a neighbor at the same refinement level at a block face
  integer, parameter :: FACE_SAME_LEVEL = 1

  !> Value indicating a neighbor at a higher refinement level at a block face
  integer, parameter :: FACE_COARSE_TO_FINE = 2

  !> Value indicating a neighbor at a lower refinement level at a block face
  integer, parameter :: FACE_FINE_TO_COARSE = 3

  !> Type to store an array of integers
  type int_array_t
     integer, allocatable :: i(:)
  end type int_array_t

  !> Type to describe a face boundary. The same data structure is defined in
  !> p4est_wrapper.c
  type, bind(c) :: bnd_face_t
     integer(c_int) :: face_type !< What kind of face (same level, ...)
     integer(c_int) :: face      !< Direction of the face
     integer(c_int) :: other_proc !< MPI rank that owns quadid(2)
     integer(c_int) :: quadid(2)  !< quadid(1) is always local, (2) can be non-local
     integer(c_int) :: offset     !< Offset for a hanging face
     integer(c_int) :: ibuf_recv  !< Index in receive buffer (not filled here)
     integer(c_int) :: ibuf_send  !< Index in send buffer (not filled here)
  end type bnd_face_t

  !> Data structure that contains the whole AMR grid and required information
  !> for ghost cell communication
  type, public :: foap4_t
     integer   :: bx(2)                           !< Block size (cells)
     integer   :: n_gc                            !< Number of ghost cells
     integer   :: max_blocks                      !< Maximum number of blocks used
     integer   :: n_vars                          !< Number of variables
     real(dp)  :: tree_length(2)                  !< Length of tree
     integer   :: ilo(2)                          !< Minimum index in a block
     integer   :: ihi(2)                          !< Maximum index in a block
     real(dp)  :: dr_level(2, 0:p4est_maxlevel-1) !< Grid spacing per level
     character(len=32), allocatable :: var_names(:) !< Names of the variables

     integer :: n_blocks              !< Number of blocks used
     integer :: gc_mesh_revision = -1 !< Revision number (global) of the mesh

     !> Level of each block
     integer, allocatable  :: block_level(:)
     !> Origin of each block
     real(dp), allocatable :: block_origin(:, :)
     !> Storage of block data uu(i, j, i_var, i_block)
     real(dp), allocatable :: uu(:, :, :, :)
     !> Refinement flag of each block. Negative means coarsen (if possible),
     !> positive means refine, and zero means keep refinement.
     integer, allocatable  :: refinement_flags(:)

     ! For communication
     type(MPI_comm)        :: mpicomm        !< MPI communicator
     integer               :: mpirank        !< MPI rank of this task
     integer               :: mpisize        !< Number of ranks in communicator
     integer, allocatable  :: recv_offset(:) !< 0:mpisize offsets for receiving
     integer, allocatable  :: send_offset(:) !< 0:mpisize offsets for sending
     real(dp), allocatable :: recv_buffer(:) !< Buffer for receiving data
     real(dp), allocatable :: send_buffer(:) !< Buffer for sending data

     !> Pointer to a structure pw_state_t defined in p4est_wrapper.c, which
     !> contains all the p4est state required for the wrapper
     type(c_ptr) :: pw

     ! For handling ghost cells on faces

     !> It is required that bx(1) == bx(2). The ghost cell data size per face is
     !> thus n_gc * bx(1)
     integer              :: gc_data_size

     !> Receive offset (per MPI rank) in recv_buffer
     integer, allocatable :: gc_recv_offset(:)
     !> Receive offset (per MPI rank) in recv_buffer for the coarse-to-fine step
     integer, allocatable :: gc_recv_offset_c2f(:)
     !> Send offset (per MPI rank) in send_buffer
     integer, allocatable :: gc_send_offset(:)
     !> Send offset (per MPI rank) in send_buffer for the coarse-to-fine step
     integer, allocatable :: gc_send_offset_c2f(:)

     ! All the arrays below are sorted by face direction

     !> gc_srl_local(1:2, n) stores the quad indices of the nth local face at
     !> the same refinement level
     integer, allocatable :: gc_srl_local(:, :)
     !> Quad index and buffer index for each srl face that is received
     integer, allocatable :: gc_srl_from_buf(:, :)
     !> Quad index and buffer index for each srl face that is sent
     integer, allocatable :: gc_srl_to_buf(:, :)
     !> gc_f2c_local(1:3, n) contains [qid fine, qid coarse, fine offset]
     integer, allocatable :: gc_f2c_local(:, :)
     !> gc_c2f_from_buf(1:3, n) contains [qid coarse, fine offset, ibuf_recv]
     integer, allocatable :: gc_c2f_from_buf(:, :)
     !> Quad index, offset and send buffer index for c2f face
     integer, allocatable :: gc_c2f_to_buf(:, :)
     !> Quad index and recv buffer index per f2c face
     integer, allocatable :: gc_f2c_from_buf(:, :)
     !> Quad index and send buffer index per f2c face
     integer, allocatable :: gc_f2c_to_buf(:, :)
     !> List of quad indices with a physical boundary
     integer, allocatable :: gc_phys(:)

     ! Since the above arrays are sorted by face direction, we can loop over
     ! them per face direction. The arrays below define the start index for each
     ! direction (with the last value being the total number of elements).

     integer :: gc_srl_local_iface(0:4)    !< Start index per face direction
     integer :: gc_srl_from_buf_iface(0:4) !< Start index per face direction
     integer :: gc_srl_to_buf_iface(0:4)   !< Start index per face direction
     integer :: gc_f2c_local_iface(0:4)    !< Start index per face direction
     integer :: gc_c2f_from_buf_iface(0:4) !< Start index per face direction
     integer :: gc_c2f_to_buf_iface(0:4)   !< Start index per face direction
     integer :: gc_f2c_from_buf_iface(0:4) !< Start index per face direction
     integer :: gc_f2c_to_buf_iface(0:4)   !< Start index per face direction
     integer :: gc_phys_iface(0:4)         !< Start index per face direction

     ! Performance information
     real(dp) :: wtime_t0 = 0.0_dp
     real(dp) :: wtime_gc_fill_round1 = 0.0_dp
     real(dp) :: wtime_gc_fill_round2 = 0.0_dp
     real(dp) :: wtime_gc_fill_buff_round1 = 0.0_dp
     real(dp) :: wtime_gc_fill_buff_round2 = 0.0_dp
     real(dp) :: wtime_adjust_ref_p4est = 0.0_dp
     real(dp) :: wtime_adjust_ref_foap4 = 0.0_dp
     real(dp) :: wtime_partition = 0.0_dp
     real(dp) :: wtime_write_grid = 0.0_dp
     real(dp) :: wtime_update_gc_pattern = 0.0_dp
     real(dp) :: wtime_exchange_buffers = 0.0_dp

  end type foap4_t

  interface
     !> Initialize p4est and MPI
     subroutine pw_initialize(pw, mpicomm, log_level) bind(c)
       import c_int, c_ptr
       type(c_ptr), intent(out)          :: pw
       integer(c_int), intent(out)       :: mpicomm
       integer(c_int), intent(in), value :: log_level
     end subroutine pw_initialize

     !> Destroy p4est data
     subroutine pw_destroy(pw) bind(c)
       import c_ptr
       type(c_ptr), intent(in), value :: pw
     end subroutine pw_destroy

     !> Finalize p4est and MPI
     subroutine pw_finalize(pw) bind(c)
       import c_ptr
       type(c_ptr), intent(in), value :: pw
     end subroutine pw_finalize

     !> Create p4est brick
     subroutine pw_set_connectivity_brick(pw, mi, ni, periodic_a, periodic_b, &
          min_level, fill_uniform, max_blocks) bind(c)
       import c_int, c_ptr
       type(c_ptr), intent(in), value    :: pw
       integer(c_int), value, intent(in) :: mi
       integer(c_int), value, intent(in) :: ni
       integer(c_int), value, intent(in) :: periodic_a
       integer(c_int), value, intent(in) :: periodic_b
       integer(c_int), value, intent(in) :: min_level
       integer(c_int), value, intent(in) :: fill_uniform
       integer(c_int), intent(in), value :: max_blocks
     end subroutine pw_set_connectivity_brick

     !> Get number of quadrants of this MPI rank
     pure function pw_get_num_local_quadrants(pw) result(n) bind(c)
       import c_int, c_ptr
       type(c_ptr), intent(in), value :: pw
       integer(c_int)                 :: n
     end function pw_get_num_local_quadrants

     !> Get number of quadrants of all MPI ranks together
     pure function pw_get_num_global_quadrants(pw) result(n) bind(c)
       import c_int, c_ptr
       type(c_ptr), intent(in), value :: pw
       integer(c_int)                 :: n
     end function pw_get_num_global_quadrants

     !> Get the global revision number of the mesh
     pure function pw_get_mesh_revision(pw) result(n) bind(c)
       import c_int, c_ptr
       type(c_ptr), intent(in), value :: pw
       integer(c_int)                 :: n
     end function pw_get_mesh_revision

     !> Get the highest refinement level of this MPI rank
     pure function pw_get_highest_local_level(pw) result(n) bind(c)
       import c_int, c_ptr
       type(c_ptr), intent(in), value :: pw
       integer(c_int)                 :: n
     end function pw_get_highest_local_level

     !> Get the coordinates and level of the quadrants
     subroutine pw_get_quadrants(pw, n_quadrants, coord, level) bind(c)
       import c_int, c_ptr, c_double
       type(c_ptr), intent(in), value     :: pw
       integer(c_int), value, intent(in)  :: n_quadrants
       real(kind=c_double), intent(inout) :: coord(*)
       integer(c_int), intent(inout)      :: level(n_quadrants)
     end subroutine pw_get_quadrants

     !> Write p4est vtk file with the MPI rank
     subroutine pw_vtk_write_file(pw, fname) bind(c)
       import c_char, c_ptr
       type(c_ptr), intent(in), value     :: pw
       character(kind=C_char), intent(in) :: fname(*)
     end subroutine pw_vtk_write_file

     !> Get all the boundaries at quadrant faces
     subroutine pw_get_all_faces(pw, n_faces, bnd_face_ptr) bind(c)
       import c_int, c_ptr
       type(c_ptr), intent(in), value :: pw
       integer(c_int), intent(out)    :: n_faces
       type(c_ptr), intent(out)       :: bnd_face_ptr
     end subroutine pw_get_all_faces

     !> Adjust the refinement according to the refinement flags
     subroutine pw_adjust_refinement(pw, n_quadrants, flags, has_changed) bind(c)
       import c_int, c_ptr
       type(c_ptr), intent(in), value    :: pw
       integer(c_int), value, intent(in) :: n_quadrants
       integer(c_int), intent(in)        :: flags(n_quadrants)
       integer(c_int), intent(out)       :: has_changed
     end subroutine pw_adjust_refinement

     !> Partition the blocks over the MPI ranks
     subroutine pw_partition(pw, n_changed_global, gfq_old, gfq_new) bind(c)
       import c_int, c_int64_t, c_ptr
       type(c_ptr), intent(in), value    :: pw
       integer(c_int), intent(out)       :: n_changed_global
       integer(c_int64_t), intent(out)   :: gfq_old(*)
       integer(c_int64_t), intent(out)   :: gfq_new(*)
     end subroutine pw_partition

     !> Transfer block data between MPI ranks
     subroutine pw_partition_transfer(pw, gfq_old, src_data, dest_data, data_size) bind(c)
       import c_int, c_int64_t, c_ptr
       type(c_ptr), intent(in), value    :: pw
       integer(c_int64_t), intent(in)    :: gfq_old(*)
       type(c_ptr), value                :: src_data
       type(c_ptr), value                :: dest_data
       integer(c_int), intent(in), value :: data_size
     end subroutine pw_partition_transfer
  end interface

  public :: f4_initialize
  public :: f4_destroy
  public :: f4_finalize
  public :: f4_construct_brick
  public :: f4_reset_wtime
  public :: f4_print_wtime
  public :: f4_write_grid
  public :: f4_get_mesh_revision
  public :: f4_get_global_highest_level
  public :: f4_get_num_local_blocks
  public :: f4_get_num_global_blocks
  public :: f4_cell_coord
  public :: f4_exchange_buffers
  public :: f4_update_ghostcells
  public :: f4_adjust_refinement
  public :: f4_partition

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

#ifdef _OPENACC
    call set_openacc_device(f4)
#endif

  end subroutine f4_initialize

#ifdef _OPENACC
  !> Set the device to be used by OpenACC. Code based on:
  !> https://docs.nvidia.com/hpc-sdk/compilers/openacc-mpi-tutorial/
  !> OpenMPI provides an environment variable, but this is not portable
  subroutine set_openacc_device(f4)
    use openacc
    type(foap4_t), intent(in)            :: f4

    interface
       ! Get a unique number to identify the host
       function gethostid() bind(C)
         import C_int
         integer (C_int) :: gethostid
       end function gethostid
    end interface

    integer :: hostids(0:f4%mpisize-1), local_procs(0:f4%mpisize-1)
    integer :: hostid, ierr, num_devices, my_device, rank, num_local_procs
    integer(acc_device_kind) :: dev_type

    dev_type = ACC_DEVICE_DEFAULT

    ! Get the hostids to determine how many processes are on this host
    hostid = gethostid()
    call MPI_Allgather(hostid, 1, MPI_INTEGER, hostids, 1, MPI_INTEGER, &
         f4%mpicomm, ierr)

    ! Determine the local MPI ranks and number them, starting at zero
    num_local_procs = 0
    local_procs     = 0

    do rank = 0, f4%mpisize-1
       if (hostid == hostids(rank)) then
          local_procs(rank) = num_local_procs
          num_local_procs = num_local_procs+1
       endif
    enddo

    num_devices = acc_get_num_devices(dev_type)
    if (num_devices < 1) error stop "No devices available on host"

    if (num_devices < num_local_procs) then
       ! Print warning only for first local process
       if (local_procs(f4%mpirank) == 0) then
          write(*, "(A,I0,A,I0,A,I0,A)") "WARNING from ", f4%mpirank, &
               ": more local processes (", num_local_procs, &
               ") than GPUs (", num_devices, ")"
       endif

       my_device = mod(local_procs(f4%mpirank), num_devices)
    else
       my_device = local_procs(f4%mpirank)
    endif

    call acc_set_device_num(my_device, dev_type)

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
  end subroutine f4_reset_wtime

  !> Print wall clock time measurements
  subroutine f4_print_wtime(f4)
    type(foap4_t), intent(inout) :: f4
    real(dp)                     :: t_total, fac

    t_total = MPI_Wtime() - f4%wtime_t0
    fac     = 1e2_dp / t_total
    write(*, "(I6,A25,F9.2,' s')") f4%mpirank, "total_time", t_total
    write(*, "(I6,A25,F9.2,' %')") f4%mpirank, "gc_fill_round1", &
         f4%wtime_gc_fill_round1 * fac
    write(*, "(I6,A25,F9.2,' %')") f4%mpirank, "gc_fill_round2", &
         f4%wtime_gc_fill_round2 * fac
    write(*, "(I6,A25,F9.2,' %')") f4%mpirank, "gc_fill_buff_round1", &
         f4%wtime_gc_fill_buff_round1 * fac
    write(*, "(I6,A25,F9.2,' %')") f4%mpirank, "gc_fill_buff_round2", &
         f4%wtime_gc_fill_buff_round2 * fac
    write(*, "(I6,A25,F9.2,' %')") f4%mpirank, "adjust_ref_p4est", &
         f4%wtime_adjust_ref_p4est * fac
    write(*, "(I6,A25,F9.2,' %')") f4%mpirank, "adjust_ref_foap4", &
         f4%wtime_adjust_ref_foap4 * fac
    write(*, "(I6,A25,F9.2,' %')") f4%mpirank, "partition", &
         f4%wtime_partition * fac
    write(*, "(I6,A25,F9.2,' %')") f4%mpirank, "write_grid", &
         f4%wtime_write_grid * fac
    write(*, "(I6,A25,F9.2,' %')") f4%mpirank, "update_gc_pattern", &
         f4%wtime_update_gc_pattern * fac
    write(*, "(I6,A25,F9.2,' %')") f4%mpirank, "exchange_buffers", &
         f4%wtime_exchange_buffers * fac
  end subroutine f4_print_wtime

  !> Destroy all data for the current mesh
  subroutine f4_destroy(f4)
    type(foap4_t), intent(inout) :: f4

    ! OpenACC - Remove data from device

    !$acc exit data delete(f4%block_level, f4%block_origin)
    !$acc exit data delete(f4%uu, f4%refinement_flags)
    !$acc exit data delete(f4%recv_buffer, f4%send_buffer)
    !$acc exit data delete(&
    !$acc &f4%gc_srl_local_iface, f4%gc_srl_from_buf_iface, f4%gc_srl_to_buf_iface, &
    !$acc &f4%gc_f2c_local_iface, f4%gc_f2c_from_buf_iface, f4%gc_f2c_to_buf_iface, &
    !$acc &f4%gc_c2f_from_buf_iface, f4%gc_c2f_to_buf_iface, f4%gc_phys_iface)
    !$acc exit data delete(f4)

    call pw_destroy(f4%pw)

    deallocate(f4%var_names)
    deallocate(f4%block_origin)
    deallocate(f4%block_level)
    deallocate(f4%refinement_flags)
    deallocate(f4%uu)
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
       n_vars, var_names, periodic, min_level, max_blocks)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: trees_per_dim(2) !< How many trees per dimension
    real(dp), intent(in)         :: tree_length(2)   !< Length of each tree
    integer, intent(in)          :: bx(2)            !< Number of cells in a block
    integer, intent(in)          :: n_gc !< Number of ghost cells
    integer, intent(in)          :: n_vars !< Number of variables
    character(len=*), intent(in) :: var_names(n_vars) !< Variable names
    logical, intent(in)          :: periodic(2) !< Periodic flag per dim
    integer, intent(in)          :: min_level   !< Refine up to this level
    integer, intent(in)          :: max_blocks  !< Maximum number of blocks
    integer                      :: i, periodic_as_int(2)

    if (bx(1) /= bx(2)) error stop "TODO: unequal bx(:) not yet supported"
    if (any(bx < 2 * n_gc)) error stop "Cannot have any(bx < 2 * n_gc)"
    if (any(iand(bx, 1) == 1)) error stop "All bx have to be even"

    f4%bx          = bx
    f4%n_gc        = n_gc
    f4%ilo         = 1 - n_gc
    f4%ihi         = bx + n_gc
    f4%n_vars      = n_vars
    f4%max_blocks  = max_blocks
    f4%tree_length = tree_length

    where (periodic)
       periodic_as_int = 1
    elsewhere
       periodic_as_int = 0
    end where

    allocate(f4%var_names(n_vars))
    do i = 1, n_vars
       f4%var_names(i) = var_names(i)
    end do

    do i = 0, P4EST_MAXLEVEL-1
       f4%dr_level(:, i) = (tree_length/bx) * 0.5**i
    end do

    call pw_set_connectivity_brick(f4%pw, &
         trees_per_dim(1), trees_per_dim(2), &
         periodic_as_int(1), periodic_as_int(2), min_level, 1, max_blocks)

    allocate(f4%block_origin(2, max_blocks))
    allocate(f4%block_level(max_blocks))
    allocate(f4%refinement_flags(max_blocks))
    allocate(f4%uu(1-n_gc:bx(1)+n_gc, 1-n_gc:bx(2)+n_gc, n_vars, max_blocks))
    f4%uu = 0.0_dp

    f4%gc_data_size = f4%bx(1) * f4%n_gc

    ! Maximum size of recv/send buffer
    i = max_blocks * 4 * f4%gc_data_size
    allocate(f4%recv_buffer(i))
    allocate(f4%send_buffer(i))
    allocate(f4%recv_offset(0:f4%mpisize))
    allocate(f4%send_offset(0:f4%mpisize))

    ! OpenACC - Copy data structure and create allocatable components
    !$acc enter data copyin(f4)
    !$acc enter data create(f4%block_level, f4%block_origin)
    !$acc enter data create(f4%uu, f4%refinement_flags)
    !$acc enter data create(f4%recv_buffer, f4%send_buffer)
    !$acc enter data create(&
    !$acc &f4%gc_srl_local_iface, f4%gc_srl_from_buf_iface, f4%gc_srl_to_buf_iface, &
    !$acc &f4%gc_f2c_local_iface, f4%gc_f2c_from_buf_iface, f4%gc_f2c_to_buf_iface, &
    !$acc &f4%gc_c2f_from_buf_iface, f4%gc_c2f_to_buf_iface, f4%gc_phys_iface)

    call f4_set_quadrants(f4)

  end subroutine f4_construct_brick

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
    use iso_fortran_env, only: error_unit
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
    !$acc update device(f4%n_blocks, f4%block_origin(:, 1:f4%n_blocks), &
    !$acc &f4%block_level(1:f4%n_blocks))
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

  !> Write the AMR grid to a file
  subroutine f4_write_grid(f4, fname, n_output, time, viewer)
    use m_xdmf_writer
    type(foap4_t), intent(inout)           :: f4
    character(len=*), intent(in)           :: fname    !< Base file name
    integer, intent(in)                    :: n_output !< Output index
    real(dp), intent(in), optional         :: time     !< Simulation time
    character(len=*), intent(in), optional :: viewer   !< Optimize for viewer
    character(len=len_trim(fname)+7)       :: full_fname
    integer                                :: n
    real(dp), allocatable                  :: dr(:, :)
    real(dp)                               :: t0, t1

    t0 = MPI_Wtime()
    ! OpenACC - get the block data from the device
    !$acc update self(f4%uu(:, :, :, 1:f4%n_blocks))

    write(full_fname, "(A,A,I06.6)") trim(fname), "_", n_output

    call pw_vtk_write_file(f4%pw, trim(full_fname) // C_null_char)

    allocate(dr(2, f4%n_blocks))

    do n = 1, f4%n_blocks
       dr(:, n) = f4%dr_level(:, f4%block_level(n))
    end do

    call xdmf_write_blocks_2DCoRect(f4%mpicomm, trim(full_fname), &
         f4%n_blocks, f4%bx+2*f4%n_gc, f4%n_vars, &
         f4%var_names, f4%n_gc, f4%block_origin(:, 1:f4%n_blocks), dr, &
         cc_data=f4%uu(:, :, :, 1:f4%n_blocks), time=time, viewer=viewer)
    t1 = MPI_Wtime()
    f4%wtime_write_grid = f4%wtime_write_grid + t1 - t0
  end subroutine f4_write_grid

  !> Return the coordinates at the center of a grid cells
  pure function f4_cell_coord(f4, i_block, i, j) result(rr)
    !$acc routine seq
    type(foap4_t), intent(in) :: f4
    integer, intent(in)       :: i_block, i, j
    real(dp)                  :: rr(2), dr(2)

    dr = f4%dr_level(:, f4%block_level(i_block))
    rr = f4%block_origin(:, i_block) + dr * [i-0.5_dp, j-0.5_dp]
  end function f4_cell_coord

  !> Update the information required to update ghost cells
  subroutine update_ghostcell_pattern(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: mesh_revision
    type(bnd_face_t), pointer    :: bnd_face(:)
    type(c_ptr)                  :: tmp
    integer                      :: n_faces

    mesh_revision = pw_get_mesh_revision(f4%pw)
    if (mesh_revision == f4%gc_mesh_revision) return

    call pw_get_all_faces(f4%pw, n_faces, tmp)

    call c_f_pointer(tmp, bnd_face, shape=[n_faces])

    call set_ghost_cell_pattern(f4, size(bnd_face), bnd_face, &
         f4%mpirank, f4%mpisize)
    f4%gc_mesh_revision = mesh_revision

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
    end if

    f4%gc_recv_offset(0) = 0
    f4%gc_send_offset(0) = 0
    f4%gc_recv_offset_c2f(0) = 0
    f4%gc_send_offset_c2f(0) = 0

    ! Determine offsets for sending and receiving with other ranks
    do rank = 1, mpisize
       if (rank-1 == mpirank) then
          f4%gc_recv_offset(rank) = f4%gc_recv_offset(rank-1)
          f4%gc_send_offset(rank) = f4%gc_send_offset(rank-1)

          f4%gc_recv_offset_c2f(rank) = f4%gc_recv_offset_c2f(rank-1)
          f4%gc_send_offset_c2f(rank) = f4%gc_send_offset_c2f(rank-1)
       else
          f4%gc_recv_offset(rank) = f4%gc_recv_offset(rank-1) + &
               f4%gc_data_size * i_same(rank-1) + &
               f4%gc_data_size/2 * i_c2f(rank-1)
          f4%gc_send_offset(rank) = f4%gc_send_offset(rank-1) + &
               f4%gc_data_size * i_same(rank-1) + &
               f4%gc_data_size/2 * i_f2c(rank-1)

          ! In a second round of communication, handle the fine side of
          ! refinement boundaries
          f4%gc_recv_offset_c2f(rank) = f4%gc_recv_offset_c2f(rank-1) + &
               f4%gc_data_size * i_f2c(rank-1)
          f4%gc_send_offset_c2f(rank) = f4%gc_send_offset_c2f(rank-1) + &
               f4%gc_data_size * i_c2f(rank-1)
       end if
    end do

    if (allocated(f4%gc_srl_local)) then
       ! OpenACC - deallocate arrays
       !$acc exit data delete(&
       !$acc &f4%gc_srl_local, f4%gc_srl_from_buf, f4%gc_srl_to_buf, &
       !$acc &f4%gc_f2c_local, f4%gc_f2c_from_buf, f4%gc_f2c_to_buf, &
       !$acc &f4%gc_c2f_from_buf, f4%gc_c2f_to_buf, f4%gc_phys)

       deallocate(f4%gc_srl_local, f4%gc_srl_from_buf, f4%gc_srl_to_buf, &
            f4%gc_f2c_local, f4%gc_f2c_from_buf, f4%gc_f2c_to_buf, &
            f4%gc_c2f_from_buf, f4%gc_c2f_to_buf, f4%gc_phys)
    end if

    ! Local ghost cell exchange at the same level
    allocate(f4%gc_srl_local(2, i_same(mpirank)))

    ! Local ghost cell exchange at refinement boundaries
    allocate(f4%gc_f2c_local(3, i_f2c(mpirank)))

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
       f4%gc_f2c_local(:, n) = [bnd_face(i)%quadid(1), &
            bnd_face(i)%quadid(2), bnd_face(i)%offset]
    end do

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
    allocate(all_f2c_from_buf%i(n))
    allocate(all_f2c_to_buf%i(n))

    ! Non-local ghost cell exchange from coarse to fine
    n = sum(i_c2f) - i_c2f(mpirank)
    allocate(f4%gc_c2f_from_buf(3, n))
    allocate(f4%gc_c2f_to_buf(3, n))
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
          i_buf_send(rank) = i_buf_send(rank) + f4%gc_data_size / 2
       end do

       ! Receiving coarse from fine
       call sort_for_recv_or_send(c2f_ix(rank), n_faces, bnd_face, recv)

       do n = 1, size(c2f_ix(rank)%i)
          i = c2f_ix(rank)%i(n) ! Index in bnd_face array
          i_c2f_from_buf = i_c2f_from_buf + 1
          all_c2f_from_buf%i(i_c2f_from_buf) = i

          bnd_face(i)%ibuf_recv = i_buf_recv(rank)
          i_buf_recv(rank) = i_buf_recv(rank) + f4%gc_data_size / 2
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

    ! OpenACC - copy/sync data to device

    !$acc enter data copyin(&
    !$acc &f4%gc_srl_local, f4%gc_srl_from_buf, f4%gc_srl_to_buf, &
    !$acc &f4%gc_f2c_local, f4%gc_f2c_from_buf, f4%gc_f2c_to_buf, &
    !$acc &f4%gc_c2f_from_buf, f4%gc_c2f_to_buf, f4%gc_phys)

    !$acc update device(&
    !$acc &f4%gc_srl_local_iface, f4%gc_srl_from_buf_iface, f4%gc_srl_to_buf_iface, &
    !$acc &f4%gc_f2c_local_iface, f4%gc_f2c_from_buf_iface, f4%gc_f2c_to_buf_iface, &
    !$acc &f4%gc_c2f_from_buf_iface, f4%gc_c2f_to_buf_iface, f4%gc_phys_iface)

  end subroutine set_ghost_cell_pattern

  !> Sort index array by face direction
  subroutine sort_by_face(ix, n_bnd_face, bnd_face, iface)
    type(int_array_t), intent(inout) :: ix !< Index array
    integer, intent(in)              :: n_bnd_face
    type(bnd_face_t), intent(in)     :: bnd_face(n_bnd_face)
    !> Start index for each face direction
    integer, intent(out)             :: iface(0:4)
    integer                          :: face_count(0:3)
    integer                          :: face_offset(0:3)
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
    do n = 1, 3
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
    integer, intent(in)  :: face_count(0:3)
    integer, intent(out) :: iface(0:4)
    integer              :: n

    iface(0) = 1
    do n = 1, 4
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

    include 'qsort_inline.f90'

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
         else
            ! For hanging faces, sort by offset
            less_than = (bnd_face(a)%offset < bnd_face(b)%offset)
         end if
      else
         ! Order by face and quadid of 'other' side
         if (bnd_face(a)%face /= bnd_face(b)%face) then
            ! Note that the face is swapped for the receiving side
            less_than = (face_swap(bnd_face(a)%face) < &
                 face_swap(bnd_face(b)%face))
         else if (bnd_face(a)%quadid(2) /= bnd_face(b)%quadid(2)) then
            less_than = (bnd_face(a)%quadid(2) < bnd_face(b)%quadid(2))
         else
            less_than = (bnd_face(a)%offset < bnd_face(b)%offset)
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
  subroutine fill_ghostcell_buffers_round_one(f4, n_vars, i_vars)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer                      :: i, j, n, ivar, iv
    integer                      :: i_f, j_f, half_bx(2)
    integer                      :: iq, i_buf, i_buf0, face

    if (maxval(f4%gc_send_offset) * n_vars > size(f4%send_buffer)) &
         error stop "send buffer too small"

    half_bx = f4%bx/2

    associate (bx => f4%bx, n_gc => f4%n_gc, uu => f4%uu)
      !$acc parallel

#:def fyp_srl_to_buf(face, jlim, ilim, i0=0, j0=0)
      !$acc loop gang private(iq, i_buf0)
      do n = f4%gc_srl_to_buf_iface(${face}$), f4%gc_srl_to_buf_iface(${face}$+1)-1
         iq = f4%gc_srl_to_buf(1, n) + 1
         i_buf0 = f4%gc_srl_to_buf(2, n) * n_vars
         !$acc loop collapse(3) private(ivar, i_buf)
         do iv = 1, n_vars
            do j = 1, ${jlim}$
               do i = 1, ${ilim}$
                  ivar = i_vars(iv)
                  i_buf = i_buf0 + (iv - 1) * ${jlim}$ * ${ilim}$ + (j - 1) * ${ilim}$ + i
                  f4%send_buffer(i_buf) = uu(${i0}$+i, ${j0}$+j, ivar, iq)
               end do
            end do
         end do
      end do
#:enddef

      @:fyp_srl_to_buf(0, bx(2), n_gc)
      @:fyp_srl_to_buf(1, bx(2), n_gc, i0=bx(1)-n_gc)
      @:fyp_srl_to_buf(2, n_gc, bx(1))
      @:fyp_srl_to_buf(3, n_gc, bx(1), j0=bx(2)-n_gc)

      ! Nonlocal fine-to-coarse boundaries, fill buffer for coarse side

#:def fyp_f2c_to_buf(face, jlim, ilim, i0=0, j0=0)
      !$acc loop gang private(iq, i_buf0)
      do n = f4%gc_f2c_to_buf_iface(${face}$), f4%gc_f2c_to_buf_iface(${face}$+1)-1
         iq = f4%gc_f2c_to_buf(1, n) + 1 ! fine block
         i_buf0 = f4%gc_f2c_to_buf(2, n) * n_vars

         !$acc loop collapse(3) private(ivar, j_f, i_f, i_buf)
         do iv = 1, n_vars
            do j = 1, ${jlim}$
               do i = 1, ${ilim}$
                  ivar = i_vars(iv)
                  j_f = ${j0}$ + 2 * j - 1
                  i_f = ${i0}$ + 2 * i - 1
                  i_buf = i_buf0 + (iv - 1) * ${ilim}$ * ${jlim}$ + (j - 1) * ${ilim}$ + i
                  f4%send_buffer(i_buf) = 0.25_dp * &
                       sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
               end do
            end do
         end do
      end do
#:enddef

      @:fyp_f2c_to_buf(0, half_bx(2), n_gc)

      @:fyp_f2c_to_buf(1, half_bx(2), n_gc, i0=bx(1) - 2*n_gc)

      @:fyp_f2c_to_buf(2, n_gc, half_bx(1))

      @:fyp_f2c_to_buf(3, n_gc, half_bx(1), j0=bx(2) - 2*n_gc)

      !$acc end parallel

      f4%recv_offset(:) = f4%gc_recv_offset * n_vars
      f4%send_offset(:) = f4%gc_send_offset * n_vars
    end associate

  end subroutine fill_ghostcell_buffers_round_one

  !> Fill buffers for the fine side of coarse-to-fine boundaries
  subroutine fill_ghostcell_buffers_round_two(f4, n_vars, i_vars)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer                      :: i, j, n, ivar, iv
    integer                      :: i_c, j_c, half_bx(2)
    integer                      :: iq, i_buf, i_buf0, offset, face
    integer                      :: half_n_gc
    logical                      :: odd_n_gc
    real(dp)                     :: fine(4)

    ! Update send/recv offsets
    f4%recv_offset(:) = f4%gc_recv_offset_c2f * n_vars
    f4%send_offset(:) = f4%gc_send_offset_c2f * n_vars

    ! If nothing to do, save time by not starting parallel region
    if (f4%gc_c2f_to_buf_iface(4) == 1) return

    if (maxval(f4%gc_send_offset_c2f) * n_vars > size(f4%send_buffer)) &
         error stop "send buffer too small"

    half_bx = f4%bx/2
    half_n_gc = f4%n_gc/2 ! Round down
    odd_n_gc  = (iand(f4%n_gc, 1) == 1)

    associate (bx => f4%bx, n_gc => f4%n_gc, uu => f4%uu)
      !$acc parallel

      face = 0
      !$acc loop gang private(iq, offset, i_buf0)
      do n = f4%gc_c2f_to_buf_iface(face), f4%gc_c2f_to_buf_iface(face+1)-1
         iq = f4%gc_c2f_to_buf(1, n) + 1 ! coarse block
         offset = f4%gc_c2f_to_buf(2, n)
         i_buf0 = f4%gc_c2f_to_buf(3, n) * n_vars

         !$acc loop collapse(3) private(ivar, j_c, i_c, i_buf)
         do iv = 1, n_vars
            do j = 1, half_bx(2)
               do i = 1, half_n_gc
                  ivar = i_vars(iv)
                  j_c = j + offset * half_bx(2)
                  i_c = i
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq), &
                       f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq), fine)

                  i_buf = i_buf0 + 4 * (((iv - 1) * half_bx(2) + (j - 1)) * half_n_gc + i - 1)
                  f4%send_buffer(i_buf+1:i_buf+4) = fine
               end do
            end do
         end do

         i_buf0 = i_buf0 + n_vars * half_bx(2) * half_n_gc * 4

         if (odd_n_gc) then
            !$acc loop collapse(2) private(ivar, j_c, i_c, i_buf)
            do iv = 1, n_vars
               do j = 1, half_bx(2)
                  ivar = i_vars(iv)
                  i_c = 1 + half_n_gc
                  j_c = j + offset * half_bx(2)
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq), &
                       f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq), fine)

                  i_buf = i_buf0 + 2 * ((iv - 1) * half_bx(2) + (j - 1))
                  f4%send_buffer(i_buf+1) = fine(1)
                  f4%send_buffer(i_buf+2) = fine(3)
               end do
            end do
         end if
      end do

      face = 1
      !$acc loop gang private(iq, offset, i_buf0)
      do n = f4%gc_c2f_to_buf_iface(face), f4%gc_c2f_to_buf_iface(face+1)-1
         iq = f4%gc_c2f_to_buf(1, n) + 1 ! coarse block
         offset = f4%gc_c2f_to_buf(2, n)
         i_buf0 = f4%gc_c2f_to_buf(3, n) * n_vars

         !$acc loop collapse(3) private(ivar, j_c, i_c, i_buf)
         do iv = 1, n_vars
            do j = 1, half_bx(2)
               do i = 1, half_n_gc
                  ivar = i_vars(iv)
                  j_c = j + offset * half_bx(2)
                  i_c = bx(1) - half_n_gc + i
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq), &
                       f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq), fine)
                  i_buf = i_buf0 + 4 * (((iv - 1) * half_bx(2) + (j - 1)) * half_n_gc + i - 1)
                  f4%send_buffer(i_buf+1:i_buf+4) = fine
               end do
            end do
         end do

         i_buf0 = i_buf0 + n_vars * half_bx(2) * half_n_gc * 4

         if (odd_n_gc) then
            !$acc loop collapse(2) private(ivar, j_c, i_c, i_buf)
            do iv = 1, n_vars
               do j = 1, half_bx(2)
                  ivar = i_vars(iv)
                  i_c = bx(1) - half_n_gc
                  j_c = j + offset * half_bx(2)
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq), &
                       f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq), fine)

                  i_buf = i_buf0 + 2 * ((iv - 1) * half_bx(2) + (j - 1))
                  f4%send_buffer(i_buf+1) = fine(2)
                  f4%send_buffer(i_buf+2) = fine(4)
               end do
            end do
         end if
      end do

      face = 2
      !$acc loop gang private(iq, offset, i_buf0)
      do n = f4%gc_c2f_to_buf_iface(face), f4%gc_c2f_to_buf_iface(face+1)-1
         iq = f4%gc_c2f_to_buf(1, n) + 1 ! coarse block
         offset = f4%gc_c2f_to_buf(2, n)
         i_buf0 = f4%gc_c2f_to_buf(3, n) * n_vars

         !$acc loop collapse(3) private(ivar, j_c, i_c, i_buf)
         do iv = 1, n_vars
            do j = 1, half_n_gc
               do i = 1, half_bx(1)
                  ivar = i_vars(iv)
                  j_c = j
                  i_c = i + offset * half_bx(1)
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq), &
                       f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq), fine)

                  i_buf = i_buf0 + 4 * (((iv - 1) * half_n_gc + (j - 1)) * half_bx(1) + i - 1)
                  f4%send_buffer(i_buf+1:i_buf+4) = fine
               end do
            end do
         end do

         i_buf0 = i_buf0 + n_vars * half_n_gc * half_bx(1) * 4

         if (odd_n_gc) then
            !$acc loop collapse(2) private(ivar, j_c, i_c, i_buf)
            do iv = 1, n_vars
               do i = 1, half_bx(1)
                  ivar = i_vars(iv)
                  j_c = 1 + half_n_gc
                  i_c = i + offset * half_bx(1)
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq), &
                       f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq), fine)

                  i_buf = i_buf0 + 2 * ((iv - 1) * half_bx(1) + (i - 1))
                  f4%send_buffer(i_buf+1) = fine(1)
                  f4%send_buffer(i_buf+2) = fine(2)
               end do
            end do
         end if
      end do

      face = 3
      !$acc loop gang private(iq, offset, i_buf0)
      do n = f4%gc_c2f_to_buf_iface(face), f4%gc_c2f_to_buf_iface(face+1)-1
         iq = f4%gc_c2f_to_buf(1, n) + 1 ! coarse block
         offset = f4%gc_c2f_to_buf(2, n)
         i_buf0 = f4%gc_c2f_to_buf(3, n) * n_vars

         !$acc loop collapse(3) private(ivar, j_c, i_c, i_buf)
         do iv = 1, n_vars
            do j = 1, half_n_gc
               do i = 1, half_bx(1)
                  ivar = i_vars(iv)
                  j_c = bx(2) - half_n_gc + j
                  i_c = i + offset * half_bx(1)
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq), &
                       f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq), fine)

                  i_buf = i_buf0 + 4 * (((iv - 1) * half_n_gc + (j - 1)) * half_bx(1) + i - 1)
                  f4%send_buffer(i_buf+1:i_buf+4) = fine
               end do
            end do
         end do

         i_buf0 = i_buf0 + n_vars * half_n_gc * half_bx(1) * 4

         if (odd_n_gc) then
            !$acc loop collapse(2) private(ivar, j_c, i_c, i_buf)
            do iv = 1, n_vars
               do i = 1, half_bx(1)
                  ivar = i_vars(iv)
                  j_c = bx(2) - half_n_gc
                  i_c = i + offset * half_bx(1)
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq), &
                       f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq), fine)

                  i_buf = i_buf0 + 2 * ((iv - 1) * half_bx(1) + (i - 1))
                  f4%send_buffer(i_buf+1) = fine(3)
                  f4%send_buffer(i_buf+2) = fine(4)
               end do
            end do
         end if
      end do

      !$acc end parallel
    end associate

  end subroutine fill_ghostcell_buffers_round_two

  !> Exchange the receive and send buffers according to the specified offsets
  !> per MPI rank
  subroutine f4_exchange_buffers(f4)
    type(foap4_t), intent(inout) :: f4
    type(MPI_Request)            :: send_req(0:f4%mpisize-1)
    type(MPI_Request)            :: recv_req(0:f4%mpisize-1)
    integer                      :: n_send, n_recv, ilo, ihi, ierr, rank
    integer, parameter           :: tag = 0

    n_send = 0
    n_recv = 0

    ! Use device pointers for device data
    !$acc host_data use_device(f4%recv_buffer, f4%send_buffer)
    do rank = 0, f4%mpisize - 1
       ilo = f4%send_offset(rank) + 1
       ihi = f4%send_offset(rank+1)

       if (ihi >= ilo) then
          n_send = n_send + 1
          call MPI_Isend(f4%send_buffer(ilo:ihi), ihi-ilo+1, &
               MPI_DOUBLE_PRECISION, rank, tag, f4%mpicomm, &
               send_req(n_send), ierr)
       end if

       ilo = f4%recv_offset(rank) + 1
       ihi = f4%recv_offset(rank+1)

       if (ihi >= ilo) then
          n_recv = n_recv + 1
          call MPI_Irecv(f4%recv_buffer(ilo:ihi), ihi-ilo+1, &
               MPI_DOUBLE_PRECISION, rank, tag, f4%mpicomm, &
               recv_req(n_recv), ierr)
       end if
    end do
    !$acc end host_data

    call MPI_Waitall(n_recv, recv_req(1:n_recv), MPI_STATUSES_IGNORE, ierr)
    call MPI_Waitall(n_send, send_req(1:n_send), MPI_STATUSES_IGNORE, ierr)

  end subroutine f4_exchange_buffers

  !> After buffers have been communicated, handle all ghost cells for "round
  !> one", which excludes coarse-to-fine interpolation
  subroutine fill_ghostcells_round_one(f4, n_vars, i_vars, bx, n_gc, &
       ilo, ihi, max_vars, max_blocks, uu)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer, intent(in)          :: bx(2)
    integer, intent(in)          :: n_gc
    integer, intent(in)          :: ilo(2)
    integer, intent(in)          :: ihi(2)
    integer, intent(in)          :: max_vars
    integer, intent(in)          :: max_blocks
    real(dp), intent(inout)      :: uu(ilo(1):ihi(1), ilo(2):ihi(2), max_vars, max_blocks)
    integer                      :: n, i, j, iq, jq, i_f, j_f, face
    integer                      :: i_buf, i_buf0, iv, ivar
    integer                      :: half_bx(2), offset

    half_bx = f4%bx/2

    !$acc parallel

    ! ----------------------------------------
    ! Fill local boundaries at the same refinement level
    ! ----------------------------------------
    face = 0
    !$acc loop gang private(iq, jq)
    do n = f4%gc_srl_local_iface(face), f4%gc_srl_local_iface(face+1)-1
       iq   = f4%gc_srl_local(1, n) + 1
       jq   = f4%gc_srl_local(2, n) + 1

       !$acc loop collapse(3) private(ivar)
       do iv = 1, n_vars
          do j = 1, bx(2)
             do i = 1, n_gc
                ivar = i_vars(iv)
                uu(-n_gc+i, j, ivar, iq) = uu(bx(1)-n_gc+i, j, ivar, jq)
                uu(bx(1)+i, j, ivar, jq) = uu(i, j, ivar, iq)
             end do
          end do
       end do
    end do

    face = 1
    !$acc loop gang private(iq, jq)
    do n = f4%gc_srl_local_iface(face), f4%gc_srl_local_iface(face+1)-1
       iq   = f4%gc_srl_local(1, n) + 1
       jq   = f4%gc_srl_local(2, n) + 1

       !$acc loop collapse(3) private(ivar)
       do iv = 1, n_vars
          do j = 1, bx(2)
             do i = 1, n_gc
                ivar = i_vars(iv)
                uu(bx(1)+i, j, ivar, iq) = uu(i, j, ivar, jq)
                uu(-n_gc+i, j, ivar, jq) = uu(bx(1)-n_gc+i, j, ivar, iq)
             end do
          end do
       end do
    end do

    face = 2
    !$acc loop gang private(iq, jq)
    do n = f4%gc_srl_local_iface(face), f4%gc_srl_local_iface(face+1)-1
       iq   = f4%gc_srl_local(1, n) + 1
       jq   = f4%gc_srl_local(2, n) + 1

       !$acc loop collapse(3) private(ivar)
       do iv = 1, n_vars
          do j = 1, n_gc
             do i = 1, bx(1)
                ivar = i_vars(iv)
                uu(i, -n_gc+j, ivar, iq) = uu(i, bx(2)-n_gc+j, ivar, jq)
                uu(i, bx(2)+j, ivar, jq) = uu(i, j, ivar, iq)
             end do
          end do
       end do
    end do

    face = 3
    !$acc loop gang private(iq, jq)
    do n = f4%gc_srl_local_iface(face), f4%gc_srl_local_iface(face+1)-1
       iq   = f4%gc_srl_local(1, n) + 1
       jq   = f4%gc_srl_local(2, n) + 1

       !$acc loop collapse(3) private(ivar)
       do iv = 1, n_vars
          do j = 1, n_gc
             do i = 1, bx(1)
                ivar = i_vars(iv)
                uu(i, bx(2)+j, ivar, iq) = uu(i, j, ivar, jq)
                uu(i, -n_gc+j, ivar, jq) = uu(i, bx(2)-n_gc+j, ivar, iq)
             end do
          end do
       end do
    end do

    ! ----------------------------------------
    ! Fill physical boundaries
    ! ----------------------------------------

    face = 0
    !$acc loop gang private(iq)
    do n = f4%gc_phys_iface(face), f4%gc_phys_iface(face+1)-1
       iq = f4%gc_phys(n) + 1

       !$acc loop collapse(3) private(ivar)
       do iv = 1, n_vars
          do j = 1, bx(2)
             do i = 1, n_gc
                ivar = i_vars(iv)
                uu(1-i, j, ivar, iq) = uu(1, j, ivar, iq) - i * &
                     (uu(2, j, ivar, iq) - uu(1, j, ivar, iq))
             end do
          end do
       end do
    end do

    face = 1
    !$acc loop gang private(iq)
    do n = f4%gc_phys_iface(face), f4%gc_phys_iface(face+1)-1
       iq = f4%gc_phys(n) + 1

       !$acc loop collapse(3) private(ivar)
       do iv = 1, n_vars
          do j = 1, bx(2)
             do i = 1, n_gc
                ivar = i_vars(iv)
                uu(bx(1)+i, j, ivar, iq) = uu(bx(1), j, ivar, iq) + i * &
                     (uu(bx(1), j, ivar, iq) - uu(bx(1)-1, j, ivar, iq))
             end do
          end do
       end do
    end do

    face = 2
    !$acc loop gang private(iq)
    do n = f4%gc_phys_iface(face), f4%gc_phys_iface(face+1)-1
       iq = f4%gc_phys(n) + 1

       !$acc loop collapse(3) private(ivar)
       do iv = 1, n_vars
          do j = 1, n_gc
             do i = 1, bx(1)
                ivar = i_vars(iv)
                uu(i, 1-j, ivar, iq) = uu(i, 1, ivar, iq) - j * &
                     (uu(i, 2, ivar, iq) - uu(i, 1, ivar, iq))
             end do
          end do
       end do
    end do

    face = 3
    !$acc loop gang private(iq)
    do n = f4%gc_phys_iface(face), f4%gc_phys_iface(face+1)-1
       iq = f4%gc_phys(n) + 1

       !$acc loop collapse(3) private(ivar)
       do iv = 1, n_vars
          do j = 1, n_gc
             do i = 1, bx(1)
                ivar = i_vars(iv)
                uu(i, bx(2)+j, ivar, iq) = uu(i, bx(2), ivar, iq) + j * &
                     (uu(i, bx(2), ivar, iq) - uu(i, bx(2)-1, ivar, iq))
             end do
          end do
       end do
    end do

    ! ----------------------------------------
    ! Fill coarse side of local fine-to-coarse refinement boundaries
    ! ----------------------------------------

    face = 0
    !$acc loop gang private(iq, jq, offset)
    do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
       iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
       jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
       offset = f4%gc_f2c_local(3, n)     ! offset

       !$acc loop collapse(3) private(i_f, j_f, ivar)
       do iv = 1, n_vars
          do j = 1, half_bx(2)
             do i = 1, n_gc
                ivar = i_vars(iv)
                j_f = 2 * j - 1
                i_f = 2 * i - 1
                uu(bx(1)+i, offset*half_bx(2)+j, ivar, jq) = 0.25_dp * &
                     sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
             end do
          end do
       end do
    end do

    face = 1
    !$acc loop gang private(iq, jq, offset)
    do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
       iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
       jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
       offset = f4%gc_f2c_local(3, n)     ! offset

       !$acc loop collapse(3) private(i_f, j_f, ivar)
       do iv = 1, n_vars
          do j = 1, half_bx(2)
             do i = 1, n_gc
                ivar = i_vars(iv)
                j_f = 2 * j - 1
                i_f = bx(1) - 2 * n_gc + 2 * i - 1
                uu(-n_gc+i, offset*half_bx(2)+j, ivar, jq) = 0.25_dp * &
                     sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
             end do
          end do
       end do
    end do

    face = 2
    !$acc loop gang private(iq, jq, offset)
    do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
       iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
       jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
       offset = f4%gc_f2c_local(3, n)     ! offset

       !$acc loop collapse(3) private(i_f, j_f, ivar)
       do iv = 1, n_vars
          do j = 1, n_gc
             do i = 1, half_bx(1)
                ivar = i_vars(iv)
                j_f = 2 * j - 1
                i_f = 2 * i - 1
                uu(offset*half_bx(1)+i, bx(2)+j, ivar, jq) = 0.25_dp * &
                     sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
             end do
          end do
       end do
    end do

    face = 3
    !$acc loop gang private(iq, jq, offset)
    do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
       iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
       jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
       offset = f4%gc_f2c_local(3, n)     ! offset

       !$acc loop collapse(3) private(i_f, j_f, ivar)
       do iv = 1, n_vars
          do j = 1, n_gc
             do i = 1, half_bx(1)
                ivar = i_vars(iv)
                j_f = bx(2) - 2 * n_gc + 2 * j - 1
                i_f = 2 * i - 1
                uu(offset*half_bx(1)+i, -n_gc+j, ivar, jq) = 0.25_dp * &
                     sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
             end do
          end do
       end do
    end do

    ! ----------------------------------------
    ! Fill ghost cells at the same refinement level from the buffer
    ! ----------------------------------------

    face = 0
    !$acc loop gang private(iq, i_buf0)
    do n = f4%gc_srl_from_buf_iface(face), f4%gc_srl_from_buf_iface(face+1)-1
       iq = f4%gc_srl_from_buf(1, n) + 1
       i_buf0 = f4%gc_srl_from_buf(2, n) * n_vars
       !$acc loop collapse(3) private(ivar, i_buf)
       do iv = 1, n_vars
          do j = 1, bx(2)
             do i = 1, n_gc
                ivar = i_vars(iv)
                i_buf = i_buf0 + ((iv - 1) * bx(2) + (j - 1)) * n_gc + i
                uu(-n_gc+i, j, ivar, iq) = f4%recv_buffer(i_buf)
             end do
          end do
       end do
    end do

    face = 1
    !$acc loop gang private(iq, i_buf0)
    do n = f4%gc_srl_from_buf_iface(face), f4%gc_srl_from_buf_iface(face+1)-1
       iq = f4%gc_srl_from_buf(1, n) + 1
       i_buf0 = f4%gc_srl_from_buf(2, n) * n_vars

       !$acc loop collapse(3) private(ivar, i_buf)
       do iv = 1, n_vars
          do j = 1, bx(2)
             do i = 1, n_gc
                ivar = i_vars(iv)
                i_buf = i_buf0 + ((iv - 1) * bx(2) + (j - 1)) * n_gc + i
                uu(bx(1)+i, j, ivar, iq) = f4%recv_buffer(i_buf)
             end do
          end do
       end do
    end do

    face = 2
    !$acc loop gang private(iq, i_buf0)
    do n = f4%gc_srl_from_buf_iface(face), f4%gc_srl_from_buf_iface(face+1)-1
       iq = f4%gc_srl_from_buf(1, n) + 1
       i_buf0 = f4%gc_srl_from_buf(2, n) * n_vars

       !$acc loop collapse(3) private(ivar, i_buf)
       do iv = 1, n_vars
          do j = 1, n_gc
             do i = 1, bx(1)
                ivar = i_vars(iv)
                i_buf = i_buf0 + ((iv - 1) * n_gc + (j - 1)) * bx(1) + i
                uu(i, -n_gc+j, ivar, iq) = f4%recv_buffer(i_buf)
             end do
          end do
       end do
    end do

    face = 3
    !$acc loop gang private(iq, i_buf0)
    do n = f4%gc_srl_from_buf_iface(face), f4%gc_srl_from_buf_iface(face+1)-1
       iq = f4%gc_srl_from_buf(1, n) + 1
       i_buf0 = f4%gc_srl_from_buf(2, n) * n_vars

       !$acc loop collapse(3) private(ivar, i_buf)
       do iv = 1, n_vars
          do j = 1, n_gc
             do i = 1, bx(1)
                ivar = i_vars(iv)
                i_buf = i_buf0 + ((iv - 1) * n_gc + (j - 1)) * bx(1) + i
                uu(i, bx(2)+j, ivar, iq) = f4%recv_buffer(i_buf)
             end do
          end do
       end do
    end do

    ! ----------------------------------------
    ! Update coarse side from buffers at coarse-to-fine buffer
    ! ----------------------------------------

    face = 0
    !$acc loop gang private(iq, offset, i_buf0)
    do n = f4%gc_c2f_from_buf_iface(face), f4%gc_c2f_from_buf_iface(face+1)-1
       iq     = f4%gc_c2f_from_buf(1, n) + 1 ! Coarse block
       offset = f4%gc_c2f_from_buf(2, n)     ! Offset
       i_buf0  = f4%gc_c2f_from_buf(3, n) * n_vars

       !$acc loop collapse(3) private(ivar, i_buf)
       do iv = 1, n_vars
          do j = 1, half_bx(2)
             do i = 1, n_gc
                ivar = i_vars(iv)
                i_buf = i_buf0 + ((iv - 1) * half_bx(2) + (j - 1)) * n_gc + i
                uu(-n_gc+i, offset*half_bx(2)+j, ivar, iq) = &
                     f4%recv_buffer(i_buf)
             end do
          end do
       end do
    end do

    face = 1
    !$acc loop gang private(iq, offset, i_buf0)
    do n = f4%gc_c2f_from_buf_iface(face), f4%gc_c2f_from_buf_iface(face+1)-1
       iq     = f4%gc_c2f_from_buf(1, n) + 1 ! Coarse block
       offset = f4%gc_c2f_from_buf(2, n)     ! Offset
       i_buf0  = f4%gc_c2f_from_buf(3, n) * n_vars

       !$acc loop collapse(3) private(ivar, i_buf)
       do iv = 1, n_vars
          do j = 1, half_bx(2)
             do i = 1, n_gc
                ivar = i_vars(iv)
                i_buf = i_buf0 + ((iv - 1) * half_bx(2) + (j - 1)) * n_gc + i
                uu(bx(1)+i, offset*half_bx(2)+j, ivar, iq) = &
                     f4%recv_buffer(i_buf)
             end do
          end do
       end do
    end do

    face = 2
    !$acc loop gang private(iq, offset, i_buf0)
    do n = f4%gc_c2f_from_buf_iface(face), f4%gc_c2f_from_buf_iface(face+1)-1
       iq     = f4%gc_c2f_from_buf(1, n) + 1 ! Coarse block
       offset = f4%gc_c2f_from_buf(2, n)     ! Offset
       i_buf0  = f4%gc_c2f_from_buf(3, n) * n_vars

       !$acc loop collapse(3) private(ivar, i_buf)
       do iv = 1, n_vars
          do j = 1, n_gc
             do i = 1, half_bx(1)
                ivar = i_vars(iv)
                i_buf = i_buf0 + ((iv - 1) * n_gc + (j - 1)) * half_bx(1) + i
                uu(offset*half_bx(1)+i, -n_gc+j, ivar, iq) = &
                     f4%recv_buffer(i_buf)
             end do
          end do
       end do
    end do

    face = 3
    !$acc loop gang private(iq, offset, i_buf0)
    do n = f4%gc_c2f_from_buf_iface(face), f4%gc_c2f_from_buf_iface(face+1)-1
       iq     = f4%gc_c2f_from_buf(1, n) + 1 ! Coarse block
       offset = f4%gc_c2f_from_buf(2, n)     ! Offset
       i_buf0  = f4%gc_c2f_from_buf(3, n) * n_vars

       !$acc loop collapse(3) private(ivar, i_buf)
       do iv = 1, n_vars
          do j = 1, n_gc
             do i = 1, half_bx(1)
                ivar = i_vars(iv)
                i_buf = i_buf0 + ((iv - 1) * n_gc + (j - 1)) * half_bx(1) + i
                uu(offset*half_bx(1)+i, bx(2)+j, ivar, iq) = &
                     f4%recv_buffer(i_buf)
             end do
          end do
       end do
    end do

    !$acc end parallel

  end subroutine fill_ghostcells_round_one

  !> Handle coarse-to-fine ghost cells on the fine side
  subroutine fill_ghostcells_round_two(f4, n_vars, i_vars, bx, n_gc, &
       ilo, ihi, max_vars, max_blocks, uu)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer, intent(in)          :: bx(2)
    integer, intent(in)          :: n_gc
    integer, intent(in)          :: ilo(2)
    integer, intent(in)          :: ihi(2)
    integer, intent(in)          :: max_vars
    integer, intent(in)          :: max_blocks
    real(dp), intent(inout)      :: uu(ilo(1):ihi(1), ilo(2):ihi(2), max_vars, max_blocks)
    integer                      :: n, i, j, iq, jq, i_c, j_c, i_f, j_f, face
    integer                      :: i_buf, i_buf0, iv, ivar
    integer                      :: half_bx(2), half_n_gc, offset
    logical                      :: odd_n_gc
    real(dp)                     :: fine(4)

    ! If nothing to do, save time by not starting parallel region
    if (f4%gc_f2c_local_iface(4) == 1 .and. &
         f4%gc_f2c_from_buf_iface(4) == 1) return

    half_bx   = f4%bx/2
    half_n_gc = f4%n_gc/2 ! Round down
    odd_n_gc  = (iand(f4%n_gc, 1) == 1)

    !$acc parallel

    ! ----------------------------------------
    ! Fill fine side of local coarse-to-fine boundaries
    ! ----------------------------------------

    face = 0
    !$acc loop gang private(iq, jq, offset)
    do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
       iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
       jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
       offset = f4%gc_f2c_local(3, n)     ! Offset

       !$acc loop collapse(3) private(ivar, j_f, j_c, i_f, i_c)
       do iv = 1, n_vars
          do j = 1, half_bx(2)
             do i = 1, half_n_gc
                ivar = i_vars(iv)
                j_f = 2 * j - 1
                j_c = j + offset * half_bx(2)
                i_c = bx(1) - half_n_gc + i
                i_f = -(2 * half_n_gc) + 2*i - 1
                call prolong_local_5point(uu(i_c, j_c, iv, jq), &
                     uu(i_c-1, j_c, iv, jq), uu(i_c+1, j_c, iv, jq), &
                     uu(i_c, j_c-1, iv, jq), uu(i_c, j_c+1, iv, jq), fine)
                uu(i_f  , j_f  , ivar, iq) = fine(1)
                uu(i_f+1, j_f  , ivar, iq) = fine(2)
                uu(i_f  , j_f+1, ivar, iq) = fine(3)
                uu(i_f+1, j_f+1, ivar, iq) = fine(4)
             end do
          end do
       end do

       if (odd_n_gc) then
          i_c = bx(1) - half_n_gc
          i_f = -n_gc + 1

          !$acc loop collapse(2) private(ivar, j_f, j_c)
          do iv = 1, n_vars
             do j = 1, half_bx(2)
                ivar = i_vars(iv)
                j_f = 2 * j - 1
                j_c = j + offset * half_bx(2)

                call prolong_local_5point(uu(i_c, j_c, iv, jq), &
                     uu(i_c-1, j_c, iv, jq), uu(i_c+1, j_c, iv, jq), &
                     uu(i_c, j_c-1, iv, jq), uu(i_c, j_c+1, iv, jq), fine)
                uu(i_f, j_f, ivar, iq) = fine(2)
                uu(i_f, j_f+1, ivar, iq) = fine(4)
             end do
          end do
       end if
    end do

    face = 1
    !$acc loop gang private(iq, jq, offset)
    do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
       iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
       jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
       offset = f4%gc_f2c_local(3, n)     ! Offset

       !$acc loop collapse(3) private(ivar, j_f, j_c, i_f, i_c)
       do iv = 1, n_vars
          do j = 1, half_bx(2)
             do i = 1, half_n_gc
                ivar = i_vars(iv)
                j_f = 2 * j - 1
                j_c = j + offset * half_bx(2)
                i_c = i
                i_f = bx(1) + 2*i - 1
                call prolong_local_5point(uu(i_c, j_c, iv, jq), &
                     uu(i_c-1, j_c, iv, jq), uu(i_c+1, j_c, iv, jq), &
                     uu(i_c, j_c-1, iv, jq), uu(i_c, j_c+1, iv, jq), fine)
                uu(i_f  , j_f  , ivar, iq) = fine(1)
                uu(i_f+1, j_f  , ivar, iq) = fine(2)
                uu(i_f  , j_f+1, ivar, iq) = fine(3)
                uu(i_f+1, j_f+1, ivar, iq) = fine(4)
             end do
          end do
       end do

       if (odd_n_gc) then
          i_c = 1 + half_n_gc
          i_f = bx(1) + n_gc

          !$acc loop collapse(2) private(ivar, j_f, j_c)
          do iv = 1, n_vars
             do j = 1, half_bx(2)
                ivar = i_vars(iv)
                j_f = 2 * j - 1
                j_c = j + offset * half_bx(2)
                call prolong_local_5point(uu(i_c, j_c, iv, jq), &
                     uu(i_c-1, j_c, iv, jq), uu(i_c+1, j_c, iv, jq), &
                     uu(i_c, j_c-1, iv, jq), uu(i_c, j_c+1, iv, jq), fine)
                uu(i_f, j_f, ivar, iq) = fine(1)
                uu(i_f, j_f+1, ivar, iq) = fine(3)
             end do
          end do
       end if
    end do

    face = 2
    !$acc loop gang private(iq, jq, offset)
    do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
       iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
       jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
       offset = f4%gc_f2c_local(3, n)     ! Offset

       !$acc loop collapse(3) private(ivar, j_f, j_c, i_f, i_c)
       do iv = 1, n_vars
          do j = 1, half_n_gc
             do i = 1, half_bx(1)
                ivar = i_vars(iv)
                j_c = bx(2) - half_n_gc + j
                j_f = -(2 * half_n_gc) + 2*j - 1
                i_f = 2 * i - 1
                i_c = i + offset * half_bx(1)
                call prolong_local_5point(uu(i_c, j_c, iv, jq), &
                     uu(i_c-1, j_c, iv, jq), uu(i_c+1, j_c, iv, jq), &
                     uu(i_c, j_c-1, iv, jq), uu(i_c, j_c+1, iv, jq), fine)
                uu(i_f  , j_f  , ivar, iq) = fine(1)
                uu(i_f+1, j_f  , ivar, iq) = fine(2)
                uu(i_f  , j_f+1, ivar, iq) = fine(3)
                uu(i_f+1, j_f+1, ivar, iq) = fine(4)
             end do
          end do
       end do

       if (odd_n_gc) then
          j_c = bx(2) - half_n_gc
          j_f = -n_gc + 1

          !$acc loop collapse(2) private(ivar, i_f, i_c)
          do iv = 1, n_vars
             do i = 1, half_bx(1)
                ivar = i_vars(iv)
                i_f = 2 * i - 1
                i_c = i + offset * half_bx(1)
                call prolong_local_5point(uu(i_c, j_c, iv, jq), &
                     uu(i_c-1, j_c, iv, jq), uu(i_c+1, j_c, iv, jq), &
                     uu(i_c, j_c-1, iv, jq), uu(i_c, j_c+1, iv, jq), fine)
                uu(i_f, j_f, ivar, iq) = fine(3)
                uu(i_f+1, j_f, ivar, iq) = fine(4)
             end do
          end do
       end if
    end do

    face = 3
    !$acc loop gang private(iq, jq, offset)
    do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
       iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
       jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
       offset = f4%gc_f2c_local(3, n)     ! Offset

       !$acc loop collapse(3) private(ivar, j_f, j_c, i_f, i_c)
       do iv = 1, n_vars
          do j = 1, half_n_gc
             do i = 1, half_bx(1)
                ivar = i_vars(iv)
                j_c = j
                j_f = bx(2) + 2*j - 1
                i_f = 2 * i - 1
                i_c = i + offset * half_bx(1)
                call prolong_local_5point(uu(i_c, j_c, iv, jq), &
                     uu(i_c-1, j_c, iv, jq), uu(i_c+1, j_c, iv, jq), &
                     uu(i_c, j_c-1, iv, jq), uu(i_c, j_c+1, iv, jq), fine)
                uu(i_f  , j_f  , ivar, iq) = fine(1)
                uu(i_f+1, j_f  , ivar, iq) = fine(2)
                uu(i_f  , j_f+1, ivar, iq) = fine(3)
                uu(i_f+1, j_f+1, ivar, iq) = fine(4)
             end do
          end do
       end do

       if (odd_n_gc) then
          j_c = 1 + half_n_gc
          j_f = bx(2) + n_gc

          !$acc loop collapse(2) private(ivar, i_f, i_c)
          do iv = 1, n_vars
             do i = 1, half_bx(1)
                ivar = i_vars(iv)
                i_f = 2 * i - 1
                i_c = i + offset * half_bx(1)
                call prolong_local_5point(uu(i_c, j_c, iv, jq), &
                     uu(i_c-1, j_c, iv, jq), uu(i_c+1, j_c, iv, jq), &
                     uu(i_c, j_c-1, iv, jq), uu(i_c, j_c+1, iv, jq), fine)
                uu(i_f, j_f, ivar, iq) = fine(1)
                uu(i_f+1, j_f, ivar, iq) = fine(2)
             end do
          end do
       end if
    end do

    ! ----------------------------------------
    ! Fill fine side of nonlocal coarse-to-fine boundaries
    ! ----------------------------------------

    face = 0
    !$acc loop gang private(iq, i_buf0)
    do n = f4%gc_f2c_from_buf_iface(face), f4%gc_f2c_from_buf_iface(face+1)-1
       iq    = f4%gc_f2c_from_buf(1, n) + 1 ! Fine block
       i_buf0 = f4%gc_f2c_from_buf(2, n) * n_vars

       !$acc loop collapse(3) private(ivar, j_f, i_f, i_buf)
       do iv = 1, n_vars
          do j = 1, half_bx(2)
             do i = 1, half_n_gc
                ivar = i_vars(iv)
                j_f = 2 * j - 1
                i_f = -(2 * half_n_gc) + 2*i - 1

                i_buf = i_buf0 + 4 * (((iv - 1) * half_bx(2) + (j - 1)) * half_n_gc + i - 1)
                uu(i_f  , j_f  , ivar, iq) = f4%recv_buffer(i_buf+1)
                uu(i_f+1, j_f  , ivar, iq) = f4%recv_buffer(i_buf+2)
                uu(i_f  , j_f+1, ivar, iq) = f4%recv_buffer(i_buf+3)
                uu(i_f+1, j_f+1, ivar, iq) = f4%recv_buffer(i_buf+4)
             end do
          end do
       end do

       if (odd_n_gc) then
          i_buf0 = i_buf0 + n_vars * half_bx(2) * half_n_gc * 4
          i_f = -n_gc + 1

          !$acc loop collapse(2) private(ivar, j_f, i_buf)
          do iv = 1, n_vars
             do j = 1, half_bx(2)
                ivar = i_vars(iv)
                j_f = 2 * j - 1

                i_buf = i_buf0 + 2 * ((iv - 1) * half_bx(2) + (j - 1))
                uu(i_f, j_f  , ivar, iq) = f4%recv_buffer(i_buf+1)
                uu(i_f, j_f+1, ivar, iq) = f4%recv_buffer(i_buf+2)
             end do
          end do
       end if
    end do

    face = 1
    !$acc loop gang private(iq, i_buf0)
    do n = f4%gc_f2c_from_buf_iface(face), f4%gc_f2c_from_buf_iface(face+1)-1
       iq    = f4%gc_f2c_from_buf(1, n) + 1 ! Fine block
       i_buf0 = f4%gc_f2c_from_buf(2, n) * n_vars

       !$acc loop collapse(3) private(ivar, j_f, i_f, i_buf)
       do iv = 1, n_vars
          do j = 1, half_bx(2)
             do i = 1, half_n_gc
                ivar = i_vars(iv)
                j_f = 2 * j - 1
                i_f = bx(1) + 2*i - 1

                i_buf = i_buf0 + 4 * (((iv - 1) * half_bx(2) + (j - 1)) * half_n_gc + i - 1)
                uu(i_f  , j_f  , ivar, iq) = f4%recv_buffer(i_buf+1)
                uu(i_f+1, j_f  , ivar, iq) = f4%recv_buffer(i_buf+2)
                uu(i_f  , j_f+1, ivar, iq) = f4%recv_buffer(i_buf+3)
                uu(i_f+1, j_f+1, ivar, iq) = f4%recv_buffer(i_buf+4)
             end do
          end do
       end do


       if (odd_n_gc) then
          i_buf0 = i_buf0 + n_vars * half_bx(2) * half_n_gc * 4
          i_f = bx(1) + n_gc

          !$acc loop collapse(2) private(ivar, j_f, i_buf)
          do iv = 1, n_vars
             do j = 1, half_bx(2)
                ivar = i_vars(iv)
                j_f = 2 * j - 1

                i_buf = i_buf0 + 2 * ((iv - 1) * half_bx(2) + (j - 1))
                uu(i_f, j_f  , ivar, iq) = f4%recv_buffer(i_buf+1)
                uu(i_f, j_f+1, ivar, iq) = f4%recv_buffer(i_buf+2)
             end do
          end do
       end if
    end do

    face = 2
    !$acc loop gang private(iq, i_buf0)
    do n = f4%gc_f2c_from_buf_iface(face), f4%gc_f2c_from_buf_iface(face+1)-1
       iq    = f4%gc_f2c_from_buf(1, n) + 1 ! Fine block
       i_buf0 = f4%gc_f2c_from_buf(2, n) * n_vars

       !$acc loop collapse(3) private(ivar, j_f, i_f, i_buf)
       do iv = 1, n_vars
          do j = 1, half_n_gc
             do i = 1, half_bx(1)
                ivar = i_vars(iv)
                j_f = -(2 * half_n_gc) + 2*j - 1
                i_f = 2 * i - 1

                i_buf = i_buf0 + 4 * (((iv - 1) * half_n_gc + (j - 1)) * half_bx(1) + i - 1)
                uu(i_f  , j_f  , ivar, iq) = f4%recv_buffer(i_buf+1)
                uu(i_f+1, j_f  , ivar, iq) = f4%recv_buffer(i_buf+2)
                uu(i_f  , j_f+1, ivar, iq) = f4%recv_buffer(i_buf+3)
                uu(i_f+1, j_f+1, ivar, iq) = f4%recv_buffer(i_buf+4)
             end do
          end do
       end do

       if (odd_n_gc) then
          i_buf0 = i_buf0 + n_vars * half_n_gc * half_bx(1) * 4
          j_f = -n_gc + 1

          !$acc loop collapse(2) private(ivar, i_f, i_buf)
          do iv = 1, n_vars
             do i = 1, half_bx(1)
                ivar = i_vars(iv)
                i_f = 2 * i - 1

                i_buf = i_buf0 + 2 * ((iv - 1) * half_bx(1) + (i - 1))
                uu(i_f  , j_f  , ivar, iq) = f4%recv_buffer(i_buf+1)
                uu(i_f+1, j_f  , ivar, iq) = f4%recv_buffer(i_buf+2)
             end do
          end do
       end if
    end do

    face = 3
    !$acc loop gang private(iq, i_buf0)
    do n = f4%gc_f2c_from_buf_iface(face), f4%gc_f2c_from_buf_iface(face+1)-1
       iq    = f4%gc_f2c_from_buf(1, n) + 1 ! Fine block
       i_buf0 = f4%gc_f2c_from_buf(2, n) * n_vars

       !$acc loop collapse(3) private(ivar, j_f, i_f, i_buf)
       do iv = 1, n_vars
          do j = 1, half_n_gc
             do i = 1, half_bx(1)
                ivar = i_vars(iv)
                j_f = bx(2) + 2*j - 1
                i_f = 2 * i - 1

                i_buf = i_buf0 + 4 * (((iv - 1) * half_n_gc + (j - 1)) * half_bx(1) + i - 1)
                uu(i_f  , j_f  , ivar, iq) = f4%recv_buffer(i_buf+1)
                uu(i_f+1, j_f  , ivar, iq) = f4%recv_buffer(i_buf+2)
                uu(i_f  , j_f+1, ivar, iq) = f4%recv_buffer(i_buf+3)
                uu(i_f+1, j_f+1, ivar, iq) = f4%recv_buffer(i_buf+4)
             end do
          end do
       end do

       if (odd_n_gc) then
          i_buf0 = i_buf0 + n_vars * half_n_gc * half_bx(1) * 4
          j_f = bx(2) + n_gc

          !$acc loop collapse(2) private(ivar, i_f, i_buf)
          do iv = 1, n_vars
             do i = 1, half_bx(1)
                ivar = i_vars(iv)
                i_f = 2 * i - 1

                i_buf = i_buf0 + 2 * ((iv - 1) * half_bx(1) + (i - 1))
                uu(i_f  , j_f  , ivar, iq) = f4%recv_buffer(i_buf+1)
                uu(i_f+1, j_f  , ivar, iq) = f4%recv_buffer(i_buf+2)
             end do
          end do
       end if
    end do

    !$acc end parallel

  end subroutine fill_ghostcells_round_two

  !> Update ghost cells for selected variables
  subroutine f4_update_ghostcells(f4, n_vars, i_vars)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    real(dp)                     :: t0, t1

    t0 = MPI_Wtime()
    call update_ghostcell_pattern(f4)
    t1 = MPI_Wtime()
    f4%wtime_update_gc_pattern = f4%wtime_update_gc_pattern + t1 - t0

    call fill_ghostcell_buffers_round_one(f4, n_vars, i_vars)
    t0 = MPI_Wtime()
    f4%wtime_gc_fill_buff_round1 = f4%wtime_gc_fill_buff_round1 + t0 - t1

    call f4_exchange_buffers(f4)
    t1 = MPI_Wtime()
    f4%wtime_exchange_buffers = f4%wtime_exchange_buffers + t1 - t0

    call fill_ghostcells_round_one(f4, n_vars, i_vars, f4%bx, f4%n_gc, &
         f4%ilo, f4%ihi, f4%n_vars, f4%max_blocks, f4%uu)

    t0 = MPI_Wtime()
    f4%wtime_gc_fill_round1 = f4%wtime_gc_fill_round1 + t0 - t1

    ! Do coarse-to-fine refinement boundaries last, so that ghost cells
    ! required for interpolation have been filled
    call fill_ghostcell_buffers_round_two(f4, n_vars, i_vars)

    t1 = MPI_Wtime()
    f4%wtime_gc_fill_buff_round2 = f4%wtime_gc_fill_buff_round2 + t1 - t0

    call f4_exchange_buffers(f4)

    t0 = MPI_Wtime()
    f4%wtime_exchange_buffers = f4%wtime_exchange_buffers + t0 - t1

    call fill_ghostcells_round_two(f4, n_vars, i_vars, f4%bx, f4%n_gc, &
         f4%ilo, f4%ihi, f4%n_vars, f4%max_blocks, f4%uu)

    t1 = MPI_Wtime()
    f4%wtime_gc_fill_round2 = f4%wtime_gc_fill_round2 + t1 - t0

  end subroutine f4_update_ghostcells

  !> Refine the mesh according to f4%refinement_flags
  subroutine f4_adjust_refinement(f4, partition_after)
    type(foap4_t), intent(inout) :: f4
    logical, intent(in)          :: partition_after
    integer                      :: n_blocks_new, n_blocks_old, iv
    integer                      :: has_changed, i_srl, i_refine, i_coarsen
    integer, allocatable         :: srl(:, :), refine(:, :), coarsen(:, :)
    integer                      :: n, i, j, i_c, j_c, i_f, j_f, n_old
    integer                      :: half_bx(2), offset_copy
    integer                      :: i_from, i_to, i_ch
    real(dp)                     :: fine(4), t0, t1

    t0 = MPI_Wtime()
    half_bx = f4%bx / 2

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
          n = n + 4
          n_old = n_old + 1
       case (-1)
          ! Block has been coarsened
          i_coarsen = i_coarsen + 1
          coarsen(:, i_coarsen) = [n_old, n]
          n = n + 1
          n_old = n_old + 4
       case default
          error stop "Refinement: difference in levels > 1"
       end select
    end do

    if (n_old /= offset_copy + n_blocks_old + 1) &
         error stop "Refinement: loops did not end simultaneously"

    t1 = MPI_Wtime()

    !$acc enter data copyin(srl, refine, coarsen, half_bx)

    ! Copy on device
    !$acc parallel loop private(i_from, i_to) async
    do n = 1, i_srl
       i_from = srl(1, n)
       i_to = srl(2, n)

       !$acc loop collapse(3)
       do iv = 1, f4%n_vars
          do j = f4%ilo(2), f4%ihi(2)
             do i = f4%ilo(1), f4%ihi(2)
                f4%uu(i, j, iv, i_to) = f4%uu(i, j, iv, i_from)
             end do
          end do
       end do
    end do

    ! Refine on device
    !$acc parallel loop private(i_from, i_to) async
    do n = 1, i_refine
       i_from = refine(1, n)
       i_to = refine(2, n)

       !$acc loop collapse(4) private(j_c, j_f, i_c, i_f, fine)
       do i_ch = 1, 4
          do iv = 1, f4%n_vars
             do j = 1, half_bx(2)
                do i = 1, half_bx(1)
                   j_c = j + child_offset(2, i_ch) * half_bx(2)
                   j_f = 2 * j - 1
                   i_c = i + child_offset(1, i_ch) * half_bx(1)
                   i_f = 2 * i - 1

                   call prolong_local_5point(f4%uu(i_c, j_c, iv, i_from), &
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

    ! Coarsen on device
    !$acc parallel loop private(i_from, i_to) async
    do n = 1, i_coarsen
       i_from = coarsen(1, n)
       i_to = coarsen(2, n)

       !$acc loop collapse(4) private(j_c, j_f, i_c, i_f)
       do i_ch = 1, 4
          do iv = 1, f4%n_vars
             do j = 1, half_bx(2)
                do i = 1, half_bx(1)
                   j_c = j + child_offset(2, i_ch) * half_bx(2)
                   j_f = 2 * j - 1
                   i_c = i + child_offset(1, i_ch) * half_bx(1)
                   i_f = 2 * i - 1

                   f4%uu(i_c, j_c, iv, i_to) = 0.25_dp * (&
                        f4%uu(i_f,   j_f, iv, i_from+i_ch-1) + &
                        f4%uu(i_f+1, j_f, iv, i_from+i_ch-1) + &
                        f4%uu(i_f,   j_f+1, iv, i_from+i_ch-1) + &
                        f4%uu(i_f+1, j_f+1, iv, i_from+i_ch-1))
                end do
             end do
          end do
       end do
    end do

    !$acc wait
    !$acc exit data delete(srl, refine, coarsen, half_bx)

    t0 = MPI_Wtime()
    f4%wtime_adjust_ref_foap4 = f4%wtime_adjust_ref_foap4 + t0 - t1

    if (partition_after) call f4_partition(f4)
  end subroutine f4_adjust_refinement

  !> Temporarily copy blocks to the end of the block array
  subroutine copy_blocks_to_end(f4, n_blocks_old, n_blocks_new, offset_copy)
    use iso_fortran_env, only: error_unit
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_blocks_old
    integer, intent(in)          :: n_blocks_new
    integer, intent(out)         :: offset_copy
    integer                      :: n, i, j, iv

    offset_copy = max(n_blocks_old, n_blocks_new)

    if (offset_copy + n_blocks_old > f4%max_blocks) then
       write(error_unit, "(A,I0,A,I0)") "ERROR: max_blocks = ", &
            f4%max_blocks, ", copy requires ", offset_copy + n_blocks_old
       error stop "Not enough block memory for copying"
    end if

    ! Copy block metadata on host
    do n = 1, n_blocks_old
       f4%block_origin(:, offset_copy+n) = f4%block_origin(:, n)
       f4%block_level(offset_copy+n)     = f4%block_level(n)
    end do

    ! Copy block solution data on device
    !$acc parallel loop
    do n = 1, n_blocks_old
       !$acc loop collapse(3)
       do iv = 1, f4%n_vars
          do j = f4%ilo(2), f4%ihi(2)
             do i = f4%ilo(1), f4%ihi(1)
                f4%uu(i, j, iv, offset_copy+n) = f4%uu(i, j, iv, n)
             end do
          end do
       end do
    end do
  end subroutine copy_blocks_to_end

  !> Method for prolongation (interpolation) of a coarse block to its children
  subroutine prolong_local_5point(center, xlo, xhi, ylo, yhi, fine)
    !$acc routine seq
    real(dp), intent(in)  :: center ! Center value
    real(dp), intent(in)  :: xlo, xhi ! x-neighbors (-1, +1)
    real(dp), intent(in)  :: ylo, yhi ! y-neighbors (-1, +1)
    real(dp), intent(out) :: fine(4)
    real(dp)              :: f(0:2), slopes_a(2), slopes_b(2)

    f(0) = center          ! Identical to coarse_y(2)
    slopes_a(1) = center - xlo
    slopes_a(2) = center - ylo
    slopes_b(1) = xhi - center
    slopes_b(2) = yhi - center
    f(1:2) = 0.25_dp * limiter_minmod(slopes_a, slopes_b)

    fine(1) = f(0) - f(1) - f(2)
    fine(2) = f(0) + f(1) - f(2)
    fine(3) = f(0) - f(1) + f(2)
    fine(4) = f(0) + f(1) + f(2)
  end subroutine prolong_local_5point

  !> Generalized minmod limiter. The parameter theta controls how dissipative
  !> the limiter is, with 1 corresponding to the minmod limiter and 2 to the MC
  !> limiter.
  elemental function limiter_gminmod(a, b, theta) result(phi)
    !$acc routine seq
    real(dp), intent(in) :: a, b, theta
    real(dp)             :: phi

    if (a * b > 0) then
       phi = sign(minval(abs([theta * a, theta * b, &
            0.5_dp * (a + b)])), a)
    else
       phi = 0.0_dp
    end if
  end function limiter_gminmod

  !> Minmod limiter
  elemental function limiter_minmod(a, b) result(phi)
    !$acc routine seq
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp)             :: phi
    phi = limiter_gminmod(a, b, 1.0_dp)
  end function limiter_minmod

  !> Partition the blocks over the MPI ranks
  !> TODO: try strided transfer, of only first N variables per block
  subroutine f4_partition(f4)
    type(foap4_t), intent(inout)   :: f4
    integer                        :: n_changed_global
    integer                        :: n_blocks_old, n_blocks_new
    integer(c_int64_t)             :: gfq_old(0:f4%mpisize)
    integer(c_int64_t)             :: gfq_new(0:f4%mpisize)
    integer                        :: dsize, offset_copy, ierr
    integer, parameter             :: tag = 0
    integer                        :: gfq_src(0:f4%mpisize)
    integer                        :: gfq_dest(0:f4%mpisize)
    integer                        :: dest_begin, dest_end
    integer                        :: src_begin, src_end
    integer                        :: gend, gbegin, n_blocks_transfer
    integer                        :: first_sender, last_sender
    integer                        :: first_receiver, last_receiver
    integer                        :: num_senders, num_receivers
    integer                        :: n_recv, n_send
    integer                        :: rank, block_ix
    type(MPI_Request), allocatable :: send_req(:), recv_req(:)
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

       num_senders = last_sender - first_sender + 1
       allocate(recv_req(num_senders))

       gend = dest_begin
       block_ix = 1

       !$acc host_data use_device(f4%uu)
       do rank = first_sender, last_sender
          gbegin = gend
          gend = min(dest_end, gfq_src(rank + 1))
          n_blocks_transfer = gend - gbegin

          if (n_blocks_transfer > 0) then
             n_recv = n_recv + 1
             call MPI_Irecv (f4%uu(:, :, :, block_ix), dsize*n_blocks_transfer, &
                  MPI_DOUBLE_PRECISION, rank, tag, f4%mpicomm, recv_req(n_recv), ierr)
             block_ix = block_ix + n_blocks_transfer
          end if
       end do
       !$acc end host_data
    end if

    if (src_end > src_begin) then
       ! Find index of first/last receiver, -1 to get zero-based index
       first_receiver = find_bracket(f4%mpisize+1, gfq_dest, src_begin) - 1
       last_receiver = find_bracket(f4%mpisize+1, gfq_dest, src_end-1) - 1

       num_receivers = last_receiver - first_receiver + 1
       allocate(send_req(num_receivers))

       gend = src_begin
       block_ix = offset_copy + 1

       !$acc host_data use_device(f4%uu)
       do rank = first_receiver, last_receiver
          gbegin = gend
          gend = min(src_end, gfq_dest(rank + 1))
          n_blocks_transfer = gend - gbegin

          if (n_blocks_transfer > 0) then
             n_send = n_send + 1
             call MPI_Isend (f4%uu(:, :, :, block_ix), dsize*n_blocks_transfer, &
                  MPI_DOUBLE_PRECISION, rank, tag, f4%mpicomm, send_req(n_send), ierr)
             block_ix = block_ix + n_blocks_transfer
          end if
       end do
       !$acc end host_data
    end if

    call MPI_Waitall(n_recv, recv_req(1:n_recv), MPI_STATUSES_IGNORE, ierr)
    call MPI_Waitall(n_send, send_req(1:n_send), MPI_STATUSES_IGNORE, ierr)

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

end module m_foap4
