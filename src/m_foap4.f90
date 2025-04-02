!> Foap4 stands for "Fortran OpenACC p4est", and implements MPI-parallel
!> quadtree/octree adaptive mesh refinement with support for OpenACC. Together
!> with p4est_wrapper.c, this module implements the required data structures and
!> methods.
!>
!> Author(s): Jannis Teunissen
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
     integer(c_int) :: other_proc !< MPI rank that owns quadid[1]
     integer(c_int) :: quadid(2)  !< quadid[0] is always local, [1] can be non-local
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
     subroutine pw_partition(pw, n_changed_global, gfq_old) bind(c)
       import c_int, c_int64_t, c_ptr
       type(c_ptr), intent(in), value    :: pw
       integer(c_int), intent(out)       :: n_changed_global
       integer(c_int64_t), intent(out)   :: gfq_old(*)
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
    integer                      :: log_level

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
  end subroutine f4_initialize

  !> Destroy all data for the current mesh
  subroutine f4_destroy(f4)
    type(foap4_t), intent(inout) :: f4

    ! Remove data from device
    !$acc exit data delete(f4)
    !$acc exit data delete(f4%block_level, f4%block_origin, f4%uu)
    !$acc exit data delete(f4%recv_buffer, f4%send_buffer)
    !$acc exit data delete(&
    !$acc &f4%gc_srl_local_iface, f4%gc_srl_from_buf_iface, f4%gc_srl_to_buf_iface, &
    !$acc &f4%gc_f2c_local_iface, f4%gc_f2c_from_buf_iface, f4%gc_f2c_to_buf_iface, &
    !$acc &f4%gc_c2f_from_buf_iface, f4%gc_c2f_to_buf_iface, f4%gc_phys_iface)

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

    integer :: i, ierr, periodic_as_int(2)

    if (bx(1) /= bx(2)) error stop "Unequal bx(:) not supported"
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

    call MPI_COMM_RANK(f4%mpicomm, f4%mpirank, ierr)
    call MPI_COMM_SIZE(f4%mpicomm, f4%mpisize, ierr)

    call pw_set_connectivity_brick(f4%pw, &
         trees_per_dim(1), trees_per_dim(2), &
         periodic_as_int(1), periodic_as_int(2), min_level, 1, max_blocks)

    allocate(f4%block_origin(2, max_blocks))
    allocate(f4%block_level(max_blocks))
    allocate(f4%refinement_flags(max_blocks))
    allocate(f4%uu(1-n_gc:bx(1)+n_gc, 1-n_gc:bx(2)+n_gc, n_vars, max_blocks))
    f4%uu(:, :, :, :) = 0.0_dp

    f4%gc_data_size = f4%bx(1) * f4%n_gc

    ! Maximum size of recv/send buffer
    i = max_blocks * 4 * f4%gc_data_size
    allocate(f4%recv_buffer(i))
    allocate(f4%send_buffer(i))
    allocate(f4%recv_offset(0:f4%mpisize))
    allocate(f4%send_offset(0:f4%mpisize))

    f4%uu = 0.0_dp

    call f4_set_quadrants(f4)

    ! Copy data structure and allocatable components to device
    !$acc enter data copyin(f4)
    !$acc enter data copyin(f4%block_level, f4%block_origin, f4%uu)
    !$acc enter data create(f4%recv_buffer, f4%send_buffer)
    !$acc enter data create(&
    !$acc &f4%gc_srl_local_iface, f4%gc_srl_from_buf_iface, f4%gc_srl_to_buf_iface, &
    !$acc &f4%gc_f2c_local_iface, f4%gc_f2c_from_buf_iface, f4%gc_f2c_to_buf_iface, &
    !$acc &f4%gc_c2f_from_buf_iface, f4%gc_c2f_to_buf_iface, f4%gc_phys_iface)

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
    type(foap4_t), intent(inout) :: f4
    integer                      :: n

    f4%n_blocks = f4_get_num_local_blocks(f4)
    if (.not. allocated(f4%block_origin)) error stop "block_origin not allocated"
    if (.not. allocated(f4%block_level)) error stop "block_level not allocated"
    if (f4%n_blocks > f4%max_blocks) error stop "n_blocks > max_blocks"

    call pw_get_quadrants(f4%pw, f4%n_blocks, &
         f4%block_origin(:, 1:f4%n_blocks), &
         f4%block_level(1:f4%n_blocks))

    do n = 1, f4%n_blocks
       f4%block_origin(:, n) = f4%block_origin(:, n) * f4%tree_length
    end do

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
    type(foap4_t), intent(in)              :: f4
    character(len=*), intent(in)           :: fname !< Base file name
    integer, intent(in)                    :: n_output !< Output index
    real(dp), intent(in), optional         :: time !< Simulation time
    character(len=*), intent(in), optional :: viewer !< Optimize for viewer
    character(len=len_trim(fname)+7)       :: full_fname
    integer                                :: n
    real(dp), allocatable                  :: dr(:, :)

    ! Get the block data from the device
    !$acc update self(f4%uu)

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
  end subroutine f4_write_grid

  !> Return the coordinates at the center of a grid cells
  pure function f4_cell_coord(f4, i_block, i, j) result(rr)
    type(foap4_t), intent(in) :: f4
    integer, intent(in)       :: i_block, i, j
    real(dp)                  :: rr(2), dr(2)

    dr = f4%dr_level(:, f4%block_level(i_block))
    rr = f4%block_origin(:, i_block) + dr * [i-0.5_dp, j-0.5_dp]
  end function f4_cell_coord

  !> Update the information required to update ghost cells
  subroutine f4_update_ghostcell_pattern(f4)
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

  end subroutine f4_update_ghostcell_pattern

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
       ! Deallocate arrays, since they probably changed size

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

    ! Copy/sync data to device

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

  !> Fill buffers for 'round 1' of the ghost cell exchange
  subroutine f4_fill_ghostcell_buffers(f4, n_vars, i_vars)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer                      :: i, j, n, ivar, iv
    integer                      :: i_f, j_f, half_bx(2)
    integer                      :: iq, i_buf, face

    if (maxval(f4%gc_send_offset) * n_vars > size(f4%send_buffer)) &
         error stop "send buffer too small"

    half_bx = f4%bx/2

    associate (bx => f4%bx, n_gc => f4%n_gc, uu => f4%uu)

      face = 0
      do n = f4%gc_srl_to_buf_iface(face), f4%gc_srl_to_buf_iface(face+1)-1
         iq = f4%gc_srl_to_buf(1, n) + 1
         i_buf = f4%gc_srl_to_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, bx(2)
               do i = 1, n_gc
                  i_buf = i_buf + 1
                  f4%send_buffer(i_buf) = uu(i, j, ivar, iq)
               end do
            end do
         end do
      end do

      face = 1
      do n = f4%gc_srl_to_buf_iface(face), f4%gc_srl_to_buf_iface(face+1)-1
         iq = f4%gc_srl_to_buf(1, n) + 1
         i_buf = f4%gc_srl_to_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, bx(2)
               do i = 1, n_gc
                  i_buf = i_buf + 1
                  f4%send_buffer(i_buf) = uu(bx(1)-n_gc+i, j, ivar, iq)
               end do
            end do
         end do
      end do

      face = 2
      do n = f4%gc_srl_to_buf_iface(face), f4%gc_srl_to_buf_iface(face+1)-1
         iq = f4%gc_srl_to_buf(1, n) + 1
         i_buf = f4%gc_srl_to_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, n_gc
               do i = 1, bx(1)
                  i_buf = i_buf + 1
                  f4%send_buffer(i_buf) = uu(i, j, ivar, iq)
               end do
            end do
         end do
      end do

      face = 3
      do n = f4%gc_srl_to_buf_iface(face), f4%gc_srl_to_buf_iface(face+1)-1
         iq = f4%gc_srl_to_buf(1, n) + 1
         i_buf = f4%gc_srl_to_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, n_gc
               do i = 1, bx(1)
                  i_buf = i_buf + 1
                  f4%send_buffer(i_buf) = uu(i, bx(2)-n_gc+j, ivar, iq)
               end do
            end do
         end do
      end do

      ! Nonlocal fine-to-coarse boundaries, fill buffer for coarse side

      face = 0
      do n = f4%gc_f2c_to_buf_iface(face), f4%gc_f2c_to_buf_iface(face+1)-1
         iq = f4%gc_f2c_to_buf(1, n) + 1 ! fine block
         i_buf = f4%gc_f2c_to_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, half_bx(2)
               j_f = 2 * j - 1
               do i = 1, n_gc
                  i_f = 2 * i - 1
                  i_buf = i_buf + 1
                  f4%send_buffer(i_buf) = 0.25_dp * &
                       sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
               end do
            end do
         end do
      end do

      face = 1
      do n = f4%gc_f2c_to_buf_iface(face), f4%gc_f2c_to_buf_iface(face+1)-1
         iq = f4%gc_f2c_to_buf(1, n) + 1 ! fine block
         i_buf = f4%gc_f2c_to_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, half_bx(2)
               j_f = 2 * j - 1
               do i = 1, n_gc
                  i_f = bx(1) - 2 * n_gc + 2 * i - 1
                  i_buf = i_buf + 1
                  f4%send_buffer(i_buf) = 0.25_dp * &
                       sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
               end do
            end do
         end do
      end do

      face = 2
      do n = f4%gc_f2c_to_buf_iface(face), f4%gc_f2c_to_buf_iface(face+1)-1
         iq = f4%gc_f2c_to_buf(1, n) + 1 ! fine block
         i_buf = f4%gc_f2c_to_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, n_gc
               j_f = 2 * j - 1
               do i = 1, half_bx(1)
                  i_f = 2 * i - 1
                  i_buf = i_buf + 1
                  f4%send_buffer(i_buf) = 0.25_dp * &
                       sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
               end do
            end do
         end do
      end do

      face = 3
      do n = f4%gc_f2c_to_buf_iface(face), f4%gc_f2c_to_buf_iface(face+1)-1
         iq = f4%gc_f2c_to_buf(1, n) + 1 ! fine block
         i_buf = f4%gc_f2c_to_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, n_gc
               j_f = bx(2) - 2 * n_gc + 2 * j - 1
               do i = 1, half_bx(1)
                  i_f = 2 * i - 1
                  i_buf = i_buf + 1
                  f4%send_buffer(i_buf) = 0.25_dp * &
                       sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
               end do
            end do
         end do
      end do

      f4%recv_offset(:) = f4%gc_recv_offset * n_vars
      f4%send_offset(:) = f4%gc_send_offset * n_vars
    end associate

  end subroutine f4_fill_ghostcell_buffers

  !> Fill buffers for 'round 2' of the ghost cell exchange
  subroutine f4_fill_ghostcell_buffers_c2f(f4, n_vars, i_vars)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer                      :: i, j, n, ivar, iv
    integer                      :: i_c, j_c, half_bx(2)
    integer                      :: iq, i_buf, offset, face
    integer                      :: half_n_gc
    logical                      :: odd_n_gc
    real(dp)                     :: fine(4)

    if (maxval(f4%gc_send_offset_c2f) * n_vars > size(f4%send_buffer)) &
         error stop "send buffer too small"

    half_bx = f4%bx/2
    half_n_gc = f4%n_gc/2 ! Round down
    odd_n_gc  = (iand(f4%n_gc, 1) == 1)

    associate (bx => f4%bx, n_gc => f4%n_gc, uu => f4%uu)

      face = 0
      do n = f4%gc_c2f_to_buf_iface(face), f4%gc_c2f_to_buf_iface(face+1)-1
         iq = f4%gc_c2f_to_buf(1, n) + 1 ! coarse block
         offset = f4%gc_c2f_to_buf(2, n)
         i_buf = f4%gc_c2f_to_buf(3, n) * n_vars

         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, half_bx(2)
               j_c = j + offset * half_bx(2)

               do i = 1, half_n_gc
                  i_c = i
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                       [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                  f4%send_buffer(i_buf+1:i_buf+4) = fine
                  i_buf = i_buf + 4
               end do

               if (odd_n_gc) then
                  i_c = 1 + half_n_gc
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                       [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                  f4%send_buffer(i_buf+1:i_buf+2) = fine([1, 3])
                  i_buf = i_buf + 2
               end if
            end do
         end do
      end do

      face = 1
      do n = f4%gc_c2f_to_buf_iface(face), f4%gc_c2f_to_buf_iface(face+1)-1
         iq = f4%gc_c2f_to_buf(1, n) + 1 ! coarse block
         offset = f4%gc_c2f_to_buf(2, n)
         i_buf = f4%gc_c2f_to_buf(3, n) * n_vars

         do iv = 1, n_vars
            ivar = i_vars(iv)

            do j = 1, half_bx(2)
               j_c = j + offset * half_bx(2)
               do i = 1, half_n_gc
                  i_c = bx(1) - half_n_gc + i
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                       [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                  f4%send_buffer(i_buf+1:i_buf+4) = fine
                  i_buf = i_buf + 4
               end do

               if (odd_n_gc) then
                  i_c = bx(1) - half_n_gc
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                       [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                  f4%send_buffer(i_buf+1:i_buf+2) = fine([2, 4])
                  i_buf = i_buf + 2
               end if
            end do
         end do
      end do

      face = 2
      do n = f4%gc_c2f_to_buf_iface(face), f4%gc_c2f_to_buf_iface(face+1)-1
         iq = f4%gc_c2f_to_buf(1, n) + 1 ! coarse block
         offset = f4%gc_c2f_to_buf(2, n)
         i_buf = f4%gc_c2f_to_buf(3, n) * n_vars

         do iv = 1, n_vars
            ivar = i_vars(iv)

            do j = 1, half_n_gc
               j_c = j
               do i = 1, half_bx(1)
                  i_c = i + offset * half_bx(1)
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                       [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                  f4%send_buffer(i_buf+1:i_buf+4) = fine
                  i_buf = i_buf + 4
               end do
            end do

            if (odd_n_gc) then
               j_c = 1 + half_n_gc
               do i = 1, half_bx(1)
                  i_c = i + offset * half_bx(1)
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                       [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                  f4%send_buffer(i_buf+1:i_buf+2) = fine([1, 2])
                  i_buf = i_buf + 2
               end do
            end if
         end do
      end do

      face = 3
      do n = f4%gc_c2f_to_buf_iface(face), f4%gc_c2f_to_buf_iface(face+1)-1
         iq = f4%gc_c2f_to_buf(1, n) + 1 ! coarse block
         offset = f4%gc_c2f_to_buf(2, n)
         i_buf = f4%gc_c2f_to_buf(3, n) * n_vars

         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, half_n_gc
               j_c = bx(2) - half_n_gc + j
               do i = 1, half_bx(1)
                  i_c = i + offset * half_bx(1)
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                       [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                  f4%send_buffer(i_buf+1:i_buf+4) = fine
                  i_buf = i_buf + 4
               end do
            end do

            if (odd_n_gc) then
               j_c = bx(2) - half_n_gc
               do i = 1, half_bx(1)
                  i_c = i + offset * half_bx(1)
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                       [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                       [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                  f4%send_buffer(i_buf+1:i_buf+2) = fine([3, 4])
                  i_buf = i_buf + 2
               end do
            end if
         end do
      end do

      f4%recv_offset(:) = f4%gc_recv_offset_c2f * n_vars
      f4%send_offset(:) = f4%gc_send_offset_c2f * n_vars
    end associate

  end subroutine f4_fill_ghostcell_buffers_c2f

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

    do rank = 0, f4%mpisize - 1
       ilo = f4%send_offset(rank) + 1
       ihi = f4%send_offset(rank+1)

       if (ihi >= ilo) then
          n_send = n_send + 1
          call mpi_isend(f4%send_buffer(ilo:ihi), ihi-ilo+1, &
               MPI_DOUBLE_PRECISION, rank, tag, f4%mpicomm, &
               send_req(n_send), ierr)
       end if

       ilo = f4%recv_offset(rank) + 1
       ihi = f4%recv_offset(rank+1)

       if (ihi >= ilo) then
          n_recv = n_recv + 1
          call mpi_irecv(f4%recv_buffer(ilo:ihi), ihi-ilo+1, MPI_DOUBLE_PRECISION, &
               rank, 0, f4%mpicomm, recv_req(n_recv), ierr)
       end if
    end do

    call mpi_waitall(n_recv, recv_req(1:n_recv), MPI_STATUSES_IGNORE, ierr)
    call mpi_waitall(n_send, send_req(1:n_send), MPI_STATUSES_IGNORE, ierr)

  end subroutine f4_exchange_buffers

  !> Update ghost cells for selected variables
  subroutine f4_update_ghostcells(f4, n_vars, i_vars)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer                      :: n, i, j, iq, jq, face
    integer                      :: i_buf, iv, ivar, i_f, j_f, i_c, j_c
    integer                      :: half_bx(2), half_n_gc, offset
    real(dp)                     :: slope, fine(4)
    logical                      :: odd_n_gc

    call f4_update_ghostcell_pattern(f4)
    call f4_fill_ghostcell_buffers(f4, n_vars, i_vars)
    call f4_exchange_buffers(f4)

    half_bx   = f4%bx/2
    half_n_gc = f4%n_gc/2 ! Round down
    odd_n_gc  = (iand(f4%n_gc, 1) == 1)

    associate (bx => f4%bx, n_gc => f4%n_gc, uu => f4%uu)

      face = 0
      do n = f4%gc_srl_from_buf_iface(face), f4%gc_srl_from_buf_iface(face+1)-1
         iq = f4%gc_srl_from_buf(1, n) + 1
         i_buf = f4%gc_srl_from_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, bx(2)
               do i = 1, n_gc
                  i_buf = i_buf + 1
                  uu(-n_gc+i, j, ivar, iq) = f4%recv_buffer(i_buf)
               end do
            end do
         end do
      end do

      face = 1
      do n = f4%gc_srl_from_buf_iface(face), f4%gc_srl_from_buf_iface(face+1)-1
         iq = f4%gc_srl_from_buf(1, n) + 1
         i_buf = f4%gc_srl_from_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, bx(2)
               do i = 1, n_gc
                  i_buf = i_buf + 1
                  uu(bx(1)+i, j, ivar, iq) = f4%recv_buffer(i_buf)
               end do
            end do
         end do
      end do

      face = 2
      do n = f4%gc_srl_from_buf_iface(face), f4%gc_srl_from_buf_iface(face+1)-1
         iq = f4%gc_srl_from_buf(1, n) + 1
         i_buf = f4%gc_srl_from_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, n_gc
               do i = 1, bx(1)
                  i_buf = i_buf + 1
                  uu(i, -n_gc+j, ivar, iq) = f4%recv_buffer(i_buf)
               end do
            end do
         end do
      end do

      face = 3
      do n = f4%gc_srl_from_buf_iface(face), f4%gc_srl_from_buf_iface(face+1)-1
         iq = f4%gc_srl_from_buf(1, n) + 1
         i_buf = f4%gc_srl_from_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, n_gc
               do i = 1, bx(1)
                  i_buf = i_buf + 1
                  uu(i, bx(2)+j, ivar, iq) = f4%recv_buffer(i_buf)
               end do
            end do
         end do
      end do

      ! Local fine-to-coarse refinement boundaries, fill coarse side
      face = 0
      do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
         iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
         jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
         offset = f4%gc_f2c_local(3, n)     ! offset
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, half_bx(2)
               j_f = 2 * j - 1
               do i = 1, n_gc
                  i_f = 2 * i - 1
                  uu(bx(1)+i, offset*half_bx(2)+j, ivar, jq) = 0.25_dp * &
                       sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
               end do
            end do
         end do
      end do

      face = 1
      do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
         iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
         jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
         offset = f4%gc_f2c_local(3, n)     ! offset
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, half_bx(2)
               j_f = 2 * j - 1
               do i = 1, n_gc
                  i_f = bx(1) - 2 * n_gc + 2 * i - 1
                  uu(-n_gc+i, offset*half_bx(2)+j, ivar, jq) = 0.25_dp * &
                       sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
               end do
            end do
         end do
      end do

      face = 2
      do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
         iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
         jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
         offset = f4%gc_f2c_local(3, n)     ! offset
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, n_gc
               j_f = 2 * j - 1
               do i = 1, half_bx(1)
                  i_f = 2 * i - 1
                  uu(offset*half_bx(1)+i, bx(2)+j, ivar, jq) = 0.25_dp * &
                       sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
               end do
            end do
         end do
      end do

      face = 3
      do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
         iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
         jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
         offset = f4%gc_f2c_local(3, n)     ! offset
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, n_gc
               j_f = bx(2) - 2 * n_gc + 2 * j - 1
               do i = 1, half_bx(1)
                  i_f = 2 * i - 1
                  uu(offset*half_bx(1)+i, -n_gc+j, ivar, jq) = 0.25_dp * &
                       sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
               end do
            end do
         end do
      end do

      ! Coarse-to-fine, update coarse side from buffers
      face = 0
      do n = f4%gc_c2f_from_buf_iface(face), f4%gc_c2f_from_buf_iface(face+1)-1
         iq     = f4%gc_c2f_from_buf(1, n) + 1 ! Coarse block
         offset = f4%gc_c2f_from_buf(2, n)     ! Offset
         i_buf  = f4%gc_c2f_from_buf(3, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, half_bx(2)
               do i = 1, n_gc
                  i_buf = i_buf + 1
                  uu(-n_gc+i, offset*half_bx(2)+j, ivar, iq) = &
                       f4%recv_buffer(i_buf)
               end do
            end do
         end do
      end do

      face = 1
      do n = f4%gc_c2f_from_buf_iface(face), f4%gc_c2f_from_buf_iface(face+1)-1
         iq     = f4%gc_c2f_from_buf(1, n) + 1 ! Coarse block
         offset = f4%gc_c2f_from_buf(2, n)     ! Offset
         i_buf  = f4%gc_c2f_from_buf(3, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, half_bx(2)
               do i = 1, n_gc
                  i_buf = i_buf + 1
                  uu(bx(1)+i, offset*half_bx(2)+j, ivar, iq) = &
                       f4%recv_buffer(i_buf)
               end do
            end do
         end do
      end do

      face = 2
      do n = f4%gc_c2f_from_buf_iface(face), f4%gc_c2f_from_buf_iface(face+1)-1
         iq     = f4%gc_c2f_from_buf(1, n) + 1 ! Coarse block
         offset = f4%gc_c2f_from_buf(2, n)     ! Offset
         i_buf  = f4%gc_c2f_from_buf(3, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, n_gc
               do i = 1, half_bx(1)
                  i_buf = i_buf + 1
                  uu(offset*half_bx(1)+i, -n_gc+j, ivar, iq) = &
                       f4%recv_buffer(i_buf)
               end do
            end do
         end do
      end do

      face = 3
      do n = f4%gc_c2f_from_buf_iface(face), f4%gc_c2f_from_buf_iface(face+1)-1
         iq     = f4%gc_c2f_from_buf(1, n) + 1 ! Coarse block
         offset = f4%gc_c2f_from_buf(2, n)     ! Offset
         i_buf  = f4%gc_c2f_from_buf(3, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, n_gc
               do i = 1, half_bx(1)
                  i_buf = i_buf + 1
                  uu(offset*half_bx(1)+i, bx(2)+j, ivar, iq) = &
                       f4%recv_buffer(i_buf)
               end do
            end do
         end do
      end do

      ! Local boundaries at the same refinement level

      face = 0
      do n = f4%gc_srl_local_iface(face), f4%gc_srl_local_iface(face+1)-1
         iq   = f4%gc_srl_local(1, n) + 1
         jq   = f4%gc_srl_local(2, n) + 1

         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, bx(2)
               do i = 1, n_gc
                  uu(-n_gc+i, j, ivar, iq) = uu(bx(1)-n_gc+i, j, ivar, jq)
                  uu(bx(1)+i, j, ivar, jq) = uu(i, j, ivar, iq)
               end do
            end do
         end do
      end do

      face = 1
      do n = f4%gc_srl_local_iface(face), f4%gc_srl_local_iface(face+1)-1
         iq   = f4%gc_srl_local(1, n) + 1
         jq   = f4%gc_srl_local(2, n) + 1

         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, bx(2)
               do i = 1, n_gc
                  uu(bx(1)+i, j, ivar, iq) = uu(i, j, ivar, jq)
                  uu(-n_gc+i, j, ivar, jq) = uu(bx(1)-n_gc+i, j, ivar, iq)
               end do
            end do
         end do
      end do

      face = 2
      do n = f4%gc_srl_local_iface(face), f4%gc_srl_local_iface(face+1)-1
         iq   = f4%gc_srl_local(1, n) + 1
         jq   = f4%gc_srl_local(2, n) + 1
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, n_gc
               do i = 1, bx(1)
                  uu(i, -n_gc+j, ivar, iq) = uu(i, bx(2)-n_gc+j, ivar, jq)
                  uu(i, bx(2)+j, ivar, jq) = uu(i, j, ivar, iq)
               end do
            end do
         end do
      end do

      face = 3
      do n = f4%gc_srl_local_iface(face), f4%gc_srl_local_iface(face+1)-1
         iq   = f4%gc_srl_local(1, n) + 1
         jq   = f4%gc_srl_local(2, n) + 1
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, n_gc
               do i = 1, bx(1)
                  uu(i, bx(2)+j, ivar, iq) = uu(i, j, ivar, jq)
                  uu(i, -n_gc+j, ivar, jq) = uu(i, bx(2)-n_gc+j, ivar, iq)
               end do
            end do
         end do
      end do

      ! Physical boundaries
      face = 0
      do n = f4%gc_phys_iface(face), f4%gc_phys_iface(face+1)-1
         iq = f4%gc_phys(n) + 1
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, bx(2)
               slope = uu(2, j, ivar, iq) - uu(1, j, ivar, iq)
               do i = 1, n_gc
                  uu(1-i, j, ivar, iq) = uu(1, j, ivar, iq) - i * slope
               end do
            end do
         end do
      end do

      face = 1
      do n = f4%gc_phys_iface(face), f4%gc_phys_iface(face+1)-1
         iq = f4%gc_phys(n) + 1
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, bx(2)
               slope = uu(bx(1), j, ivar, iq) - uu(bx(1)-1, j, ivar, iq)
               do i = 1, n_gc
                  uu(bx(1)+i, j, ivar, iq) = uu(bx(1), j, ivar, iq) + i * slope
               end do
            end do
         end do
      end do

      face = 2
      do n = f4%gc_phys_iface(face), f4%gc_phys_iface(face+1)-1
         iq = f4%gc_phys(n) + 1
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, n_gc
               do i = 1, bx(1)
                  slope = uu(i, 2, ivar, iq) - uu(i, 1, ivar, iq)
                  uu(i, 1-j, ivar, iq) = uu(i, 1, ivar, iq) - j * slope
               end do
            end do
         end do
      end do

      face = 3
      do n = f4%gc_phys_iface(face), f4%gc_phys_iface(face+1)-1
         iq = f4%gc_phys(n) + 1
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, n_gc
               do i = 1, bx(1)
                  slope = uu(i, bx(2), ivar, iq) - uu(i, bx(2)-1, ivar, iq)
                  uu(i, bx(2)+j, ivar, iq) = uu(i, bx(2), ivar, iq) + j * slope
               end do
            end do
         end do
      end do

      ! Do coarse-to-fine refinement boundaries last, so that ghost cells
      ! required for interpolation have been filled
      call f4_fill_ghostcell_buffers_c2f(f4, n_vars, i_vars)
      call f4_exchange_buffers(f4)

      ! Local coarse-to-fine boundaries, fill fine side
      face = 0
      do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
         iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
         jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
         offset = f4%gc_f2c_local(3, n)     ! Offset

         do iv = 1, n_vars
            ivar = i_vars(iv)

            do j = 1, half_bx(2)
               j_f = 2 * j - 1
               j_c = j + offset * half_bx(2)
               do i = 1, half_n_gc
                  i_c = bx(1) - half_n_gc + i
                  i_f = -(2 * half_n_gc) + 2*i - 1
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, jq), &
                       [f4%uu(i_c-1, j_c, iv, jq), f4%uu(i_c+1, j_c, iv, jq)], &
                       [f4%uu(i_c, j_c-1, iv, jq), f4%uu(i_c, j_c+1, iv, jq)], fine)
                  f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, iq) = reshape(fine, [2, 2])
               end do

               if (odd_n_gc) then
                  i_c = bx(1) - half_n_gc
                  i_f = -n_gc + 1
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, jq), &
                       [f4%uu(i_c-1, j_c, iv, jq), f4%uu(i_c+1, j_c, iv, jq)], &
                       [f4%uu(i_c, j_c-1, iv, jq), f4%uu(i_c, j_c+1, iv, jq)], fine)
                  f4%uu(i_f, j_f:j_f+1, ivar, iq) = fine([2, 4])
               end if
            end do
         end do
      end do

      face = 1
      do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
         iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
         jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
         offset = f4%gc_f2c_local(3, n)     ! Offset

         do iv = 1, n_vars
            ivar = i_vars(iv)

            do j = 1, half_bx(2)
               j_f = 2 * j - 1
               j_c = j + offset * half_bx(2)
               do i = 1, half_n_gc
                  i_c = i
                  i_f = bx(1) + 2*i - 1
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, jq), &
                       [f4%uu(i_c-1, j_c, iv, jq), f4%uu(i_c+1, j_c, iv, jq)], &
                       [f4%uu(i_c, j_c-1, iv, jq), f4%uu(i_c, j_c+1, iv, jq)], fine)
                  f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, iq) = reshape(fine, [2, 2])
               end do

               if (odd_n_gc) then
                  i_c = 1 + half_n_gc
                  i_f = bx(1) + n_gc
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, jq), &
                       [f4%uu(i_c-1, j_c, iv, jq), f4%uu(i_c+1, j_c, iv, jq)], &
                       [f4%uu(i_c, j_c-1, iv, jq), f4%uu(i_c, j_c+1, iv, jq)], fine)
                  f4%uu(i_f, j_f:j_f+1, ivar, iq) = fine([1, 3])
               end if
            end do
         end do
      end do

      face = 2
      do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
         iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
         jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
         offset = f4%gc_f2c_local(3, n)     ! Offset

         do iv = 1, n_vars
            ivar = i_vars(iv)

            do j = 1, half_n_gc
               j_c = bx(2) - half_n_gc + j
               j_f = -(2 * half_n_gc) + 2*j - 1
               do i = 1, half_bx(1)
                  i_f = 2 * i - 1
                  i_c = i + offset * half_bx(1)
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, jq), &
                       [f4%uu(i_c-1, j_c, iv, jq), f4%uu(i_c+1, j_c, iv, jq)], &
                       [f4%uu(i_c, j_c-1, iv, jq), f4%uu(i_c, j_c+1, iv, jq)], fine)
                  f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, iq) = reshape(fine, [2, 2])
               end do
            end do

            if (odd_n_gc) then
               j_c = bx(2) - half_n_gc
               j_f = -n_gc + 1
               do i = 1, half_bx(1)
                  i_f = 2 * i - 1
                  i_c = i + offset * half_bx(1)
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, jq), &
                       [f4%uu(i_c-1, j_c, iv, jq), f4%uu(i_c+1, j_c, iv, jq)], &
                       [f4%uu(i_c, j_c-1, iv, jq), f4%uu(i_c, j_c+1, iv, jq)], fine)
                  f4%uu(i_f:i_f+1, j_f, ivar, iq) = fine([3, 4])
               end do
            end if
         end do
      end do

      face = 3
      do n = f4%gc_f2c_local_iface(face), f4%gc_f2c_local_iface(face+1)-1
         iq     = f4%gc_f2c_local(1, n) + 1 ! Fine block
         jq     = f4%gc_f2c_local(2, n) + 1 ! coarse block
         offset = f4%gc_f2c_local(3, n)     ! Offset

         do iv = 1, n_vars
            ivar = i_vars(iv)

            do j = 1, half_n_gc
               j_c = j
               j_f = bx(2) + 2*j - 1
               do i = 1, half_bx(1)
                  i_f = 2 * i - 1
                  i_c = i + offset * half_bx(1)
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, jq), &
                       [f4%uu(i_c-1, j_c, iv, jq), f4%uu(i_c+1, j_c, iv, jq)], &
                       [f4%uu(i_c, j_c-1, iv, jq), f4%uu(i_c, j_c+1, iv, jq)], fine)
                  f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, iq) = reshape(fine, [2, 2])
               end do
            end do

            if (odd_n_gc) then
               j_c = 1 + half_n_gc
               j_f = bx(2) + n_gc
               do i = 1, half_bx(1)
                  i_f = 2 * i - 1
                  i_c = i + offset * half_bx(1)
                  call prolong_local_5point(f4%uu(i_c, j_c, iv, jq), &
                       [f4%uu(i_c-1, j_c, iv, jq), f4%uu(i_c+1, j_c, iv, jq)], &
                       [f4%uu(i_c, j_c-1, iv, jq), f4%uu(i_c, j_c+1, iv, jq)], fine)
                  f4%uu(i_f:i_f+1, j_f, ivar, iq) = fine([1, 2])
               end do
            end if
         end do

      end do

      ! Nonlocal coarse-to-fine boundaries, fill fine side

      face = 0
      do n = f4%gc_f2c_from_buf_iface(face), f4%gc_f2c_from_buf_iface(face+1)-1
         iq    = f4%gc_f2c_from_buf(1, n) + 1 ! Fine block
         i_buf = f4%gc_f2c_from_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, half_bx(2)
               j_f = 2 * j - 1
               do i = 1, half_n_gc
                  i_f = -(2 * half_n_gc) + 2*i - 1
                  f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, iq) = &
                       reshape(f4%recv_buffer(i_buf+1:i_buf+4), [2, 2])
                  i_buf = i_buf + 4
               end do

               if (odd_n_gc) then
                  i_f = -n_gc + 1
                  f4%uu(i_f, j_f:j_f+1, ivar, iq) = &
                       f4%recv_buffer(i_buf+1:i_buf+2)
                  i_buf = i_buf + 2
               end if
            end do
         end do
      end do

      face = 1
      do n = f4%gc_f2c_from_buf_iface(face), f4%gc_f2c_from_buf_iface(face+1)-1
         iq    = f4%gc_f2c_from_buf(1, n) + 1 ! Fine block
         i_buf = f4%gc_f2c_from_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, half_bx(2)
               j_f = 2 * j - 1
               do i = 1, half_n_gc
                  i_f = bx(1) + 2*i - 1
                  f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, iq) = &
                       reshape(f4%recv_buffer(i_buf+1:i_buf+4), [2, 2])
                  i_buf = i_buf + 4
               end do

               if (odd_n_gc) then
                  i_f = bx(1) + n_gc
                  f4%uu(i_f, j_f:j_f+1, ivar, iq) = &
                       f4%recv_buffer(i_buf+1:i_buf+2)
                  i_buf = i_buf + 2
               end if
            end do
         end do
      end do

      face = 2
      do n = f4%gc_f2c_from_buf_iface(face), f4%gc_f2c_from_buf_iface(face+1)-1
         iq    = f4%gc_f2c_from_buf(1, n) + 1 ! Fine block
         i_buf = f4%gc_f2c_from_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, half_n_gc
               j_f = -(2 * half_n_gc) + 2*j - 1
               do i = 1, half_bx(1)
                  i_f = 2 * i - 1
                  f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, iq) = &
                       reshape(f4%recv_buffer(i_buf+1:i_buf+4), [2, 2])
                  i_buf = i_buf + 4
               end do
            end do

            if (odd_n_gc) then
               j_f = -n_gc + 1
               do i = 1, half_bx(1)
                  i_f = 2 * i - 1
                  f4%uu(i_f:i_f+1, j_f, ivar, iq) = &
                       f4%recv_buffer(i_buf+1:i_buf+2)
                  i_buf = i_buf + 2
               end do
            end if
         end do
      end do

      face = 3
      do n = f4%gc_f2c_from_buf_iface(face), f4%gc_f2c_from_buf_iface(face+1)-1
         iq    = f4%gc_f2c_from_buf(1, n) + 1 ! Fine block
         i_buf = f4%gc_f2c_from_buf(2, n) * n_vars
         do iv = 1, n_vars
            ivar = i_vars(iv)
            do j = 1, half_n_gc
               j_f = bx(2) + 2*j - 1
               do i = 1, half_bx(1)
                  i_f = 2 * i - 1
                  f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, iq) = &
                       reshape(f4%recv_buffer(i_buf+1:i_buf+4), [2, 2])
                  i_buf = i_buf + 4
               end do
            end do

            if (odd_n_gc) then
               j_f = bx(2) + n_gc
               do i = 1, half_bx(1)
                  i_f = 2 * i - 1
                  f4%uu(i_f:i_f+1, j_f, ivar, iq) = &
                       f4%recv_buffer(i_buf+1:i_buf+2)
                  i_buf = i_buf + 2
               end do
            end if
         end do
      end do
    end associate

  end subroutine f4_update_ghostcells

  !> Refine the mesh according to f4%refinement_flags
  subroutine f4_adjust_refinement(f4, partition)
    type(foap4_t), intent(inout) :: f4
    logical, intent(in)          :: partition
    integer                      :: n, n_blocks_new, n_blocks_old, k, iv
    integer                      :: has_changed

    n_blocks_old = f4%n_blocks
    call pw_adjust_refinement(f4%pw, f4%n_blocks, &
         f4%refinement_flags(1:f4%n_blocks), has_changed)

    if (has_changed == 0) return

    n_blocks_new = pw_get_num_local_quadrants(f4%pw)

    if (n_blocks_new + f4%n_blocks > f4%max_blocks) &
         error stop "Not enough block memory during refinement"

    call copy_blocks_to_end(f4, n_blocks_old, n_blocks_new)

    call f4_set_quadrants(f4)

    k = n_blocks_new + 1
    n = 1
    do while (n <= n_blocks_new)
       select case (f4%block_level(n) - f4%block_level(k))
       case (0)
          ! Same refinement level
          call copy_block(f4, k, n)
          n = n + 1
          k = k + 1
       case (1)
          ! Block has been refined
          do iv = 1, f4%n_vars
             call prolong_to_blocks(f4, k, n, iv)
          end do
          n = n + 4
          k = k + 1
       case (-1)
          ! Block has been coarsened
          do iv = 1, f4%n_vars
             call coarsen_from_blocks(f4, k, n, iv)
          end do
          n = n + 1
          k = k + 4
       case default
          error stop "Refinement: difference in levels > 1"
       end select
    end do

    if (k /= n_blocks_new + n_blocks_old + 1) &
         error stop "Refinement: loops did not end simultaneously"

    if (partition) call f4_partition(f4)
  end subroutine f4_adjust_refinement

  ! TODO: openacc
  subroutine copy_block(f4, i_from, i_to)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: i_from, i_to
    f4%block_origin(:, i_to) = f4%block_origin(:, i_from)
    f4%block_level(i_to)     = f4%block_level(i_from)
    f4%uu(:, :, :, i_to)     = f4%uu(:, :, :, i_from)
  end subroutine copy_block

  subroutine copy_blocks_to_end(f4, n_blocks_old, n_blocks_new)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_blocks_old
    integer, intent(in)          :: n_blocks_new
    integer                      :: n

    ! Copy blocks to end of list. The backward order ensures old data is not
    ! overwritten before it is copied.
    do n = n_blocks_old, 1, -1
       call copy_block(f4, n, n_blocks_new+n)
    end do
  end subroutine copy_blocks_to_end

  !> Coarsen a family of child blocks to their parent
  ! TODO: openacc
  subroutine coarsen_from_blocks(f4, i_from, i_to, iv)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: i_from
    integer, intent(in)          :: i_to
    integer, intent(in)          :: iv
    integer                      :: n, i, j, i_c, j_c, i_f, j_f
    integer                      :: half_bx(2), ix_offset(2)

    half_bx = f4%bx / 2

    do n = 1, 4
       ix_offset = child_offset(:, n) * half_bx

       do j = 1, half_bx(2)
          j_c = j + ix_offset(2)
          j_f = 2 * j - 1
          do i = 1, half_bx(1)
             i_c = i + ix_offset(1)
             i_f = 2 * i - 1

             f4%uu(i_c, j_c, iv, i_to) = 0.25_dp * (&
                  f4%uu(i_f,   j_f, iv, i_from+n-1) + &
                  f4%uu(i_f+1, j_f, iv, i_from+n-1) + &
                  f4%uu(i_f,   j_f+1, iv, i_from+n-1) + &
                  f4%uu(i_f+1, j_f+1, iv, i_from+n-1))
          end do
       end do
    end do
  end subroutine coarsen_from_blocks

  !> Prolong to family of child blocks
  ! TODO: openacc
  subroutine prolong_to_blocks(f4, i_from, i_to, iv)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: i_from
    integer, intent(in)          :: i_to
    integer, intent(in)          :: iv
    integer                      :: n, i, j, i_c, j_c, i_f, j_f
    integer                      :: half_bx(2), ix_offset(2)
    real(dp)                     :: fine(4)

    half_bx = f4%bx / 2

    do n = 1, 4
       ix_offset = child_offset(:, n) * half_bx

       do j = 1, half_bx(2)
          j_c = j + ix_offset(2)
          j_f = 2 * j - 1
          do i = 1, half_bx(1)
             i_c = i + ix_offset(1)
             i_f = 2 * i - 1

             call prolong_local_5point(f4%uu(i_c, j_c, iv, i_from), &
                  [f4%uu(i_c-1, j_c, iv, i_from), f4%uu(i_c+1, j_c, iv, i_from)], &
                  [f4%uu(i_c, j_c-1, iv, i_from), f4%uu(i_c, j_c+1, iv, i_from)], &
                  fine)

             f4%uu(i_f,   j_f, iv, i_to+n-1)   = fine(1)
             f4%uu(i_f+1, j_f, iv, i_to+n-1)   = fine(2)
             f4%uu(i_f,   j_f+1, iv, i_to+n-1) = fine(3)
             f4%uu(i_f+1, j_f+1, iv, i_to+n-1) = fine(4)
          end do
       end do
    end do
  end subroutine prolong_to_blocks

  !> Method for prolongation (interpolation) of a coarse block to its children
  subroutine prolong_local_5point(coarse_c, coarse_x, coarse_y, fine)
    real(dp), intent(in)  :: coarse_c    ! Center value
    real(dp), intent(in)  :: coarse_x(2) ! x-neighbors (-1, +1)
    real(dp), intent(in)  :: coarse_y(2) ! y-neighbors (-1, +1)
    real(dp), intent(out) :: fine(4)
    real(dp)              :: f(0:2), slopes_a(2), slopes_b(2)

    f(0) = coarse_c          ! Identical to coarse_y(2)
    slopes_a = [coarse_c - coarse_x(1), coarse_c - coarse_y(1)]
    slopes_b = [coarse_x(2) - coarse_c, coarse_y(2) - coarse_c]
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
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp)             :: phi
    phi = limiter_gminmod(a, b, 1.0_dp)
  end function limiter_minmod

  !> Partition the blocks over the MPI ranks
  subroutine f4_partition(f4)
    type(foap4_t), intent(inout), target :: f4
    integer                              :: n_changed_global
    integer                              :: n_blocks_old, n_blocks_new
    integer(c_int64_t)                   :: gfq_old(0:f4%mpisize)
    integer                              :: dsize

    n_blocks_old = f4%n_blocks
    call pw_partition(f4%pw, n_changed_global, gfq_old)

    ! No need to do anything if the global number of blocks shipped is zero
    if (n_changed_global == 0) return

    n_blocks_new = pw_get_num_local_quadrants(f4%pw)

    if (n_blocks_new + n_blocks_old > f4%max_blocks) &
         error stop "Not enough block memory during partitioning"

    call copy_blocks_to_end(f4, n_blocks_old, n_blocks_new)

    ! Size of a block
    dsize = product(f4%bx + 2 * f4%n_gc) * f4%n_vars * storage_size(1.0_dp)/8

    ! Call p4est routine to transfer the block data
    call pw_partition_transfer(f4%pw, gfq_old, &
         c_loc(f4%uu(:, :, :, n_blocks_new+1)), c_loc(f4%uu), dsize)

    call f4_set_quadrants(f4)
  end subroutine f4_partition

end module m_foap4
