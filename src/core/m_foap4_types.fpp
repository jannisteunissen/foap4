#:include 'definitions.fpp'
module m_foap4_types_${NDIM}$d
  use, intrinsic :: iso_c_binding
  use mpi_f08

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  integer, parameter, public :: ndim = ${NDIM}$

  !> Maximum refinement level in p4est
  integer, parameter, public :: P4EST_MAXLEVEL = 30

  !> The opposite of the faces 0-3
#:if NDIM == 2
  integer, parameter, public :: face_swap(0:2*ndim-1) = [1, 0, 3, 2]
#:elif NDIM == 3
  integer, parameter, public :: face_swap(0:2*ndim-1) = [1, 0, 3, 2, 5, 4]
#:endif

  !> The offset of the children
#:if NDIM == 2
  integer, parameter, public :: f4_child_offset(ndim, 4) = reshape(&
       [0,0, 1,0, 0,1, 1,1], [2,4])
#:elif NDIM == 3
  integer, parameter, public :: f4_child_offset(ndim, 8) = reshape( &
       [0,0,0, 1,0,0, 0,1,0, 1,1,0, &
       0,0,1, 1,0,1, 0,1,1, 1,1,1], [3,8])
#:endif
  !$acc declare create(f4_child_offset)

  !> The dimension of faces
#:if NDIM == 2
  integer, parameter, public :: f4_face_dim(0:2*ndim-1) = [1, 1, 2, 2]
#:elif NDIM == 3
  integer, parameter, public :: f4_face_dim(0:2*ndim-1) = [1, 1, 2, 2, 3, 3]
#:endif
  !$acc declare create(f4_face_dim)

  !> Index of the face pointing to -x
  integer, parameter, public :: f4_face_xlo = 0

  !> Index of the face pointing to +x
  integer, parameter, public :: f4_face_xhi = 1

  !> Index of the face pointing to -y
  integer, parameter, public :: f4_face_ylo = 2

  !> Index of the face pointing to +y
  integer, parameter, public :: f4_face_yhi = 3

#:if NDIM == 3
  !> Index of the face pointing to -z
  integer, parameter, public :: f4_face_zlo = 4

  !> Index of the face pointing to +z
  integer, parameter, public :: f4_face_zhi = 5
#:endif

  !> Value indicating a physical boundary at a block face
  integer, parameter, public :: FACE_BOUNDARY = 0

  !> Value indicating a neighbor at the same refinement level at a block face
  integer, parameter, public :: FACE_SAME_LEVEL = 1

  !> Value indicating a neighbor at a higher refinement level at a block face
  integer, parameter, public :: FACE_COARSE_TO_FINE = 2

  !> Value indicating a neighbor at a lower refinement level at a block face
  integer, parameter, public :: FACE_FINE_TO_COARSE = 3

  !> Dirichlet boundary condition
  integer, parameter, public :: f4_bc_dirichlet = 0

  !> Neumann boundary condition
  integer, parameter, public :: f4_bc_neumann = 1

  !> extrapolation boundary condition
  integer, parameter, public :: f4_bc_linear_extrap = 2

  !> values will be set in ghost cells directly
  integer, parameter, public :: f4_bc_fixed_value = 3

  !> Type to store an array of integers
  type, public :: int_array_t
     integer, allocatable :: i(:)
  end type int_array_t

  !> Type to describe a face boundary. The same data structure is defined in
  !> p4est_wrapper.c
  type, bind(c), public :: bnd_face_t
     integer(c_int) :: face_type !< What kind of face (same level, ...)
     integer(c_int) :: face      !< Direction of the face
     integer(c_int) :: other_proc !< MPI rank that owns quadid(2)
     integer(c_int) :: quadid(2)  !< quadid(1) is always local, (2) can be non-local
     integer(c_int) :: offset(NDIM-1) !< Offset for a hanging face
     integer(c_int) :: ibuf_recv  !< Index in receive buffer (not filled here)
     integer(c_int) :: ibuf_send  !< Index in send buffer (not filled here)
  end type bnd_face_t

  !> Data structure that contains the whole AMR grid and required information
  !> for ghost cell communication
  type, public :: foap4_t
     integer   :: bx(ndim)           !< Block size (cells)
     integer   :: hbx(ndim)          !< Half of block size (cells)
     integer   :: n_gc               !< Number of ghost cells
     integer   :: max_blocks         !< Maximum number of blocks used
     integer   :: n_vars_nontemporal !< Number of non-temporal variables
     integer   :: n_vars_temporal    !< Number of temporal variables
     integer   :: n_vars             !< Total number of variables (excl. temporal copies)
     integer   :: n_vars_all         !< Total number of variables (incl. temporal copies)
     integer   :: n_temporal_states  !< Number of temporal states
     real(dp)  :: tree_length(ndim)                  !< Length of tree
     !> Coefficient for generalized minmod limiter for prolongation
     real(dp)  :: gminmod_theta_prolong = 1.0_dp
     integer   :: ilo(ndim)                          !< Minimum index in a block
     integer   :: ihi(ndim)                          !< Maximum index in a block
     real(dp)  :: dr_level(ndim, 0:p4est_maxlevel-1) !< Grid spacing per level
     character(len=32), allocatable :: var_names(:) !< Names of the variables

     integer :: n_blocks              !< Number of blocks used
     integer :: gc_mesh_revision = -1 !< Revision number (global) of the mesh

     !> Array storing simple physical boundary conditions bc_simple_type(ivar, iface)
     integer, allocatable :: bc_simple_type(:, :)
     !> Array storing value of boundary conditions bc_simple(ivar, iface)
     real(dp), allocatable :: bc_simple(:, :)

     !> Level of each block
     integer, allocatable  :: block_level(:)
     !> Origin of each block
     real(dp), allocatable :: block_origin(:, :)
     !> Refinement flag of each block. Negative means coarsen (if possible),
     !> positive means refine, and zero means keep refinement.
     integer, allocatable  :: refinement_flags(:)
     !> Index into bc_data array (only used if > 0), shaped bc_data_ix(face, i_block)
     integer, allocatable  :: bc_data_ix(:, :)
     !> Index into bflux array (only used if > 0), shaped bflux_ix(face, i_block)
     integer, allocatable  :: bflux_ix(:, :)

     !> Storage of block data uu(i, j, [k,] i_var, i_block)
     real(dp), allocatable :: uu(@{DTIMES(:)}@, :, :)
     !> Storage of boundary flux * dt, as bflux(i, [j,] i_var, i_face)
     real(dp), allocatable :: bflux(@{DTIMES(:)}@, :)
     !> Storage of boundary condition data bc_data(i, [j,] i_var, i_face)
     real(dp), allocatable :: bc_data(@{DTIMES(:)}@, :)
     !> Storage of boundary condition type bc_data_type(i, [j,] i_var, i_face)
     integer, allocatable :: bc_data_type(@{DTIMES(:)}@, :)

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

     !> Number of ghost cells on a SRL face boundary
     integer              :: gc_data_size
     !> Number of ghost cells on a coarse-to-fine face boundary
     integer              :: gc_data_size_c2f
     !> Number of cells transmitted for flux fixing
     integer              :: gc_data_size_fluxfix

     !> Receive offset (per MPI rank) in recv_buffer
     integer, allocatable :: gc_recv_offset(:)
     !> Receive offset (per MPI rank) in recv_buffer for the coarse-to-fine step
     integer, allocatable :: gc_recv_offset_c2f(:)
     !> Receive offset (per MPI rank) in recv_buffer for flux fixing
     integer, allocatable :: gc_recv_offset_fluxfix(:)
     !> Send offset (per MPI rank) in send_buffer
     integer, allocatable :: gc_send_offset(:)
     !> Send offset (per MPI rank) in send_buffer for the coarse-to-fine step
     integer, allocatable :: gc_send_offset_c2f(:)
     !> Send offset (per MPI rank) in send_buffer for flux fixing
     integer, allocatable :: gc_send_offset_fluxfix(:)

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

     !> Same as gc_f2c_to_buf but for flux fixing (different buffer offsets)
     integer, allocatable :: gc_f2c_to_buf_fluxfix(:, :)
     !> Same as gc_c2f_from_buf but for flux fixing (different buffer offsets)
     integer, allocatable :: gc_c2f_from_buf_fluxfix(:, :)

     ! Since the above arrays are sorted by face direction, we can loop over
     ! them per face direction. The arrays below define the start index for each
     ! direction (with the last value being the total number of elements).

     integer :: gc_srl_local_iface(0:2*ndim)
     integer :: gc_srl_from_buf_iface(0:2*ndim)
     integer :: gc_srl_to_buf_iface(0:2*ndim)
     integer :: gc_f2c_local_iface(0:2*ndim)
     integer :: gc_c2f_from_buf_iface(0:2*ndim)
     integer :: gc_c2f_to_buf_iface(0:2*ndim)
     integer :: gc_f2c_from_buf_iface(0:2*ndim)
     integer :: gc_f2c_to_buf_iface(0:2*ndim)
     integer :: gc_phys_iface(0:2*ndim)

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
     real(dp) :: wtime_flux_fix = 0.0_dp

     !> Optional procedure to set boundary conditions after changing the mesh.
     !> Should be set before the mesh is created.
     procedure(bc_callback_t), pointer, nopass :: bc_callback => null()
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
     subroutine pw_set_connectivity_brick(pw, mi, periodic, &
          min_level, fill_uniform, max_blocks) bind(c)
       import c_int, c_ptr, ndim
       type(c_ptr), intent(in), value    :: pw
       integer(c_int), intent(in) :: mi(ndim)
       integer(c_int), intent(in) :: periodic(ndim)
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
  end interface

  abstract interface
     !> Get cell-centered data for a grid block
     subroutine subr_cc_data_2D(n, cc_data)
       import
       integer, intent(in)     :: n !< Index of block
       real(dp), intent(inout) :: cc_data(:, :, :)
     end subroutine subr_cc_data_2D

     !> Get cell-centered data for a grid block
     subroutine subr_cc_data_3D(n, cc_data)
       import
       integer, intent(in)     :: n !< Index of block
       real(dp), intent(inout) :: cc_data(:, :, :, :)
     end subroutine subr_cc_data_3D

     subroutine bc_callback_t(f4)
       import
       type(foap4_t), intent(inout) :: f4
     end subroutine bc_callback_t
  end interface

  ! Public interfaces for subroutines
  public :: pw_initialize
  public :: pw_destroy
  public :: pw_finalize
  public :: pw_set_connectivity_brick
  public :: pw_get_quadrants
  public :: pw_vtk_write_file
  public :: pw_get_all_faces
  public :: pw_adjust_refinement
  public :: pw_partition
  public :: subr_cc_data_2D
  public :: subr_cc_data_3D
  public :: bc_callback_t

  ! Public interfaces for functions
  public :: pw_get_num_local_quadrants
  public :: pw_get_num_global_quadrants
  public :: pw_get_mesh_revision
  public :: pw_get_highest_local_level

end module m_foap4_types_${NDIM}$d
