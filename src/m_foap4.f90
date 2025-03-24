! Foap4 stands for "Fortran OpenACC p4est". Together with p4est_wrapper.c,
! this module implements the required data structures and methods.
!
! Author(s): Jannis Teunissen
module m_foap4
  use, intrinsic :: iso_c_binding
  use mpi_f08

  implicit none
  private

  integer, parameter, private :: dp = kind(0.0d0)

  integer, parameter :: P4EST_MAXLEVEL = 30
  integer, parameter :: face_swap(0:3) = [1, 0, 3, 2]
  integer, parameter :: child_offset(2, 4) = reshape([0,0,1,0,0,1,1,1], [2,4])

  integer, parameter :: FACE_BOUNDARY = 0
  integer, parameter :: FACE_SAME_LEVEL = 1
  integer, parameter :: FACE_COARSE_TO_FINE = 2
  integer, parameter :: FACE_FINE_TO_COARSE = 3

  type int_array_t
     integer, allocatable :: i(:)
  end type int_array_t

  type, bind(c) :: bnd_face_t
     integer(c_int) :: face_type
     integer(c_int) :: face
     integer(c_int) :: other_proc
     integer(c_int) :: quadid(2)
     integer(c_int) :: extra
  end type bnd_face_t

  type face_gc_t
     integer              :: mesh_revision = -1
     integer              :: data_size
     integer, allocatable :: recv_offset(:)
     integer, allocatable :: recv_offset_c2f(:)
     integer, allocatable :: send_offset(:)
     integer, allocatable :: send_offset_c2f(:)
     integer, allocatable :: same_local(:, :)
     integer, allocatable :: same_from_buf(:, :)
     integer, allocatable :: same_to_buf(:, :)
     integer, allocatable :: c2f_local(:, :)
     integer, allocatable :: c2f_from_buf(:, :)
     integer, allocatable :: c2f_to_buf(:, :)
     integer, allocatable :: f2c_from_buf(:, :)
     integer, allocatable :: f2c_to_buf(:, :)
     integer, allocatable :: phys(:, :)
  end type face_gc_t

  type, public :: foap4_t
     integer   :: bx(2)                         ! Block size (cells)
     integer   :: n_gc                          ! Number of ghost cells
     integer   :: n_blocks                      ! Number of blocks used
     integer   :: max_blocks                    ! Maximum number of blocks used
     integer   :: n_vars                        ! Number of variables
     real(dp)  :: tree_length(2)                ! Length of tree
     real(dp)  :: dr_lvl(2, 0:p4est_maxlevel-1) ! Grid spacing per level
     character(len=32), allocatable :: var_names(:) ! Names of the variables

     ! The data per block
     integer, allocatable  :: block_level(:)      ! Level of each block
     real(dp), allocatable :: block_origin(:, :)  ! Origin of each block
     real(dp), allocatable :: uu(:, :, :, :)      ! Block storage
     integer, allocatable  :: refinement_flags(:) ! Refinement flags

     ! For communication
     integer, allocatable :: recv_offset(:) ! 0:mpisize offsets for receiving
     integer, allocatable :: send_offset(:) ! 0:mpisize offsets for sending
     real(dp), allocatable :: recv_buffer(:)
     real(dp), allocatable :: send_buffer(:)

     ! Information about ghost cells on faces
     type(face_gc_t) :: face_gc

     ! p4est state
     type(c_ptr) :: pw
  end type foap4_t

  interface
     subroutine pw_initialize(pw, mpicomm, max_blocks) bind(c)
       import c_int, c_ptr
       type(c_ptr), intent(out)          :: pw
       integer(c_int), intent(out)       :: mpicomm
       integer(c_int), intent(in), value :: max_blocks
     end subroutine pw_initialize

     subroutine pw_set_connectivity_brick(pw, mi, ni, periodic_a, periodic_b, &
          min_level, fill_uniform) bind(c)
       import c_int, c_ptr
       type(c_ptr), intent(in), value    :: pw
       integer(c_int), value, intent(in) :: mi
       integer(c_int), value, intent(in) :: ni
       integer(c_int), value, intent(in) :: periodic_a
       integer(c_int), value, intent(in) :: periodic_b
       integer(c_int), value, intent(in) :: min_level
       integer(c_int), value, intent(in) :: fill_uniform
     end subroutine pw_set_connectivity_brick

     pure function pw_get_num_local_quadrants(pw) result(n) bind(c)
       import c_int, c_ptr
       type(c_ptr), intent(in), value :: pw
       integer(c_int)                 :: n
     end function pw_get_num_local_quadrants

     pure function pw_get_mesh_revision(pw) result(n) bind(c)
       import c_int, c_ptr
       type(c_ptr), intent(in), value :: pw
       integer(c_int)                 :: n
     end function pw_get_mesh_revision

     subroutine pw_get_quadrants(pw, n_quadrants, coord, level) bind(c)
       import c_int, c_ptr, c_double
       type(c_ptr), intent(in), value     :: pw
       integer(c_int), value, intent(in)  :: n_quadrants
       real(kind=c_double), intent(inout) :: coord(*)
       integer(c_int), intent(inout)      :: level(n_quadrants)
     end subroutine pw_get_quadrants

     subroutine pw_vtk_write_file(pw, fname) bind(c)
       import c_char, c_ptr
       type(c_ptr), intent(in), value     :: pw
       character(kind=C_char), intent(in) :: fname(*)
     end subroutine pw_vtk_write_file

     subroutine pw_finalize_mpi_and_p4est(pw) bind(c)
       import c_ptr
       type(c_ptr), intent(in), value :: pw
     end subroutine pw_finalize_mpi_and_p4est

     subroutine pw_get_all_faces(pw, n_faces, bnd_face_ptr) bind(c)
       import c_int, c_ptr
       type(c_ptr), intent(in), value :: pw
       integer(c_int), intent(out)    :: n_faces
       type(c_ptr), intent(out)       :: bnd_face_ptr
     end subroutine pw_get_all_faces

     subroutine pw_adjust_refinement(pw, n_quadrants, flags, has_changed) bind(c)
       import c_int, c_ptr
       type(c_ptr), intent(in), value    :: pw
       integer(c_int), value, intent(in) :: n_quadrants
       integer(c_int), intent(in)        :: flags(n_quadrants)
       integer(c_int), intent(out)       :: has_changed
     end subroutine pw_adjust_refinement
  end interface

  type(MPI_comm), protected  :: mpicomm
  integer, public, protected :: mpirank, mpisize

  public :: f4_initialize_grid
  public :: f4_finalize_grid
  public :: f4_write_grid
  public :: f4_get_num_local_blocks
  public :: f4_cell_coord
  public :: f4_exchange_buffers
  public :: f4_update_ghostcells
  public :: f4_adjust_refinement

contains

  subroutine f4_initialize_grid(n_blocks_per_dim, tree_length, bx, n_gc, &
       n_vars, var_names, periodic, min_level, max_blocks, f4)
    integer, intent(in)          :: n_blocks_per_dim(2)
    real(dp), intent(in)         :: tree_length(2)
    integer, intent(in)          :: bx(2)
    integer, intent(in)          :: n_gc
    integer, intent(in)          :: n_vars
    character(len=*), intent(in) :: var_names(n_vars)
    logical, intent(in)          :: periodic(2)
    integer, intent(in)          :: min_level
    integer, intent(in)          :: max_blocks
    type(foap4_t), intent(out)   :: f4

    integer :: i, ierr, periodic_as_int(2)

    if (bx(1) /= bx(2)) error stop "Unequal bx(:) not supported"
    if (any(bx < 2 * n_gc)) error stop "Cannot have any(bx < 2 * n_gc)"
    if (any(iand(bx, 1) == 1)) error stop "All bx have to be even"

    f4%bx         = bx
    f4%n_gc       = n_gc
    f4%n_vars     = n_vars
    f4%max_blocks = max_blocks
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
       f4%dr_lvl(:, i) = (tree_length/bx) * 0.5**i
    end do

    call pw_initialize(f4%pw, mpicomm%MPI_VAL, max_blocks)

    call MPI_COMM_RANK(mpicomm, mpirank, ierr)
    call MPI_COMM_SIZE(mpicomm, mpisize, ierr)

    call pw_set_connectivity_brick(f4%pw, &
         n_blocks_per_dim(1), n_blocks_per_dim(2), &
         periodic_as_int(1), periodic_as_int(2), min_level, 1)

    allocate(f4%block_origin(2, max_blocks))
    allocate(f4%block_level(max_blocks))
    allocate(f4%refinement_flags(max_blocks))
    allocate(f4%uu(1-n_gc:bx(1)+n_gc, 1-n_gc:bx(2)+n_gc, n_vars, max_blocks))
    f4%uu(:, :, :, :) = 0.0_dp

    f4%face_gc%data_size = f4%bx(1) * f4%n_gc

    ! Maximum size of recv/send buffer
    i = max_blocks * 4 * f4%face_gc%data_size
    allocate(f4%recv_buffer(i))
    allocate(f4%send_buffer(i))
    allocate(f4%recv_offset(0:mpisize))
    allocate(f4%send_offset(0:mpisize))

    f4%uu = 0.0_dp

    call f4_get_quadrants(f4)

  end subroutine f4_initialize_grid

  subroutine f4_finalize_grid(f4)
    type(foap4_t), intent(inout) :: f4

    ! TODO: deallocate

    call pw_finalize_mpi_and_p4est(f4%pw)
  end subroutine f4_finalize_grid

  pure integer function f4_get_num_local_blocks(f4)
    type(foap4_t), intent(in) :: f4
    f4_get_num_local_blocks = pw_get_num_local_quadrants(f4%pw)
  end function f4_get_num_local_blocks

  subroutine f4_get_quadrants(f4)
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
  end subroutine f4_get_quadrants

  subroutine f4_write_grid(f4, fname, time, viewer)
    use m_xdmf_writer
    type(foap4_t), intent(in)              :: f4
    character(len=*), intent(in)           :: fname
    real(dp), intent(in), optional         :: time
    character(len=*), intent(in), optional :: viewer
    integer                                :: n
    real(dp), allocatable                  :: dr(:, :)

    call pw_vtk_write_file(f4%pw, trim(fname) // C_null_char)

    allocate(dr(2, f4%n_blocks))

    do n = 1, f4%n_blocks
       dr(:, n) = f4%dr_lvl(:, f4%block_level(n))
    end do

    call xdmf_write_blocks_2DCoRect(mpicomm, trim(fname), &
         f4%n_blocks, f4%bx+2*f4%n_gc, f4%n_vars, &
         f4%var_names, f4%n_gc, f4%block_origin(:, 1:f4%n_blocks), dr, &
         cc_data=f4%uu(:, :, :, 1:f4%n_blocks), time=time, viewer=viewer)
  end subroutine f4_write_grid

  pure function f4_cell_coord(f4, i_block, i, j) result(rr)
    type(foap4_t), intent(in) :: f4
    integer, intent(in)       :: i_block, i, j
    real(dp)                  :: rr(2), dr(2)

    dr = f4%dr_lvl(:, f4%block_level(i_block))
    rr = f4%block_origin(:, i_block) + dr * [i-0.5_dp, j-0.5_dp]
  end function f4_cell_coord

  subroutine f4_update_ghostcell_pattern(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: mesh_revision
    type(bnd_face_t), pointer    :: bnd_face(:)
    type(c_ptr)                  :: tmp
    integer                      :: n_faces

    mesh_revision = pw_get_mesh_revision(f4%pw)
    if (mesh_revision == f4%face_gc%mesh_revision) return

    call pw_get_all_faces(f4%pw, n_faces, tmp)

    call c_f_pointer(tmp, bnd_face, shape=[n_faces])

    call f4_get_ghost_cell_pattern(f4%face_gc, size(bnd_face), bnd_face, &
         mpirank, mpisize)
    f4%face_gc%mesh_revision = mesh_revision
  end subroutine f4_update_ghostcell_pattern

  subroutine f4_get_ghost_cell_pattern(face_gc, n_faces, bnd_face, mpirank, mpisize)
    type(face_gc_t), intent(inout)  :: face_gc
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
    integer :: i_same_nonlocal
    integer :: i_f2c_to_buf, i_f2c_from_buf
    integer :: i_c2f_to_buf, i_c2f_from_buf

    type(int_array_t) :: same_ix(0:mpisize-1)
    type(int_array_t) :: c2f_ix(0:mpisize-1)
    type(int_array_t) :: f2c_ix(0:mpisize-1)
    type(int_array_t) :: ix_send

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
    if (.not. allocated(face_gc%recv_offset)) then
       allocate(face_gc%recv_offset(0:mpisize))
       allocate(face_gc%send_offset(0:mpisize))
       allocate(face_gc%recv_offset_c2f(0:mpisize))
       allocate(face_gc%send_offset_c2f(0:mpisize))
    end if

    face_gc%recv_offset(0) = 0
    face_gc%send_offset(0) = 0
    face_gc%recv_offset_c2f(0) = 0
    face_gc%send_offset_c2f(0) = 0

    ! Determine offsets for sending and receiving with other ranks
    do rank = 1, mpisize
       if (rank-1 == mpirank) then
          face_gc%recv_offset(rank) = face_gc%recv_offset(rank-1)
          face_gc%send_offset(rank) = face_gc%send_offset(rank-1)

          face_gc%recv_offset_c2f(rank) = face_gc%recv_offset_c2f(rank-1)
          face_gc%send_offset_c2f(rank) = face_gc%send_offset_c2f(rank-1)
       else
          face_gc%recv_offset(rank) = face_gc%recv_offset(rank-1) + &
               face_gc%data_size * (i_same(rank-1) + i_c2f(rank-1))
          face_gc%send_offset(rank) = face_gc%send_offset(rank-1) + &
               face_gc%data_size * (i_same(rank-1) + i_f2c(rank-1))

          ! In a second round of communication, handle the fine side of
          ! refinement boundaries
          face_gc%recv_offset_c2f(rank) = face_gc%recv_offset_c2f(rank-1) + &
               face_gc%data_size * 2 * i_f2c(rank-1)
          face_gc%send_offset_c2f(rank) = face_gc%send_offset_c2f(rank-1) + &
               face_gc%data_size * 2 * i_c2f(rank-1)
       end if
    end do

    if (allocated(face_gc%same_local)) then
       ! Deallocate, since they probably changed size
       deallocate(face_gc%same_local, &
            face_gc%phys, &
            face_gc%same_from_buf, &
            face_gc%same_to_buf, &
            face_gc%c2f_local, &
            face_gc%c2f_from_buf, &
            face_gc%c2f_to_buf, &
            face_gc%f2c_from_buf, &
            face_gc%f2c_to_buf)
    end if

    ! Local ghost cell exchange at the same level
    allocate(face_gc%same_local(3, i_same(mpirank)))

    ! Local ghost cell exchange at refinement boundaries
    allocate(face_gc%c2f_local(4, i_c2f(mpirank)))

    ! Physical boundaries
    allocate(face_gc%phys(2, i_phys))

    ! To store indices for faces requiring communication
    do rank = 0, mpisize - 1
       if (rank /= mpirank) then
          allocate(same_ix(rank)%i(i_same(rank)))
          allocate(c2f_ix(rank)%i(i_c2f(rank)))
          allocate(f2c_ix(rank)%i(i_f2c(rank)))
       end if
    end do

    i_same = 0
    i_phys = 0
    i_c2f  = 0
    i_f2c  = 0

    ! Store information needed for ghost cell exchange
    do n = 1, n_faces
       rank = bnd_face(n)%other_proc

       if (bnd_face(n)%face_type == FACE_SAME_LEVEL) then
          i_same(rank) = i_same(rank) + 1

          if (rank == mpirank) then
             face_gc%same_local(:, i_same(rank)) = [bnd_face(n)%quadid(1), &
                  bnd_face(n)%quadid(2), bnd_face(n)%face]
          else
             same_ix(rank)%i(i_same(rank)) = n
          end if
       else if (bnd_face(n)%face_type == FACE_BOUNDARY) then
          i_phys = i_phys + 1
          face_gc%phys(:, i_phys) = [bnd_face(n)%quadid(1), bnd_face(n)%face]
       else if (bnd_face(n)%face_type == FACE_COARSE_TO_FINE) then
          i_c2f(rank) = i_c2f(rank) + 1

          if (rank == mpirank) then
             face_gc%c2f_local(:, i_c2f(rank)) = [bnd_face(n)%quadid(1), &
                  bnd_face(n)%quadid(2), bnd_face(n)%extra, bnd_face(n)%face]
          else
             c2f_ix(rank)%i(i_c2f(rank)) = n
          end if
       else if (bnd_face(n)%face_type == FACE_FINE_TO_COARSE) then
          i_f2c(rank) = i_f2c(rank) + 1

          ! Always non-local; local ones are stored as FACE_COARSE_TO_FINE
          f2c_ix(rank)%i(i_f2c(rank)) = n
       else
          error stop "Unknown face type"
       end if
    end do

    ! Non-local ghost cell exchange at the same level
    n = sum(i_same) - i_same(mpirank)
    allocate(face_gc%same_from_buf(3, n))
    allocate(face_gc%same_to_buf(3, n))

    ! Non-local ghost cell exchange from fine to coarse
    n = sum(i_f2c) - i_f2c(mpirank)
    allocate(face_gc%f2c_from_buf(4, n))
    allocate(face_gc%f2c_to_buf(4, n))

    ! Non-local ghost cell exchange from coarse to fine
    n = sum(i_c2f) - i_c2f(mpirank)
    allocate(face_gc%c2f_from_buf(3, n))
    allocate(face_gc%c2f_to_buf(3, n))

    i_same_nonlocal = 0
    i_f2c_to_buf = 0
    i_f2c_from_buf = 0
    i_c2f_to_buf = 0
    i_c2f_from_buf = 0

    do rank = 0, mpisize - 1
       if (rank == mpirank) cycle

       ! Incremental offsets for sending and receiving data
       i_buf_recv(rank) = face_gc%recv_offset(rank)
       i_buf_send(rank) = face_gc%send_offset(rank)

       ! Boundaries at the same level
       ix_send = same_ix(rank)
       call sort_for_recv_or_send(same_ix(rank), n_faces, bnd_face, recv)
       call sort_for_recv_or_send(ix_send, n_faces, bnd_face, send)

       do n = 1, size(same_ix(rank)%i)
          i = same_ix(rank)%i(n) ! Index in bnd_face array
          j = ix_send%i(n)       ! Index in bnd_face array

          i_same_nonlocal = i_same_nonlocal + 1

          face_gc%same_from_buf(:, i_same_nonlocal) = &
               [bnd_face(i)%quadid(1), i_buf_recv(rank), bnd_face(i)%face]
          face_gc%same_to_buf(:, i_same_nonlocal) = &
               [bnd_face(j)%quadid(1), i_buf_send(rank), bnd_face(j)%face]

          i_buf_recv(rank) = i_buf_recv(rank) + face_gc%data_size
          i_buf_send(rank) = i_buf_send(rank) + face_gc%data_size
       end do

       ! Sending from fine to coarse
       call sort_for_recv_or_send(f2c_ix(rank), n_faces, bnd_face, send)

       do n = 1, size(f2c_ix(rank)%i)
          i = f2c_ix(rank)%i(n) ! Index in bnd_face array
          i_f2c_to_buf = i_f2c_to_buf + 1

          face_gc%f2c_to_buf(:, i_f2c_to_buf) = &
               [bnd_face(i)%quadid(1), bnd_face(i)%extra, i_buf_send(rank), &
               bnd_face(i)%face]
          i_buf_send(rank) = i_buf_send(rank) + face_gc%data_size
       end do

       ! Receiving coarse from fine
       call sort_for_recv_or_send(c2f_ix(rank), n_faces, bnd_face, recv)

       do n = 1, size(c2f_ix(rank)%i)
          i = c2f_ix(rank)%i(n) ! Index in bnd_face array
          i_c2f_from_buf = i_c2f_from_buf + 1

          face_gc%c2f_from_buf(:, i_c2f_from_buf) = &
               [bnd_face(i)%quadid(1), i_buf_recv(rank), bnd_face(i)%face]
          i_buf_recv(rank) = i_buf_recv(rank) + face_gc%data_size
       end do

       ! After the above ghost cells have been update, we can handle the fine
       ! side of refinement boundaries. This involves a new round of
       ! communication, so reset buffer offsets.
       i_buf_recv(rank) = face_gc%recv_offset_c2f(rank)
       i_buf_send(rank) = face_gc%send_offset_c2f(rank)

       ! Sending from coarse to fine
       call sort_for_recv_or_send(c2f_ix(rank), n_faces, bnd_face, send)

       do n = 1, size(c2f_ix(rank)%i)
          i = c2f_ix(rank)%i(n) ! Index in bnd_face array
          i_c2f_to_buf = i_c2f_to_buf + 1

          face_gc%c2f_to_buf(:, i_c2f_to_buf) = &
               [bnd_face(i)%quadid(1), i_buf_send(rank), bnd_face(i)%face]
          i_buf_send(rank) = i_buf_send(rank) + 2 * face_gc%data_size
       end do

       ! Receiving fine from coarse
       call sort_for_recv_or_send(f2c_ix(rank), n_faces, bnd_face, recv)

       do n = 1, size(f2c_ix(rank)%i)
          i = f2c_ix(rank)%i(n) ! Index in bnd_face array
          i_f2c_from_buf = i_f2c_from_buf + 1

          face_gc%f2c_from_buf(:, i_f2c_from_buf) = &
               [bnd_face(i)%quadid(1), bnd_face(i)%extra, i_buf_recv(rank), &
               bnd_face(i)%face]
          i_buf_recv(rank) = i_buf_recv(rank) + 2 * face_gc%data_size
       end do

    end do

  end subroutine f4_get_ghost_cell_pattern

  ! Sort face boundaries for receiving or sending
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
         ! Order by quadid and face of 'this' side
         if (bnd_face(a)%quadid(1) /= bnd_face(b)%quadid(1)) then
            less_than = (bnd_face(a)%quadid(1) < bnd_face(b)%quadid(1))
         else
            less_than = (bnd_face(a)%face < bnd_face(b)%face)
         end if
      else
         ! Order by quadid and face of 'other' side
         if (bnd_face(a)%quadid(2) /= bnd_face(b)%quadid(2)) then
            less_than = (bnd_face(a)%quadid(2) < bnd_face(b)%quadid(2))
         else
            ! Note that the face is swapped for the receiving side
            less_than = (face_swap(bnd_face(a)%face) < face_swap(bnd_face(b)%face))
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

  subroutine f4_fill_ghostcell_buffers(f4, n_vars, i_vars)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer                      :: i, j, n, ivar, iv
    integer                      :: i_f, j_f, half_bx(2)
    integer                      :: iq, jq, i_buf, face

    if (maxval(f4%face_gc%send_offset) * n_vars > size(f4%send_buffer)) &
         error stop "send buffer too small"

    half_bx = f4%bx/2

    associate (bx => f4%bx, n_gc => f4%n_gc, uu => f4%uu, face_gc => f4%face_gc)
      do n = 1, size(face_gc%same_to_buf, 2)
         ! +1 to account for index offset between C and Fortran
         iq = face_gc%same_to_buf(1, n) + 1
         i_buf = face_gc%same_to_buf(2, n) * n_vars
         face = face_gc%same_to_buf(3, n)

         select case (face)
         case (0)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, bx(2)
                  do i = 1, n_gc
                     i_buf = i_buf + 1
                     f4%send_buffer(i_buf) = uu(i, j, ivar, iq)
                  end do
               end do
            end do
         case (1)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, bx(2)
                  do i = 1, n_gc
                     i_buf = i_buf + 1
                     f4%send_buffer(i_buf) = uu(bx(1)-n_gc+i, j, ivar, iq)
                  end do
               end do
            end do
         case (2)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, n_gc
                  do i = 1, bx(1)
                     i_buf = i_buf + 1
                     f4%send_buffer(i_buf) = uu(i, j, ivar, iq)
                  end do
               end do
            end do
         case (3)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, n_gc
                  do i = 1, bx(1)
                     i_buf = i_buf + 1
                     f4%send_buffer(i_buf) = uu(i, bx(2)-n_gc+j, ivar, iq)
                  end do
               end do
            end do
         end select
      end do

      do n = 1, size(face_gc%f2c_to_buf, 2)
         iq = face_gc%f2c_to_buf(1, n) + 1 ! fine block
         jq = face_gc%f2c_to_buf(2, n) + 1 ! fine block
         i_buf = face_gc%f2c_to_buf(3, n) * n_vars
         face = face_gc%f2c_to_buf(4, n)

         select case (face)
         case (0)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               ! First fine block
               do j = 1, half_bx(2)
                  j_f = 2 * j - 1
                  do i = 1, n_gc
                     i_f = 2 * i - 1
                     i_buf = i_buf + 1
                     f4%send_buffer(i_buf) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
                  end do
               end do

               ! Second fine block
               do j = 1, half_bx(2)
                  j_f = 2 * j - 1
                  do i = 1, n_gc
                     i_f = 2 * i - 1
                     i_buf = i_buf + 1
                     f4%send_buffer(i_buf) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, jq))
                  end do
               end do
            end do

         case (1)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               ! First fine block
               do j = 1, half_bx(2)
                  j_f = 2 * j - 1
                  do i = 1, n_gc
                     i_f = bx(1) - 2 * n_gc + 2 * i - 1
                     i_buf = i_buf + 1
                     f4%send_buffer(i_buf) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
                  end do
               end do

               ! Second fine block
               do j = 1, half_bx(2)
                  j_f = 2 * j - 1
                  do i = 1, n_gc
                     i_f = bx(1) - 2 * n_gc + 2 * i - 1
                     i_buf = i_buf + 1
                     f4%send_buffer(i_buf) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, jq))
                  end do
               end do
            end do

         case (2)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               ! First fine block
               do j = 1, n_gc
                  j_f = 2 * j - 1
                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1
                     i_buf = i_buf + 1
                     f4%send_buffer(i_buf) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
                  end do
               end do

               ! Second fine block
               do j = 1, n_gc
                  j_f = 2 * j - 1
                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1
                     i_buf = i_buf + 1
                     f4%send_buffer(i_buf) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, jq))
                  end do
               end do
            end do

         case (3)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               ! First fine block
               do j = 1, n_gc
                  j_f = bx(2) - 2 * n_gc + 2 * j - 1
                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1
                     i_buf = i_buf + 1
                     f4%send_buffer(i_buf) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, iq))
                  end do
               end do

               ! Second fine block
               do j = 1, n_gc
                  j_f = bx(2) - 2 * n_gc + 2 * j - 1
                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1
                     i_buf = i_buf + 1
                     f4%send_buffer(i_buf) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, jq))
                  end do
               end do
            end do
         end select
      end do

      f4%recv_offset = face_gc%recv_offset * n_vars
      f4%send_offset = face_gc%send_offset * n_vars
    end associate

  end subroutine f4_fill_ghostcell_buffers

  subroutine f4_fill_ghostcell_buffers_c2f(f4, n_vars, i_vars)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer                      :: i, j, n, ivar, iv
    integer                      :: i_c, j_c, half_bx(2)
    integer                      :: iq, i_buf, face
    integer                      :: half_n_gc
    logical                      :: odd_n_gc
    real(dp)                     :: fine(4)

    if (maxval(f4%face_gc%send_offset_c2f) * n_vars > size(f4%send_buffer)) &
         error stop "send buffer too small"

    half_bx = f4%bx/2
    half_n_gc = f4%n_gc/2 ! Round down
    odd_n_gc  = (iand(f4%n_gc, 1) == 1)

    associate (bx => f4%bx, n_gc => f4%n_gc, uu => f4%uu, face_gc => f4%face_gc)

      do n = 1, size(face_gc%c2f_to_buf, 2)
         iq = face_gc%c2f_to_buf(1, n) + 1 ! coarse block
         i_buf = face_gc%c2f_to_buf(2, n) * n_vars
         face = face_gc%c2f_to_buf(3, n)

         select case (face)
         case (0)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, half_bx(2)
                  do i = 1, half_n_gc
                     i_c = i

                     ! First fine block
                     j_c = j
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%send_buffer(i_buf+1:i_buf+4) = fine
                     i_buf = i_buf + 4

                     ! Second fine block
                     j_c = j_c + half_bx(2)
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%send_buffer(i_buf+1:i_buf+4) = fine
                     i_buf = i_buf + 4
                  end do

                  if (odd_n_gc) then
                     i_c = 1 + half_n_gc

                     ! First fine block
                     j_c = j
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%send_buffer(i_buf+1:i_buf+2) = fine([1, 3])
                     i_buf = i_buf + 2

                     ! Second fine block
                     j_c = j_c + half_bx(2)
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%send_buffer(i_buf+1:i_buf+2) = fine([1, 3])
                     i_buf = i_buf + 2
                  end if
               end do
            end do

         case (1)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               do j = 1, half_bx(2)
                  do i = 1, half_n_gc
                     i_c = bx(1) - half_n_gc + i

                     ! First fine block
                     j_c = j
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%send_buffer(i_buf+1:i_buf+4) = fine
                     i_buf = i_buf + 4

                     ! Second fine block
                     j_c = j_c + half_bx(2)
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%send_buffer(i_buf+1:i_buf+4) = fine
                     i_buf = i_buf + 4
                  end do

                  if (odd_n_gc) then
                     i_c = bx(1) - half_n_gc

                     ! First fine block
                     j_c = j
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%send_buffer(i_buf+1:i_buf+2) = fine([2, 4])
                     i_buf = i_buf + 2

                     ! Second fine block
                     j_c = j_c + half_bx(2)
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%send_buffer(i_buf+1:i_buf+2) = fine([2, 4])
                     i_buf = i_buf + 2
                  end if
               end do
            end do

         case (2)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               do j = 1, half_n_gc
                  j_c = j
                  do i = 1, half_bx(1)
                     ! First fine block
                     i_c = i
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%send_buffer(i_buf+1:i_buf+4) = fine
                     i_buf = i_buf + 4

                     ! Second fine block
                     i_c = i_c + half_bx(1)
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
                     ! First fine block
                     i_c = i
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%send_buffer(i_buf+1:i_buf+2) = fine([1, 2])
                     i_buf = i_buf + 2

                     ! Second fine block
                     i_c = i_c + half_bx(2)
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%send_buffer(i_buf+1:i_buf+2) = fine([1, 2])
                     i_buf = i_buf + 2
                  end do
               end if
            end do

         case (3)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, half_n_gc
                  j_c = bx(2) - half_n_gc + j
                  do i = 1, half_bx(1)
                     ! First fine block
                     i_c = i
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%send_buffer(i_buf+1:i_buf+4) = fine
                     i_buf = i_buf + 4

                     ! Second fine block
                     i_c = i_c + half_bx(1)
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
                     ! First fine block
                     i_c = i
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%send_buffer(i_buf+1:i_buf+2) = fine([3, 4])
                     i_buf = i_buf + 2

                     ! Second fine block
                     i_c = i_c + half_bx(2)
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%send_buffer(i_buf+1:i_buf+2) = fine([3, 4])
                     i_buf = i_buf + 2
                  end do
               end if
            end do
         end select
      end do

      f4%recv_offset = face_gc%recv_offset_c2f * n_vars
      f4%send_offset = face_gc%send_offset_c2f * n_vars
    end associate

  end subroutine f4_fill_ghostcell_buffers_c2f

  ! Exchange the receive and send buffers according to the specified offsets
  ! per MPI rank
  subroutine f4_exchange_buffers(f4)
    type(foap4_t), intent(inout) :: f4
    type(MPI_Request)            :: send_req(0:mpisize-1)
    type(MPI_Request)            :: recv_req(0:mpisize-1)
    integer                      :: n_send, n_recv, ilo, ihi, ierr, rank
    integer, parameter           :: tag = 0

    n_send = 0
    n_recv = 0

    do rank = 0, mpisize - 1
       ilo = f4%send_offset(rank) + 1
       ihi = f4%send_offset(rank+1)

       if (ihi >= ilo) then
          n_send = n_send + 1
          call mpi_isend(f4%send_buffer(ilo:ihi), ihi-ilo+1, &
               MPI_DOUBLE_PRECISION, rank, tag, mpicomm, &
               send_req(n_send), ierr)
       end if

       ilo = f4%recv_offset(rank) + 1
       ihi = f4%recv_offset(rank+1)

       if (ihi >= ilo) then
          n_recv = n_recv + 1
          call mpi_irecv(f4%recv_buffer(ilo:ihi), ihi-ilo+1, MPI_DOUBLE_PRECISION, &
               rank, 0, mpicomm, recv_req(n_recv), ierr)
       end if
    end do

    call mpi_waitall(n_recv, recv_req(1:n_recv), MPI_STATUSES_IGNORE, ierr)
    call mpi_waitall(n_send, send_req(1:n_send), MPI_STATUSES_IGNORE, ierr)

  end subroutine f4_exchange_buffers

  subroutine f4_update_ghostcells(f4, n_vars, i_vars)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_vars
    integer, intent(in)          :: i_vars(n_vars)
    integer                      :: n, i, j, iq, jq, kq, face
    integer                      :: i_buf, iv, ivar, i_f, j_f, i_c, j_c
    integer                      :: half_bx(2), half_n_gc
    real(dp)                     :: slope, fine(4)
    logical                      :: odd_n_gc

    call f4_update_ghostcell_pattern(f4)
    call f4_fill_ghostcell_buffers(f4, n_vars, i_vars)
    call f4_exchange_buffers(f4)

    half_bx   = f4%bx/2
    half_n_gc = f4%n_gc/2 ! Round down
    odd_n_gc  = (iand(f4%n_gc, 1) == 1)

    associate (bx => f4%bx, n_gc => f4%n_gc, uu => f4%uu, face_gc => f4%face_gc)
      do n = 1, size(face_gc%same_from_buf, 2)
         ! +1 to account for index offset between C and Fortran
         iq = face_gc%same_from_buf(1, n) + 1
         i_buf = face_gc%same_from_buf(2, n) * n_vars
         face = face_gc%same_from_buf(3, n)

         select case (face)
         case (0)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, bx(2)
                  do i = 1, n_gc
                     i_buf = i_buf + 1
                     uu(-n_gc+1, j, ivar, iq) = f4%recv_buffer(i_buf)
                  end do
               end do
            end do
         case (1)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, bx(2)
                  do i = 1, n_gc
                     i_buf = i_buf + 1
                     uu(bx(1)+i, j, ivar, iq) = f4%recv_buffer(i_buf)
                  end do
               end do
            end do
         case (2)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, n_gc
                  do i = 1, bx(1)
                     i_buf = i_buf + 1
                     uu(i, -n_gc+j, ivar, iq) = f4%recv_buffer(i_buf)
                  end do
               end do
            end do
         case (3)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, n_gc
                  do i = 1, bx(1)
                     i_buf = i_buf + 1
                     uu(i, bx(2)+j, ivar, iq) = f4%recv_buffer(i_buf)
                  end do
               end do
            end do
         end select
      end do

      ! Local coarse-to-fine refinement boundaries, fill coarse side
      do n = 1, size(face_gc%c2f_local, 2)
         iq   = face_gc%c2f_local(1, n) + 1 ! Coarse block
         jq   = face_gc%c2f_local(2, n) + 1 ! First fine block
         kq   = face_gc%c2f_local(3, n) + 1 ! Second fine block
         face = face_gc%c2f_local(4, n)

         select case (face)
         case (0)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               ! First fine block
               do j = 1, half_bx(2)
                  j_f = 2 * j - 1
                  do i = 1, n_gc
                     i_f = bx(1) - 2 * n_gc + 2 * i - 1
                     uu(-n_gc+i, j, ivar, iq) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, jq))
                  end do
               end do

               ! Second fine block
               do j = 1, half_bx(2)
                  j_f = 2 * j - 1
                  do i = 1, n_gc
                     i_f = bx(1) - 2 * n_gc + 2 * i - 1
                     uu(-n_gc+i, half_bx(2)+j, ivar, iq) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, kq))
                  end do
               end do
            end do
         case (1)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               ! First fine block
               do j = 1, half_bx(2)
                  j_f = 2 * j - 1
                  do i = 1, n_gc
                     i_f = 2 * i - 1
                     uu(bx(1)+i, j, ivar, iq) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, jq))
                  end do
               end do

               ! Second fine block
               do j = 1, half_bx(2)
                  j_f = 2 * j - 1
                  do i = 1, n_gc
                     i_f = 2 * i - 1
                     uu(bx(1)+i, half_bx(2)+j, ivar, iq) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, kq))
                  end do
               end do
            end do
         case (2)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               ! First fine block
               do j = 1, n_gc
                  j_f = bx(2) - 2 * n_gc + 2 * j - 1
                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1
                     uu(i, -n_gc+j, ivar, iq) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, jq))
                  end do
               end do

               ! Second fine block
               do j = 1, n_gc
                  j_f = bx(2) - 2 * n_gc + 2 * j - 1
                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1
                     uu(half_bx(1)+i, -n_gc+j, ivar, iq) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, kq))
                  end do
               end do
            end do
         case (3)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               ! First fine block
               do j = 1, n_gc
                  j_f = 2 * j - 1
                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1
                     uu(i, bx(2)+j, ivar, iq) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, jq))
                  end do
               end do

               ! Second fine block
               do j = 1, n_gc
                  j_f = 2 * j - 1
                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1
                     uu(half_bx(1)+i, bx(2)+j, ivar, iq) = 0.25_dp * &
                          sum(uu(i_f:i_f+1, j_f:j_f+1, ivar, kq))
                  end do
               end do
            end do
         end select
      end do

      ! Coarse-to-fine, update coarse side from buffers
      do n = 1, size(face_gc%c2f_from_buf, 2)
         iq    = face_gc%c2f_from_buf(1, n) + 1 ! Coarse block
         i_buf = face_gc%c2f_from_buf(2, n) * n_vars
         face  = face_gc%c2f_from_buf(3, n)

         select case (face)
         case (0)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               ! First fine block
               do j = 1, half_bx(2)
                  do i = 1, n_gc
                     i_buf = i_buf + 1
                     uu(-n_gc+i, j, ivar, iq) = f4%recv_buffer(i_buf)
                  end do
               end do

               ! Second fine block
               do j = 1, half_bx(2)
                  do i = 1, n_gc
                     i_buf = i_buf + 1
                     uu(-n_gc+i, half_bx(2)+j, ivar, iq) = f4%recv_buffer(i_buf)
                  end do
               end do
            end do
         case (1)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               ! First fine block
               do j = 1, half_bx(2)
                  do i = 1, n_gc
                     i_buf = i_buf + 1
                     uu(bx(1)+i, j, ivar, iq) = f4%recv_buffer(i_buf)
                  end do
               end do

               ! Second fine block
               do j = 1, half_bx(2)
                  do i = 1, n_gc
                     i_buf = i_buf + 1
                     uu(bx(1)+i, half_bx(2)+j, ivar, iq) = f4%recv_buffer(i_buf)
                  end do
               end do
            end do
         case (2)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               ! First fine block
               do j = 1, n_gc
                  do i = 1, half_bx(1)
                     i_buf = i_buf + 1
                     uu(i, -n_gc+j, ivar, iq) = f4%recv_buffer(i_buf)
                  end do
               end do

               ! Second fine block
               do j = 1, n_gc
                  do i = 1, half_bx(1)
                     i_buf = i_buf + 1
                     uu(half_bx(1)+i, -n_gc+j, ivar, iq) = f4%recv_buffer(i_buf)
                  end do
               end do
            end do
         case (3)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               ! First fine block
               do j = 1, n_gc
                  do i = 1, half_bx(1)
                     i_buf = i_buf + 1
                     uu(i, bx(2)+j, ivar, iq) = f4%recv_buffer(i_buf)
                  end do
               end do

               ! Second fine block
               do j = 1, n_gc
                  do i = 1, half_bx(1)
                     i_buf = i_buf + 1
                     uu(half_bx(1)+i, bx(2)+j, ivar, iq) = f4%recv_buffer(i_buf)
                  end do
               end do
            end do
         end select
      end do

      ! Local boundaries at the same refinement level
      do n = 1, size(face_gc%same_local, 2)
         iq   = face_gc%same_local(1, n) + 1
         jq   = face_gc%same_local(2, n) + 1
         face = face_gc%same_local(3, n)

         if (face == 1) then
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, bx(2)
                  do i = 1, n_gc
                     uu(bx(1)+i, j, ivar, iq) = uu(i, j, ivar, jq)
                     uu(-n_gc+i, j, ivar, jq) = uu(bx(1)-n_gc+i, j, ivar, iq)
                  end do
               end do
            end do
         else if (face == 3) then
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, n_gc
                  do i = 1, bx(1)
                     uu(i, bx(2)+j, ivar, iq) = uu(i, j, ivar, jq)
                     uu(i, -n_gc+j, ivar, jq) = uu(i, bx(2)-n_gc+j, ivar, iq)
                  end do
               end do
            end do
         else
            error stop "face not implemented"
         end if
      end do

      ! Physical boundaries
      do n = 1, size(face_gc%phys, 2)
         iq = face_gc%phys(1, n) + 1
         face = face_gc%phys(2, n)

         select case (face)
         case (0)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, bx(2)
                  slope = uu(2, j, ivar, iq) - uu(1, j, ivar, iq)
                  do i = 1, n_gc
                     uu(1-i, j, ivar, iq) = uu(1, j, ivar, iq) - i * slope
                  end do
               end do
            end do
         case (1)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, bx(2)
                  slope = uu(bx(1), j, ivar, iq) - uu(bx(1)-1, j, ivar, iq)
                  do i = 1, n_gc
                     uu(bx(1)+i, j, ivar, iq) = uu(bx(1), j, ivar, iq) + i * slope
                  end do
               end do
            end do
         case (2)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, n_gc
                  do i = 1, bx(1)
                     slope = uu(i, 2, ivar, iq) - uu(i, 1, ivar, iq)
                     uu(i, 1-j, ivar, iq) = uu(i, 1, ivar, iq) - j * slope
                  end do
               end do
            end do
         case (3)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, n_gc
                  do i = 1, bx(1)
                     slope = uu(i, bx(2), ivar, iq) - uu(i, bx(2)-1, ivar, iq)
                     uu(i, bx(2)+j, ivar, iq) = uu(i, bx(2), ivar, iq) + j * slope
                  end do
               end do
            end do
         end select
      end do

      ! Do coarse-to-fine refinement boundaries last, so that ghost cells
      ! required for interpolation have been filled
      call f4_fill_ghostcell_buffers_c2f(f4, n_vars, i_vars)
      call f4_exchange_buffers(f4)

      ! Local coarse-to-fine boundaries, fill fine side
      do n = 1, size(face_gc%c2f_local, 2)
         iq   = face_gc%c2f_local(1, n) + 1 ! Coarse block
         jq   = face_gc%c2f_local(2, n) + 1 ! First fine block
         kq   = face_gc%c2f_local(3, n) + 1 ! Second fine block
         face = face_gc%c2f_local(4, n)

         select case (face)
         case (0)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               do j = 1, half_bx(2)
                  j_f = 2 * j - 1
                  do i = 1, half_n_gc
                     i_c = i
                     i_f = bx(1) + 2*i - 1

                     ! First fine block
                     j_c = j
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, jq) = reshape(fine, [2, 2])

                     ! Second fine block
                     j_c = j_c + half_bx(2)
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, kq) = reshape(fine, [2, 2])
                  end do

                  if (odd_n_gc) then
                     i_c = 1 + half_n_gc
                     i_f = bx(1) + n_gc

                     ! First fine block
                     j_c = j
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f, j_f:j_f+1, ivar, jq) = fine([1, 3])

                     ! Second fine block
                     j_c = j_c + half_bx(2)
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f, j_f:j_f+1, ivar, kq) = fine([1, 3])
                  end if
               end do
            end do

         case (1)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               do j = 1, half_bx(2)
                  j_f = 2 * j - 1
                  do i = 1, half_n_gc
                     i_c = bx(1) - half_n_gc + i
                     i_f = -(2 * half_n_gc) + 2*i - 1

                     ! First fine block
                     j_c = j
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, jq) = reshape(fine, [2, 2])

                     ! Second fine block
                     j_c = j_c + half_bx(2)
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, kq) = reshape(fine, [2, 2])
                  end do

                  if (odd_n_gc) then
                     i_c = bx(1) - half_n_gc
                     i_f = -n_gc + 1

                     ! First fine block
                     j_c = j
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f, j_f:j_f+1, ivar, jq) = fine([2, 4])

                     ! Second fine block
                     j_c = j_c + half_bx(2)
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f, j_f:j_f+1, ivar, kq) = fine([2, 4])
                  end if
               end do
            end do

         case (2)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               do j = 1, half_n_gc
                  j_c = j
                  j_f = bx(2) + 2*j - 1
                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1

                     ! First fine block
                     i_c = i
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, jq) = reshape(fine, [2, 2])

                     ! Second fine block
                     i_c = i_c + half_bx(1)
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, kq) = reshape(fine, [2, 2])
                  end do
               end do

               if (odd_n_gc) then
                  j_c = 1 + half_n_gc
                  j_f = bx(2) + n_gc

                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1

                     ! First fine block
                     i_c = i
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f:i_f+1, j_f, ivar, jq) = fine([1, 2])

                     ! Second fine block
                     i_c = i_c + half_bx(2)
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f:i_f+1, j_f, ivar, kq) = fine([1, 2])
                  end do
               end if
            end do

         case (3)
            do iv = 1, n_vars
               ivar = i_vars(iv)

               do j = 1, half_n_gc
                  j_c = bx(2) - half_n_gc + j
                  j_f = -(2 * half_n_gc) + 2*j - 1
                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1

                     ! First fine block
                     i_c = i
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, jq) = reshape(fine, [2, 2])

                     ! Second fine block
                     i_c = i_c + half_bx(1)
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, kq) = reshape(fine, [2, 2])
                  end do
               end do

               if (odd_n_gc) then
                  j_c = bx(2) - half_n_gc
                  j_f = -n_gc + 1

                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1

                     ! First fine block
                     i_c = i
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f:i_f+1, j_f, ivar, jq) = fine([3, 4])

                     ! Second fine block
                     i_c = i_c + half_bx(2)
                     call prolong_local_5point(f4%uu(i_c, j_c, iv, iq), &
                          [f4%uu(i_c-1, j_c, iv, iq), f4%uu(i_c+1, j_c, iv, iq)], &
                          [f4%uu(i_c, j_c-1, iv, iq), f4%uu(i_c, j_c+1, iv, iq)], fine)
                     f4%uu(i_f:i_f+1, j_f, ivar, kq) = fine([3, 4])
                  end do
               end if
            end do
         end select
      end do

      ! Nonlocal coarse-to-fine boundaries, fill fine side
      do n = 1, size(face_gc%f2c_from_buf, 2)
         iq    = face_gc%f2c_from_buf(1, n) + 1 ! Fine block
         jq    = face_gc%f2c_from_buf(2, n) + 1 ! Fine block
         i_buf = face_gc%f2c_from_buf(3, n) * n_vars
         face  = face_gc%f2c_from_buf(4, n)

         select case (face)
         case (1)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, half_bx(2)
                  j_f = 2 * j - 1
                  do i = 1, half_n_gc
                     i_f = bx(1) + 2*i - 1
                     ! First fine block
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, iq) = &
                          reshape(f4%recv_buffer(i_buf+1:i_buf+4), [2, 2])
                     i_buf = i_buf + 4

                     ! Second fine block
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, jq) = &
                          reshape(f4%recv_buffer(i_buf+1:i_buf+4), [2, 2])
                     i_buf = i_buf + 4
                  end do

                  if (odd_n_gc) then
                     i_f = bx(1) + n_gc

                     ! First fine block
                     f4%uu(i_f, j_f:j_f+1, ivar, iq) = &
                          f4%recv_buffer(i_buf+1:i_buf+2)
                     i_buf = i_buf + 2

                     ! Second fine block
                     f4%uu(i_f, j_f:j_f+1, ivar, jq) = &
                          f4%recv_buffer(i_buf+1:i_buf+2)
                     i_buf = i_buf + 2
                  end if
               end do
            end do

         case (0)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, half_bx(2)
                  j_f = 2 * j - 1
                  do i = 1, half_n_gc
                     i_f = -(2 * half_n_gc) + 2*i - 1

                     ! First fine block
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, iq) = &
                          reshape(f4%recv_buffer(i_buf+1:i_buf+4), [2, 2])
                     i_buf = i_buf + 4

                     ! Second fine block
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, jq) = &
                          reshape(f4%recv_buffer(i_buf+1:i_buf+4), [2, 2])
                     i_buf = i_buf + 4
                  end do

                  if (odd_n_gc) then
                     i_f = -n_gc + 1

                     ! First fine block
                     f4%uu(i_f, j_f:j_f+1, ivar, iq) = &
                          f4%recv_buffer(i_buf+1:i_buf+2)
                     i_buf = i_buf + 2

                     ! Second fine block
                     f4%uu(i_f, j_f:j_f+1, ivar, jq) = &
                          f4%recv_buffer(i_buf+1:i_buf+2)
                     i_buf = i_buf + 2
                  end if
               end do
            end do

         case (3)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, half_n_gc
                  j_f = bx(2) + 2*j - 1
                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1

                     ! First fine block
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, iq) = &
                          reshape(f4%recv_buffer(i_buf+1:i_buf+4), [2, 2])
                     i_buf = i_buf + 4

                     ! Second fine block
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, jq) = &
                          reshape(f4%recv_buffer(i_buf+1:i_buf+4), [2, 2])
                     i_buf = i_buf + 4
                  end do
               end do

               if (odd_n_gc) then
                  j_f = bx(2) + n_gc
                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1

                     ! First fine block
                     f4%uu(i_f:i_f+1, j_f, ivar, iq) = &
                          f4%recv_buffer(i_buf+1:i_buf+2)
                     i_buf = i_buf + 2

                     ! Second fine block
                     f4%uu(i_f:i_f+1, j_f, ivar, jq) = &
                          f4%recv_buffer(i_buf+1:i_buf+2)
                     i_buf = i_buf + 2
                  end do
               end if
            end do

         case (2)
            do iv = 1, n_vars
               ivar = i_vars(iv)
               do j = 1, half_n_gc
                  j_f = -(2 * half_n_gc) + 2*j - 1
                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1

                     ! First fine block
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, iq) = &
                          reshape(f4%recv_buffer(i_buf+1:i_buf+4), [2, 2])
                     i_buf = i_buf + 4

                     ! Second fine block
                     f4%uu(i_f:i_f+1, j_f:j_f+1, ivar, jq) = &
                          reshape(f4%recv_buffer(i_buf+1:i_buf+4), [2, 2])
                     i_buf = i_buf + 4
                  end do
               end do

               if (odd_n_gc) then
                  j_f = -n_gc + 1
                  do i = 1, half_bx(1)
                     i_f = 2 * i - 1

                     ! First fine block
                     f4%uu(i_f:i_f+1, j_f, ivar, iq) = &
                          f4%recv_buffer(i_buf+1:i_buf+2)
                     i_buf = i_buf + 2

                     ! Second fine block
                     f4%uu(i_f:i_f+1, j_f, ivar, jq) = &
                          f4%recv_buffer(i_buf+1:i_buf+2)
                     i_buf = i_buf + 2
                  end do
               end if
            end do
         end select
      end do
    end associate

  end subroutine f4_update_ghostcells

  subroutine f4_adjust_refinement(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, n_blocks_new, n_blocks_old, k, iv
    integer                      :: has_changed

    f4%refinement_flags(1:f4%n_blocks) = 0
    if (mpirank == 0) f4%refinement_flags(1) = 1

    n_blocks_old = f4%n_blocks
    call pw_adjust_refinement(f4%pw, f4%n_blocks, &
         f4%refinement_flags(1:f4%n_blocks), has_changed)

    if (has_changed == 0) return

    n_blocks_new = pw_get_num_local_quadrants(f4%pw)

    if (n_blocks_new + f4%n_blocks > f4%max_blocks) &
         error stop "Not enough memory for tree copy during refinement"

    ! Copy blocks to end of list. The backward order ensures old data is not
    ! overwritten before it is copied.
    do n = n_blocks_old, 1, -1
       call copy_block(f4, n, n_blocks_new+n)
    end do

    call f4_get_quadrants(f4)

    k = n_blocks_new + 1
    n = 1
    do
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

       if (n == n_blocks_new + 1) then
          if (k /= n_blocks_new + n_blocks_old + 1) &
               error stop "Refinement: loops do not end simultaneously"
          exit
       end if
    end do
  end subroutine f4_adjust_refinement

  ! Copy block data
  subroutine copy_block(f4, i_from, i_to)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: i_from, i_to
    f4%block_origin(:, i_to) = f4%block_origin(:, i_from)
    f4%block_level(i_to)     = f4%block_level(i_from)
    f4%uu(:, :, :, i_to)     = f4%uu(:, :, :, i_from)
  end subroutine copy_block

  ! Coarsen a family of child blocks to their parent
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

  ! Prolong to family of child blocks
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

  subroutine prolong_local_5point(coarse_c, coarse_x, coarse_y, fine)
    real(dp), intent(in)  :: coarse_c    ! Center value
    real(dp), intent(in)  :: coarse_x(2) ! x-neighbors (-1, +1)
    real(dp), intent(in)  :: coarse_y(2) ! y-neighbors (-1, +1)
    real(dp), intent(out) :: fine(4)
    real(dp)              :: f(0:2), slopes_a(2), slopes_b(2)

    f(0) = coarse_c          ! Identical to coarse_y(2)
    slopes_a = [coarse_c - coarse_x(1), coarse_c - coarse_y(1)]
    slopes_b = [coarse_x(2) - coarse_c, coarse_y(2) - coarse_c]
    f(1:2) = 0.25_dp * af_limiter_minmod(slopes_a, slopes_b)

    fine(1) = f(0) - f(1) - f(2)
    fine(2) = f(0) + f(1) - f(2)
    fine(3) = f(0) - f(1) + f(2)
    fine(4) = f(0) + f(1) + f(2)
  end subroutine prolong_local_5point

  !> Generalized minmod limiter. The parameter theta controls how dissipative
  !> the limiter is, with 1 corresponding to the minmod limiter and 2 to the MC
  !> limiter.
  elemental function af_limiter_gminmod(a, b, theta) result(phi)
    real(dp), intent(in) :: a, b, theta
    real(dp)             :: phi

    if (a * b > 0) then
       phi = sign(minval(abs([theta * a, theta * b, &
            0.5_dp * (a + b)])), a)
    else
       phi = 0.0_dp
    end if
  end function af_limiter_gminmod

  elemental function af_limiter_minmod(a, b) result(phi)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp)             :: phi
    phi = af_limiter_gminmod(a, b, 1.0_dp)
  end function af_limiter_minmod

end module m_foap4
