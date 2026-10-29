!> Module for writing output
!>
!> The XDMF part was inspired by the implementation in the bhac code, see
!> src/amrvacio/amrio.t in
!> https://gitlab.itp.uni-frankfurt.de/BHAC-release/bhac
#:include 'definitions_ndim.fpp'
#:include 'definitions_parallel.fpp'
module m_io_${NDIM}$d
  use mpi_f08
  use, intrinsic :: iso_c_binding
  use iso_fortran_env, only: int64
  use m_foap4_types_${NDIM}$d

  implicit none
  private

  type :: int_list_${NDIM}$d
     integer, allocatable :: id(:)
  end type int_list_${NDIM}$d

  ! For hashing
  type :: key_t
     integer :: x(NDIM+1)
  end type key_t

  ! Minimal version of https://github.com/jannisteunissen/ffhash
  type :: ffh_t
     integer                   :: n_buckets     = 0
     integer                   :: n_keys_stored = 0
     integer                   :: hash_mask     = 0
     !> 0 = empty, 1 = filled
     integer, allocatable      :: flags(:)
     type(key_t), allocatable  :: keys(:)
     integer, allocatable      :: vals(:)
   contains
     procedure :: store_value => ffh_store_value
     procedure :: get_value   => ffh_get_value
  end type ffh_t

  !> Number of bytes per stored (cell-centered) value
  integer, parameter :: bytes_per_val = storage_size(1.0_fp) / 8

  !> Relative tolerance used to detect domain boundaries
  real(dp), parameter :: bnd_rel_tol = 1e-13_dp

  ! Public methods
  public :: io_write_grid
  public :: io_xdmf_write_blocks_${NDIM}$DCoRect

contains

  !> Write the AMR grid to a file
  subroutine io_write_grid(f4, fname, n_output, viewer, write_p4vtu, n_gc_out)
    type(foap4_t), intent(inout)           :: f4
    character(len=*), intent(in)           :: fname    !< Base file name
    integer, intent(in)                    :: n_output !< Output index
    character(len=*), intent(in), optional :: viewer   !< Optimize for viewer
    !> Also write p4est vtu file (default: false)
    logical, intent(in), optional          :: write_p4vtu
    !> Number of ghost cells to include in output
    integer, intent(in), optional          :: n_gc_out
    character(len=len_trim(fname)+7)       :: full_fname
    integer                                :: n, out_gc
    real(dp), allocatable                  :: dr(:, :)
    real(dp)                               :: t0, t1

    t0 = MPI_Wtime()

    out_gc = 0; if (present(n_gc_out)) out_gc = n_gc_out
    if (out_gc < 0 .or. out_gc > f4%n_gc) error stop "Invalid n_gc_out"

    ! Get the block data from the device
    $:UPDATE_SELF('f4%uu(' + DTIMES(':') + ', :, 1:f4%n_blocks)')

    write(full_fname, "(A,A,I06.6)") trim(fname), "_", n_output

    if (present(write_p4vtu)) then
       if (write_p4vtu) then
          call pw_vtk_write_file(f4%pw, trim(full_fname) // C_null_char)
       end if
    end if

    allocate(dr(NDIM, f4%n_blocks))

    do n = 1, f4%n_blocks
       dr(:, n) = f4%dr_level(:, f4%block_level(n))
    end do

    call io_xdmf_write_blocks_${NDIM}$DCoRect(f4%mpicomm, trim(full_fname), &
         f4%n_blocks, f4%bx, f4%n_vars, &
         f4%var_names(1:f4%n_vars), f4%n_gc, out_gc, &
         f4%block_origin(:, 1:f4%n_blocks), dr, f4%block_level(1:f4%n_blocks),&
         f4%r_min, f4%r_max, &
         get_block_cc_data=get_block_data, time=f4%time, viewer=viewer)
    t1 = MPI_Wtime()
    f4%wtimes(f4_timer_write_grid) = f4%wtimes(f4_timer_write_grid) + t1 - t0

  contains

    subroutine get_block_data(i_block, cc_data)
      integer, intent(in)     :: i_block
      real(fp), intent(inout) :: cc_data(@{DTIMES(:)}@, :)
      cc_data = f4%uu(@{DTIMES(:)}@, 1:f4%n_vars, i_block)
    end subroutine get_block_data

  end subroutine io_write_grid

  !> Write block data to binary files (one per task) and a single .xdmf header file
  subroutine io_xdmf_write_blocks_${NDIM}$DCoRect(mpicomm, filename, n_blocks, nx, n_cc, &
       cc_names, in_gc, out_gc, origin, dr, level, r_min, r_max, get_block_cc_data, &
       time, viewer, fill_diagonals)
    type(MPI_comm), intent(in)   :: mpicomm            !< MPI communicator
    character(len=*), intent(in) :: filename           !< File name without extension
    integer, intent(in)          :: n_blocks           !< Number of blocks
    integer, intent(in)          :: nx(NDIM)           !< Size of the blocks (excl. ghost cells)
    integer, intent(in)          :: n_cc               !< Number of variables
    character(len=*), intent(in) :: cc_names(n_cc)     !< Names of variables
    integer, intent(in)          :: in_gc              !< Number of ghost cells in input
    integer, intent(in)          :: out_gc             !< Number of ghost cells to write
    !> Origin of each block
    real(dp), intent(in)         :: origin(NDIM, n_blocks)
    real(dp), intent(in)         :: dr(NDIM, n_blocks) !< Grid spacing of each block
    integer, intent(in)          :: level(n_blocks)    !< Level of each block
    real(dp)                     :: r_min(NDIM)        !< Min. coordinate of domain
    real(dp)                     :: r_max(NDIM)        !< Max. coordinate of domain
    !> Method to get cell-centered data
    procedure(subr_cc_data_${NDIM}$D) :: get_block_cc_data
    !> Simulation time
    real(dp), intent(in), optional :: time
    !> Which viewer (visit, paraview) will be used
    character(len=*), intent(in), optional :: viewer
    !> Whether to fill diagional ghost cells through extrapolation (default: true)
    logical, intent(in), optional :: fill_diagonals

    integer               :: mpirank, mpisize, ierr, n_super
    integer               :: coord_ix(NDIM)
    logical               :: extrap_diag
    character(len=20)     :: for_viewer

    ! Per super-block metadata (local to this rank)
    integer, allocatable  :: super_lvl(:)
    integer, allocatable  :: super_lo(:, :), super_hi(:, :)
    integer, allocatable  :: super_ghost_lo(:, :), super_ghost_hi(:, :)
    real(dp), allocatable :: super_dr(:, :), super_origin(:, :)
    integer(int64), allocatable :: super_offset(:)

    ! Hash table mapping block index -> local block id
    type(ffh_t)           :: h

    type(int_list_${NDIM}$d), allocatable :: super_block_ids(:)

    for_viewer = "visit"; if (present(viewer)) for_viewer = viewer
    extrap_diag = .true.; if (present(fill_diagonals)) extrap_diag = fill_diagonals

    call get_coord_permutation(for_viewer, coord_ix)

    call MPI_COMM_RANK(mpicomm, mpirank, ierr)
    call MPI_COMM_SIZE(mpicomm, mpisize, ierr)

    ! Build the super-blocks (rectangles of same-level blocks)
    call build_superblocks(h, n_blocks, nx, origin, dr, level, r_min, r_max, &
         out_gc, for_viewer, n_super, super_lvl, super_lo, super_hi, &
         super_ghost_lo, super_ghost_hi, super_dr, super_origin, super_block_ids)

    ! Write this rank's binary file and record byte offsets
    call write_superblock_data(filename, mpirank, nx, in_gc, out_gc, n_cc, &
         n_super, super_lo, super_hi, super_ghost_lo, super_ghost_hi, &
         super_block_ids, get_block_cc_data, super_offset)

    ! Gather metadata on rank 0 and write the XDMF header
    call write_xdmf_header(mpicomm, mpirank, mpisize, filename, nx, n_cc, &
         cc_names, coord_ix, n_super, super_lo, super_hi, super_ghost_lo, &
         super_ghost_hi, super_dr, super_origin, super_offset, time)

    call MPI_Barrier(mpicomm, ierr)

  end subroutine io_xdmf_write_blocks_${NDIM}$DCoRect

  !> Determine coordinate permutation for the given viewer
  subroutine get_coord_permutation(for_viewer, coord_ix)
    character(len=*), intent(in) :: for_viewer
    integer, intent(out)         :: coord_ix(NDIM)
    integer                      :: i

    select case (for_viewer)
    case ("visit")
       coord_ix(:) = [(i, i = 1, NDIM)]
    case ("paraview")
       coord_ix(:) = [(i, i = NDIM, 1, -1)]
    case default
       error stop "viewer can be: visit, paraview"
    end select
  end subroutine get_coord_permutation

  !> Group same-level blocks into rectangular super-blocks and determine their
  !> geometry and ghost-cell layout.
  subroutine build_superblocks(h, n_blocks, nx, origin, dr, level, r_min, &
       r_max, out_gc, for_viewer, n_super, super_lvl, super_lo, super_hi, &
       super_ghost_lo, super_ghost_hi, super_dr, super_origin, super_block_ids)
    type(ffh_t), intent(inout)         :: h
    integer, intent(in)                :: n_blocks, nx(NDIM), out_gc
    integer, intent(in)                :: level(n_blocks)
    real(dp), intent(in)               :: origin(NDIM, n_blocks)
    real(dp), intent(in)               :: dr(NDIM, n_blocks)
    real(dp), intent(in)               :: r_min(NDIM), r_max(NDIM)
    character(len=*), intent(in)       :: for_viewer
    integer, intent(out)               :: n_super
    integer, allocatable, intent(out)  :: super_lvl(:)
    integer, allocatable, intent(out)  :: super_lo(:, :), super_hi(:, :)
    integer, allocatable, intent(out)  :: super_ghost_lo(:, :), super_ghost_hi(:, :)
    real(dp), allocatable, intent(out) :: super_dr(:, :), super_origin(:, :)
    !> For each super-block, the list of block indices it contains, in the
    !> same KJI loop order used when writing the data.
    type(int_list_${NDIM}$d), allocatable, intent(out) :: super_block_ids(:)

    integer, allocatable :: super_id(:), blk_ix(:, :)
    integer, allocatable :: seed_id(:), tmp_lo(:, :), tmp_hi(:, :)
    integer              :: n, dim, seed, lvl, i_block, ix, status, cnt
    integer              :: lo(NDIM), hi(NDIM), test_lo(NDIM), test_hi(NDIM)
    integer              :: bx(NDIM), ghost_lo(NDIM), ghost_hi(NDIM), ${IJK}$
    real(dp)             :: r0(NDIM), r1(NDIM)
    logical              :: bnd_lo(NDIM), bnd_hi(NDIM)
    type(key_t)          :: key

    allocate(super_id(n_blocks))
    allocate(blk_ix(NDIM, n_blocks))
    allocate(seed_id(n_blocks))
    allocate(tmp_lo(NDIM, n_blocks))
    allocate(tmp_hi(NDIM, n_blocks))
    super_id = 0

    do n = 1, n_blocks
       do dim = 1, NDIM
          blk_ix(dim, n) = nint((origin(dim, n) - r_min(dim)) / &
               (dr(dim, n) * nx(dim)))
       end do
       key%x = [level(n), blk_ix(:, n)]
       call h%store_value(key, n, ix, existing_key_is_error=.true.)
    end do

    ! First pass: grow super-blocks (rectangles of same-level blocks)
    n_super = 0
    do seed = 1, n_blocks
       if (super_id(seed) /= 0) cycle

       n_super = n_super + 1
       lvl = level(seed)
       lo  = blk_ix(:, seed)
       hi  = blk_ix(:, seed)

       ! Try to extend the rectangle one dimension at a time
       do dim = 1, NDIM
          ! Extend in +dim direction
          do
             test_lo = lo; test_hi = hi
             test_lo(dim) = hi(dim) + 1
             test_hi(dim) = hi(dim) + 1
             if (.not. slab_available(h, lvl, test_lo, test_hi, super_id)) exit
             hi(dim) = test_hi(dim)
          end do
          ! Extend in -dim direction
          do
             test_lo = lo; test_hi = hi
             test_lo(dim) = lo(dim) - 1
             test_hi(dim) = lo(dim) - 1
             if (.not. slab_available(h, lvl, test_lo, test_hi, super_id)) exit
             lo(dim) = test_lo(dim)
          end do
       end do

       call mark_slab(h, lvl, lo, hi, n_super, super_id)
       seed_id(n_super) = seed
       tmp_lo(:, n_super) = lo
       tmp_hi(:, n_super) = hi
    end do

    deallocate(super_id, blk_ix)

    ! Now n_super is known: allocate outputs
    allocate(super_lvl(n_super))
    allocate(super_origin(NDIM, n_super))
    allocate(super_lo(NDIM, n_super))
    allocate(super_hi(NDIM, n_super))
    allocate(super_ghost_lo(NDIM, n_super))
    allocate(super_ghost_hi(NDIM, n_super))
    allocate(super_dr(NDIM, n_super))
    allocate(super_block_ids(n_super))

    ! Second pass: fill geometry, ghost-cell layout and block id lists
    do n = 1, n_super
       super_lo(:, n) = tmp_lo(:, n)
       super_hi(:, n) = tmp_hi(:, n)
       super_dr(:, n) = dr(:, seed_id(n))
       super_lvl(n)   = level(seed_id(n))

       key%x = [level(seed_id(n)), tmp_lo(:, n)]
       call h%get_value(key, i_block, status)
       super_origin(:, n) = origin(:, i_block)

       ! Store the contained block ids in the same order used when writing
       lo = tmp_lo(:, n)
       hi = tmp_hi(:, n)
       allocate(super_block_ids(n)%id(product(hi - lo + 1)))
       cnt = 0
       do @{KJI_LOOP_array_to_array(lo, hi)}@
          key%x = [super_lvl(n), ${IJK}$]
          call h%get_value(key, i_block, status)
          cnt = cnt + 1
          super_block_ids(n)%id(cnt) = i_block
       end do; ${KJI_CLOSE_LOOP}$

       ! Determine number of ghost cells on each side
       bx = (tmp_hi(:, n) - tmp_lo(:, n) + 1) * nx
       r0 = super_origin(:, n)
       r1 = r0 + super_dr(:, n) * bx
       call check_boundary(r0, r1, r_min, r_max, bnd_lo, bnd_hi)

       if (for_viewer == "visit") then
          ghost_lo = out_gc
          ghost_hi = out_gc
          where (bnd_lo) ghost_lo = 0
          where (bnd_hi) ghost_hi = 0
       else
          ghost_lo = 0
          ghost_hi = 0
       end if

       super_ghost_lo(:, n) = ghost_lo
       super_ghost_hi(:, n) = ghost_hi
    end do

  end subroutine build_superblocks

  !> Assemble super-block data and write it to this rank's binary file. The
  !> starting byte offset of each super-block is recorded in super_offset.
  subroutine write_superblock_data(filename, mpirank, nx, in_gc, out_gc, &
       n_cc, n_super, super_lo, super_hi, super_ghost_lo, &
       super_ghost_hi, super_block_ids, get_block_cc_data, super_offset)
    character(len=*), intent(in)       :: filename
    integer, intent(in)                :: mpirank, nx(NDIM), in_gc, out_gc, n_cc
    integer, intent(in)                :: n_super
    integer, intent(in)                :: super_lo(:, :), super_hi(:, :)
    integer, intent(in)                :: super_ghost_lo(:, :), super_ghost_hi(:, :)
    type(int_list_${NDIM}$d), intent(in) :: super_block_ids(:)
    procedure(subr_cc_data_${NDIM}$D)  :: get_block_cc_data
    integer(int64), allocatable, intent(out) :: super_offset(:)

    character(len=len_trim(filename)+20) :: binary_fname
    integer                :: my_unit, n, i_block, cnt
    integer                :: ilo(NDIM), ihi(NDIM), jlo(NDIM), jhi(NDIM)
    integer                :: lo(NDIM), hi(NDIM), bx(NDIM)
    integer                :: ghost_lo(NDIM), ghost_hi(NDIM), ${IJK}$
    integer(int64)         :: cur_offset
#:if NDIM == 2
    real(fp), allocatable  :: cc_block(:, :, :), cc_super(:, :, :)
#:elif NDIM == 3
    real(fp), allocatable  :: cc_block(:, :, :, :), cc_super(:, :, :, :)
#:endif

    allocate(super_offset(n_super))

    call get_fname_rank(trim(filename), '.bin', mpirank, binary_fname)
    open(newunit=my_unit, file=trim(binary_fname), form='unformatted', &
         access='stream', status='replace')

#:if NDIM == 2
    allocate(cc_block(-in_gc+1:nx(1)+in_gc, -in_gc+1:nx(2)+in_gc, n_cc))
#:elif NDIM == 3
    allocate(cc_block(-in_gc+1:nx(1)+in_gc, -in_gc+1:nx(2)+in_gc, &
         -in_gc+1:nx(3)+in_gc, n_cc))
#:endif

    cur_offset = 0

    do n = 1, n_super
       ghost_lo = super_ghost_lo(:, n)
       ghost_hi = super_ghost_hi(:, n)
       bx = (super_hi(:, n) - super_lo(:, n) + 1) * nx

#:if NDIM == 2
       allocate(cc_super(-ghost_lo(1)+1:bx(1)+ghost_hi(1), &
            -ghost_lo(2)+1:bx(2)+ghost_hi(2), n_cc))
#:elif NDIM == 3
       allocate(cc_super(-ghost_lo(1)+1:bx(1)+ghost_hi(1), &
            -ghost_lo(2)+1:bx(2)+ghost_hi(2), &
            -ghost_lo(3)+1:bx(3)+ghost_hi(3), n_cc))
#:endif
       lo = super_lo(:, n)
       hi = super_hi(:, n)

       cnt = 0
       do @{KJI_LOOP_array_to_array(lo, hi)}@
          cnt = cnt + 1
          i_block = super_block_ids(n)%id(cnt)
          call get_block_cc_data(i_block, cc_block)

          ! A bit inefficient, but simpler than doing it on the superblock
          call fill_diagonal_gc(nx, in_gc, out_gc, n_cc, cc_block)

          ! Index on super-block
          ilo = ([${IJK}$] - lo) * nx + 1
          ihi = ilo + nx - 1

          ! Index on block
          jlo = 1
          jhi = nx

          where ([${IJK}$] == lo)
             ilo = ilo - ghost_lo
             jlo = jlo - ghost_lo
          end where

          where ([${IJK}$] == hi)
             ihi = ihi + ghost_hi
             jhi = jhi + ghost_hi
          end where

#:if NDIM == 2
          cc_super(ilo(1):ihi(1), ilo(2):ihi(2), :) = &
               cc_block(jlo(1):jhi(1), jlo(2):jhi(2), :)
#:elif NDIM == 3
          cc_super(ilo(1):ihi(1), ilo(2):ihi(2), ilo(3):ihi(3), :) = &
               cc_block(jlo(1):jhi(1), jlo(2):jhi(2), jlo(3):jhi(3), :)
#:endif
       end do; ${KJI_CLOSE_LOOP}$

       super_offset(n) = cur_offset
       write(my_unit) cc_super
       cur_offset = cur_offset + &
            size(cc_super, kind=int64) * int(bytes_per_val, int64)
       deallocate(cc_super)
    end do

    close(my_unit)

  end subroutine write_superblock_data

  !> Gather super-block metadata onto rank 0 and write the XDMF header file.
  subroutine write_xdmf_header(mpicomm, mpirank, mpisize, filename, nx, n_cc, &
       cc_names, coord_ix, n_super, super_lo, super_hi, super_ghost_lo, &
       super_ghost_hi, super_dr, super_origin, super_offset, time)
    type(MPI_comm), intent(in)   :: mpicomm
    integer, intent(in)          :: mpirank, mpisize, nx(NDIM), n_cc
    character(len=*), intent(in) :: filename, cc_names(n_cc)
    integer, intent(in)          :: coord_ix(NDIM), n_super
    integer, intent(in)          :: super_lo(:, :), super_hi(:, :)
    integer, intent(in)          :: super_ghost_lo(:, :), super_ghost_hi(:, :)
    real(dp), intent(in)         :: super_dr(:, :), super_origin(:, :)
    integer(int64), intent(in)   :: super_offset(:)
    real(dp), intent(in), optional :: time

    integer                :: my_unit, ierr, rank, n
    integer                :: n_total, idx
    integer, allocatable   :: blocks_per_rank(:), displ(:)

    ! Gathered metadata, shaped for convenient indexing:
    !   g_bx(NDIM, n), g_ghost(NDIM, 2, n): (:,1,:)=lo (:,2,:)=hi
    real(dp), allocatable  :: g_origin(:, :), g_dr(:, :)
    integer, allocatable   :: g_bx(:, :), g_ghost(:, :, :)
    integer(int64), allocatable :: g_offset(:)

    ! Local send arrays (this rank's super-blocks)
    integer, allocatable   :: bx_send(:, :), ghost_send(:, :, :)
    character(len=len_trim(filename)+20) :: binary_fname
    character(len=len_trim(filename)+20) :: binary_basename

    allocate(blocks_per_rank(0:mpisize-1))
    call MPI_ALLGATHER(n_super, 1, MPI_INTEGER, blocks_per_rank, 1, &
         MPI_INTEGER, mpicomm, ierr)

    n_total = sum(blocks_per_rank)

    ! Displacements (in units of super-blocks) for Gatherv
    allocate(displ(0:mpisize-1))
    displ(0) = 0
    do rank = 1, mpisize-1
       displ(rank) = displ(rank-1) + blocks_per_rank(rank-1)
    end do

    ! Prepare local send arrays with the "reshaped" layout
    allocate(bx_send(NDIM, n_super))
    allocate(ghost_send(NDIM, 2, n_super))
    bx_send(:, 1:n_super)      = super_hi(:, 1:n_super) - super_lo(:, 1:n_super) + 1
    ghost_send(:, 1, 1:n_super) = super_ghost_lo(:, 1:n_super)
    ghost_send(:, 2, 1:n_super) = super_ghost_hi(:, 1:n_super)

    ! Allocate receive buffers on rank 0
    if (mpirank == 0) then
       allocate(g_origin(NDIM, n_total))
       allocate(g_dr(NDIM, n_total))
       allocate(g_bx(NDIM, n_total))
       allocate(g_ghost(NDIM, 2, n_total))
       allocate(g_offset(n_total))
    else
       allocate(g_origin(NDIM, 0), g_dr(NDIM, 0), g_bx(NDIM, 0))
       allocate(g_ghost(NDIM, 2, 0), g_offset(0))
    end if

    ! Gather each field. Counts/displacements are per super-block times the
    ! number of values per super-block for that field.
    call MPI_Gatherv(super_origin, NDIM*n_super, MPI_DOUBLE_PRECISION, &
         g_origin, NDIM*blocks_per_rank, NDIM*displ, MPI_DOUBLE_PRECISION, &
         0, mpicomm, ierr)
    call MPI_Gatherv(super_dr, NDIM*n_super, MPI_DOUBLE_PRECISION, &
         g_dr, NDIM*blocks_per_rank, NDIM*displ, MPI_DOUBLE_PRECISION, &
         0, mpicomm, ierr)
    call MPI_Gatherv(bx_send, NDIM*n_super, MPI_INTEGER, &
         g_bx, NDIM*blocks_per_rank, NDIM*displ, MPI_INTEGER, &
         0, mpicomm, ierr)
    call MPI_Gatherv(ghost_send, 2*NDIM*n_super, MPI_INTEGER, &
         g_ghost, 2*NDIM*blocks_per_rank, 2*NDIM*displ, MPI_INTEGER, &
         0, mpicomm, ierr)
    call MPI_Gatherv(super_offset, n_super, MPI_INTEGER8, &
         g_offset, blocks_per_rank, displ, MPI_INTEGER8, &
         0, mpicomm, ierr)

    if (mpirank == 0) then
       open(newunit=my_unit, file=trim(filename) // '.xdmf', action="write", &
            status='replace')

       write(my_unit, "(a)") '<?xml version="1.0" encoding="US-ASCII"?>'
       write(my_unit, "(a)") '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
       write(my_unit, "(a)") '<Xdmf Version="2.0">'
       write(my_unit, "(a)") '<Domain>'
       write(my_unit, "(a)") '<Grid Name="Mesh" GridType="Collection">'

       if (present(time)) then
          write(my_unit, *) '  <Time Value="', time, '" />'
       end if

       ! Write one <Grid> per super-block, grouped by rank so we know the
       ! corresponding binary file name.
       do rank = 0, mpisize-1
          call get_fname_rank(trim(filename), '.bin', rank, binary_fname)
          call get_basename(binary_fname, binary_basename)

          do n = 1, blocks_per_rank(rank)
             idx = displ(rank) + n
             call write_xdmf_grid_block(my_unit, trim(binary_basename), &
                  nx, n_cc, cc_names, coord_ix, idx, &
                  g_origin(:, idx), g_dr(:, idx), &
                  g_bx(:, idx), g_ghost(:, 1, idx), &
                  g_ghost(:, 2, idx), g_offset(idx))
          end do
       end do

       write(my_unit, "(a)") '</Grid>'
       write(my_unit, "(a)") '</Domain>'
       write(my_unit, "(a)") '</Xdmf>'
       close(my_unit)

       print *, "Wrote ", trim(filename) // '.xdmf'
    end if

  end subroutine write_xdmf_header

  !> Write the XDMF <Grid> element for a single super-block.
  subroutine write_xdmf_grid_block(my_unit, binary_basename, nx, n_cc, &
       cc_names, coord_ix, block_id, origin, dr, bx, ghost_lo, ghost_hi, offset)
    integer, intent(in)          :: my_unit
    character(len=*), intent(in) :: binary_basename
    integer, intent(in)          :: nx(NDIM), n_cc, coord_ix(NDIM), block_id
    character(len=*), intent(in) :: cc_names(n_cc)
    real(dp), intent(in)         :: origin(NDIM), dr(NDIM)
    integer, intent(in)          :: bx(NDIM), ghost_lo(NDIM), ghost_hi(NDIM)
    integer(int64), intent(in)   :: offset

    integer        :: iv
    integer        :: n_cells(NDIM)
    integer(int64) :: byte_off
    real(dp)       :: r0(NDIM)

    ! Number of cells to use for rendering (including ghost cells)
    n_cells = nx * bx + ghost_lo + ghost_hi

    ! Adjust origin to include ghost cells
    r0 = origin - ghost_lo * dr

    write(my_unit, "(a,I0,a)") &
         '  <Grid Name="MeshBlock', block_id, '" GridType="Uniform">'
#:if NDIM == 2
    write(my_unit, "(a,I0,' ',I0,a)") &
         '    <Topology TopologyType="2DCoRectMesh" Dimensions="', &
         n_cells(2)+1, n_cells(1)+1, '"/>'
    write(my_unit, "(a,4(I0,' ')a)") &
         '    <Information Name="GhostOffsets" Value="', &
         ghost_lo(2), ghost_hi(2), ghost_lo(1), ghost_hi(1), '"/>'
    write(my_unit, "(a)") &
         '    <Geometry GeometryType="ORIGIN_DXDY">'
#:elif NDIM == 3
    write(my_unit, "(a,I0,a,I0,' ',I0,' ',I0,a)") &
         '    <Topology TopologyType="', NDIM, 'DCoRectMesh" Dimensions="', &
         n_cells(3)+1, n_cells(2)+1, n_cells(1)+1, '"/>'
    write(my_unit, "(a,6(I0,' ')a)") &
         '    <Information Name="GhostOffsets" Value="', &
         ghost_lo(3), ghost_hi(3), ghost_lo(2), ghost_hi(2), &
         ghost_lo(1), ghost_hi(1), '"/>'
    write(my_unit, "(a)") &
         '    <Geometry GeometryType="ORIGIN_DXDYDZ">'
#:endif
    write(my_unit, "(a,I0,a)") '      <DataItem Dimensions="', NDIM, '">'
    write(my_unit, "(*(ES24.17))") r0(coord_ix)
    write(my_unit, *) '     </DataItem>'
    write(my_unit, "(a,I0,a)") '      <DataItem Dimensions="', NDIM, '">'
    write(my_unit, "(*(ES24.17))") dr(coord_ix)
    write(my_unit, *) '     </DataItem>'
    write(my_unit, "(a)") '    </Geometry>'

    ! Write cell-centered data
    do iv = 1, n_cc
       ! byte offset of variable iv within this super-block
       byte_off = offset + int(iv-1, int64) * &
            product(int(n_cells, int64)) * bytes_per_val

       write(my_unit, "(a,a,a)") '    <Attribute Name="', &
            trim(cc_names(iv)), '" Center="Cell">'
#:if NDIM == 2
       write(my_unit, "(a,I0,' ',I0,a,I0,a,I0,a)") &
            '      <DataItem Dimensions="', n_cells(2), n_cells(1), &
            '" Format="Binary" NumberType="Float" Precision="', &
            bytes_per_val, '" Seek="', byte_off, '">'
#:elif NDIM == 3
       write(my_unit, "(a,I0,' ',I0,' ',I0,a,I0,a,I0,a)") &
            '      <DataItem Dimensions="', n_cells(3), n_cells(2), n_cells(1), &
            '" Format="Binary" NumberType="Float" Precision="', &
            bytes_per_val, '" Seek="', byte_off, '">'
#:endif
       write(my_unit, "(a)") trim(binary_basename)
       write(my_unit, "(a)") '      </DataItem>'
       write(my_unit, "(a)") '    </Attribute>'
    end do

    write(my_unit, "(a)") '  </Grid>'

  end subroutine write_xdmf_grid_block

  !> Fill diagonal/edge/corner ghost cells by local linear extrapolation.
  !> The face (side) ghost cells must already be filled.
  subroutine fill_diagonal_gc(nx, in_gc, out_gc, n_cc, cc)
    integer, intent(in)     :: nx(NDIM) !< Number of interior cells per dim
    integer, intent(in)     :: in_gc    !< Number of ghost layers in input
    integer, intent(in)     :: out_gc   !< Number of ghost layers to write
    integer, intent(in)     :: n_cc     !< Number of variables
#:if NDIM == 2
    real(fp), intent(inout) :: cc(-in_gc+1:nx(1)+in_gc, &
         -in_gc+1:nx(2)+in_gc, n_cc)
    integer                 :: c1, c2
#:elif NDIM == 3
    real(fp), intent(inout) :: cc(-in_gc+1:nx(1)+in_gc, &
         -in_gc+1:nx(2)+in_gc, -in_gc+1:nx(3)+in_gc, n_cc)
    integer                 :: c1, c2, c3, dim, n, o_dims(2)
    integer                 :: ia(NDIM), ib(NDIM), ic(NDIM)
#:endif
    integer                 :: g, ix(NDIM), di(NDIM)

    if (out_gc > in_gc) error stop "out_gc must be <= in_gc"

    do g = 1, out_gc

#:if NDIM == 2
       ! 2D: fill the four corner ghost regions
       do c2 = 0, 1
          do c1 = 0, 1
             ! Corner index for this ghost layer g
             ix(1) = merge(1 - g, nx(1) + g, c1 == 0)
             ix(2) = merge(1 - g, nx(2) + g, c2 == 0)
             ! Direction pointing back into the domain
             di(1) = merge(1, -1, c1 == 0)
             di(2) = merge(1, -1, c2 == 0)

             cc(ix(1), ix(2), :) = &
                  cc(ix(1)+di(1), ix(2), :) + &
                  cc(ix(1),       ix(2)+di(2), :) - &
                  cc(ix(1)+di(1), ix(2)+di(2), :)
          end do
       end do

#:elif NDIM == 3
       ! 3D: first fill the 12 edges, then the 8 corners

       ! Edges parallel to dimension dim
       do dim = 1, NDIM
          o_dims = [1 + mod(dim, NDIM), 1 + mod(dim + 1, NDIM)]

          do c2 = 0, 1
             do c1 = 0, 1
                di = 0
                ix = 0

                di(o_dims(1)) = merge(1, -1, c1 == 0)
                di(o_dims(2)) = merge(1, -1, c2 == 0)

                ix(o_dims(1)) = merge(1 - g, nx(o_dims(1)) + g, c1 == 0)
                ix(o_dims(2)) = merge(1 - g, nx(o_dims(2)) + g, c2 == 0)

                ia = ix; ia(o_dims(1)) = ia(o_dims(1)) + di(o_dims(1))
                ib = ix; ib(o_dims(2)) = ib(o_dims(2)) + di(o_dims(2))
                ic = ix + di

                do n = 1, nx(dim)
                   ix(dim) = n; ia(dim) = n; ib(dim) = n; ic(dim) = n
                   cc(ix(1), ix(2), ix(3), :) = &
                        cc(ia(1), ia(2), ia(3), :) + &
                        cc(ib(1), ib(2), ib(3), :) - &
                        cc(ic(1), ic(2), ic(3), :)
                end do
             end do
          end do
       end do

       ! Corners
       do c3 = 0, 1
          do c2 = 0, 1
             do c1 = 0, 1
                ! Corner index for this ghost layer g
                ix(1) = merge(1 - g, nx(1) + g, c1 == 0)
                ix(2) = merge(1 - g, nx(2) + g, c2 == 0)
                ix(3) = merge(1 - g, nx(3) + g, c3 == 0)
                ! Direction pointing back into the domain
                di(1) = merge(1, -1, c1 == 0)
                di(2) = merge(1, -1, c2 == 0)
                di(3) = merge(1, -1, c3 == 0)

                cc(ix(1), ix(2), ix(3), :) = &
                     cc(ix(1),       ix(2)+di(2), ix(3)+di(3), :) + &
                     cc(ix(1)+di(1), ix(2),       ix(3)+di(3), :) + &
                     cc(ix(1)+di(1), ix(2)+di(2), ix(3), :) - &
                     2 * cc(ix(1)+di(1), ix(2)+di(2), ix(3)+di(3), :)
             end do
          end do
       end do
#:endif
    end do

  end subroutine fill_diagonal_gc

  !> Check that the whole slab between lo and hi exists and is not yet
  !> assigned to a super-block
  logical function slab_available(h, lvl, lo, hi, ids)
    type(ffh_t), intent(inout) :: h
    integer, intent(in)        :: lvl, lo(NDIM), hi(NDIM)
    integer, intent(in)        :: ids(:)
    integer                    :: ${IJK}$, i_block, status
    type(key_t)                :: key

    slab_available = .false.
    do @{KJI_LOOP_array_to_array(lo, hi)}@
       key%x = [lvl, ${IJK}$]
       call h%get_value(key, i_block, status)
       if (status == -1) then
          return
       else if (ids(i_block) /= 0) then
          return
       end if
    end do; ${KJI_CLOSE_LOOP}$
    slab_available = .true.
  end function slab_available

  !> Mark all blocks in the slab from lo to hi as belonging to the current
  !> super-block.
  subroutine mark_slab(h, lvl, lo, hi, n_super, ids)
    type(ffh_t), intent(inout) :: h
    integer, intent(in)        :: lvl, lo(NDIM), hi(NDIM)
    integer, intent(in)        :: n_super
    integer, intent(inout)     :: ids(:)
    integer                    :: ${IJK}$, i_block, status
    type(key_t)                :: key

    do @{KJI_LOOP_array_to_array(lo, hi)}@
       key%x = [lvl, ${IJK}$]
       call h%get_value(key, i_block, status)
       ids(i_block) = n_super
    end do; ${KJI_CLOSE_LOOP}$
  end subroutine mark_slab

  !> Get basename of fullpath
  subroutine get_basename(fullpath, out_basename)
    character(len=*), intent(in)  :: fullpath
    character(len=*), intent(out) :: out_basename
    integer :: p
    p = scan(trim(fullpath), '/', back=.true.)
    out_basename = fullpath(p+1:len_trim(fullpath))
  end subroutine get_basename

  !> Return the name of the binary file for mpirank
  subroutine get_fname_rank(filename, extension, mpirank, fname)
    character(len=*), intent(in)    :: filename
    character(len=*), intent(in)    :: extension
    integer, intent(in)             :: mpirank
    character(len=*), intent(inout) :: fname

    write(fname, '(A,A,I6.6,A)') trim(filename), "_", mpirank, trim(extension)
  end subroutine get_fname_rank

  !> Check whether a block face coincides with a domain boundary. This is
  !> important for rendering ghost cells in Visit.
  subroutine check_boundary(r0, r1, r0_domain, r1_domain, bnd_lo, bnd_hi)
    real(dp), intent(in) :: r0(NDIM), r1(NDIM)
    real(dp), intent(in) :: r0_domain(NDIM), r1_domain(NDIM)
    logical, intent(out) :: bnd_lo(NDIM), bnd_hi(NDIM)
    real(dp)             :: dmax(NDIM)

    dmax = bnd_rel_tol * maxval(r1_domain - r0_domain)
    bnd_lo = abs(r0 - r0_domain) < dmax
    bnd_hi = abs(r1 - r1_domain) < dmax
  end subroutine check_boundary

  pure logical function keys_equal(a, b)
    type(key_t), intent(in) :: a, b
    keys_equal = all(a%x == b%x)
  end function keys_equal

  !> Simple polynomial hash of the integer key array
  pure integer function ffh_hash(key, mask) result(i)
    type(key_t), intent(in) :: key
    integer, intent(in)     :: mask
    integer                 :: k
    integer(int64)          :: h
    integer(int64), parameter :: prime = 1000003_int64

    h = 1469598903_int64
    do k = 1, size(key%x)
       h = ieor(h, int(key%x(k), int64))
       h = h * prime
    end do
    i = iand(int(iand(h, int(huge(1), int64)), kind=4), mask)
  end function ffh_hash

  !> Resize (and rehash) the table to n_new buckets (power of two)
  subroutine ffh_resize(h, n_new)
    type(ffh_t), intent(inout) :: h
    integer, intent(in)        :: n_new
    integer, allocatable       :: old_flags(:), old_vals(:)
    type(key_t), allocatable   :: old_keys(:)
    integer                    :: i, j, step, old_n

    old_n = h%n_buckets
    call move_alloc(h%flags, old_flags)
    call move_alloc(h%keys,  old_keys)
    call move_alloc(h%vals,  old_vals)

    allocate(h%flags(0:n_new-1), h%keys(0:n_new-1), h%vals(0:n_new-1))
    h%flags        = 0
    h%n_buckets    = n_new
    h%hash_mask    = n_new - 1
    h%n_keys_stored = 0

    do j = 0, old_n-1
       if (old_flags(j) == 1) then
          i = ffh_hash(old_keys(j), h%hash_mask)
          do step = 1, h%n_buckets
             if (h%flags(i) == 0) exit
             i = iand(i + step, h%hash_mask)
          end do
          h%flags(i)     = 1
          h%keys(i)      = old_keys(j)
          h%vals(i)      = old_vals(j)
          h%n_keys_stored = h%n_keys_stored + 1
       end if
    end do
  end subroutine ffh_resize

  !> Store a key-value pair. existing_key_is_error triggers an error stop.
  subroutine ffh_store_value(h, key, val, ix, existing_key_is_error)
    class(ffh_t), intent(inout)   :: h
    type(key_t), intent(in)       :: key
    integer, intent(in)           :: val
    integer, intent(out)          :: ix
    logical, intent(in), optional :: existing_key_is_error
    integer                       :: i, step
    logical                       :: err_if_exists

    err_if_exists = .false.
    if (present(existing_key_is_error)) err_if_exists = existing_key_is_error

    ! Grow if load factor > 0.7
    if (h%n_buckets == 0) then
       call ffh_resize(h, 64)
    else if (h%n_keys_stored*10 >= h%n_buckets*7) then
       call ffh_resize(h, 2*h%n_buckets)
    end if

    i = ffh_hash(key, h%hash_mask)
    do step = 1, h%n_buckets
       if (h%flags(i) == 0) exit
       if (keys_equal(h%keys(i), key)) then
          if (err_if_exists) error stop "ffh: key already present"
          ix = i
          return
       end if
       i = iand(i + step, h%hash_mask)
    end do

    h%flags(i)      = 1
    h%keys(i)       = key
    h%vals(i)       = val
    h%n_keys_stored = h%n_keys_stored + 1
    ix = i
  end subroutine ffh_store_value

  !> Get value for a key. status = -1 if not found, else the index (>= 0).
  pure subroutine ffh_get_value(h, key, val, status)
    class(ffh_t), intent(in) :: h
    type(key_t), intent(in)  :: key
    integer, intent(inout)   :: val
    integer, intent(out)     :: status
    integer                  :: i, step

    status = -1
    if (h%n_buckets == 0) return

    i = ffh_hash(key, h%hash_mask)
    do step = 1, h%n_buckets
       if (h%flags(i) == 0) return           ! empty -> not found
       if (keys_equal(h%keys(i), key)) then
          val    = h%vals(i)
          status = i
          return
       end if
       i = iand(i + step, h%hash_mask)
    end do
  end subroutine ffh_get_value

end module m_io_${NDIM}$d
