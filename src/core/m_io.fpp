!> Module for writing output
!>
!> The XDMF part was inspired by the implementation in the bhac code, see
!> src/amrvacio/amrio.t in
!> https://gitlab.itp.uni-frankfurt.de/BHAC-release/bhac
#:include 'definitions_ndim.fpp'
#:include 'definitions_parallel.fpp'
module m_io_${NDIM}$d
  use mpi_f08
  use m_io_hash_${NDIM}$d
  use, intrinsic :: iso_c_binding
  use iso_fortran_env, only: int64
  use m_foap4_types_${NDIM}$d

  implicit none
  private

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
    f4%wtime_write_grid = f4%wtime_write_grid + t1 - t0

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
    integer, parameter           :: NDIM = ${NDIM}$
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

    integer                              :: my_unit, n, m, iv, mpirank, mpisize, ierr
    integer                              :: rank, n_prev_blocks, coord_ix(NDIM)
    integer, allocatable                 :: blocks_per_rank(:)
#:if NDIM == 2
    real(fp), allocatable                :: cc_block(:, :, :), cc_super(:, :, :)
#:elif NDIM == 3
    real(fp), allocatable                :: cc_block(:, :, :, :), cc_super(:, :, :, :)
#:endif
    character(len=len_trim(filename)+20) :: binary_fname, binary_basename
    character(len=20)                    :: for_viewer
    logical                              :: extrap_diag
    integer                              :: tag, status
    real(dp), allocatable                :: dr_recvbuf(:), dr_sendbuf(:)
    real(dp), allocatable                :: origin_recvbuf(:), origin_sendbuf(:)
    integer, allocatable                 :: bx_recvbuf(:), bx_sendbuf(:)
    integer, allocatable                 :: ghost_recvbuf(:), ghost_sendbuf(:)
    integer(int64), allocatable          :: offset_sendbuf(:), offset_recvbuf(:)
    integer(int64)                       :: byte_off
    type(mpi_request)                    :: requests(5)
    real(dp)                             :: r0(NDIM), r1(NDIM)
    logical                              :: bnd_lo(NDIM), bnd_hi(NDIM)
    integer                              :: ix_lo(NDIM), n_cells(NDIM)
    integer                              :: ghost_lo(NDIM), ghost_hi(NDIM)

    integer, allocatable :: blk_ix(:, :)   ! integer index (NDIM) of each block
    integer, allocatable :: super_id(:)    ! super-block id of each block (0 = none)
    integer, allocatable :: super_lvl(:)   ! level of super-block
    integer, allocatable :: super_lo(:, :) ! lower index of super-block
    integer, allocatable :: super_hi(:, :) ! upper index of super-block
    integer, allocatable :: super_ghost_lo(:, :) ! Num. ghost cells on lower side
    integer, allocatable :: super_ghost_hi(:, :) ! Num. ghost cells on upper side
    real(dp), allocatable :: super_dr(:, :) ! grid spacing of super-block
    real(dp), allocatable :: super_origin(:, :) ! origin of super-block
    integer(int64), allocatable :: super_offset(:)  ! byte offset per super-block
    integer(int64) :: cur_offset
    integer              :: n_super, ix, ilo(NDIM), ihi(NDIM)
    integer              :: jlo(NDIM), jhi(NDIM)
    integer              :: lvl, dim, seed, bx(NDIM), i_block, ${IJK}$
    integer              :: lo(NDIM), hi(NDIM), test_lo(NDIM), test_hi(NDIM)
    type(ffh_t)          :: h
    type(key_t)          :: key

    for_viewer = "visit"; if (present(viewer)) for_viewer = viewer
    extrap_diag = .true.; if (present(fill_diagonals)) extrap_diag = fill_diagonals

    select case (for_viewer)
    case ("visit")
       coord_ix(:) = [(i, i = 1, NDIM)]
    case ("paraview")
       coord_ix(:) = [(i, i = NDIM, 1, -1)]
    case default
       error stop "viewer can be: visit, paraview"
    end select

    call MPI_COMM_RANK(mpicomm, mpirank, ierr)
    call MPI_COMM_SIZE(mpicomm, mpisize, ierr)

    allocate(super_id(n_blocks))
    allocate(super_lvl(n_blocks))
    allocate(super_origin(NDIM, n_blocks))
    allocate(super_lo(NDIM, n_blocks))
    allocate(super_hi(NDIM, n_blocks))
    allocate(super_ghost_lo(NDIM, n_blocks))
    allocate(super_ghost_hi(NDIM, n_blocks))
    allocate(super_dr(NDIM, n_blocks))
    allocate(blk_ix(NDIM, n_blocks))
    cur_offset = 0
    super_id = 0

    do n = 1, n_blocks
       ! Determine integer index along each dimension
       do dim = 1, NDIM
          blk_ix(dim, n) = nint((origin(dim, n) - r_min(dim)) / &
               (dr(dim, n) * nx(dim)))
       end do
       key%x = [level(n), blk_ix(:, n)]
       call h%store_value(key, n, ix, existing_key_is_error=.true.)
    end do

    ! Grow super-blocks (rectangles of same-level blocks)
    n_super = 0
    do seed = 1, n_blocks
       if (super_id(seed) /= 0) cycle

       n_super = n_super + 1
       lvl = level(seed)
       lo  = blk_ix(:, seed)
       hi  = blk_ix(:, seed)

       ! Try to extend the rectangle one dimension at a time
       do dim = 1, NDIM
          test_lo = lo
          test_hi = hi

          ! Extend in +dim direction
          do
             test_lo(dim) = hi(dim) + 1
             test_hi(dim) = hi(dim) + 1
             if (.not. slab_available(h, lvl, test_lo, test_hi, &
                  n_blocks, super_id)) exit
             hi(dim) = test_hi(dim)
          end do

          ! Extend in -dim direction
          do
             test_lo(dim) = lo(dim) - 1
             test_hi(dim) = lo(dim) - 1
             if (.not. slab_available(h, lvl, test_lo, test_hi, &
                  n_blocks, super_id)) exit
             lo(dim) = test_lo(dim)
          end do
       end do

       call mark_slab(h, lvl, lo, hi, n_super, n_blocks, super_id)
       super_lo(:, n_super) = lo
       super_hi(:, n_super) = hi
       super_dr(:, n_super) = dr(:, seed)
       super_lvl(n_super) = level(seed)

       ! Get origin
       key%x = [lvl, lo]
       call h%get_value(key, i_block, status)
       super_origin(:, n_super) = origin(:, i_block)

       ! Determine number of ghost cells on each side
       bx = (super_hi(:, n_super) - super_lo(:, n_super) + 1) * nx
       r0 = super_origin(:, n_super)
       r1 = r0 + super_dr(:, n_super) * bx
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
       super_ghost_lo(:, n_super) = ghost_lo
       super_ghost_hi(:, n_super) = ghost_hi
    end do

    allocate(blocks_per_rank(0:mpisize-1))

    call MPI_ALLGATHER(n_super, 1, MPI_INTEGER, blocks_per_rank, 1, &
         MPI_INTEGER, mpicomm, ierr)

    ! Write binary file
    call get_fname_rank(trim(filename), '.bin', mpirank, binary_fname)
    open(newunit=my_unit, file=trim(binary_fname), form='unformatted', &
         access='stream', status='replace')

#:if NDIM == 2
    allocate(cc_block(-in_gc+1:nx(1)+in_gc, -in_gc+1:nx(2)+in_gc, n_cc))
#:elif NDIM == 3
    allocate(cc_block(-in_gc+1:nx(1)+in_gc, -in_gc+1:nx(2)+in_gc, &
         -in_gc+1:nx(3)+in_gc, n_cc))
#:endif

    allocate(super_offset(n_super))

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

       do @{KJI_LOOP_array_to_array(lo, hi)}@
          key%x = [super_lvl(n), ${IJK}$]
          call h%get_value(key, i_block, status)
          call get_block_cc_data(i_block, cc_block)
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
       cur_offset = cur_offset + size(cc_super, kind=int64) * &
            (storage_size(1.0_fp)/8)
       deallocate(cc_super)
    end do

    close(my_unit)

    if (mpirank == 0) then
       ! Write header
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
    end if

    tag = 0
    n_prev_blocks = 0

    do rank = 0, mpisize-1

       if (mpirank == rank) then
          allocate(origin_sendbuf(NDIM*n_super))
          allocate(dr_sendbuf(NDIM*n_super))
          allocate(bx_sendbuf(NDIM*n_super))
          allocate(ghost_sendbuf(2*NDIM*n_super))
          allocate(offset_sendbuf(n_super))
          origin_sendbuf(:) = pack(super_origin(:, 1:n_super), .true.)
          dr_sendbuf(:) = pack(super_dr(:, 1:n_super), .true.)
          bx_sendbuf(:) = pack((super_hi(:, 1:n_super) - &
               super_lo(:, 1:n_super) + 1), .true.)
          ghost_sendbuf(1:NDIM*n_super) = pack(super_ghost_lo(:, 1:n_super), .true.)
          ghost_sendbuf(NDIM*n_super+1:) = pack(super_ghost_hi(:, 1:n_super), .true.)
          offset_sendbuf = super_offset(1:n_super)

          call MPI_Isend(origin_sendbuf, NDIM*n_super, &
               MPI_DOUBLE_PRECISION, 0, tag, mpicomm, requests(1), ierr)
          call MPI_Isend(dr_sendbuf, NDIM*n_super, &
               MPI_DOUBLE_PRECISION, 0, tag+1, mpicomm, requests(2), ierr)
          call MPI_Isend(bx_sendbuf, NDIM*n_super, &
               MPI_INTEGER, 0, tag+2, mpicomm, requests(3), ierr)
          call MPI_Isend(ghost_sendbuf, 2*NDIM*n_super, &
               MPI_INTEGER, 0, tag+3, mpicomm, requests(4), ierr)
          call MPI_Isend(offset_sendbuf, n_super, &
               MPI_INTEGER8, 0, tag+4, mpicomm, requests(5), ierr)
          call MPI_Waitall(5, requests, MPI_STATUSES_IGNORE, ierr)
          deallocate(origin_sendbuf, dr_sendbuf, bx_sendbuf, &
               ghost_sendbuf, offset_sendbuf)
       end if

       if (mpirank == 0) then
          allocate(origin_recvbuf(NDIM*blocks_per_rank(rank)))
          allocate(dr_recvbuf(NDIM*blocks_per_rank(rank)))
          allocate(bx_recvbuf(NDIM*blocks_per_rank(rank)))
          allocate(ghost_recvbuf(2*NDIM*blocks_per_rank(rank)))
          allocate(offset_recvbuf(blocks_per_rank(rank)))

          call MPI_Recv(origin_recvbuf, NDIM*blocks_per_rank(rank), &
               MPI_DOUBLE_PRECISION, rank, tag, mpicomm, MPI_STATUS_IGNORE, ierr)
          call MPI_Recv(dr_recvbuf, NDIM*blocks_per_rank(rank), &
               MPI_DOUBLE_PRECISION, rank, tag+1, mpicomm, MPI_STATUS_IGNORE, ierr)
          call MPI_Recv(bx_recvbuf, NDIM*blocks_per_rank(rank), &
               MPI_INTEGER, rank, tag+2, mpicomm, MPI_STATUS_IGNORE, ierr)
          call MPI_Recv(ghost_recvbuf, 2*NDIM*blocks_per_rank(rank), &
               MPI_INTEGER, rank, tag+3, mpicomm, MPI_STATUS_IGNORE, ierr)
          call MPI_Recv(offset_recvbuf, blocks_per_rank(rank), &
               MPI_INTEGER8, rank, tag+4, mpicomm, MPI_STATUS_IGNORE, ierr)

          ! Get name corresponding to this rank
          call get_fname_rank(trim(filename), '.bin', rank, binary_fname)
          call get_basename(binary_fname, binary_basename)

          do n = 1, blocks_per_rank(rank)
             ghost_lo = ghost_recvbuf((n-1)*NDIM+1:n*NDIM)
             m = n + blocks_per_rank(rank)
             ghost_hi = ghost_recvbuf((m-1)*NDIM+1:m*NDIM)

             ix_lo = 0
             ! Number of cells to use for rendering
             n_cells = nx * bx_recvbuf((n-1)*NDIM+1:n*NDIM) + &
                  ghost_lo + ghost_hi

             ! Adjust origin to include ghost cells
             r0 = origin_recvbuf((n-1)*NDIM+1:n*NDIM)
             r0 = r0 - ghost_lo * dr_recvbuf((n-1)*NDIM+1:n*NDIM)

             write(my_unit, "(a,I0,a)") &
                  '  <Grid Name="MeshBlock', n + n_prev_blocks, &
                  '" GridType="Uniform">'
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
             write(my_unit, "(*(ES24.17))") dr_recvbuf((n-1)*NDIM + coord_ix)
             write(my_unit, *) '     </DataItem>'
             write(my_unit, "(a)") '    </Geometry>'

             ! Write cell-centered data
             do iv = 1, n_cc
                ! byte offset of variable iv within this super-block
                byte_off = offset_recvbuf(n) + int(iv-1, int64) * &
                     product(int(n_cells, int64)) * (storage_size(1.0_fp)/8)

                write(my_unit, "(a,a,a)") '    <Attribute Name="', &
                     trim(cc_names(iv)), '" Center="Cell">'
#:if NDIM == 2
                write(my_unit, "(a,I0,' ',I0,a,I0,a,I0,a)") &
                     '      <DataItem Dimensions="', n_cells(2), n_cells(1), &
                     '" Format="Binary" NumberType="Float" Precision="', &
                     storage_size(1.0_fp)/8, '" Seek="', byte_off, '">'
#:elif NDIM == 3
                write(my_unit, "(a,I0,' ',I0,' ',I0,a,I0,a,I0,a)") &
                     '      <DataItem Dimensions="', n_cells(3), n_cells(2), n_cells(1), &
                     '" Format="Binary" NumberType="Float" Precision="', &
                     storage_size(1.0_fp)/8, '" Seek="', byte_off, '">'
#:endif
                write(my_unit, "(a)") trim(binary_basename)
                write(my_unit, "(a)") '      </DataItem>'
                write(my_unit, "(a)") '    </Attribute>'
             end do

             write(my_unit, "(a)") '  </Grid>'
          end do

          n_prev_blocks = n_prev_blocks + blocks_per_rank(rank)
          deallocate(origin_recvbuf)
          deallocate(dr_recvbuf)
          deallocate(bx_recvbuf)
          deallocate(ghost_recvbuf)
          deallocate(offset_recvbuf)
       end if
    end do

    if (mpirank == 0) then
       ! Complete header
       write(my_unit, "(a)") '</Grid>'
       write(my_unit, "(a)") '</Domain>'
       write(my_unit, "(a)") '</Xdmf>'
       close(my_unit)

       print *, "Wrote ", trim(filename) // '.xdmf'
    end if

    call MPI_Barrier(mpicomm, ierr)

  end subroutine io_xdmf_write_blocks_${NDIM}$DCoRect

  !> Return the name of the binary file for mpirank
  subroutine get_fname_rank(filename, extension, mpirank, fname)
    character(len=*), intent(in)    :: filename
    character(len=*), intent(in)    :: extension
    integer, intent(in)             :: mpirank
    character(len=20)               :: suffix
    character(len=*), intent(inout) :: fname

    write(suffix, '(A,I06.6)') "_", mpirank
    fname = trim(filename) // trim(suffix) // trim(extension)
  end subroutine get_fname_rank

  !> Get basename of fullpath
  subroutine get_basename(fullpath, out_basename)
    character(len=*), intent(in) :: fullpath
    character(len=*), intent(out) :: out_basename
    integer :: i, last_slash, len_path

    len_path = len_trim(fullpath)
    last_slash = 0

    ! Find the position of the last slash '/'
    do i = len_path, 1, -1
       if (fullpath(i:i) == '/') then
          last_slash = i
          exit
       end if
    end do

    if (last_slash == 0) then
       ! No slash found, entire string is the basename
       out_basename = fullpath(1:len_path)
    else
       ! Extract substring after the last '/'
       out_basename = fullpath(last_slash+1:len_path)
    end if
  end subroutine get_basename

  !> Check whether a block face coincides with a domain boundary. This is
  !> important for rendering ghost cells in Visit.
  subroutine check_boundary(r0, r1, r0_domain, r1_domain, bnd_lo, bnd_hi)
    real(dp), intent(in) :: r0(NDIM), r1(NDIM)
    real(dp), intent(in) :: r0_domain(NDIM), r1_domain(NDIM)
    logical, intent(out) :: bnd_lo(NDIM), bnd_hi(NDIM)
    real(dp)             :: dmax
    integer              :: idim

    dmax = 1e-13_dp * maxval(r1_domain - r0_domain)

    do idim = 1, NDIM
       bnd_lo(idim) = abs(r0(idim) - r0_domain(idim)) < dmax
       bnd_hi(idim) = abs(r1(idim) - r1_domain(idim)) < dmax
    end do
  end subroutine check_boundary

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


  !> Check that the whole slab at position pos (perpendicular to dim,
  !> spanning lo:hi in the other dims) exists and is unassigned.
  logical function slab_available(h, lvl, lo, hi, n_blocks, ids)
    type(ffh_t), intent(inout) :: h
    integer, intent(in)        :: lvl, lo(NDIM), hi(NDIM)
    integer, intent(in)        :: n_blocks
    integer, intent(in)        :: ids(n_blocks)
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

  !> Mark all blocks in the slab as belonging to the current super-block.
  subroutine mark_slab(h, lvl, lo, hi, n_super, n_blocks, ids)
    type(ffh_t), intent(inout) :: h
    integer, intent(in)        :: lvl, lo(NDIM), hi(NDIM)
    integer, intent(in)        :: n_super
    integer, intent(in)        :: n_blocks
    integer, intent(inout)     :: ids(n_blocks)
    integer                    :: ${IJK}$, i_block, status
    type(key_t)                :: key

    do @{KJI_LOOP_array_to_array(lo, hi)}@
       key%x = [lvl, ${IJK}$]
       call h%get_value(key, i_block, status)
       ids(i_block) = n_super
    end do; ${KJI_CLOSE_LOOP}$
  end subroutine mark_slab

end module m_io_${NDIM}$d
