!> Module for writing output
!>
!> The XDMF part was inspired by the implementation in the bhac code, see
!> src/amrvacio/amrio.t in
!> https://gitlab.itp.uni-frankfurt.de/BHAC-release/bhac
#:include 'definitions.fpp'
module m_io_${NDIM}$d
  use mpi_f08
  use, intrinsic :: iso_c_binding
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

    ! OpenACC - get the block data from the device
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
         f4%block_origin(:, 1:f4%n_blocks), dr, &
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
       cc_names, in_gc, out_gc, origin, dr, cc_data, get_block_cc_data, time, viewer)
    integer, parameter             :: NDIM = ${NDIM}$
    type(MPI_comm), intent(in)     :: mpicomm !< MPI communicator
    character(len=*), intent(in)   :: filename !< File name without extension
    integer, intent(in)            :: n_blocks !< Number of blocks
    integer, intent(in)            :: nx(NDIM) !< Size of the blocks (excl. ghost cells)
    integer, intent(in)            :: n_cc !< Number of variables
    character(len=*), intent(in)   :: cc_names(n_cc) !< Names of variables
    integer, intent(in)            :: in_gc !< Number of ghost cells in input
    integer, intent(in)            :: out_gc !< Number of ghost cells to write
    !> Origin of each block (incl. ghost cells)
    real(dp), intent(in)           :: origin(NDIM, n_blocks)
    real(dp), intent(in)           :: dr(NDIM, n_blocks) !< Grid spacing of each block
    !> Cell-centered data
#:if NDIM == 2
    real(fp), intent(in), optional :: cc_data(-in_gc+1:nx(1)+in_gc, &
         -in_gc+1:nx(2)+in_gc, n_cc, n_blocks)
#:elif NDIM == 3
    real(fp), intent(in), optional :: cc_data(-in_gc+1:nx(1)+in_gc, &
         -in_gc+1:nx(2)+in_gc, -in_gc+1:nx(3)+in_gc, n_cc, n_blocks)
#:endif
    !> Method to get cell-centered data
    procedure(subr_cc_data_${NDIM}$D), optional :: get_block_cc_data
    !> Simulation time
    real(dp), intent(in), optional :: time
    !> Which viewer (visit, paraview) will be used
    character(len=*), intent(in), optional :: viewer

    integer                              :: my_unit, n, iv, mpirank, mpisize, ierr
    integer                              :: rank, n_prev_blocks, coord_ix(NDIM)
    integer, allocatable                 :: blocks_per_rank(:)
#:if NDIM == 2
    real(fp), allocatable                :: cc_block(:, :, :)
#:elif NDIM == 3
    real(fp), allocatable                :: cc_block(:, :, :, :)
#:endif
    character(len=len_trim(filename)+20) :: binary_fname, binary_basename
    character(len=20)                    :: for_viewer
    integer                              :: i, tag
    real(dp), allocatable                :: dr_recvbuf(:), origin_recvbuf(:)
    real(dp), allocatable                :: dr_sendbuf(:), origin_sendbuf(:)
    type(mpi_request)                    :: requests(2)

    for_viewer = "visit"; if (present(viewer)) for_viewer = viewer

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
    allocate(blocks_per_rank(0:mpisize-1))

    blocks_per_rank = 0
    blocks_per_rank(mpirank) = n_blocks
    call MPI_ALLGATHER(n_blocks, 1, MPI_INTEGER, blocks_per_rank, 1, &
         MPI_INTEGER, mpicomm, ierr)

    ! Write binary file
    call get_fname_rank(trim(filename), '.bin', mpirank, binary_fname)
    open(newunit=my_unit, file=trim(binary_fname), form='unformatted', &
         access='stream', status='replace')

    if (present(cc_data)) then
#:if NDIM == 2
       write(my_unit) cc_data(-out_gc+1:nx(1)+out_gc, &
         -out_gc+1:nx(2)+out_gc, :, :)
#:elif NDIM == 3
       write(my_unit) cc_data(-out_gc+1:nx(1)+out_gc, &
         -out_gc+1:nx(2)+out_gc, -out_gc+1:nx(3)+out_gc, :, :)
#:endif
    else if (present(get_block_cc_data)) then
#:if NDIM == 2
       allocate(cc_block(-in_gc+1:nx(1)+in_gc, -in_gc+1:nx(2)+in_gc, n_cc))
#:elif NDIM == 3
       allocate(cc_block(-in_gc+1:nx(1)+in_gc, -in_gc+1:nx(2)+in_gc, &
            -in_gc+1:nx(2)+in_gc, n_cc))
#:endif

       do n = 1, n_blocks
          call get_block_cc_data(n, cc_block)
#:if NDIM == 2
          write(my_unit) cc_block(-out_gc+1:nx(1)+out_gc, &
               -out_gc+1:nx(2)+out_gc, :)
#:elif NDIM == 3
          write(my_unit) cc_block(-out_gc+1:nx(1)+out_gc, &
               -out_gc+1:nx(2)+out_gc, -out_gc+1:nx(3)+out_gc, :)
#:endif
       end do
    else
       error stop "Either cc_data or get_block_cc_data should be given"
    end if

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

       ! if (n_gc > 0) then
       !    write(my_unit, "(a,i0,a)") '  <GhostZones Value="', n_gc, '" />'
       ! end if
    end if

    tag = 0
    n_prev_blocks = 0

    do rank = 0, mpisize-1

       if (mpirank == rank) then
          allocate(origin_sendbuf(NDIM*n_blocks))
          allocate(dr_sendbuf(NDIM*n_blocks))
          origin_sendbuf(:) = pack(origin, .true.)
          dr_sendbuf(:) = pack(dr, .true.)

          call MPI_Isend(origin_sendbuf, NDIM*n_blocks, &
               MPI_DOUBLE_PRECISION, 0, tag, mpicomm, requests(1), ierr)
          call MPI_Isend(dr_sendbuf, NDIM*n_blocks, &
               MPI_DOUBLE_PRECISION, 0, tag, mpicomm, requests(2), ierr)
       end if

       if (mpirank == 0) then
          allocate(origin_recvbuf(NDIM*blocks_per_rank(rank)))
          allocate(dr_recvbuf(NDIM*blocks_per_rank(rank)))

          call MPI_Recv(origin_recvbuf, NDIM*blocks_per_rank(rank), &
               MPI_DOUBLE_PRECISION, rank, tag, mpicomm, MPI_STATUS_IGNORE, ierr)
          call MPI_Recv(dr_recvbuf, NDIM*blocks_per_rank(rank), &
               MPI_DOUBLE_PRECISION, rank, tag, mpicomm, MPI_STATUS_IGNORE, ierr)

          ! Get name corresponding to this rank
          call get_fname_rank(trim(filename), '.bin', rank, binary_fname)
          call get_basename(binary_fname, binary_basename)

          do n = 1, blocks_per_rank(rank)
             write(my_unit, "(a,I0,a)") &
                  '  <Grid Name="MeshBlock', n + n_prev_blocks, &
                  '" GridType="Uniform">'
#:if NDIM == 2
             write(my_unit, "(a,I0,a,I0,' ',I0,a)") &
                  '    <Topology TopologyType="', NDIM, 'DCoRectMesh" Dimensions="', &
                  nx(2)+1, nx(1)+1, '"/>'
             write(my_unit, "(a)") &
                  '    <Geometry GeometryType="ORIGIN_DXDY">'
#:elif NDIM == 3
             write(my_unit, "(a,I0,a,I0,' ',I0,' ',I0,a)") &
                  '    <Topology TopologyType="', NDIM, 'DCoRectMesh" Dimensions="', &
                  nx(3)+1, nx(2)+1, nx(1)+1, '"/>'
             write(my_unit, "(a)") &
                  '    <Geometry GeometryType="ORIGIN_DXDYDZ">'
#:endif
             write(my_unit, "(a,I0,a)") '      <DataItem Dimensions="', NDIM, '">'

#:if NDIM == 2
             write(my_unit, "(2ES24.17)") origin_recvbuf((n-1)*NDIM + coord_ix)
#:elif NDIM == 3
             write(my_unit, "(3ES24.17)") origin_recvbuf((n-1)*NDIM + coord_ix)
#:endif
             write(my_unit, *) '      </DataItem>'
             write(my_unit, "(a,I0,a)") '      <DataItem Dimensions="', NDIM, '">'
#:if NDIM == 2
             write(my_unit, "(2ES24.17)") dr_recvbuf((n-1)*NDIM + coord_ix)
#:elif NDIM == 3
             write(my_unit, "(3ES24.17)") dr_recvbuf((n-1)*NDIM + coord_ix)
#:endif
             write(my_unit, *) '      </DataItem>'
             write(my_unit, "(a)") '    </Geometry>'

             ! Write cell-centered data
             do iv = 1, n_cc
                write(my_unit, "(a,a,a)") '    <Attribute Name="', &
                     trim(cc_names(iv)), '" Center="Cell">'
#:if NDIM == 2
                write(my_unit, "(a,I0,a,I0,a)") &
                     '      <DataItem ItemType="HyperSlab" Dimensions="',&
                     nx(2), ' ', nx(1), '">'
                write(my_unit, "(a, 12(I0,' '),a)") &
                     '        <DataItem Dimensions="3 4"> ', &
                     n-1, iv-1, out_gc, out_gc, &            ! start
                     1, 1, 1, 1, &                       ! stride
                     1, 1, nx(2), nx(1), & ! count
                     '</DataItem>'
                write(my_unit, "(a, 4(I0,' '),a,I0,a)") &
                     '        <DataItem Dimensions="', n_blocks, n_cc, &
                     nx(2) + 2*out_gc, nx(1) + 2*out_gc, &
                     '" Format="Binary" NumberType="Float" Precision="', &
                storage_size(1.0_fp)/8, '">'
#:elif NDIM == 3
                write(my_unit, "(a,I0,a,I0,a,I0,a)") &
                     '      <DataItem ItemType="HyperSlab" Dimensions="',&
                     nx(3), ' ', nx(2), ' ', nx(1), '">'
                write(my_unit, "(a, 15(I0,' '),a)") &
                     '        <DataItem Dimensions="3 5"> ', &
                     n-1, iv-1, out_gc, out_gc, out_gc, & ! start
                     1, 1, 1, 1, 1, &               ! stride
                     1, 1, nx(3), nx(2), nx(1), & ! count
                     '</DataItem>'

                write(my_unit, "(a, 5(I0,' '),a,I0,a)") &
                     '        <DataItem Dimensions="', n_blocks, n_cc, &
                     nx(3) + 2*out_gc, nx(2) + 2*out_gc, nx(1) + 2*out_gc, &
                     '" Format="Binary" NumberType="Float" Precision="', &
                     storage_size(1.0_fp)/8, '">'
#:endif
                write(my_unit, "(a)") trim(binary_basename)
                write(my_unit, "(a)") '        </DataItem>'
                write(my_unit, "(a)") '      </DataItem>'
                write(my_unit, "(a)") '    </Attribute>'
             end do
             write(my_unit, "(a)") '  </Grid>'
          end do

          n_prev_blocks = n_prev_blocks + blocks_per_rank(rank)
          deallocate(origin_recvbuf)
          deallocate(dr_recvbuf)
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

end module m_io_${NDIM}$d
