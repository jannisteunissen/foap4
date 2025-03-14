module m_xdmf_writer

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  abstract interface
     !> Get cell-centered data for a grid block
     subroutine subr_cc_data(n, cc_data)
       import
       integer, intent(in)     :: n !< Index of block
       real(dp), intent(inout) :: cc_data(:, :, :)
     end subroutine subr_cc_data
  end interface

  ! Public methods
  public :: xdmf_write_blocks_2DCoRect

contains

  subroutine xdmf_write_blocks_2DCoRect(mpicomm, filename, n_blocks, nx, n_cc, &
       cc_names, n_gc, origin, dr, cc_data, get_block_cc_data, time, viewer)
    use mpi_f08
    type(MPI_comm), intent(in)  :: mpicomm !< MPI communicator
    character(len=*), intent(in)   :: filename            !< File name without extension
    integer, intent(in)            :: n_blocks            !< Number of blocks
    integer, intent(in)            :: nx(2)               !< Size of the blocks (incl. ghost cells)
    integer, intent(in)            :: n_cc               !< Number of variables
    character(len=*), intent(in)   :: cc_names(n_cc)    !< Names of variables
    integer, intent(in)            :: n_gc                !< Number of ghost cells
    real(dp), intent(in)           :: origin(2, n_blocks) !< Origin of each block (incl. ghost cells)
    real(dp), intent(in)           :: dr(2, n_blocks)     !< Grid spacing of each block
    real(dp), intent(in), optional :: cc_data(nx(1), nx(2), n_cc, n_blocks) !< Cell-centered data
    procedure(subr_cc_data), optional :: get_block_cc_data   !< Method to get cell-centered data
    real(dp), intent(in), optional :: time                !< Simulation time
    character(len=*), intent(in), optional :: viewer !< Which viewer (visit, paraview) will be used

    integer                              :: my_unit, n, iv, mpirank, mpisize, ierr
    integer                              :: rank, n_prev_blocks, coord_ix(2)
    integer, allocatable                 :: blocks_per_rank(:)
    real(dp), allocatable                :: cc_block(:, :, :)
    character(len=len_trim(filename)+10) :: binary_fname
    character(len=20)                    :: suffix, for_viewer

    for_viewer = "visit"; if (present(viewer)) for_viewer = viewer
    select case (for_viewer)
       case ("visit")
          coord_ix = [1, 2]
       case ("paraview")
          coord_ix = [2, 1]
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

    write(suffix, '(A,I04.4)') "_", mpirank

    ! Write binary file
    binary_fname = trim(filename) // trim(suffix) // '.bin'
    open(newunit=my_unit, file=trim(binary_fname), form='unformatted', &
         access='stream', status='replace')

    if (present(cc_data)) then
       write(my_unit) cc_data
    else if (present(get_block_cc_data)) then
       allocate(cc_block(nx(1), nx(2), n_cc))

       do n = 1, n_blocks
          call get_block_cc_data(n, cc_block)
          write(my_unit) cc_block
       end do
    else
       error stop "Either cc_data or get_block_cc_data should be given"
    end if

    close(my_unit)

    if (mpirank == 0) then
       ! Write header
       open(newunit=my_unit, file=trim(filename) // '.xdmf', action="write")

       write(my_unit, "(a)") '<?xml version="1.0" encoding="US-ASCII"?>'
       write(my_unit, "(a)") '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
       write(my_unit, "(a)") '<Xdmf Version="2.0">'
       write(my_unit, "(a)") '<Domain>'
       write(my_unit, "(a)") '<Grid Name="Mesh" GridType="Collection">'

       if (present(time)) then
          write(my_unit, "(a,E16.8,a)") '  <Time Value="', time, '" />'
       end if

       ! if (n_gc > 0) then
       !    write(my_unit, "(a,i0,a)") '  <GhostZones Value="', n_gc, '" />'
       ! end if

       close(my_unit)
    end if

    n_prev_blocks = 0

    do rank = 0, mpisize-1
       if (mpirank == rank) then
          open(newunit=my_unit, file=trim(filename) // '.xdmf', action="write", &
               access="append")

          do n = 1, n_blocks
             write(my_unit, "(a,I0,a)") &
                  '  <Grid Name="MeshBlock', n + n_prev_blocks, '" GridType="Uniform">'
             write(my_unit, "(a,I0,' ',I0,a)") &
                  '    <Topology TopologyType="2DCoRectMesh" Dimensions="', &
                  nx(2)+1-2*n_gc, nx(1)+1-2*n_gc, '"/>'
             write(my_unit, "(a)") &
                  '    <Geometry GeometryType="ORIGIN_DXDY">'
             write(my_unit, "(a,I0,a)") '      <DataItem Dimensions="', 2, '">'
             write(my_unit, *) origin(coord_ix, n)
             write(my_unit, *) '      </DataItem>'
             write(my_unit, "(a,I0,a)") '      <DataItem Dimensions="', 2, '">'
             write(my_unit, *) dr(coord_ix, n)
             write(my_unit, *) '      </DataItem>'
             write(my_unit, "(a)") '    </Geometry>'

             ! Write cell-centered data
             do iv = 1, n_cc
                write(my_unit, "(a,a,a)") '    <Attribute Name="', trim(cc_names(iv)), &
                     '" Center="Cell">'
                write(my_unit, "(a,I0,a,I0,a)") '      <DataItem ItemType="HyperSlab" Dimensions="',&
                     nx(2)-2*n_gc, ' ', nx(1)-2*n_gc, '">'
                write(my_unit, "(a, 12(I0,' '),a)") '        <DataItem Dimensions="3 4"> ', &
                     n-1, iv-1, n_gc, n_gc, &        ! start
                     1, 1, 1, 1, &             ! stride
                     1, 1, nx(2)-2*n_gc, nx(1)-2*n_gc, &     ! count
                     '</DataItem>'

                write(my_unit, "(a, 4(I0,' '),a,a,a)") &
                     '        <DataItem Dimensions="', n_blocks, n_cc, nx(2), nx(1), &
                     '" Format="Binary" NumberType="Float" Precision="8">'
                write(my_unit, "(a)") trim(binary_fname)
                write(my_unit, "(a)") '        </DataItem>'
                write(my_unit, "(a)") '      </DataItem>'
                write(my_unit, "(a)") '    </Attribute>'
             end do
             write(my_unit, "(a)") '  </Grid>'
          end do

          close(my_unit)
       end if

       call MPI_BARRIER(mpicomm, ierr)
       n_prev_blocks = n_prev_blocks + blocks_per_rank(rank)
    end do

    if (mpirank == 0) then
       ! Complete header
       open(newunit=my_unit, file=trim(filename) // '.xdmf', action="write", &
            access="append")
       write(my_unit, "(a)") '</Grid>'
       write(my_unit, "(a)") '</Domain>'
       write(my_unit, "(a)") '</Xdmf>'
       close(my_unit)
    end if

  end subroutine xdmf_write_blocks_2DCoRect

end module m_xdmf_writer
