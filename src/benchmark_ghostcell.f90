program benchmark_gc

  use mpi_f08
  use m_foap4

  implicit none
  integer, parameter :: dp = kind(0.0d0)

  type(foap4_t), target :: f4
  integer               :: n_output, min_level, n_refine_steps
  integer               :: n_iterations, max_blocks, bx(2)
  logical               :: write_grid

  call f4_initialize(f4, "error")

  n_output       = 0
  n_iterations   = 1000
  min_level      = 6
  n_refine_steps = 0
  max_blocks     = 40000
  write_grid     = .false.
  bx(:)          = 32

  call benchmark_ghostcell(f4, bx, n_iterations, [1e-2_dp, 1e-2_dp], &
       min_level, n_refine_steps, max_blocks, write_grid, &
       "output/benchmark_gc", n_output)

  if (f4%mpirank == 0) call f4_print_wtime(f4)
  call f4_finalize(f4)

contains

  subroutine benchmark_ghostcell(f4, bx, n_iterations, refine_location, &
       min_level, n_refine_steps, max_blocks, write_grid, base_name, n_output)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: bx(2)
    integer, intent(in)          :: n_iterations
    real(dp), intent(in)         :: refine_location(2)
    integer, intent(in)          :: min_level
    integer, intent(in)          :: n_refine_steps
    integer, intent(in)          :: max_blocks
    logical, intent(in)          :: write_grid
    character(len=*), intent(in) :: base_name
    integer, intent(inout)       :: n_output
    integer, parameter           :: n_blocks_per_dim(2) = [1, 1]
    real(dp), parameter          :: block_length(2)     = [1.0_dp, 1.0_dp]
    integer, parameter           :: n_gc                = 1
    integer, parameter           :: n_vars              = 2
    character(len=20)            :: var_names(n_vars)   = ['rho', 'phi']
    logical, parameter           :: periodic(2)         = [.false., .false.]
    integer                      :: n, n_ghostcells, n_blocks_global
    real(dp)                     :: t0, t1
    real(dp)                     :: t_total, ghostcells_per_ns

    call f4_construct_brick(f4, n_blocks_per_dim, block_length, bx, n_gc, &
         n_vars, var_names, periodic, min_level, max_blocks)

    call set_init_cond(f4)
    call f4_update_ghostcells(f4, 2, [1, 2])

    do n = 1, n_refine_steps
       call set_refinement_flag(f4, refine_location)
       call f4_adjust_refinement(f4, .true.)
       call f4_update_ghostcells(f4, 2, [1, 2])
    end do

    if (write_grid) then
       n_output = n_output + 1
       call f4_write_grid(f4, base_name, n_output)
    end if

    t0 = MPI_Wtime()
    do n = 1, n_iterations
       call f4_update_ghostcells(f4, 2, [1, 2])
    end do
    t1 = MPI_Wtime()

    t_total = (t1 - t0)
    n_blocks_global = f4_get_num_global_blocks(f4)
    n_ghostcells = n_blocks_global * n_vars * 2 * n_gc * sum(f4%bx)
    ghostcells_per_ns = 1e-9_dp * n_iterations/t_total * n_ghostcells

    if (f4%mpirank == 0) then
       write(*, "(A,F14.3)") " Ghostcells/ns:        ", ghostcells_per_ns
       write(*, "(A,I14)")   " n_blocks_global:      ", n_blocks_global
       write(*, "(A,F14.3)") " Global mesh size (MB):", n_blocks_global * &
            n_vars * 0.5_dp**20 * product(f4%bx + 2 * f4%n_gc)
    end if

    call f4_destroy(f4)
  end subroutine benchmark_ghostcell

  pure real(dp) function rho_init(x, y)
    real(dp), intent(in) :: x, y
    rho_init = x + y
  end function rho_init

  pure real(dp) function phi_init(x, y)
    real(dp), intent(in) :: x, y
    phi_init = x - y
  end function phi_init

  subroutine set_init_cond(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, i, j
    real(dp)                     :: rr(2)

    !$acc parallel loop
    do n = 1, f4%n_blocks
       !$acc loop collapse(2) private(rr)
       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)
             rr = f4_cell_coord(f4, n, i, j)
             f4%uu(i, j, 1, n) = rho_init(rr(1), rr(2))
             f4%uu(i, j, 2, n) = phi_init(rr(1), rr(2))
          end do
       end do
    end do
  end subroutine set_init_cond

  subroutine set_refinement_flag(f4, r0)
    type(foap4_t), intent(inout) :: f4
    real(dp), intent(in)         :: r0(2)
    integer                      :: n, lvl
    real(dp)                     :: rmin(2), rmax(2)

    do n = 1, f4%n_blocks
       lvl = f4%block_level(n)
       rmin = f4%block_origin(:, n)
       rmax = rmin + f4%bx * f4%dr_level(:, lvl)

       if (all(r0 >= rmin .and. r0 <= rmax)) then
          f4%refinement_flags(n) = 1
       else
          f4%refinement_flags(n) = 0
       end if
    end do
  end subroutine set_refinement_flag

end program benchmark_gc
