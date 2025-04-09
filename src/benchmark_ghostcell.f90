program benchmark_gc

  use m_foap4

  implicit none
  integer, parameter :: dp = kind(0.0d0)

  type(foap4_t), target :: f4
  integer :: n_output

  call f4_initialize(f4, "error")

  n_output = 0
  call benchmark_ghostcell(f4, 10000, [1e-2_dp, 1e-2_dp], &
       "output/benchmark_gc", n_output)

  if (f4%mpirank == 0) call f4_print_wtime(f4)
  call f4_finalize(f4)

contains

  subroutine benchmark_ghostcell(f4, n_update_gc, refine_location, &
       base_name, n_output)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_update_gc
    real(dp), intent(in)         :: refine_location(2)
    character(len=*), intent(in) :: base_name
    integer, intent(inout)       :: n_output
    integer, parameter           :: n_blocks_per_dim(2) = [1, 1]
    real(dp), parameter          :: block_length(2)     = [1.0_dp, 1.0_dp]
    integer, parameter           :: bx(2)               = [16, 16]
    integer, parameter           :: n_gc                = 1
    integer, parameter           :: n_vars              = 2
    character(len=20)            :: var_names(n_vars)   = ['rho', 'phi']
    logical, parameter           :: periodic(2)         = [.false., .false.]
    integer, parameter           :: min_level           = 5
    integer, parameter           :: max_blocks          = 8000
    integer, parameter           :: n_refine_steps      = 2
    integer                      :: n, n_ghostcells
    integer                      :: t_start, t_end, count_rate
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

    n_output = n_output + 1
    call f4_write_grid(f4, base_name, n_output)

    call system_clock(t_start, count_rate)
    do n = 1, n_update_gc
       call f4_update_ghostcells(f4, 2, [1, 2])
    end do
    call system_clock(t_end, count_rate)

    t_total = (t_end - t_start)/real(count_rate, dp)
    n_ghostcells = f4_get_num_global_blocks(f4) * n_vars * &
         2 * n_gc * sum(f4%bx)
    ghostcells_per_ns = 1e-9_dp * n_update_gc/t_total * n_ghostcells

    if (f4%mpirank == 0) then
       write(*, "(A,F12.6)") " Ghostcells/ns:", ghostcells_per_ns
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

    do n = 1, f4%n_blocks
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
