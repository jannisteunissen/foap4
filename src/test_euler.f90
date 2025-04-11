program euler
  use iso_fortran_env, only: int64
  use mpi_f08
  use m_foap4
  use m_config
  use m_euler

  implicit none
  integer, parameter :: dp = kind(0.0d0)

  type(foap4_t) :: f4
  type(CFG_t)   :: cfg

  integer  :: min_level   = 2
  integer  :: max_blocks  = 1000
  integer  :: bx(2)       = [32, 32]
  integer  :: num_outputs = 40
  integer  :: n_iter      = 100
  logical  :: periodic(2) = .true.
  integer  :: n_gc        = 2
  real(dp) :: dt          = 1e-3_dp

  call f4_initialize(f4, "error")

  call CFG_update_from_arguments(cfg)
  call CFG_add_get(cfg, 'num_outputs', num_outputs, 'Write this many output files')
  call CFG_add_get(cfg, 'min_level', min_level, 'Initial refinement level')
  call CFG_add_get(cfg, 'max_blocks', max_blocks, 'Max. number of blocks')
  call CFG_add_get(cfg, 'n_iter', n_iter, 'Number of iterations')
  call CFG_add_get(cfg, 'periodic', periodic, 'Whether the domain is periodic')
  call CFG_add_get(cfg, 'bx', bx, 'Size of grid blocks')
  call CFG_add_get(cfg, 'n_gc', n_gc, 'Number of ghost cells')
  call CFG_add_get(cfg, 'dt', dt, 'Time step')
  call CFG_check(cfg)

  if (n_gc < 2) error stop "n_gc < 2"

  call test_euler(f4, bx, min_level, max_blocks, &
       num_outputs, "output/test_euler")

  if (f4%mpirank == 0) call f4_print_wtime(f4)
  call f4_finalize(f4)

contains

  subroutine test_euler(f4, bx, min_level, max_blocks, num_outputs, base_name)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: bx(2)
    integer, intent(in)          :: min_level
    integer, intent(in)          :: max_blocks
    integer, intent(in)          :: num_outputs
    character(len=*), intent(in) :: base_name
    integer, parameter           :: n_blocks_per_dim(2) = [1, 1]
    real(dp), parameter          :: block_length(2)     = [1.0_dp, 1.0_dp]
    integer, parameter           :: n_gc                = 2
    logical, parameter           :: periodic(2)         = [.true., .true.]
    real(dp), parameter          :: cfl_number          = 0.5_dp
    integer                      :: n_output !n, prev_mesh_revision
    integer                      :: n_iterations, ierr
    integer(int64)               :: sum_local_blocks, sum_global_blocks
    logical                      :: write_this_step
    real(dp)                     :: dt, dt_lim, dt_output
    real(dp)                     :: time, end_time, t0, t1

    time = 0.0_dp
    end_time = 1.0_dp
    dt_output = end_time / max(real(num_outputs, dp), 1e-100_dp)
    n_output = 0
    n_iterations = 0
    sum_local_blocks = 0
    dt_lim = 0.0_dp

    call f4_construct_brick(f4, n_blocks_per_dim, block_length, bx, n_gc, &
         n_variables, var_names, periodic, min_level, max_blocks)

    call set_initial_conditions(f4, "sod")

    ! if (do_refinement) then
    !    do n = 1, 10
    !       prev_mesh_revision = f4_get_mesh_revision(f4)
    !       call f4_update_ghostcells(f4, 1, [i_rho])
    !       call set_refinement_flag(f4)
    !       call f4_adjust_refinement(f4, .true.)
    !       call set_initial_conditions(f4)

    !       if (f4_get_mesh_revision(f4) == prev_mesh_revision) exit
    !    end do
    ! end if

    if (dt_output < end_time) call f4_write_grid(f4, base_name, n_output, time)
    n_output = n_output + 1

    t0 = MPI_Wtime()

    do while (time < end_time)
       n_iterations = n_iterations + 1
       dt = dt_lim

       write_this_step = (time + dt > n_output * dt_output)
       if (write_this_step) dt = n_output * dt_output - time

       call advance_heuns_method(f4, dt, dt_lim)
       time = time + dt

       call MPI_Allreduce(MPI_IN_PLACE, dt_lim, 1, MPI_DOUBLE_PRECISION, &
            MPI_MIN, f4%mpicomm, ierr)

       if (write_this_step) then
          call f4_write_grid(f4, base_name, n_output, time)
          n_output = n_output + 1
       end if

       ! if (do_refinement) then
       !    call f4_update_ghostcells(f4, 1, [i_rho])
       !    call set_refinement_flag(f4)
       !    call f4_adjust_refinement(f4, .true.)

       !    call f4_get_global_highest_level(f4, highest_level)
       ! TODO: reduce dt
       ! end if

       sum_local_blocks = sum_local_blocks + f4_get_num_local_blocks(f4)
    end do

    t1 = MPI_Wtime()

    call MPI_Reduce(sum_local_blocks, sum_global_blocks, 1, MPI_INTEGER8, &
         MPI_SUM, 0, f4%mpicomm, ierr)

    if (f4%mpirank == 0) then
       print *, "n_iterations:    ", n_iterations
       print *, "n_blocks_global: ", sum_global_blocks/n_iterations
       print *, "block size:      ", bx
       write(*, "(A,F14.3)") " unknowns/ns:     ", 1e-9_dp * &
            sum_global_blocks * (product(f4%bx) * 2 / (t1 - t0))
    end if

    call f4_destroy(f4)
  end subroutine test_euler

  pure subroutine to_primitive(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(n_vars_euler)

    u(ix_momx) = u(ix_momx)/u(ix_rho)
    u(ix_momy) = u(ix_momy)/u(ix_rho)

    u(ix_e) = (euler_gamma-1.0_dp) * (u(ix_e) - &
         0.5_dp * u(ix_rho)* (u(ix_momx)**2 + u(ix_momy)**2))

  end subroutine to_primitive

  pure subroutine to_conservative(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(n_vars_euler)

    ! Compute energy from pressure and kinetic energy
    u(ix_e) = u(ix_e) * inv_gamma_m1 + &
         0.5_dp * u(ix_rho) * (u(ix_momx)**2 + u(ix_momy)**2)

    ! Compute momentum from density and velocity components
    u(ix_momx) = u(ix_rho) * u(ix_momx)
    u(ix_momy) = u(ix_rho) * u(ix_momy)
  end subroutine to_conservative

  subroutine muscl_flux_euler_prim(u, flux_dim, flux, wmax)
    !$acc routine seq
    real(dp), intent(in)  :: u(5, n_vars_euler)
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(n_vars_euler, 2)
    real(dp), intent(out) :: wmax
    real(dp)              :: wmax1, wmax2
    real(dp)              :: uL(n_vars_euler), uR(n_vars_euler), wL, wR
    real(dp)              :: flux_l(n_vars_euler), flux_r(n_vars_euler)

    ! Construct uL, uR for first cell face
    uL = u(2, :) + 0.5_dp * vanleer(u(2, :) - u(1, :), u(3, :) - u(2, :))
    uR = u(3, :) - 0.5_dp * vanleer(u(3, :) - u(2, :), u(4, :) - u(3, :))

    call euler_flux(uL, flux_dim, flux_l, wL)
    call euler_flux(uR, flux_dim, flux_r, wR)
    wmax1 = max(wL, wR)

    call to_conservative(uL)
    call to_conservative(uR)
    flux(:, 1) = 0.5_dp * (flux_l + flux_r - wmax1 * (uR - uL))

    ! Construct uL, uR for second cell face
    uL = u(3, :) + 0.5_dp * vanleer(u(3, :) - u(2, :), u(4, :) - u(3, :))
    uR = u(4, :) - 0.5_dp * vanleer(u(4, :) - u(3, :), u(5, :) - u(4, :))

    call euler_flux(uL, flux_dim, flux_l, wL)
    call euler_flux(uR, flux_dim, flux_r, wR)
    wmax2 = max(wL, wR)

    call to_conservative(uL)
    call to_conservative(uR)
    flux(:, 2) = 0.5_dp * (flux_l + flux_r - wmax2 * (uR - uL))

    wmax = max(wmax1, wmax2)
  end subroutine muscl_flux_euler_prim

  subroutine euler_flux(u, flux_dim, flux, w)
    !$acc routine seq
    real(dp), intent(in)  :: u(n_vars_euler)
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(n_vars_euler)
    real(dp), intent(out) :: w

    ! Density flux
    flux(ix_rho) = u(ix_rho) * u(ix_mom(flux_dim))

    ! Momentum flux with pressure term
    flux(ix_momx) = u(ix_rho) * u(ix_momx) * u(ix_mom(flux_dim))
    flux(ix_momy) = u(ix_rho) * u(ix_momy) * u(ix_mom(flux_dim))
    flux(ix_mom(flux_dim)) = flux(ix_mom(flux_dim)) + u(ix_e)

    ! Energy flux
    flux(ix_e) = u(ix_mom(flux_dim)) * (u(ix_e) * inv_gamma_m1 + &
         0.5_dp * u(ix_rho) * (u(ix_momx)**2 + u(ix_momy)**2) + u(ix_e))

    w = sqrt(euler_gamma * u(ix_e) / u(ix_rho)) + abs(u(ix_mom(flux_dim)))
  end subroutine euler_flux

  subroutine advance_heuns_method(f4, dt, dt_lim)
    type(foap4_t), intent(inout) :: f4
    real(dp), intent(in)         :: dt
    real(dp), intent(out)        :: dt_lim

    call forward_euler(f4, f4%bx, f4%ilo, f4%ihi, f4%n_vars, f4%n_blocks, &
         dt, dt_lim, f4%uu, 0, 1, [0], [1.0_dp], n_vars_euler)
    call forward_euler(f4, f4%bx, f4%ilo, f4%ihi, f4%n_vars, f4%n_blocks, &
         0.5_dp*dt, dt_lim, f4%uu, n_vars_euler, 2, [0, n_vars_euler], &
         [0.5_dp, 0.5_dp], 0)
  end subroutine advance_heuns_method

  subroutine set_initial_conditions(f4, test_case)
    type(foap4_t), intent(inout) :: f4
    character(len=*), intent(in) :: test_case
    integer                      :: n, i, j
    real(dp)                     :: rr(2)
    real(dp)                     :: u0(n_vars_euler, 4)

    select case (test_case)
    case ("first")
       u0(ix_e, :)      = [1.0_dp, 0.4_dp, 0.0439_dp, 0.15_dp]
       u0(ix_rho, :)    = [1.0_dp, 0.5197_dp, 0.1072_dp, 0.2579_dp]
       u0(ix_momx, :) = [0.0_dp, -0.7259_dp, -0.7259_dp, 0.0_dp]
       u0(ix_momy, :) = [0.0_dp, 0.0_dp, -1.4045_dp, -1.4045_dp]
    case ("sixth")
       u0(ix_e, :)      = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
       u0(ix_rho, :)    = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
       u0(ix_momx, :) = [0.75_dp, 0.75_dp, -0.75_dp, -0.75_dp]
       u0(ix_momy, :) = [-0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp]
    case ("sod")
       ! 1D Sod shock test case
       u0(ix_rho, :)    = [0.125_dp, 1.0_dp, 1.0_dp, 0.125_dp]
       u0(ix_e, :)      = [0.1_dp, 1.0_dp, 1.0_dp, 0.1_dp]
       u0(ix_momx, :) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
       u0(ix_momy, :) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    case default
       error stop "Unknown test case"
    end select

    do n = 1, 4
       call to_conservative(u0(:, n))
    end do

    !$acc parallel loop
    do n = 1, f4%n_blocks
       !$acc loop collapse(2) private(rr)
       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)
             rr = f4_cell_coord(f4, n, i, j)

             if (rr(1) > 0.5_dp .and. rr(2) > 0.5_dp) then
                f4%uu(i, j, i_vars_grid, n) = u0(:, 1)
             elseif (rr(1) <= 0.5_dp .and. rr(2) >= 0.5_dp) then
                f4%uu(i, j, i_vars_grid, n) = u0(:, 2)
             elseif (rr(1) <= 0.5_dp .and. rr(2) <= 0.5_dp) then
                f4%uu(i, j, i_vars_grid, n) = u0(:, 3)
             else
                f4%uu(i, j, i_vars_grid, n) = u0(:, 4)
             end if
          end do
       end do
    end do
  end subroutine set_initial_conditions

  subroutine forward_euler(f4, bx, lo, hi, n_vars, n_blocks, dt, dt_lim, uu, &
       s_deriv, n_prev, s_prev, w_prev, s_out)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_blocks, bx(2), lo(2), hi(2), n_vars
    real(dp), intent(in)         :: dt
    real(dp), intent(out)        :: dt_lim
    real(dp), intent(inout)      :: uu(lo(1):hi(1), lo(2):hi(2), n_vars, n_blocks)
    integer, intent(in)          :: s_deriv        !< State to compute derivatives from
    integer, intent(in)          :: n_prev         !< Number of previous states
    integer, intent(in)          :: s_prev(n_prev) !< Previous states
    real(dp), intent(in)         :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)          :: s_out          !< Output state
    integer                      :: n, i, j, m, iv, level
    real(dp)                     :: inv_dr(2), wmax(2), max_cfl
    real(dp)                     :: fx(n_vars_euler, 2), fy(n_vars_euler, 2)
    real(dp)                     :: tmp(5, n_vars_euler), u(n_vars_euler)
    real(dp)                     :: uprim(lo(1):hi(1), lo(2):hi(2), n_vars_euler)
    real(dp)                     :: dvar(bx(1), bx(2), n_vars_euler)

    call f4_update_ghostcells(f4, n_vars_euler, i_vars_grid+s_deriv)

    dt_lim = 1e100_dp

    !$acc parallel loop private(level, inv_dr, max_cfl, uprim, dvar) &
    !$acc &reduction(min:dt_lim)
    do n = 1, n_blocks
       level = f4%block_level(n)
       inv_dr = 1/f4%dr_level(:, level)
       max_cfl = 0.0_dp

       !$acc loop collapse(2) private(u)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ! Convert to primitive
             u = uu(i, j, i_vars_grid+s_deriv, n)
             call to_primitive(u)
             uprim(i, j, :) = u
          end do
       end do

       !$acc loop collapse(2) private(tmp, fx, fy, wmax) reduction(max:max_cfl)
       do j = 1, bx(2)
          do i = 1, bx(1)
             ! Compute x and y fluxes
             tmp = uprim(i-2:i+2, j, :)
             call muscl_flux_euler_prim(tmp, 1, fx, wmax(1))

             tmp = uprim(i, j-2:j+2, :)
             call muscl_flux_euler_prim(tmp, 2, fy, wmax(2))

             ! Keep track of changes in variables
             dvar(i, j, :) = dt * &
                  ((fx(:, 1) - fx(:, 2)) * inv_dr(1) + &
                  (fy(:, 1) - fy(:, 2)) * inv_dr(2))

             max_cfl = max(max_cfl, sum(wmax * inv_dr))
          end do
       end do

       dt_lim = min(dt_lim, 1/max_cfl)

       ! Set output state after computations are done, since s_out can be
       ! equal to s_deriv and s_prev
       !$acc loop collapse(3) private(m)
       do iv = 1, n_vars_euler
          do j = 1, bx(2)
             do i = 1, bx(1)
                do m = 1, n_prev
                   ! Add weighted previous states
                   dvar(i, j, iv) = dvar(i, j, iv) + &
                        uu(i, j, i_vars_grid(iv)+s_prev(m), n) * w_prev(m)
                end do
                uu(i, j, i_vars_grid(iv)+s_out, n) = dvar(i, j, iv)
             end do
          end do
       end do
    end do

  end subroutine forward_euler

  elemental pure real(dp) function minmod(a, b)
    !$acc routine seq
    real(dp), intent(in) :: a, b

    if (a * b <= 0) then
       minmod = 0.0_dp
    else if (abs(a) < abs(b)) then
       minmod = a
    else
       minmod = b
    end if
  end function minmod

  elemental pure real(dp) function vanleer(a, b) result(phi)
    !$acc routine seq
    real(dp), intent(in) :: a, b
    real(dp)             :: ab

    ab = a * b
    if (ab > 0) then
       phi = 2 * ab / (a + b)
    else
       phi = 0
    end if
  end function vanleer

end program euler
