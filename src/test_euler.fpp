program euler
  use iso_fortran_env, only: int64
  use mpi_f08
  use m_foap4_2d
  use m_config
  use m_euler

  implicit none
  integer, parameter :: dp = kind(0.0d0)

  type(foap4_t) :: f4
  type(CFG_t)   :: cfg

  integer  :: min_level   = 2
  integer  :: max_level   = 4
  integer  :: max_blocks  = 1000
  integer  :: bx(2)       = [32, 32]
  integer  :: num_outputs = 40
  logical  :: periodic(2) = .true.
  logical  :: do_refinement = .false.
  integer  :: n_gc        = 2
  real(dp) :: end_time    = 2.0_dp
  real(dp) :: c_refine = 0.8_dp
  real(dp) :: c_derefine = 0.2_dp
  real(dp) :: c_eps = 0.01_dp
  integer :: n_steps_refinement = 4
  character(len=10) :: test_case = "sod"
  character(len=40) :: integrator_name = "heuns_method"

  call f4_initialize(f4, "error")

  call CFG_update_from_arguments(cfg)
  call CFG_add_get(cfg, 'num_outputs', num_outputs, 'Write this many output files')
  call CFG_add_get(cfg, 'min_level', min_level, 'Minimum refinement level')
  call CFG_add_get(cfg, 'max_level', max_level, 'Maximum refinement level')
  call CFG_add_get(cfg, 'do_refinement', do_refinement, 'Perform refinement')
  call CFG_add_get(cfg, 'c_refine', c_refine, 'Coefficient for refinement')
  call CFG_add_get(cfg, 'c_derefine', c_refine, 'Coefficient for derefinement')
  call CFG_add_get(cfg, 'c_eps', c_refine, 'Used in refinement criterion')
  call CFG_add_get(cfg, 'n_steps_refinement', n_steps_refinement, &
       'Perform refinement every N steps')
  call CFG_add_get(cfg, 'max_blocks', max_blocks, 'Max. number of blocks')
  call CFG_add_get(cfg, 'periodic', periodic, 'Whether the domain is periodic')
  call CFG_add_get(cfg, 'bx', bx, 'Size of grid blocks')
  call CFG_add_get(cfg, 'n_gc', n_gc, 'Number of ghost cells')
  call CFG_add_get(cfg, 'end_time', end_time, 'End time')
  call CFG_add_get(cfg, 'test_case', test_case, 'Which test case to run')
  call CFG_add_get(cfg, 'time_integrator', integrator_name, 'Time integrator')
  call CFG_check(cfg)

  if (n_gc < 2) error stop "n_gc < 2"

  call test_euler(f4, bx, min_level, max_blocks, &
       num_outputs, "output/test_euler", test_case, end_time, integrator_name)

  if (f4%mpirank == 0) call f4_print_wtime(f4)
  call f4_finalize(f4)

contains

  subroutine test_euler(f4, bx, min_level, max_blocks, num_outputs, base_name, &
       test_case, end_time, integrator_name)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: bx(2)
    integer, intent(in)          :: min_level
    integer, intent(in)          :: max_blocks
    integer, intent(in)          :: num_outputs
    character(len=*), intent(in) :: base_name
    character(len=*), intent(in) :: test_case
    real(dp), intent(in)         :: end_time
    character(len=40), intent(in) :: integrator_name
    integer, parameter           :: n_blocks_per_dim(2) = [1, 1]
    real(dp), parameter          :: block_length(2)     = [1.0_dp, 1.0_dp]
    integer, parameter           :: n_gc                = 2
    logical                      :: periodic(2)         = [.true., .true.]
    real(dp), parameter          :: cfl_number          = 0.5_dp
    integer                      :: n_output !n, prev_mesh_revision
    integer                      :: n, n_iterations, ierr, prev_mesh_revision
    integer(int64)               :: sum_local_blocks, sum_global_blocks
    logical                      :: write_this_step, temporal(n_vars_euler)
    integer                      :: integrator, n_time_states
    integer                      :: highest_level, prev_highest_level
    real(dp)                     :: dt, dt_lim, dt_output
    real(dp)                     :: time, t0, t1, rho_sum, rho_initial_sum

    time = 0.0_dp
    dt_output = end_time / max(real(num_outputs, dp), 1e-100_dp)
    n_output = 0
    n_iterations = 0
    sum_local_blocks = 0
    dt_lim = 0.0_dp
    temporal(:) = .true.

    integrator = f4_get_time_integrator(trim(integrator_name))
    n_time_states = f4_advance_num_copies(integrator)

    if (test_case == "rt") periodic(2) = .false.

    call f4_construct_brick(f4, n_blocks_per_dim, block_length, bx, n_gc, &
         n_vars_euler, var_names, temporal, n_time_states, periodic, &
         min_level, max_blocks, f4_bc_neumann, 0.0_dp)

    if (test_case == "rt") then
       call f4_set_physical_boundary(f4, i_momy, 2, f4_bc_dirichlet, 0.0_dp)
       call f4_set_physical_boundary(f4, i_momy, 3, f4_bc_dirichlet, 0.0_dp)
    end if

    call set_initial_conditions(f4, test_case)

    if (do_refinement) then
       do n = 1, 10
          prev_mesh_revision = f4_get_mesh_revision(f4)
          call f4_update_ghostcells(f4, n_vars_euler, i_vars_grid)
          call f4_set_refinement_flags_diff2(f4, min_level, max_level, &
               i_rho, c_refine, c_derefine, c_eps)
          call f4_adjust_refinement(f4, .true.)
          call set_initial_conditions(f4, test_case)

          if (f4_get_mesh_revision(f4) == prev_mesh_revision) exit
       end do
    end if

    call f4_compute_sum(f4, i_rho, rho_initial_sum)

    if (dt_output < end_time) call f4_write_grid(f4, base_name, n_output, time)
    n_output = n_output + 1

    t0 = MPI_Wtime()

    do while (time < end_time)
       n_iterations = n_iterations + 1
       dt = dt_lim

       write_this_step = (time + dt > n_output * dt_output)
       if (write_this_step) dt = n_output * dt_output - time

       call f4_advance(f4, dt, dt_lim, time, integrator, forward_euler)

       call MPI_Allreduce(MPI_IN_PLACE, dt_lim, 1, MPI_DOUBLE_PRECISION, &
            MPI_MIN, f4%mpicomm, ierr)

       if (write_this_step) then
          call f4_write_grid(f4, base_name, n_output, time)
          call f4_compute_sum(f4, i_rho, rho_sum)
          if (f4%mpirank == 0) then
             write(*, "(A,E12.4)") " Conservation error: ", &
                  rho_sum - rho_initial_sum
          end if
          n_output = n_output + 1
       end if

       if (do_refinement .and. &
            mod(n_iterations, n_steps_refinement) == 0) then
          call f4_get_global_highest_level(f4, prev_highest_level)
          call f4_update_ghostcells(f4, n_vars_euler, i_vars_grid)
          call f4_set_refinement_flags_diff2(f4, min_level, max_level, &
               i_rho, c_refine, c_derefine, c_eps)
          call f4_adjust_refinement(f4, .true.)

          call f4_get_global_highest_level(f4, highest_level)

          if (highest_level > prev_highest_level) then
             dt_lim = 0.5_dp * dt_lim
          end if
       end if

       sum_local_blocks = sum_local_blocks + f4_get_num_local_blocks(f4)
    end do

    t1 = MPI_Wtime()

    call MPI_Reduce(sum_local_blocks, sum_global_blocks, 1, MPI_INTEGER8, &
         MPI_SUM, 0, f4%mpicomm, ierr)

    if (f4%mpirank == 0) then
       print *, "n_iterations:    ", n_iterations
       print *, "n_blocks_global: ", sum_global_blocks/n_iterations
       print *, "block size:      ", bx
       write(*, "(A,F14.3)") " cell updates/ns:     ", 1e-9_dp * &
            sum_global_blocks * (product(f4%bx) * f4%n_temporal_states / (t1 - t0))
    end if

    call f4_destroy(f4)
  end subroutine test_euler

  pure subroutine to_primitive(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(n_vars_euler)
    real(dp) :: inv_rho

    inv_rho = 1/u(i_rho)
    u(i_momx) = u(i_momx) * inv_rho
    u(i_momy) = u(i_momy) * inv_rho

    u(i_e) = (euler_gamma-1.0_dp) * (u(i_e) - &
         0.5_dp * u(i_rho)* (u(i_momx)**2 + u(i_momy)**2))

  end subroutine to_primitive

  pure subroutine to_conservative(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(n_vars_euler)

    ! Compute energy from pressure and kinetic energy
    u(i_e) = u(i_e) * inv_gamma_m1 + &
         0.5_dp * u(i_rho) * (u(i_momx)**2 + u(i_momy)**2)

    ! Compute momentum from density and velocity components
    u(i_momx) = u(i_rho) * u(i_momx)
    u(i_momy) = u(i_rho) * u(i_momy)
  end subroutine to_conservative

  subroutine muscl_flux_euler_prim(u, flux_dim, flux, wmax)
    !$acc routine seq
    real(dp), intent(in)  :: u(5, n_vars_euler)
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(n_vars_euler, 2)
    real(dp), intent(out) :: wmax
    real(dp)              :: wmax1, wmax2

    call muscl_flux_euler_prim_side(u, 2, flux_dim, flux(:, 1), wmax1)
    call muscl_flux_euler_prim_side(u, 3, flux_dim, flux(:, 2), wmax2)
    wmax = max(wmax1, wmax2)
  end subroutine muscl_flux_euler_prim

  subroutine muscl_flux_euler_prim_side(u, i0, flux_dim, flux, wmax)
    !$acc routine seq
    real(dp), intent(in)  :: u(5, n_vars_euler)
    integer, intent(in)   :: i0
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(n_vars_euler)
    real(dp), intent(out) :: wmax
    real(dp)              :: uL(n_vars_euler), uR(n_vars_euler), wL, wR
    real(dp)              :: flux_l(n_vars_euler), flux_r(n_vars_euler)

    ! Construct uL, uR for cell face
    uL = u(i0, :) + 0.5_dp * vanleer(u(i0, :) - u(i0-1, :), &
         u(i0+1, :) - u(i0, :))
    uR = u(i0+1, :) - 0.5_dp * vanleer(u(i0+1, :) - u(i0, :), &
         u(i0+2, :) - u(i0+1, :))

    call euler_flux(uL, flux_dim, flux_l, wL)
    call euler_flux(uR, flux_dim, flux_r, wR)
    wmax = max(wL, wR)

    call to_conservative(uL)
    call to_conservative(uR)
    flux(:) = 0.5_dp * (flux_l + flux_r - wmax * (uR - uL))

  end subroutine muscl_flux_euler_prim_side

  subroutine euler_flux(u, flux_dim, flux, w)
    !$acc routine seq
    real(dp), intent(in)  :: u(n_vars_euler)
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(n_vars_euler)
    real(dp), intent(out) :: w

    ! Density flux
    flux(i_rho) = u(i_rho) * u(i_mom(flux_dim))

    ! Momentum flux with pressure term
    flux(i_momx) = u(i_rho) * u(i_momx) * u(i_mom(flux_dim))
    flux(i_momy) = u(i_rho) * u(i_momy) * u(i_mom(flux_dim))
    flux(i_mom(flux_dim)) = flux(i_mom(flux_dim)) + u(i_e)

    ! Energy flux
    flux(i_e) = u(i_mom(flux_dim)) * (u(i_e) * inv_gamma_m1 + &
         0.5_dp * u(i_rho) * (u(i_momx)**2 + u(i_momy)**2) + u(i_e))

    w = sqrt(euler_gamma * u(i_e) / u(i_rho)) + abs(u(i_mom(flux_dim)))
  end subroutine euler_flux

  subroutine set_initial_conditions(f4, test_case)
    type(foap4_t), intent(inout) :: f4
    character(len=*), intent(in) :: test_case
    integer                      :: n
    real(dp)                     :: u0(n_vars_euler, 4)

    select case (test_case)
    case ("first")
       u0(i_e, :)      = [1.0_dp, 0.4_dp, 0.0439_dp, 0.15_dp]
       u0(i_rho, :)    = [1.0_dp, 0.5197_dp, 0.1072_dp, 0.2579_dp]
       u0(i_momx, :) = [0.0_dp, -0.7259_dp, -0.7259_dp, 0.0_dp]
       u0(i_momy, :) = [0.0_dp, 0.0_dp, -1.4045_dp, -1.4045_dp]

       do n = 1, 4
          call to_conservative(u0(:, n))
       end do
       call set_initial_quadrants(f4, u0)
    case ("sixth")
       u0(i_e, :)      = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
       u0(i_rho, :)    = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
       u0(i_momx, :) = [0.75_dp, 0.75_dp, -0.75_dp, -0.75_dp]
       u0(i_momy, :) = [-0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp]

       do n = 1, 4
          call to_conservative(u0(:, n))
       end do
       call set_initial_quadrants(f4, u0)
    case ("sod")
       ! 1D Sod shock test case
       u0(i_rho, :)    = [0.125_dp, 1.0_dp, 1.0_dp, 0.125_dp]
       u0(i_e, :)      = [0.1_dp, 1.0_dp, 1.0_dp, 0.1_dp]
       u0(i_momx, :) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
       u0(i_momy, :) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]

       do n = 1, 4
          call to_conservative(u0(:, n))
       end do
       call set_initial_quadrants(f4, u0)
    case ("rt")
       call set_rayleigh_taylor(f4)
    case default
       error stop "Unknown test case, options: rt, sod, first, sixth"
    end select

  end subroutine set_initial_conditions

  subroutine set_initial_quadrants(f4, u0)
    type(foap4_t), intent(inout) :: f4
    real(dp), intent(in)         :: u0(n_vars_euler, 4)
    integer                      :: n, i, j
    real(dp)                     :: rr(2)

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
  end subroutine set_initial_quadrants

  subroutine set_rayleigh_taylor(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, i, j
    real(dp)                     :: y0, rho_high, rho_low, width
    real(dp)                     :: p_interface, kx, rr(2)
    real(dp), parameter          :: pi = acos(-1.0_dp)

    ! The location of interface
    y0 = 0.8d0

    ! Width of the sinusoidal fluctuations on the interface
    width = 0.05d0

    ! High and low density
    rho_high = 1.0_dp
    rho_low  = 0.1_dp

    ! Pressure at interface
    p_interface = 1.0_dp

    ! Wavelength
    kx = 2 * pi

    !$acc parallel loop
    do n = 1, f4%n_blocks
       !$acc loop collapse(2) private(rr)
       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)
             rr = f4_cell_coord(f4, n, i, j)

             f4%uu(i, j, i_momx, n) = 0.0_dp
             f4%uu(i, j, i_momy, n) = 0.0_dp

             if (rr(2) > y0 + width * sin(kx * rr(1))) then
                f4%uu(i, j, i_rho, n) = rho_high
             else
                f4%uu(i, j, i_rho, n) = rho_low
             end if

             f4%uu(i, j, i_e, n) = inv_gamma_m1 * (p_interface - 1.0_dp * &
                  f4%uu(i, j, i_rho, n) * (rr(2) - y0))
          end do
       end do
    end do
  end subroutine set_rayleigh_taylor

  subroutine forward_euler(f4, dt, dt_lim, time, s_deriv, n_prev, s_prev, w_prev, &
       s_out, i_step, n_steps)
    type(foap4_t), intent(inout) :: f4
    real(dp), intent(in)         :: dt
    real(dp), intent(inout)      :: dt_lim
    real(dp), intent(in)         :: time
    integer, intent(in)          :: s_deriv        !< State to compute derivatives from
    integer, intent(in)          :: n_prev         !< Number of previous states
    integer, intent(in)          :: s_prev(n_prev) !< Previous states
    real(dp), intent(in)         :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)          :: s_out          !< Output state
    integer, intent(in)          :: i_step         !< Step of the integrator
    integer, intent(in)          :: n_steps        !< Total number of steps
    integer                      :: n, i, j, m, iv, level
    real(dp)                     :: inv_dr(2), wmax(2), max_cfl
    real(dp)                     :: fx(n_vars_euler, 2), fy(n_vars_euler, 2)
    real(dp)                     :: tmp(5, n_vars_euler), u(n_vars_euler)
    real(dp)                     :: uprim(f4%ilo(1):f4%ihi(1), f4%ilo(2):f4%ihi(2), n_vars_euler)
    real(dp)                     :: dvar(n_vars_euler)

    call f4_update_ghostcells(f4, n_vars_euler, i_vars_grid+s_deriv)

    max_cfl = 0.0_dp

    !$acc parallel loop private(level, inv_dr, max_cfl, uprim) &
    !$acc &reduction(max:max_cfl)
    do n = 1, f4%n_blocks
       level = f4%block_level(n)
       inv_dr = 1/f4%dr_level(:, level)

       !$acc loop collapse(2) private(u)
       do j = f4%ilo(2), f4%ihi(2)
          do i = f4%ilo(1), f4%ihi(1)
             if ((j >= 1 .and. j <= f4%bx(2)) .or. (i >= 1 .and. i <= f4%bx(1))) then
                ! Convert to primitive
                u = f4%uu(i, j, i_vars_grid+s_deriv, n)
                call to_primitive(u)
                uprim(i, j, :) = u
             end if
          end do
       end do

       !$acc loop collapse(2) private(tmp, fx, fy, dvar, wmax, iv, m) &
       !$acc &reduction(max:max_cfl)
       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)
             ! Compute x and y fluxes
             tmp = uprim(i-2:i+2, j, :)
             call muscl_flux_euler_prim(tmp, 1, fx, wmax(1))

             tmp = uprim(i, j-2:j+2, :)
             call muscl_flux_euler_prim(tmp, 2, fy, wmax(2))

             max_cfl = max(max_cfl, sum(wmax * inv_dr))

             ! Change due to fluxes
             dvar(:) = dt * ((fx(:, 1) - fx(:, 2)) * inv_dr(1) + &
                  (fy(:, 1) - fy(:, 2)) * inv_dr(2))

             ! Store boundary fluxes
             if (i == 1) f4%bflux(j, 0, i_rho:i_rho+n_vars_euler-1, n) = dt * fx(:, 1)
             if (i == bx(1)) f4%bflux(j, 1, i_rho:i_rho+n_vars_euler-1, n) = dt * fx(:, 2)
             if (j == 1) f4%bflux(i, 2, i_rho:i_rho+n_vars_euler-1, n) = dt * fy(:, 1)
             if (j == bx(2)) f4%bflux(i, 3, i_rho:i_rho+n_vars_euler-1, n) = dt * fy(:, 2)

             ! Change due to source terms
             ! TODO: enable/disable this source term in a nice way
             ! TODO: use variable for source term
             dvar(i_momy) = dvar(i_momy) + dt * (-1.0_dp) * f4%uu(i, j, i_rho+s_deriv, n)
             dvar(i_e) = dvar(i_e) + dt * (-1.0_dp) * f4%uu(i, j, i_momy+s_deriv, n)

             ! Set output state
             do iv = 1, n_vars_euler
                do m = 1, n_prev
                   ! Add weighted previous states
                   dvar(iv) = dvar(iv) + &
                        f4%uu(i, j, i_vars_grid(iv)+s_prev(m), n) * w_prev(m)
                end do
                f4%uu(i, j, i_vars_grid(iv)+s_out, n) = dvar(iv)
             end do
          end do
       end do
    end do

    dt_lim = 1/max_cfl

    call f4_fix_c2f_flux(f4, n_vars_euler, i_vars_grid, s_out)

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

end program
