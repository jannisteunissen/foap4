#:include '../core/definitions.fpp'
#:set LIMITER = getvar('USE_LIMITER', 'vanleer')
#:set FLUX_SCHEME = getvar('USE_FLUX_SCHEME', 'hll')

program euler
  use iso_fortran_env, only: int64
  use mpi_f08
  use m_foap4_${NDIM}$d
  use m_physics_euler_${NDIM}$d
  use m_rk_${NDIM}$d
  use m_io_${NDIM}$d
  use m_amr_flags_${NDIM}$d
  use m_config

  implicit none

  include 'limiter_${LIMITER}$_definitions.f90'
  integer, parameter :: n_gc = limiter_num_ghostcells

  type(foap4_t) :: f4
  type(CFG_t)   :: cfg

  integer           :: min_level          = 2
  integer           :: max_level          = 4
  integer           :: max_blocks         = 1000
  integer           :: blocks_per_dim(NDIM) = 1
  integer           :: bx(NDIM)           = 32
  integer           :: num_outputs        = 40
  logical           :: do_refinement      = .false.
  real(dp)          :: end_time           = 2.0_dp
  real(dp)          :: c_refine           = 0.8_dp
  real(dp)          :: c_derefine         = 0.2_dp
  real(dp)          :: c_eps              = 0.01_dp
  real(dp)          :: c_abs              = 1e-10_dp
  integer           :: n_steps_refinement = 4
  real(dp)          :: cfl_number         = 0.5_dp
  real(dp)          :: load_imbalance_threshold = 1.1_dp
  character(len=10) :: test_case          = "sod"
  character(len=40) :: integrator_name    = "heuns_method"
  character(len=200) :: output_prefix    = "output/test_euler_${NDIM}$d"
  logical           :: write_vtu = .false.

  call f4_initialize(f4, "error")

  call CFG_update_from_arguments(cfg)
  call CFG_add_get(cfg, 'num_outputs', num_outputs, 'Write this many output files')
  call CFG_add_get(cfg, 'output_prefix', output_prefix, 'Prefix for output files')
  call CFG_add_get(cfg, 'write_vtu', write_vtu, 'Also write p4est vtu files')
  call CFG_add_get(cfg, 'min_level', min_level, 'Minimum refinement level')
  call CFG_add_get(cfg, 'max_level', max_level, 'Maximum refinement level')
  call CFG_add_get(cfg, 'do_refinement', do_refinement, 'Perform refinement')
  call CFG_add_get(cfg, 'load_imbalance_threshold', load_imbalance_threshold, &
       'Threshold for partitioning')
  call CFG_add_get(cfg, 'c_refine', c_refine, 'Coefficient for refinement')
  call CFG_add_get(cfg, 'c_derefine', c_derefine, 'Coefficient for derefinement')
  call CFG_add_get(cfg, 'c_eps', c_eps, 'Filter coefficient for AMR')
  call CFG_add_get(cfg, 'c_abs', c_abs, 'Density threshold for AMR')
  call CFG_add_get(cfg, 'n_steps_refinement', n_steps_refinement, &
       'Perform refinement every N steps')
  call CFG_add_get(cfg, 'max_blocks', max_blocks, 'Max. number of blocks')
  call CFG_add_get(cfg, 'bx', bx, 'Size of grid blocks')
  call CFG_add_get(cfg, 'blocks_per_dim', blocks_per_dim, &
       'Number of blocks (per dimension) on coarse grid')
  call CFG_add_get(cfg, 'end_time', end_time, 'End time')
  call CFG_add_get(cfg, 'test_case', test_case, 'Which test case to run')
  call CFG_add_get(cfg, 'time_integrator', integrator_name, 'Time integrator')
  call CFG_add_get(cfg, 'cfl_number', cfl_number, 'CFL number')
  call CFG_check(cfg)

  call test_euler(f4, bx, min_level, max_blocks, &
       num_outputs, trim(output_prefix), test_case, &
       end_time, integrator_name)

  call f4_print_wtime(f4)
  call f4_finalize(f4)

contains

  subroutine test_euler(f4, bx, min_level, max_blocks, num_outputs, base_name, &
       test_case, end_time, integrator_name)
    type(foap4_t), intent(inout)  :: f4
    integer, intent(in)           :: bx(NDIM)
    integer, intent(in)           :: min_level
    integer, intent(in)           :: max_blocks
    integer, intent(in)           :: num_outputs
    character(len=*), intent(in)  :: base_name
    character(len=*), intent(in)  :: test_case
    real(dp), intent(in)          :: end_time
    character(len=40), intent(in) :: integrator_name
    real(dp)                      :: block_length(NDIM), domain_length(NDIM)
    logical                       :: periodic(NDIM) = .true.
    integer                       :: n_output
    integer                       :: n, n_iterations, ierr, prev_mesh_revision
    integer(int64)                :: sum_local_blocks, sum_global_blocks
    logical                       :: write_this_step
    integer                       :: integrator, n_time_states
    integer                       :: highest_level, prev_highest_level
    real(dp)                      :: dt, dt_lim, dt_output
    real(dp)                      :: t0, t1, rho_sum, rho_initial_sum

    f4%time = 0.0_dp
    dt_output = end_time / max(real(num_outputs, dp), 1e-100_dp)
    n_output = 0
    n_iterations = 0
    sum_local_blocks = 0
    dt_lim = 0.0_dp

    integrator = rk_get_integrator_by_name(trim(integrator_name))
    n_time_states = rk_advance_num_copies(integrator)

    select case (test_case)
    case ("rt")
       ! Rayleigh-Taylor test case, as also used in e.g. MPI-AMRVAC
       periodic(NDIM) = .false.
       domain_length = 1.0_dp
       block_length = domain_length / blocks_per_dim

       call f4_construct_brick(f4, blocks_per_dim, block_length, bx, n_gc, &
            n_vars_all, var_names, var_temporal, n_time_states, periodic, &
            min_level, max_blocks, f4_bc_neumann, 0.0_dp)

       call f4_set_bc_scalar(f4, i_mom0+NDIM, 2*(NDIM-1), &
            f4_bc_dirichlet, 0.0_dp)
       call f4_set_bc_scalar(f4, i_mom0+NDIM, 2*(NDIM-1)+1, &
            f4_bc_dirichlet, 0.0_dp)

       call euler_initialize(5/3.0_dp, -1.0_dp, 1e-12_dp, 1e-12_dp)
    case ("sod")
       ! Standard 1D Sod shock tube test
       domain_length = 1.0_dp
       block_length = domain_length / blocks_per_dim

       call f4_construct_brick(f4, blocks_per_dim, block_length, bx, n_gc, &
            n_vars_all, var_names, var_temporal, n_time_states, periodic, &
            min_level, max_blocks, f4_bc_neumann, 0.0_dp)

       call euler_initialize(5/3.0_dp, 0.0_dp, 1e-12_dp, 1e-12_dp)
    case ("vortex")
       ! Isentropic vortex from Shu "Essentially non-oscillatory and
       ! weighted essentially non-oscillatory schemes for hyperbolic
       ! conservation laws" (1998)
       domain_length = 10.0_dp
       block_length = domain_length / blocks_per_dim

       call f4_construct_brick(f4, blocks_per_dim, block_length, bx, n_gc, &
            n_vars_all, var_names, var_temporal, n_time_states, periodic, &
            min_level, max_blocks, f4_bc_neumann, 0.0_dp)

       call euler_initialize(1.4_dp, 0.0_dp, 1e-12_dp, 1e-12_dp)
    case default
       error stop "Unknown test case, options: rt, sod, vortex"
    end select

    call set_initial_conditions(f4, test_case)

    if (do_refinement) then
       do n = 1, 10
          prev_mesh_revision = f4_get_mesh_revision(f4)
          call f4_update_ghostcells(f4, n_tvars, i_tvars, 0)
          call amr_flags_diff2(f4, min_level, max_level, &
               i_rho, c_refine, c_derefine, c_eps, c_abs)
          call f4_adjust_refinement(f4, load_imbalance_threshold)
          call set_initial_conditions(f4, test_case)

          if (f4_get_mesh_revision(f4) == prev_mesh_revision) exit
       end do
    end if

    call f4_compute_sum(f4, i_rho, rho_initial_sum)

    if (dt_output <= end_time) then
       call io_write_grid(f4, base_name, n_output, write_p4vtu=write_vtu)
    end if
    n_output = n_output + 1

    ! Only report timings for main time loop
    call f4_reset_wtime(f4)

    t0 = MPI_Wtime()

    do while (f4%time < end_time)
       n_iterations = n_iterations + 1
       dt = cfl_number * dt_lim

       write_this_step = (f4%time + dt > n_output * dt_output)
       if (write_this_step) dt = n_output * dt_output - f4%time
       call rk_advance(f4, dt, dt_lim, integrator, feuler_finite_volume)

       if (write_this_step) then
          call io_write_grid(f4, base_name, n_output, write_p4vtu=write_vtu)
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
          call f4_update_ghostcells(f4, n_tvars, i_tvars, 0)
          call amr_flags_diff2(f4, min_level, max_level, &
               i_rho, c_refine, c_derefine, c_eps, c_abs)
          call f4_adjust_refinement(f4, load_imbalance_threshold)

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

  subroutine set_initial_conditions(f4, test_case)
    type(foap4_t), intent(inout) :: f4
    character(len=*), intent(in) :: test_case

    select case (test_case)
    case ("sod")
       call set_initial_conditions_sod(f4)
    case ("rt")
       call set_initial_conditions_rt(f4)
    case ("vortex")
       call set_initial_conditions_vortex(f4)
    case default
       error stop "Unknown test case, options: rt, sod, vortex"
    end select
  end subroutine set_initial_conditions

  subroutine set_initial_conditions_sod(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, ${IJK}$
    real(dp)                     :: rr(NDIM)
    real(fp)                     :: u0(n_tvars, 2)

    ! 1D Sod shock test case
    u0(i_rho, :) = [1.0_fp, 0.125_fp]
    u0(i_e, :)   = [1.0_fp, 0.1_fp]
    u0(i_mom0+1:i_mom0+NDIM, :) = 0.0_fp

    call to_conservative(u0(:, 1))
    call to_conservative(u0(:, 2))

    !$acc parallel loop default(present) copyin(u0)
    do n = 1, f4%n_blocks
       !$acc loop collapse(NDIM) private(rr)
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          rr = f4_cell_coord(f4, n, ${IJK}$)

          if (rr(1) < 0.5_dp) then
             f4%uu(${IJK}$, i_tvars, n) = u0(:, 1)
          else
             f4%uu(${IJK}$, i_tvars, n) = u0(:, 2)
          end if
       end do; ${KJI_CLOSE_LOOP}$
    end do
  end subroutine set_initial_conditions_sod

  subroutine set_initial_conditions_rt(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, ${IJK}$
    real(fp)                     :: rho_high, rho_low
    real(dp)                     :: p_interface, k_vec(NDIM-1)
    real(dp)                     :: rr(NDIM), h0, dh
    real(fp), parameter          :: pi = acos(-1.0_fp)

    ! The location of interface
    h0 = 0.8_fp

    ! Width of the sinusoidal fluctuations on the interface
    dh = 0.05_fp

    ! High and low density
    rho_high = 1.0_fp
    rho_low  = 0.1_fp

    ! Pressure at interface
    p_interface = 1.0_fp

    ! Wavelength
    k_vec(1) = 2 * pi
#:if NDIM == 3
    k_vec(2) = 2 * pi
#:endif

    !$acc parallel loop default(present) copyin(k_vec)
    do n = 1, f4%n_blocks
       !$acc loop collapse(NDIM) private(rr)
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          rr = f4_cell_coord(f4, n, ${IJK}$)

          f4%uu(${IJK}$, i_mom0+1:i_mom0+NDIM, n) = 0.0_fp

          if (rr(NDIM) > h0 + dh * product(sin(k_vec * rr(1:NDIM-1)))) then
             f4%uu(${IJK}$, i_rho, n) = rho_high
          else
             f4%uu(${IJK}$, i_rho, n) = rho_low
          end if

          f4%uu(${IJK}$, i_e, n) = real(euler_inv_gamma_m1 * (p_interface + &
               f4%uu(${IJK}$, i_rho, n) * euler_gravity * (rr(NDIM) - h0)), fp)
       end do; ${KJI_CLOSE_LOOP}$
    end do
  end subroutine set_initial_conditions_rt

  subroutine set_initial_conditions_vortex(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, ${IJK}$
    real(fp)                     :: rr(NDIM), rc(NDIM), r2
    real(fp)                     :: rho, p, T
    real(fp)                     :: v(NDIM), dv(NDIM), dT
    real(fp), parameter          :: pi = acos(-1.0_fp)

    ! Isentropic vortex parameters (Shu 1998)
    real(fp), parameter :: beta    = 5.0_fp ! Vortex strength
    real(fp), parameter :: rho_inf = 1.0_fp ! Free-stream density
    real(fp), parameter :: p_inf   = 1.0_fp ! Free-stream pressure
    real(fp), parameter :: L       = 5.0_fp ! Domain size
    real(fp)            :: v_inf(NDIM), T_inf, c_inf

    T_inf      = p_inf / rho_inf
    c_inf      = sqrt(euler_gamma * p_inf / rho_inf)
    v_inf(:)   = 0.0_fp
    v_inf(1:2) = 1.0_fp

    ! Vortex center (middle of domain)
    rc(:) = 5.0_fp

    !$acc parallel loop default(present) copyin(rc, v_inf)
    do n = 1, f4%n_blocks
       !$acc loop collapse(NDIM) private(rr, r2, dv, dT, T, rho, p, v)
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          rr = real(f4_cell_coord(f4, n, ${IJK}$), fp)

          ! Distance squared from vortex center
          r2 = (rr(1) - rc(1))**2 + (rr(2) - rc(2))**2

          ! Velocity perturbations
          dv(1) = -beta / (2.0_fp * pi) * (rr(2) - rc(2)) * exp(0.5_fp * (1.0_fp - r2))
          dv(2) = beta / (2.0_fp * pi) * (rr(1) - rc(1)) * exp(0.5_fp * (1.0_fp - r2))
#:if NDIM == 3
          dv(3) = 0.0_fp
#:endif

          ! Temperature perturbation
          dT = -(euler_gamma - 1.0_fp) * beta**2 / &
               (8.0_fp * euler_gamma * pi**2) * exp(1.0_fp - r2)

          ! Primitive variables
          T   = T_inf + dT
          rho = rho_inf * (T / T_inf)**euler_inv_gamma_m1
          p   = rho * T
          v  = v_inf + dv

          ! Set conservative variables
          f4%uu(${IJK}$, i_rho, n) = rho
          f4%uu(${IJK}$, i_mom0+1, n) = rho * v(1)
          f4%uu(${IJK}$, i_mom0+2, n) = rho * v(2)
#:if NDIM == 3
          f4%uu(${IJK}$, i_mom0+3, n) = rho * v(3)
#:endif
          ! Total energy: e = p/(gamma-1) + 0.5*rho*v^2
          f4%uu(${IJK}$, i_e, n) = p * euler_inv_gamma_m1 + &
               0.5_fp * rho * (v(1)**2 + v(2)**2)
       end do; ${KJI_CLOSE_LOOP}$
    end do
  end subroutine set_initial_conditions_vortex

  subroutine source_term(u_prim, source)
    real(fp), intent(in)    :: u_prim(n_tvars)
    real(fp), intent(inout) :: source(n_tvars)
    source(i_mom0+NDIM) = euler_gravity * u_prim(i_rho)
    source(i_e)         = euler_gravity * u_prim(i_rho) * u_prim(i_mom0+NDIM)
  end subroutine source_term

  #:include 'physics_euler.fpp'

  #:include 'flux_finite_volume.fpp'

  include 'flux_${FLUX_SCHEME}$.f90'

  include 'limiter_${LIMITER}$.f90'

end program
