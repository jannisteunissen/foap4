#:include 'definitions.fpp'
#:set LIMITER = 'vanleer'
#:set FLUX_SCHEME = 'tvdlf'

program euler
  use iso_fortran_env, only: int64
  use mpi_f08
  use m_foap4_${NDIM}$d
  use m_config
  use m_physics_euler_${NDIM}$d

  implicit none

  include 'limiter_${LIMITER}$_definitions.f90'
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: NDIM = ${NDIM}$
  integer, parameter :: n_gc = limiter_num_ghostcells

  type(foap4_t) :: f4
  type(CFG_t)   :: cfg

  integer           :: min_level          = 2
  integer           :: max_level          = 4
  integer           :: max_blocks         = 1000
  integer           :: blocks_per_dim(NDIM) = 1
  integer           :: bx(NDIM)           = 32
  integer           :: num_outputs        = 40
  logical           :: periodic(NDIM)     = .true.
  logical           :: do_refinement      = .false.
  real(dp)          :: end_time           = 2.0_dp
  real(dp)          :: c_refine           = 0.8_dp
  real(dp)          :: c_derefine         = 0.2_dp
  real(dp)          :: c_eps              = 0.01_dp
  integer           :: n_steps_refinement = 4
  real(dp)          :: cfl_number         = 0.5_dp
  character(len=10) :: test_case          = "sod"
  character(len=40) :: integrator_name    = "heuns_method"

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
  call CFG_add_get(cfg, 'blocks_per_dim', blocks_per_dim, &
       'Number of blocks (per dimension) on coarse grid')
  call CFG_add_get(cfg, 'end_time', end_time, 'End time')
  call CFG_add_get(cfg, 'test_case', test_case, 'Which test case to run')
  call CFG_add_get(cfg, 'time_integrator', integrator_name, 'Time integrator')
  call CFG_add_get(cfg, 'cfl_number', cfl_number, 'CFL number')
  call CFG_check(cfg)

  call test_euler(f4, bx, min_level, max_blocks, &
       num_outputs, "output/test_euler_${NDIM}$d", test_case, &
       end_time, integrator_name)

  if (f4%mpirank == 0) call f4_print_wtime(f4)
  call f4_finalize(f4)

contains

  subroutine test_euler(f4, bx, min_level, max_blocks, num_outputs, base_name, &
       test_case, end_time, integrator_name)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: bx(NDIM)
    integer, intent(in)          :: min_level
    integer, intent(in)          :: max_blocks
    integer, intent(in)          :: num_outputs
    character(len=*), intent(in) :: base_name
    character(len=*), intent(in) :: test_case
    real(dp), intent(in)         :: end_time
    character(len=40), intent(in) :: integrator_name
    real(dp), parameter          :: block_length(NDIM) = 1.0_dp
    logical                      :: periodic(NDIM) = .true.
    integer                      :: n_output
    integer                      :: n, n_iterations, ierr, prev_mesh_revision
    integer(int64)               :: sum_local_blocks, sum_global_blocks
    logical                      :: write_this_step, temporal(n_vars)
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

    if (test_case == "rt") periodic(NDIM) = .false.

    call f4_construct_brick(f4, blocks_per_dim, block_length, bx, n_gc, &
         n_vars, var_names, temporal, n_time_states, periodic, &
         min_level, max_blocks, f4_bc_neumann, 0.0_dp)

    if (test_case == "rt") then
       call f4_set_bc_scalar(f4, i_mom0+NDIM, 2*(NDIM-1), &
            f4_bc_dirichlet, 0.0_dp)
       call f4_set_bc_scalar(f4, i_mom0+NDIM, 2*(NDIM-1)+1, &
            f4_bc_dirichlet, 0.0_dp)
       gravity_constant = 1.0_dp
    else
       gravity_constant = 0.0_dp
    end if

    !$acc update device(gravity_constant)

    call set_initial_conditions(f4, test_case)

    if (do_refinement) then
       do n = 1, 10
          prev_mesh_revision = f4_get_mesh_revision(f4)
          call f4_update_ghostcells(f4, n_vars, i_vars)
          call f4_set_refinement_flags_diff2(f4, min_level, max_level, &
               i_rho, c_refine, c_derefine, c_eps)
          call f4_adjust_refinement(f4, .true.)
          call set_initial_conditions(f4, test_case)

          if (f4_get_mesh_revision(f4) == prev_mesh_revision) exit
       end do
    end if

    call f4_compute_sum(f4, i_rho, rho_initial_sum)

    if (dt_output <= end_time) call f4_write_grid(f4, base_name, n_output, time)
    n_output = n_output + 1

    t0 = MPI_Wtime()

    do while (time < end_time)
       n_iterations = n_iterations + 1
       dt = cfl_number * dt_lim

       write_this_step = (time + dt > n_output * dt_output)
       if (write_this_step) dt = n_output * dt_output - time
       call f4_advance(f4, dt, dt_lim, time, integrator, feuler_finite_volume)

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
          call f4_update_ghostcells(f4, n_vars, i_vars)
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

  subroutine set_initial_conditions(f4, test_case)
    type(foap4_t), intent(inout) :: f4
    character(len=*), intent(in) :: test_case

    select case (test_case)
    case ("sod")
       call set_initial_conditions_sod(f4)
    case ("rt")
       call set_initial_conditions_rt(f4)
    case default
       error stop "Unknown test case, options: rt, sod"
    end select
  end subroutine set_initial_conditions

  subroutine set_initial_conditions_sod(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, ${IJK}$
    real(dp)                     :: rr(NDIM)
    real(dp)                     :: u0(n_vars, 2)

    ! 1D Sod shock test case
    u0(i_rho, :) = [1.0_dp, 0.125_dp]
    u0(i_e, :)   = [1.0_dp, 0.1_dp]
    u0(i_mom0+1:i_mom0+NDIM, :) = 0.0_dp

    call to_conservative(u0(:, 1))
    call to_conservative(u0(:, 2))

    !$acc parallel loop
    do n = 1, f4%n_blocks
       !$acc loop collapse(NDIM) private(rr)
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          rr = f4_cell_coord(f4, n, ${IJK}$)

          if (rr(1) < 0.5_dp) then
             f4%uu(${IJK}$, i_vars, n) = u0(:, 1)
          else
             f4%uu(${IJK}$, i_vars, n) = u0(:, 2)
          end if
       end do; ${KJI_CLOSE_LOOP}$
    end do
  end subroutine set_initial_conditions_sod

  subroutine set_initial_conditions_rt(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, ${IJK}$
    real(dp)                     :: h0, dh, rho_high, rho_low
    real(dp)                     :: p_interface, k_vec(NDIM-1), rr(NDIM)
    real(dp), parameter          :: pi = acos(-1.0_dp)

    ! The location of interface
    h0 = 0.8d0

    ! Width of the sinusoidal fluctuations on the interface
    dh = 0.05d0

    ! High and low density
    rho_high = 1.0_dp
    rho_low  = 0.1_dp

    ! Pressure at interface
    p_interface = 1.0_dp

    ! Wavelength
    k_vec(1) = 2 * pi
#:if NDIM == 3
    k_vec(2) = 4 * pi
#:endif

    !$acc parallel loop
    do n = 1, f4%n_blocks
       !$acc loop collapse(NDIM) private(rr)
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          rr = f4_cell_coord(f4, n, ${IJK}$)

          f4%uu(${IJK}$, i_mom0+1:i_mom0+NDIM, n) = 0.0_dp

          if (rr(NDIM) > h0 + dh * product(sin(k_vec * rr(1:NDIM-1)))) then
             f4%uu(${IJK}$, i_rho, n) = rho_high
          else
             f4%uu(${IJK}$, i_rho, n) = rho_low
          end if

          f4%uu(${IJK}$, i_e, n) = inv_gamma_m1 * (p_interface - 1.0_dp * &
               f4%uu(${IJK}$, i_rho, n) * (rr(NDIM) - h0))
       end do; ${KJI_CLOSE_LOOP}$
    end do
  end subroutine set_initial_conditions_rt

  subroutine source_term(u_prim, source)
    real(dp), intent(in)    :: u_prim(n_vars)
    real(dp), intent(inout) :: source(n_vars)
    source(i_mom0+NDIM) = -gravity_constant * u_prim(i_rho)
    source(i_e)         = -gravity_constant * source(i_mom0+NDIM)
  end subroutine source_term

  #:include 'physics_euler.fpp'

  #:include 'flux_finite_volume.fpp'

  include 'flux_${FLUX_SCHEME}$.f90'

  include 'limiter_${LIMITER}$.f90'

end program
