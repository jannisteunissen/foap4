#:include 'definitions.fpp'
#:set LIMITER = 'vanleer'
#:set FLUX_SCHEME = 'tvdlf'

program euler
  use iso_fortran_env, only: int64
  use mpi_f08
  use m_foap4_2d
  use m_config
  use m_euler

  implicit none

  include 'limiters/${LIMITER}$_definitions.f90'
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: NDIM = ${NDIM}$
  integer, parameter :: n_gc = limiter_num_ghostcells

  type(foap4_t) :: f4
  type(CFG_t)   :: cfg

  integer  :: min_level   = 2
  integer  :: max_level   = 4
  integer  :: max_blocks  = 1000
  integer  :: bx(2)       = [32, 32]
  integer  :: num_outputs = 40
  logical  :: periodic(2) = .true.
  logical  :: do_refinement = .false.
  real(dp) :: end_time    = 2.0_dp
  real(dp) :: c_refine = 0.8_dp
  real(dp) :: c_derefine = 0.2_dp
  real(dp) :: c_eps = 0.01_dp
  integer :: n_steps_refinement = 4
  real(dp) :: cfl_number       = 0.5_dp
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
  call CFG_add_get(cfg, 'end_time', end_time, 'End time')
  call CFG_add_get(cfg, 'test_case', test_case, 'Which test case to run')
  call CFG_add_get(cfg, 'time_integrator', integrator_name, 'Time integrator')
  call CFG_add_get(cfg, 'cfl_number', cfl_number, 'CFL number')
  call CFG_check(cfg)

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

    if (test_case == "rt") periodic(2) = .false.

    call f4_construct_brick(f4, n_blocks_per_dim, block_length, bx, n_gc, &
         n_vars, var_names, temporal, n_time_states, periodic, &
         min_level, max_blocks, f4_bc_neumann, 0.0_dp)

    if (test_case == "rt") then
       call f4_set_physical_boundary(f4, i_mom(2), 2, f4_bc_dirichlet, 0.0_dp)
       call f4_set_physical_boundary(f4, i_mom(2), 3, f4_bc_dirichlet, 0.0_dp)
    end if

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

    if (dt_output < end_time) call f4_write_grid(f4, base_name, n_output, time)
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
    integer                      :: n
    real(dp)                     :: u0(n_vars, 4)

    select case (test_case)
    case ("first")
       u0(i_e, :)      = [1.0_dp, 0.4_dp, 0.0439_dp, 0.15_dp]
       u0(i_rho, :)    = [1.0_dp, 0.5197_dp, 0.1072_dp, 0.2579_dp]
       u0(i_mom(1), :) = [0.0_dp, -0.7259_dp, -0.7259_dp, 0.0_dp]
       u0(i_mom(2), :) = [0.0_dp, 0.0_dp, -1.4045_dp, -1.4045_dp]

       do n = 1, 4
          call to_conservative(u0(:, n))
       end do
       call set_initial_quadrants(f4, u0)
    case ("sixth")
       u0(i_e, :)      = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
       u0(i_rho, :)    = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
       u0(i_mom(1), :) = [0.75_dp, 0.75_dp, -0.75_dp, -0.75_dp]
       u0(i_mom(2), :) = [-0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp]

       do n = 1, 4
          call to_conservative(u0(:, n))
       end do
       call set_initial_quadrants(f4, u0)
    case ("sod")
       ! 1D Sod shock test case
       u0(i_rho, :)    = [0.125_dp, 1.0_dp, 1.0_dp, 0.125_dp]
       u0(i_e, :)      = [0.1_dp, 1.0_dp, 1.0_dp, 0.1_dp]
       u0(i_mom(1), :) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
       u0(i_mom(2), :) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]

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
    real(dp), intent(in)         :: u0(n_vars, 4)
    integer                      :: n, i, j
    real(dp)                     :: rr(2)

    !$acc parallel loop
    do n = 1, f4%n_blocks
       !$acc loop collapse(2) private(rr)
       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)
             rr = f4_cell_coord(f4, n, i, j)

             if (rr(1) > 0.5_dp .and. rr(2) > 0.5_dp) then
                f4%uu(i, j, i_vars, n) = u0(:, 1)
             elseif (rr(1) <= 0.5_dp .and. rr(2) >= 0.5_dp) then
                f4%uu(i, j, i_vars, n) = u0(:, 2)
             elseif (rr(1) <= 0.5_dp .and. rr(2) <= 0.5_dp) then
                f4%uu(i, j, i_vars, n) = u0(:, 3)
             else
                f4%uu(i, j, i_vars, n) = u0(:, 4)
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

             f4%uu(i, j, i_mom(1), n) = 0.0_dp
             f4%uu(i, j, i_mom(2), n) = 0.0_dp

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

  #:include 'physics/euler.fpp'

  #:include 'flux_schemes/finite_volume.fpp'

  include 'flux_schemes/${FLUX_SCHEME}$.f90'

  include 'limiters/${LIMITER}$.f90'

end program
