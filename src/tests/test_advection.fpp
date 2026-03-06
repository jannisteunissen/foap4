#:include 'definitions.fpp'
#:set LIMITER = getvar('USE_LIMITER', 'weno5')
#:set FLUX_SCHEME = getvar('USE_FLUX_SCHEME', 'hll')

program test_adv
  use iso_fortran_env, only: int64
  use mpi_f08
  use m_foap4_${NDIM}$d
  use m_physics_advection_${NDIM}$d
  use m_rk_${NDIM}$d
  use m_io_${NDIM}$d
  use m_amr_flags_${NDIM}$d
  use m_config

  implicit none

  include 'limiter_${LIMITER}$_definitions.f90'
  integer, parameter :: dp   = kind(0.0d0)
  integer, parameter :: n_gc = limiter_num_ghostcells

  real(dp) :: velocity(NDIM) = 1.0_dp
  real(dp) :: cfl_number     = 0.5_dp
  real(dp) :: end_time       = 1.0_dp
  real(dp) :: c_refine       = 0.1_dp
  real(dp) :: c_derefine     = 0.0125_dp
  real(dp) :: c_eps          = 0.01_dp
  real(dp) :: c_abs          = 1e-6_dp

  logical           :: do_refinement        = .true.
  integer           :: max_refinement_level = 3
  integer           :: min_refinement_level = 1
  integer           :: n_steps_refinement   = 4
  integer           :: max_blocks           = 2000
  integer           :: blocks_per_dim(NDIM) = 1
  integer           :: bx(NDIM)             = 32
  integer           :: num_outputs          = 40
  real(dp)          :: load_imbalance_threshold = 1.1_dp
  character(len=40) :: integrator_name      = "heuns_method"
  logical           :: write_vtu = .false.

  type(foap4_t) :: f4
  type(CFG_t) :: cfg

  call f4_initialize(f4, "error")

  call CFG_update_from_arguments(cfg)
  call CFG_add_get(cfg, 'num_outputs', num_outputs, 'Write this many output files')
  call CFG_add_get(cfg, 'write_vtu', write_vtu, 'Also write p4est vtu files')
  call CFG_add_get(cfg, 'do_refinement', do_refinement, 'Perform refinement')
  call CFG_add_get(cfg, 'load_imbalance_threshold', load_imbalance_threshold, &
       'Threshold for partitioning')
  call CFG_add_get(cfg, 'min_level', min_refinement_level, &
       'Minimum refinement level in the domain')
  call CFG_add_get(cfg, 'max_level', max_refinement_level, &
       'Maximum refinement level in the domain')
  call CFG_add_get(cfg, 'c_refine', c_refine, 'Coefficient for refinement')
  call CFG_add_get(cfg, 'c_derefine', c_derefine, 'Coefficient for derefinement')
  call CFG_add_get(cfg, 'c_eps', c_eps, 'Filter coefficient for AMR')
  call CFG_add_get(cfg, 'c_abs', c_abs, 'Density threshold for AMR')
  call CFG_add_get(cfg, 'n_steps_refinement', n_steps_refinement, &
       'Perform refinement every N steps')
  call CFG_add_get(cfg, 'bx', bx, 'Size of grid blocks')
  call CFG_add_get(cfg, 'blocks_per_dim', blocks_per_dim, &
       'Number of blocks (per dimension) on coarse grid')
  call CFG_add_get(cfg, 'max_blocks', max_blocks, 'Max. number of blocks')
  call CFG_add_get(cfg, 'velocity', velocity, 'Velocity')
  call CFG_add_get(cfg, 'end_time', end_time, 'End time')
  call CFG_add_get(cfg, 'time_integrator', integrator_name, 'Time integrator')
  call CFG_add_get(cfg, 'cfl_number', cfl_number, 'CFL number')
  call CFG_check(cfg)

  if (max_refinement_level < min_refinement_level) &
       error stop "max_refinement_level < min_refinement_level"

  call test_advection(f4, bx, do_refinement, max_blocks, &
       num_outputs, "output/test_adv_${NDIM}$d", end_time, integrator_name)

  if (f4%mpirank == 0) call f4_print_wtime(f4)
  call f4_finalize(f4)

contains

  subroutine test_advection(f4, bx, do_refinement, &
       max_blocks, num_outputs, base_name, end_time, integrator_name)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: bx(NDIM)
    logical, intent(in)          :: do_refinement
    integer, intent(in)          :: max_blocks
    integer, intent(in)          :: num_outputs
    character(len=*), intent(in) :: base_name
    real(dp), intent(in)         :: end_time
    character(len=40), intent(in) :: integrator_name
    real(dp), parameter          :: block_length(NDIM) = 1.0_dp
    logical, parameter           :: periodic(NDIM) = .true.
    real(dp), parameter          :: cfl_number = 0.5_dp
    integer                      :: n, prev_mesh_revision, n_output
    integer                      :: highest_level, prev_highest_level, n_iterations, ierr
    integer(int64)               :: sum_local_blocks, sum_global_blocks
    logical                      :: write_this_step
    integer                      :: integrator, n_time_states
    real(dp)                     :: dt, dt_lim, dt_output
    real(dp)                     :: t0, t1
    real(dp)                     :: rho_initial_sum, rho_sum

    call advection_initialize(velocity)

    f4%time = 0.0_dp
    dt_lim = 0.0_dp
    dt_output = end_time / max(real(num_outputs, dp), 1e-100_dp)
    n_output = 0
    n_iterations = 0
    sum_local_blocks = 0

    integrator = rk_get_integrator_by_name(trim(integrator_name))
    n_time_states = rk_advance_num_copies(integrator)

    call f4_construct_brick(f4, blocks_per_dim, block_length, bx, n_gc, &
         n_vars_all, var_names, var_temporal, n_time_states, periodic, &
         min_refinement_level, max_blocks, f4_bc_dirichlet, 0.0_dp)

    call set_init_cond(f4)

    if (do_refinement) then
       do n = 1, 10
          prev_mesh_revision = f4_get_mesh_revision(f4)
          call f4_update_ghostcells(f4, n_tvars, i_tvars)
          call amr_flags_diff2(f4, min_refinement_level, max_refinement_level, &
               i_rho, c_refine, c_derefine, c_eps, c_abs)
          call f4_adjust_refinement(f4, load_imbalance_threshold)
          call set_init_cond(f4)

          if (f4_get_mesh_revision(f4) == prev_mesh_revision) exit
       end do
    end if

    call f4_compute_sum(f4, i_rho, rho_initial_sum)
    call f4_get_global_highest_level(f4, prev_highest_level)

    if (dt_output <= end_time) then
       call io_write_grid(f4, base_name, n_output, write_p4vtu=write_vtu)
    end if
    n_output = n_output + 1

    t0 = MPI_Wtime()

    do while (f4%time <= end_time)
       n_iterations = n_iterations + 1
       dt = cfl_number * dt_lim

       write_this_step = (f4%time + dt >= n_output * dt_output)
       if (write_this_step) dt = n_output * dt_output - f4%time

       call rk_advance(f4, dt, dt_lim, integrator, feuler_finite_volume)

       if (write_this_step) then
          call set_error(f4)
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
          call f4_update_ghostcells(f4, n_tvars, i_tvars)
          call amr_flags_diff2(f4, min_refinement_level, max_refinement_level, &
               i_rho, c_refine, c_derefine, c_eps, c_abs)
          call f4_adjust_refinement(f4, load_imbalance_threshold)
          call f4_get_global_highest_level(f4, highest_level)

          if (highest_level > prev_highest_level) then
             dt_lim = 0.5_dp * dt_lim
          end if
          prev_highest_level = highest_level
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
       write(*, "(A,F14.3)") " cell updates/ns: ", 1e-9_dp * &
            sum_global_blocks * (product(f4%bx) * 2 / (t1 - t0))
    end if

    call f4_destroy(f4)
  end subroutine test_advection

  subroutine set_init_cond(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, ${IJK}$
    real(dp)                     :: rr(NDIM)

    !$acc parallel loop default(present)
    do n = 1, f4%n_blocks
       !$acc loop collapse(NDIM) private(rr)
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          rr = f4_cell_coord(f4, n, ${IJK}$)
#:if NDIM == 2
          f4%uu(${IJK}$, i_rho, n) = rho_solution(rr(1), rr(2), 0.0_dp)
#:elif NDIM == 3
          f4%uu(${IJK}$, i_rho, n) = rho_solution(rr(1), rr(2), rr(3), 0.0_dp)
#:endif
       end do; ${KJI_CLOSE_LOOP}$
    end do
  end subroutine set_init_cond

  subroutine set_error(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, ${IJK}$
    real(dp)                     :: rr(NDIM), sol

    !$acc parallel loop default(present)
    do n = 1, f4%n_blocks
       !$acc loop collapse(NDIM) private(rr, sol)
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          rr = f4_cell_coord(f4, n, ${IJK}$)
#:if NDIM == 2
          sol = rho_solution(rr(1), rr(2), f4%time)
#:elif NDIM == 3
          sol = rho_solution(rr(1), rr(2), rr(3), f4%time)
#:endif
          f4%uu(${IJK}$, i_error, n) = f4%uu(${IJK}$, i_rho, n) - sol
       end do; ${KJI_CLOSE_LOOP}$
    end do
  end subroutine set_error

  pure real(dp) function rho_solution(${XYZ}$, t)
    !$acc routine seq
    real(dp), intent(in) :: ${XYZ}$
    real(dp), intent(in) :: t
    real(dp)             :: distance, q
    real(dp), parameter  :: radius = 0.1_dp
    real(dp), parameter  :: border = 0.05_dp
#:if NDIM == 2
    real(dp) :: x0, y0
#:elif NDIM == 3
    real(dp) :: x0, y0, z0
#:endif

#:if NDIM == 2
    x0 = modulo(x - advection_velocity(1) * t, 1.0_dp)
    y0 = modulo(y - advection_velocity(2) * t, 1.0_dp)
    distance = sqrt((x0 - 0.5_dp)**2 + (y0 - 0.5_dp)**2)
#:elif NDIM == 3
    x0 = modulo(x - advection_velocity(1) * t, 1.0_dp)
    y0 = modulo(y - advection_velocity(2) * t, 1.0_dp)
    z0 = modulo(z - advection_velocity(3) * t, 1.0_dp)
    distance = sqrt((x - 0.5_dp)**2 + (y - 0.5_dp)**2 + (z - 0.5_dp)**2)
#:endif

    if (distance < radius - border) then
       rho_solution = 1.0_dp
    else if (distance < radius) then
       ! cubic smoothstep: 1 - 3 q^2 + 2 q^3, with q in [0,1]
       q = (distance - radius + border)/border
       rho_solution = 1.0_dp - (3.0_dp * q**2 - 2.0_dp * q**3)
    else
       rho_solution = 0.0_dp
    end if
  end function rho_solution

#:include 'physics_advection.fpp'

#:include 'flux_finite_volume.fpp'

  include 'flux_${FLUX_SCHEME}$.f90'

  include 'limiter_${LIMITER}$.f90'

end program
