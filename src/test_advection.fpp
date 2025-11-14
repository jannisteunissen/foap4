#:include 'definitions.fpp'

#:set LIMITER = getvar('LIMITER', 'vanleer')
#:set FLUX_SCHEME = getvar('FLUX_SCHEME', 'tvdlf')

program test_adv
  use iso_fortran_env, only: int64
  use mpi_f08
  use m_foap4_${NDIM}$d
  use m_config
  use m_advection

  implicit none

  include 'limiters/${LIMITER}$_definitions.f90'
  integer, parameter :: dp   = kind(0.0d0)
  integer, parameter :: NDIM = ${NDIM}$
  integer, parameter :: n_gc = limiter_num_ghostcells

  logical, parameter :: temporal(n_vars) = .true.
  real(dp)           :: cfl_number       = 0.5_dp
  real(dp)           :: end_time         = 1.0_dp
  real(dp)           :: c_refine         = 0.8_dp
  real(dp)           :: c_derefine       = 0.2_dp
  real(dp)           :: c_eps            = 0.01_dp

  logical           :: do_refinement        = .true.
  integer           :: max_refinement_level = 5
  integer           :: min_refinement_level = 2
  integer           :: max_blocks           = 2000
  integer           :: bx(NDIM)             = 32
  integer           :: num_outputs          = 40
  character(len=40) :: integrator_name      = "heuns_method"

  type(foap4_t) :: f4
  type(CFG_t) :: cfg

  call f4_initialize(f4, "error")

  call CFG_update_from_arguments(cfg)
  call CFG_add_get(cfg, 'num_outputs', num_outputs, 'Write this many output files')
  call CFG_add_get(cfg, 'do_refinement', do_refinement, 'Perform refinement')
  call CFG_add_get(cfg, 'min_refinement_level', min_refinement_level, &
       'Minimum refinement level in the domain')
  call CFG_add_get(cfg, 'max_refinement_level', max_refinement_level, &
       'Maximum refinement level in the domain')
  call CFG_add_get(cfg, 'c_refine', c_refine, 'Coefficient for refinement')
  call CFG_add_get(cfg, 'c_derefine', c_refine, 'Coefficient for derefinement')
  call CFG_add_get(cfg, 'c_eps', c_refine, 'Used in refinement criterion')
  call CFG_add_get(cfg, 'bx', bx, 'Size of grid blocks')
  call CFG_add_get(cfg, 'max_blocks', max_blocks, 'Max. number of blocks')
  call CFG_add_get(cfg, 'velocity', velocity(1:NDIM), 'Velocity')
  !$acc update device(velocity(1:NDIM))
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
    integer, parameter           :: n_blocks_per_dim(NDIM) = 1
    real(dp), parameter          :: block_length(NDIM) = 1.0_dp
    logical, parameter           :: periodic(NDIM) = .true.
    real(dp), parameter          :: cfl_number = 0.5_dp
    integer                      :: n, prev_mesh_revision, n_output
    integer                      :: highest_level, prev_highest_level, n_iterations, ierr
    integer(int64)               :: sum_local_blocks, sum_global_blocks
    logical                      :: write_this_step
    integer                      :: integrator, n_time_states
    real(dp)                     :: dt, dt_lim, dt_output
    real(dp)                     :: time, t0, t1
    real(dp)                     :: rho_initial_sum, rho_sum

    time = 0.0_dp
    dt_lim = 0.0_dp
    dt_output = end_time / max(real(num_outputs, dp), 1e-100_dp)
    n_output = 0
    n_iterations = 0
    sum_local_blocks = 0

    integrator = f4_get_time_integrator(trim(integrator_name))
    n_time_states = f4_advance_num_copies(integrator)

    call f4_construct_brick(f4, n_blocks_per_dim, block_length, bx, n_gc, &
         n_vars, var_names, temporal, n_time_states, periodic, &
         min_refinement_level, max_blocks, f4_bc_dirichlet, 0.0_dp)

    call set_init_cond(f4)

    if (do_refinement) then
       do n = 1, 10
          prev_mesh_revision = f4_get_mesh_revision(f4)
          call f4_update_ghostcells(f4, 1, i_vars)
          call f4_set_refinement_flags_diff2(f4, min_refinement_level, &
               max_refinement_level, i_rho, c_refine, c_derefine, c_eps)
          call f4_adjust_refinement(f4, .true.)
          call set_init_cond(f4)

          if (f4_get_mesh_revision(f4) == prev_mesh_revision) exit
       end do
    end if

    call f4_compute_sum(f4, i_rho, rho_initial_sum)
    call f4_get_global_highest_level(f4, prev_highest_level)

    if (dt_output <= end_time) call f4_write_grid(f4, base_name, n_output, time)
    n_output = n_output + 1

    t0 = MPI_Wtime()

    do while (time <= end_time)
       n_iterations = n_iterations + 1
       dt = cfl_number * dt_lim

       write_this_step = (time + dt >= n_output * dt_output)
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

       if (do_refinement) then
          call f4_update_ghostcells(f4, 1, i_vars)
          call f4_set_refinement_flags_diff2(f4, min_refinement_level, &
               max_refinement_level, i_rho, c_refine, c_derefine, c_eps)
          call f4_adjust_refinement(f4, .true.)
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

    !$acc parallel loop
    do n = 1, f4%n_blocks
       !$acc loop collapse(NDIM) private(rr)
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          rr = f4_cell_coord(f4, n, ${IJK}$)
          f4%uu(${IJK}$, i_rho, n) = rho_init(@{DINDEX(rr)}@)
       end do; ${KJI_CLOSE_LOOP}$
    end do
  end subroutine set_init_cond

#:if NDIM == 2
  pure real(dp) function rho_init(x, y)
    !$acc routine seq
    real(dp), intent(in) :: x, y

    if (sqrt((x-0.5_dp)**2 + (y-0.5_dp)**2) < 0.1_dp) then
       rho_init = 1.0_dp
    else
       rho_init = 0.0_dp
    end if
  end function rho_init
#:elif NDIM == 3
  pure real(dp) function rho_init(x, y, z)
    !$acc routine seq
    real(dp), intent(in) :: x, y, z

    if (sqrt((x-0.5_dp)**2 + (y-0.5_dp)**2 + (z-0.5_dp)**2) < 0.1_dp) then
       rho_init = 1.0_dp
    else
       rho_init = 0.0_dp
    end if
  end function rho_init
#:endif

  pure subroutine get_flux(flux_dim, u, flux)
    !$acc routine seq
    integer, intent(in)   :: flux_dim
    real(dp), intent(in)  :: u(n_vars)
    real(dp), intent(out) :: flux(n_vars)
    flux(1) = velocity(flux_dim) * u(1)
  end subroutine get_flux

  pure subroutine to_primitive(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(n_vars)
  end subroutine to_primitive

  pure subroutine get_max_wavespeed(flux_dim, u, cmax)
    !$acc routine seq
    integer, intent(in)   :: flux_dim
    real(dp), intent(in)  :: u(n_vars)
    real(dp), intent(out) :: cmax
    cmax = abs(velocity(flux_dim))
  end subroutine get_max_wavespeed

  pure subroutine get_min_max_wavespeed(flux_dim, n_vars, u, cmin, cmax)
    !$acc routine seq
    integer, intent(in)   :: flux_dim
    integer, intent(in)   :: n_vars
    real(dp), intent(in)  :: u(n_vars)
    real(dp), intent(out) :: cmin
    real(dp), intent(out) :: cmax
    cmin = min(velocity(flux_dim), 0.0_dp)
    cmax = max(velocity(flux_dim), 0.0_dp)
  end subroutine get_min_max_wavespeed

#:include 'flux_schemes/finite_volume.fpp'

  include 'flux_schemes/${FLUX_SCHEME}$.f90'

  include 'limiters/${LIMITER}$.f90'

end program
