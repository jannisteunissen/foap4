#:include 'definitions_ndim.fpp'
#:include 'definitions_parallel.fpp'
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
  integer, parameter :: n_gc = limiter_num_ghostcells

  real(dp) :: velocity(NDIM) = 1.0_dp
  real(dp) :: r0(NDIM)       = 0.5_dp
  real(dp) :: cfl_number     = 0.5_dp
  real(dp) :: dt_max         = 5e-2_dp
  real(dp) :: end_time       = 1.0_dp
  real(dp) :: c_refine       = 0.1_dp
  real(dp) :: c_derefine     = 0.0125_dp
  real(dp) :: c_eps          = 0.01_dp
  real(dp) :: c_abs          = 1e-6_dp

  logical           :: do_refinement            = .true.
  integer           :: max_refinement_level     = ndim - 1
  integer           :: min_refinement_level     = 1
  integer           :: n_steps_refinement       = 4
  integer           :: max_blocks               = 2000
  integer           :: blocks_per_dim(NDIM)     = 1
  integer           :: bx(NDIM)                 = 16
  integer           :: num_outputs              = 4
  integer           :: n_gc_out                 = 1
  integer           :: velocity_type            = 1
  real(dp)          :: load_imbalance_threshold = 1.1_dp
  character(len=40) :: integrator_name          = "heuns_method"
  character(len=40) :: viewer                   = "visit"
  logical           :: write_vtu                = .false.
  logical           :: use_gaussian             = .false.

  type(foap4_t) :: f4
  type(CFG_t) :: cfg

  call f4_initialize(f4, "error")

  call CFG_update_from_arguments(cfg)
  call CFG_add_get(cfg, 'num_outputs', num_outputs, 'Write this many output files')
  call CFG_add_get(cfg, 'write_vtu', write_vtu, 'Also write p4est vtu files')
  call CFG_add_get(cfg, 'use_gaussian', use_gaussian, 'Use Gaussian solution')
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
  call CFG_add_get(cfg, 'velocity_type', velocity_type, &
       'Velocity type (1: uniform, 2: rotation, 3: deform)')
  call CFG_add_get(cfg, 'velocity', velocity, 'Velocity (for uniform profile)')
  call CFG_add_get(cfg, 'r0', r0, 'Center of initial solution')
  call CFG_add_get(cfg, 'end_time', end_time, 'End time')
  call CFG_add_get(cfg, 'time_integrator', integrator_name, 'Time integrator')
  call CFG_add_get(cfg, 'cfl_number', cfl_number, 'CFL number')
  call CFG_add_get(cfg, 'dt_max', dt_max, 'Max. dt (important for deform)')
  call CFG_add_get(cfg, 'viewer', viewer, &
       'Write XDMF output for this viewer (visit or paraview)')
  call CFG_check(cfg)

  if (max_refinement_level < min_refinement_level) &
       error stop "max_refinement_level < min_refinement_level"
  if (velocity_type < 1 .or. velocity_type > 3) &
       error stop "velocity type should be between 1 and 3"

  call test_advection(f4, bx, do_refinement, max_blocks, &
       num_outputs, "output/test_adv_${NDIM}$d", end_time, integrator_name)

  call f4_print_wtime(f4)
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
    integer                      :: n, prev_mesh_revision, n_output
    integer                      :: highest_level, prev_highest_level, n_iterations, ierr
    integer(int64)               :: sum_local_blocks, sum_global_blocks
    logical                      :: write_this_step
    integer                      :: integrator, n_time_states
    real(dp)                     :: dt, dt_lim, dt_output
    real(dp)                     :: t0, t1
    real(dp)                     :: rho_initial_sum, rho_sum, l1_err, l2_err

    call advection_initialize(velocity, use_gaussian, velocity_type, r0)

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
         min_refinement_level, max_blocks, f4_bc_dirichlet, 0.0_dp, .true.)

    call set_init_cond(f4)

    if (do_refinement) then
       do n = 1, 10
          prev_mesh_revision = f4_get_mesh_revision(f4)
          call f4_update_ghostcells(f4, n_tvars, i_tvars, 0)
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
       call io_write_grid(f4, base_name, n_output, write_p4vtu=write_vtu, &
            viewer=viewer, n_gc_out=n_gc_out)
    end if
    n_output = n_output + 1

    t0 = MPI_Wtime()

    do while (f4%time < end_time)
       n_iterations = n_iterations + 1
       dt = min(cfl_number * dt_lim, dt_max)
       write_this_step = (f4%time + dt >= n_output * dt_output)
       if (write_this_step) dt = n_output * dt_output - f4%time

       call rk_advance(f4, dt, dt_lim, integrator, feuler_finite_volume)

       if (write_this_step) then
          call set_error(f4)
          call compute_error_norms(f4, i_error, l1_err, l2_err)
          call io_write_grid(f4, base_name, n_output, write_p4vtu=write_vtu, &
               viewer=viewer, n_gc_out=n_gc_out)
          call f4_compute_sum(f4, i_rho, rho_sum)
          if (f4%mpirank == 0) then
             write(*, "(A,E12.4)") " Conservation error: ", &
                  rho_sum - rho_initial_sum
             write(*, "(A,E12.4)") " L1 error: ", l1_err
             write(*, "(A,E12.4)") " L2 error: ", l2_err
          end if
          n_output = n_output + 1
       end if

       if (do_refinement .and. &
            mod(n_iterations, n_steps_refinement) == 0) then
          call f4_update_ghostcells(f4, n_tvars, i_tvars, 0)
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

    ${PARALLEL_LOOP_FLAT('collapse(NDIM+1) private(rr)')}$ ${DEFAULT_PRESENT()}$
    do n = 1, f4%n_blocks
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          call f4_cell_coord(f4, n, ${IJK}$, rr)
          call rho_solution(@{DINDEX(rr)}@, 0.0_dp, f4%uu(${IJK}$, i_rho, n))
       end do; ${KJI_CLOSE_LOOP}$
    end do
  end subroutine set_init_cond

  subroutine set_error(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, ${IJK}$
    real(dp)                     :: rr(NDIM)
    real(fp)                     :: sol

    ${PARALLEL_LOOP_FLAT('collapse(NDIM+1) private(rr, sol)')}$ ${DEFAULT_PRESENT()}$
    do n = 1, f4%n_blocks
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          call f4_cell_coord(f4, n, ${IJK}$, rr)
          call rho_solution(@{DINDEX(rr)}@, f4%time, sol)
          f4%uu(${IJK}$, i_error, n) = f4%uu(${IJK}$, i_rho, n) - sol
       end do; ${KJI_CLOSE_LOOP}$
    end do
  end subroutine set_error

  subroutine compute_error_norms(f4, i_err, l1_err, l2_err)
    type(foap4_t), intent(in) :: f4
    integer, intent(in)       :: i_err
    real(dp), intent(out)     :: l1_err, l2_err
    integer                   :: level, ${IJK}$, n, ierror
    real(dp)                  :: dvol, my_l1, my_l2

    my_l1 = 0.0_dp
    my_l2 = 0.0_dp

    ${PARALLEL_LOOP_FLAT('collapse(ndim+1) private(level, dvol) reduction(+:my_l1, my_l2)')}$ ${DEFAULT_PRESENT()}$
    do n = 1, f4%n_blocks
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          level = f4%block_level(n)
#:if NDIM == 2
          dvol = f4%dr_level(1, level) * f4%dr_level(2, level)
#:elif NDIM == 3
          dvol = f4%dr_level(1, level) * f4%dr_level(2, level) * &
               f4%dr_level(3, level)
#:endif

          my_l1 = my_l1 + abs(f4%uu(${IJK}$, i_err, n)) * dvol
          my_l2 = my_l2 + f4%uu(${IJK}$, i_err, n)**2 * dvol
       end do; ${KJI_CLOSE_LOOP}$
    end do

    call MPI_Allreduce(MPI_IN_PLACE, my_l1, 1, MPI_DOUBLE_PRECISION, &
         MPI_SUM, f4%mpicomm, ierror)
    call MPI_Allreduce(MPI_IN_PLACE, my_l2, 1, MPI_DOUBLE_PRECISION, &
         MPI_SUM, f4%mpicomm, ierror)

    l1_err = my_l1
    l2_err = sqrt(my_l2)
  end subroutine compute_error_norms

  subroutine rho_solution(${XYZ}$, t, rho)
    ${ROUTINE_SEQ()}$
    real(dp), intent(in)  :: ${XYZ}$
    real(dp), intent(in)  :: t
    real(fp), intent(out) :: rho
    real(dp)              :: distance, q
    real(dp), parameter   :: radius = 0.1_dp
    real(dp), parameter   :: border = 0.05_dp
    real(dp)              :: x0(NDIM), dx(NDIM), cos_t, sin_t

    x0(1) = x
    x0(2) = y
#:if NDIM == 3
    x0(3) = z
#:endif

    select case (advection_velocity_type)
    case (1)
       ! Uniform velocity field with periodic boundaries
       x0 = x0 - advection_velocity * t

    case (2)
       ! Angular rotation around domain center with 1 rad/s
       dx = x0 - 0.5_dp
       cos_t = cos(t)
       sin_t = sin(t)

       x0(1) = 0.5_dp + cos_t * dx(1) - sin_t * dx(2)
       x0(2) = 0.5_dp + sin_t * dx(1) + cos_t * dx(2)
    case (3)
       ! No simple analytic solution is available, except for t = T, where T is
       ! the period of the deformation
    case default
    end select

    ! Keep the point inside the periodic domain
    x0 = modulo(x0, 1.0_dp)

    ! Avoid temporary array
    dx = x0 - advection_r0
    distance = sqrt(dot_product(dx, dx))

    if (advection_use_gaussian) then
       rho = real(exp(-(distance/radius)**2), fp)
    else
       if (distance < radius - border) then
          rho = 1.0_fp
       else if (distance < radius) then
          ! cubic smoothstep: 1 - 3 q^2 + 2 q^3, with q in [0,1]
          q = (distance - radius + border)/border
          rho = real(1.0_dp - (3.0_dp * q**2 - 2.0_dp * q**3), fp)
       else
          rho = 0.0_fp
       end if
    end if
  end subroutine rho_solution

  pure subroutine get_velocity(flux_dim, v, i0, n, ${IJK}$, f4)
    ${ROUTINE_SEQ()}$
    integer, intent(in)       :: flux_dim
    real(fp), intent(out)     :: v
    integer, intent(in)       :: i0      ! 0 if lower face, 1 if upper face
    integer, intent(in)       :: n       ! block index
    integer, intent(in)       :: ${IJK}$ ! i, j, k
    type(foap4_t), intent(in) :: f4
    real(dp)                  :: dr(ndim)
    real(fp)                  :: vel(ndim), rr(ndim), cost
    real(fp), parameter       :: pi = acos(-1.0_fp)
    real(fp), parameter       :: inv_T = 1/2.0_dp ! Inverse period

    select case (advection_velocity_type)
    case (1)
       v = advection_velocity(flux_dim)
    case (2, 3)
       dr = f4%dr_level(:, f4%block_level(n))
       rr(1) = real(f4%block_origin(1, n) + dr(1) * (i - 0.5_dp), fp)
       rr(2) = real(f4%block_origin(2, n) + dr(2) * (j - 0.5_dp), fp)
#:if NDIM == 3
       rr(3) = real(f4%block_origin(3, n) + dr(3) * (k - 0.5_dp), fp)
#:endif
       rr(flux_dim) = rr(flux_dim) + real((i0 - 0.5_dp) * dr(flux_dim), fp)

       if (advection_velocity_type == 2) then
          ! Clockwise solid-body rotation
          vel(1) = rr(2) - 0.5_fp
          vel(2) = -(rr(1) - 0.5_fp)
#:if NDIM == 3
          vel(3) = 0.0_fp
#:endif
       else
          cost = cos(inv_T * pi * real(f4%time, fp))
#:if NDIM == 2
          ! Deforming deformation, see eq. (9.5) in doi:10.1137/0733033
          vel(1) = -sin(pi * rr(1))**2 * sin(2 * pi * rr(2)) * cost
          vel(2) =  sin(2 * pi * rr(1)) * sin(pi * rr(2))**2 * cost
#:elif NDIM == 3
          ! Deforming deformation in 3D, see eq. (11.2) in doi:10.1137/0733033
          vel(1) = 2 * sin(pi * rr(1))**2 * sin(2 * pi * rr(2)) * &
               sin(2*pi*rr(3)) * cost
          vel(2) = -sin(2 * pi * rr(1)) * sin(pi * rr(2))**2 * &
               sin(2*pi*rr(3)) * cost
          vel(3) = -sin(2 * pi * rr(1)) * sin(2 * pi * rr(2)) * &
               sin(pi*rr(3))**2 * cost
#:endif
       end if
       v = vel(flux_dim)
    case default
       v = 0.0_fp
    end select
  end subroutine get_velocity

#:include 'physics_advection.fpp'

#:include 'flux_finite_volume.fpp'

  include 'flux_scheme_${FLUX_SCHEME}$_${NDIM}$d.f90'

  include 'limiter_${LIMITER}$.f90'

end program
