!> Test case double mach reflection test case (2D).
!>
!> Test case from: Woodward, P., Colella, P. "The numerical simulation of
!> two-dimensional fluid flow with strong shocks"
!> J. Comput. Phys. 54, 115–173 (1984)
#:assert NDIM == 2
#:include '../core/definitions.fpp'
#:set LIMITER = 'vanleer'
#:set FLUX_SCHEME = 'hll'
program test_dmr
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
  real(dp), parameter :: pi = acos(-1.0_dp)

  type(foap4_t) :: f4
  type(CFG_t)   :: cfg

  integer           :: min_level          = 2
  integer           :: max_level          = 4
  integer           :: max_blocks         = 1000
  integer           :: bx(NDIM)           = 24
  integer           :: num_outputs        = 40
  logical           :: periodic(NDIM)     = .false.
  logical           :: do_refinement      = .false.
  real(dp)          :: end_time           = 0.2_dp
  real(dp)          :: c_refine           = 0.8_dp
  real(dp)          :: c_derefine         = 0.2_dp
  real(dp)          :: c_eps              = 0.01_dp
  real(dp)          :: c_abs              = 1e-10_dp
  integer           :: i_var_refinement   = 1
  integer           :: n_steps_refinement = 4
  real(dp)          :: cfl_number         = 0.5_dp
  real(dp)          :: load_imbalance_threshold = 1.1_dp
  character(len=40) :: integrator_name    = "heuns_method"

  integer  :: blocks_per_dim(NDIM)
  real(dp) :: domain_length(NDIM)
  real(dp) :: block_length(NDIM)
  real(dp) :: dmr_xr, dmr_shock_angle, dmr_vx
  real(dp) :: dmr_mach, dmr_tan_angle
  real(fp) :: dmr_uL(n_tvars), dmr_uR(n_tvars)

  call euler_initialize(1.4_dp, 0.0_dp, 1e-12_dp, 1e-12_dp)

  blocks_per_dim(:) = [4, 1]
  domain_length(:) = [4.0_dp, 1.0_dp]
  block_length = domain_length / blocks_per_dim

  ! Initial bottom coordinate of shock
  dmr_xr = 1/6.0_dp

  ! mach number
  dmr_mach = 10

  ! Shock angle to the x-axis
  dmr_shock_angle = 60.0_dp * (pi/180)

  ! tan(shock_angle)
  dmr_tan_angle = tan(dmr_shock_angle)

  ! Shock velocity in x-direction
  dmr_vx = dmr_mach / cos(0.5_dp * pi - dmr_shock_angle)

  ! Post-shock
  dmr_uL = real([8.0_dp, 4.125_dp * sqrt(3.0_dp), -4.125_dp, 116.5_dp], fp)
  call to_conservative(dmr_uL)
  ${ENTER_DATA_COPYIN('dmr_uL')}$

  ! Pre-shock
  dmr_uR = real([1.4_dp, 0.0_dp, 0.0_dp, 1.0_dp], fp)
  call to_conservative(dmr_uR)
  ${ENTER_DATA_COPYIN('dmr_uR')}$

  call f4_initialize(f4, "error")

  call CFG_update_from_arguments(cfg)
  call CFG_add_get(cfg, 'num_outputs', num_outputs, 'Write this many output files')
  call CFG_add_get(cfg, 'min_level', min_level, 'Minimum refinement level')
  call CFG_add_get(cfg, 'max_level', max_level, 'Maximum refinement level')
  call CFG_add_get(cfg, 'do_refinement', do_refinement, 'Perform refinement')
  call CFG_add_get(cfg, 'c_refine', c_refine, 'Coefficient for refinement')
  call CFG_add_get(cfg, 'c_derefine', c_derefine, 'Coefficient for derefinement')
  call CFG_add_get(cfg, 'c_eps', c_eps, 'Filter coefficient for AMR')
  call CFG_add_get(cfg, 'c_abs', c_abs, 'Density threshold for AMR')
  call CFG_add_get(cfg, 'n_steps_refinement', n_steps_refinement, &
       'Perform refinement every N steps')
  call CFG_add_get(cfg, 'max_blocks', max_blocks, 'Max. number of blocks')
  call CFG_add_get(cfg, 'periodic', periodic, 'Whether the domain is periodic')
  call CFG_add_get(cfg, 'bx', bx, 'Size of grid blocks')
  call CFG_add_get(cfg, 'blocks_per_dim', blocks_per_dim, &
       'Number of blocks (per dimension) on coarse grid')
  call CFG_add_get(cfg, 'end_time', end_time, 'End time')
  call CFG_add_get(cfg, 'time_integrator', integrator_name, 'Time integrator')
  call CFG_add_get(cfg, 'cfl_number', cfl_number, 'CFL number')
  call CFG_check(cfg)

  call euler_dmr(f4, bx, min_level, max_blocks, &
       num_outputs, "output/test_euler_dmr_${NDIM}$d", &
       end_time, integrator_name)

  call f4_print_wtime(f4)
  call f4_finalize(f4)

contains

  subroutine euler_dmr(f4, bx, min_level, max_blocks, num_outputs, &
       base_name, end_time, integrator_name)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: bx(NDIM)
    integer, intent(in)          :: min_level
    integer, intent(in)          :: max_blocks
    integer, intent(in)          :: num_outputs
    character(len=*), intent(in) :: base_name
    real(dp), intent(in)         :: end_time
    character(len=40), intent(in) :: integrator_name
    integer                      :: n_output
    integer                      :: n, n_iterations, ierr, prev_mesh_revision
    integer(int64)               :: sum_local_blocks, sum_global_blocks
    logical                      :: write_this_step
    integer                      :: integrator, n_time_states
    integer                      :: highest_level, prev_highest_level
    real(dp)                     :: dt, dt_lim, dt_output
    real(dp)                     :: t0, t1

    f4%time = 0.0_dp
    dt_output = end_time / max(real(num_outputs, dp), 1e-100_dp)
    n_output = 0
    n_iterations = 0
    sum_local_blocks = 0
    dt_lim = 0.0_dp

    integrator = rk_get_integrator_by_name(trim(integrator_name))
    n_time_states = rk_advance_num_copies(integrator)

    f4%bc_callback => bc_callback

    call f4_construct_brick(f4, blocks_per_dim, block_length, bx, n_gc, &
         n_vars_all, var_names, var_temporal, n_time_states, periodic, &
         min_level, max_blocks, f4_bc_neumann, 0.0_dp)

    call set_initial_conditions(f4)

    if (do_refinement) then
       do n = 1, 10
          prev_mesh_revision = f4_get_mesh_revision(f4)
          call f4_update_ghostcells(f4, n_tvars, i_tvars, 0)
          call amr_flags_diff2(f4, min_level, max_level, &
               i_var_refinement, c_refine, c_derefine, c_eps, c_abs)
          call f4_adjust_refinement(f4, load_imbalance_threshold)
          call set_initial_conditions(f4)

          if (f4_get_mesh_revision(f4) == prev_mesh_revision) exit
       end do
    end if

    if (dt_output <= end_time) call io_write_grid(f4, base_name, n_output)
    n_output = n_output + 1

    t0 = MPI_Wtime()

    do while (f4%time < end_time)
       n_iterations = n_iterations + 1
       dt = cfl_number * dt_lim

       call bc_callback(f4)

       write_this_step = (f4%time + dt > n_output * dt_output)
       if (write_this_step) dt = n_output * dt_output - f4%time
       call rk_advance(f4, dt, dt_lim, integrator, feuler_finite_volume)

       if (write_this_step) then
          call io_write_grid(f4, base_name, n_output)
          n_output = n_output + 1
       end if

       if (do_refinement .and. &
            mod(n_iterations, n_steps_refinement) == 0) then
          call f4_get_global_highest_level(f4, prev_highest_level)
          call f4_update_ghostcells(f4, n_tvars, i_tvars, 0)
          call amr_flags_diff2(f4, min_level, max_level, &
               i_var_refinement, c_refine, c_derefine, c_eps, c_abs)
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
  end subroutine euler_dmr

  subroutine set_initial_conditions(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, ${IJK}$
    real(dp)                     :: rr(NDIM), x_boundary

    ${PARALLEL_LOOP()}$ ${DEFAULT_PRESENT()}$
    do n = 1, f4%n_blocks
       ${LOOP('collapse(NDIM) private(rr, x_boundary)')}$
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          rr = f4_cell_coord(f4, n, ${IJK}$)

          x_boundary = dmr_xr + rr(2) / dmr_tan_angle

          if (rr(1) < x_boundary) then
             f4%uu(${IJK}$, i_tvars0+1:i_tvars0+n_tvars, n) = dmr_uL
          else
             f4%uu(${IJK}$, i_tvars0+1:i_tvars0+n_tvars, n) = dmr_uR
          end if

       end do; ${KJI_CLOSE_LOOP}$
    end do
  end subroutine set_initial_conditions

  subroutine bc_callback(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: face, i_block, ix, n, i
    real(dp)                     :: rr(NDIM), x_boundary
#:if NDIM == 3
    integer                      :: j
#:endif

    x_boundary = dmr_xr + domain_length(2) / dmr_tan_angle + &
         dmr_vx * f4%time

    face = f4_face_ylo

    ${PARALLEL_LOOP('private(i_block, ix)')}$
    do n = f4%gc_phys_iface(face), f4%gc_phys_iface(face+1)-1
       i_block = f4%gc_phys(n) + 1
       f4%bc_data_ix(face, i_block) = abs(f4%bc_data_ix(face, i_block))
       ix = f4%bc_data_ix(face, i_block)

       ${LOOP('private(rr)')}$
       do i = 1, f4%bx(1)
          rr = f4_block_face_coord(f4, i_block, face, i)

          if (rr(1) < dmr_xr) then
             f4%bc_data(i, i_tvars0+1:i_tvars0+n_tvars, ix) = dmr_uL
             f4%bc_data_type(i, i_tvars0+1:i_tvars0+n_tvars, ix) = f4_bc_fixed_value
          else
             ! Neumann zero for all except y-velocity, which is Dirichlet zero
             f4%bc_data(i, i_tvars0+1:i_tvars0+n_tvars, ix) = 0.0_dp
             f4%bc_data_type(i, i_tvars0+1:i_tvars0+n_tvars, ix) = f4_bc_neumann
             f4%bc_data_type(i, i_mom0+2, ix) = f4_bc_dirichlet
          end if
       end do
    end do

    face = f4_face_yhi

    ${PARALLEL_LOOP('private(i_block, ix)')}$
    do n = f4%gc_phys_iface(face), f4%gc_phys_iface(face+1)-1
       i_block = f4%gc_phys(n) + 1
       f4%bc_data_ix(face, i_block) = abs(f4%bc_data_ix(face, i_block))
       ix = f4%bc_data_ix(face, i_block)

       ${LOOP('private(rr)')}$
       do i = 1, f4%bx(1)
          rr = f4_block_face_coord(f4, i_block, face, i)

          if (rr(1) < x_boundary) then
             f4%bc_data(i, i_tvars0+1:i_tvars0+n_tvars, ix) = dmr_uL
          else
             f4%bc_data(i, i_tvars0+1:i_tvars0+n_tvars, ix) = dmr_uR
          end if

          f4%bc_data_type(i, i_tvars0+1:i_tvars0+n_tvars, ix) = f4_bc_fixed_value
       end do
    end do

    face = f4_face_xlo

    ${PARALLEL_LOOP('private(i_block, ix)')}$
    do n = f4%gc_phys_iface(face), f4%gc_phys_iface(face+1)-1
       i_block = f4%gc_phys(n) + 1
       f4%bc_data_ix(face, i_block) = abs(f4%bc_data_ix(face, i_block))
       ix = f4%bc_data_ix(face, i_block)

       ${LOOP('private(rr)')}$
       do i = 1, f4%bx(1)
          rr = f4_block_face_coord(f4, i_block, face, i)

          f4%bc_data(i, i_tvars0+1:i_tvars0+n_tvars, ix) = dmr_uL
          f4%bc_data_type(i, i_tvars0+1:i_tvars0+n_tvars, ix) = f4_bc_fixed_value
       end do
    end do

  end subroutine bc_callback

  subroutine source_term(u_prim, source)
    real(fp), intent(in)    :: u_prim(n_tvars)
    real(fp), intent(inout) :: source(n_tvars)
  end subroutine source_term

  #:include 'physics_euler.fpp'

  #:include 'flux_finite_volume.fpp'

  include 'flux_${FLUX_SCHEME}$.f90'

  include 'limiter_${LIMITER}$.f90'

end program test_dmr
