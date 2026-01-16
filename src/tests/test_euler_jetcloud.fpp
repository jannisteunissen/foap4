!> Test case for jet-cloud problem, based on tests/hd/jet_cloud in MPI-AMRVAC
#:include '../core/definitions.fpp'
#:set LIMITER = 'minmod'
#:set FLUX_SCHEME = 'tvdlf'
program test_euler_jetcloud
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
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: n_gc = limiter_num_ghostcells

  type(foap4_t) :: f4
  type(CFG_t)   :: cfg

  integer           :: min_level          = 2
  integer           :: max_level          = 4
  integer           :: max_blocks         = 1000
  integer           :: bx(NDIM)           = 32
  integer           :: num_outputs        = 40
  logical           :: periodic(NDIM)     = .false.
  logical           :: do_refinement      = .false.
  real(dp)          :: end_time           = 5.0_dp
  real(dp)          :: c_refine           = 0.8_dp
  real(dp)          :: c_derefine         = 0.2_dp
  real(dp)          :: c_eps              = 0.01_dp
  real(dp)          :: c_abs              = 1e-10_dp
  integer           :: n_steps_refinement = 4
  real(dp)          :: cfl_number         = 0.5_dp
  character(len=10) :: test_case          = "sod"
  character(len=40) :: integrator_name    = "heuns_method"

  integer  :: blocks_per_dim(NDIM)
  real(dp) :: domain_length(NDIM)
  real(dp) :: block_length(NDIM)
  real(dp) :: beta, eta_jet, ca, mach, rc
  real(dp) :: r_cloud(NDIM), cloud_sigma
  real(dp) :: rho_jet

  blocks_per_dim(1) = 2
  blocks_per_dim(2:) = 1
  domain_length(1) = 32.0_dp
  domain_length(2:) = 16.0_dp
  block_length = domain_length / blocks_per_dim

  ! jet to cloud density ratio parameter
  beta = 0.04d0
  ! jet to ambient density
  eta_jet = 3.00d0
  ! ambient sound speed
  ca      = 1.00d0
  ! jet Mach number
  mach    = 10.0d0
  ! cloud to jet radii ratio
  rc      = 1.5d0
  ! Cloud location
  r_cloud = 0.5_dp * domain_length
  r_cloud(2) = r_cloud(2) + 0.075_dp * domain_length(2)
  ! Cloud width
  cloud_sigma = 0.75d0 * rc
  ! Jet density
  rho_jet = 1.0_dp

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
  call CFG_add_get(cfg, 'test_case', test_case, 'Which test case to run')
  call CFG_add_get(cfg, 'time_integrator', integrator_name, 'Time integrator')
  call CFG_add_get(cfg, 'cfl_number', cfl_number, 'CFL number')
  call CFG_check(cfg)

  call euler_jetcloud(f4, bx, min_level, max_blocks, &
       num_outputs, "output/test_euler_jetcloud_${NDIM}$d", &
       end_time, integrator_name)

  if (f4%mpirank == 0) call f4_print_wtime(f4)
  call f4_finalize(f4)

contains

  subroutine euler_jetcloud(f4, bx, min_level, max_blocks, num_outputs, &
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
    real(dp)                     :: time, t0, t1

    time = 0.0_dp
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

    ! Set boundary conditions
    do n = 1, n_tvars
       ! The values will be set by bc_callback
       call f4_set_bc_scalar(f4, i_tvars(n), f4_face_xlo, f4_bc_fixed_value, 0.0_dp)
    end do

    call f4_set_bc_scalar(f4, i_mom0+2, f4_face_ylo, f4_bc_dirichlet, 0.0_dp)
    call f4_set_bc_scalar(f4, i_mom0+2, f4_face_yhi, f4_bc_dirichlet, 0.0_dp)
#:if NDIM == 3
    call f4_set_bc_scalar(f4, i_mom0+3, f4_face_zlo, f4_bc_dirichlet, 0.0_dp)
    call f4_set_bc_scalar(f4, i_mom0+3, f4_face_zhi, f4_bc_dirichlet, 0.0_dp)
#:endif

    call set_initial_conditions(f4)

    if (do_refinement) then
       do n = 1, 10
          prev_mesh_revision = f4_get_mesh_revision(f4)
          call f4_update_ghostcells(f4, n_tvars, i_tvars)
          call amr_flags_diff2(f4, min_level, max_level, &
               i_rho, c_refine, c_derefine, c_eps, c_abs)
          call f4_adjust_refinement(f4, .true.)
          call set_initial_conditions(f4)

          if (f4_get_mesh_revision(f4) == prev_mesh_revision) exit
       end do
    end if

    if (dt_output <= end_time) call io_write_grid(f4, base_name, n_output, time)
    n_output = n_output + 1

    t0 = MPI_Wtime()

    do while (time < end_time)
       n_iterations = n_iterations + 1
       dt = cfl_number * dt_lim

       write_this_step = (time + dt > n_output * dt_output)
       if (write_this_step) dt = n_output * dt_output - time
       call rk_advance(f4, dt, dt_lim, time, integrator, feuler_finite_volume)

       if (write_this_step) then
          call io_write_grid(f4, base_name, n_output, time)
          n_output = n_output + 1
       end if

       if (do_refinement .and. &
            mod(n_iterations, n_steps_refinement) == 0) then
          call f4_get_global_highest_level(f4, prev_highest_level)
          call f4_update_ghostcells(f4, n_tvars, i_tvars)
          call amr_flags_diff2(f4, min_level, max_level, &
               i_rho, c_refine, c_derefine, c_eps, c_abs)
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
  end subroutine euler_jetcloud

  subroutine set_initial_conditions(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, ${IJK}$
    real(dp)                     :: rr(NDIM), d_inlet, d2_cloud
    real(dp)                     :: u(n_tvars)

    !$acc parallel loop present(f4%uu)
    do n = 1, f4%n_blocks
       !$acc loop collapse(NDIM) private(rr)
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          rr = f4_cell_coord(f4, n, ${IJK}$)

          d_inlet = get_inlet_distance(rr, domain_length)
          d2_cloud = sum((rr - r_cloud)**2)

          ! Ambient medium
          u(i_rho) = 1/eta_jet + (1/beta**2) * exp(-d2_cloud/cloud_sigma**2)
          u(i_mom0+1:i_mom0+NDIM) = 0.0_dp
          u(i_e) = ca**2 / (euler_gamma * eta_jet)

          if (d_inlet < 1.0_dp .and. rr(1) < 2.5_dp) then
             u(i_rho) = rho_jet
             u(i_mom0+1) = mach * ca
          end if

          call to_conservative(u)

          f4%uu(${IJK}$, i_tvars0+1:i_tvars0+n_tvars, n) = u
       end do; ${KJI_CLOSE_LOOP}$
    end do
  end subroutine set_initial_conditions

  pure real(dp) function get_inlet_distance(r, domain_length) result(d_inlet)
    real(dp), intent(in) :: r(NDIM)
    real(dp), intent(in) :: domain_length(NDIM)
    real(dp)             :: domain_center(NDIM)

    domain_center = 0.5_dp * domain_length
#:if NDIM == 2
    d_inlet = abs(r(2) - domain_center(2))
#:elif NDIM == 3
    d_inlet = norm2(r(2:3) - domain_center(2:3))
#:endif
  end function get_inlet_distance

  subroutine bc_callback(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: face, i_block, ix, n, i
    real(dp)                     :: rr(NDIM), d_inlet, u(n_tvars)
#:if NDIM == 3
    integer                      :: j
#:endif

    face = f4_face_xlo

    !$acc parallel loop private(i_block, ix) independent
    do n = f4%gc_phys_iface(face), f4%gc_phys_iface(face+1)-1
       i_block = f4%gc_phys(n) + 1
       f4%bc_data_ix(face, i_block) = abs(f4%bc_data_ix(face, i_block))
       ix = f4%bc_data_ix(face, i_block)

#:if NDIM == 2
       !$acc loop private(rr, d_inlet, u)
       do i = 1, f4%bx(1)
          rr = f4_block_face_coord(f4, i_block, face, i)
          d_inlet = get_inlet_distance(rr, domain_length)

          u(i_rho) = 1/eta_jet
          u(i_mom0+1:i_mom0+ndim) = 0.0_dp
          u(i_e) = ca**2 / (euler_gamma * eta_jet)

          if (d_inlet < 1) then
             u(i_rho) = rho_jet
             u(i_mom0+1) = mach * ca
          end if

          call to_conservative(u)
          f4%bc_data(i, i_tvars0+1:i_tvars0+n_tvars, ix) = u
       end do
#:elif NDIM == 3
       !$acc loop collapse(2) private(rr, d_inlet, u)
       do j = 1, f4%bx(1)
          do i = 1, f4%bx(1)
             rr = f4_block_face_coord(f4, i_block, face, i, j)
             d_inlet = get_inlet_distance(rr, domain_length)
             f4%bc_data(i, j, i_phi, ix) = phi_init(rr(1), rr(2), rr(3))

             if (d_inlet < 1) then
                u(i_rho) = rho_jet
                u(i_mom0+1) = mach * ca
                u(i_mom0+2:i_mom0+ndim) = 0.0_dp
                u(i_e) = ca**2 / (euler_gamma * eta_jet)
             else
                u(i_rho) = 1/eta_jet
                u(i_mom0+1:i_mom0+ndim) = 0.0_dp
                u(i_e) = ca**2 / (euler_gamma * eta_jet)
             end if

             call to_conservative(u)
             f4%bc_data(i_tvars0+1:i_tvars0+n_tvars, ix) = u
          end do
       end do
#:endif
    end do

  end subroutine bc_callback

  subroutine source_term(u_prim, source)
    real(dp), intent(in)    :: u_prim(n_tvars)
    real(dp), intent(inout) :: source(n_tvars)
  end subroutine source_term

  #:include 'physics_euler.fpp'

  #:include 'flux_finite_volume.fpp'

  include 'flux_${FLUX_SCHEME}$.f90'

  include 'limiter_${LIMITER}$.f90'

end program test_euler_jetcloud
