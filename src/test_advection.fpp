#:include 'definitions.fpp'
program test_adv
  use iso_fortran_env, only: int64
  use mpi_f08
  use m_foap4_${NDIM}$d
  use m_config

  implicit none
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: NDIM = ${NDIM}$

  integer, parameter :: n_vars            = 1
  integer, parameter :: i_rho             = 1
  character(len=20)  :: var_names(n_vars) = ['rho']
  logical, parameter :: temporal(n_vars)  = .true.
  real(dp)           :: velocity(NDIM)    = 1.0_dp
  real(dp)           :: end_time          = 1.0_dp
  real(dp)           :: c_refine          = 0.8_dp
  real(dp)           :: c_derefine        = 0.2_dp
  real(dp)           :: c_eps             = 0.01_dp

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
  call CFG_add_get(cfg, 'velocity', velocity, 'Velocity')
  call CFG_add_get(cfg, 'end_time', end_time, 'End time')
  call CFG_add_get(cfg, 'time_integrator', integrator_name, 'Time integrator')
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
    integer, parameter           :: n_gc = 2
    logical, parameter           :: periodic(NDIM) = .true.
    real(dp), parameter          :: cfl_number = 0.5_dp
    integer                      :: n, prev_mesh_revision, n_output
    integer                      :: highest_level, n_iterations, ierr
    integer(int64)               :: sum_local_blocks, sum_global_blocks
    logical                      :: write_this_step
    integer                      :: integrator, n_time_states
    real(dp)                     :: dt, dt_lim, dt_output, min_dr(NDIM)
    real(dp)                     :: time, t0, t1
    real(dp)                     :: rho_initial_sum, rho_sum

    time = 0.0_dp
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
          call f4_update_ghostcells(f4, 1, [i_rho])
          call f4_set_refinement_flags_diff2(f4, min_refinement_level, &
               max_refinement_level, 1, [i_rho], c_refine, c_derefine, c_eps)
          call f4_adjust_refinement(f4, .true.)
          call set_init_cond(f4)

          if (f4_get_mesh_revision(f4) == prev_mesh_revision) exit
       end do
    end if

    call f4_compute_sum(f4, i_rho, rho_initial_sum)

    if (dt_output < end_time) call f4_write_grid(f4, base_name, n_output, time)
    n_output = n_output + 1

    call f4_get_global_highest_level(f4, highest_level)
    min_dr = f4%dr_level(:, highest_level)

    t0 = MPI_Wtime()

    do while (time < end_time)
       n_iterations = n_iterations + 1

       dt = cfl_number / (sum(abs(velocity)/min_dr) + epsilon(1.0_dp))
       write_this_step = (time + dt > n_output * dt_output)
       if (write_this_step) dt = n_output * dt_output - time

       call f4_advance(f4, dt, dt_lim, time, integrator, forward_euler)

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
          call f4_update_ghostcells(f4, 1, [i_rho])
          call f4_set_refinement_flags_diff2(f4, min_refinement_level, &
               max_refinement_level, 1, [i_rho], c_refine, c_derefine, c_eps)
          call f4_adjust_refinement(f4, .true.)

          call f4_get_global_highest_level(f4, highest_level)
          min_dr = f4%dr_level(:, highest_level)
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

  subroutine forward_euler(f4, dt, dt_lim, time, s_deriv, n_prev, s_prev, w_prev, &
       s_out, i_step, n_steps)
    type(foap4_t), intent(inout) :: f4
    real(dp), intent(in)         :: dt
    real(dp), intent(inout)      :: dt_lim         !< Time step limit
    real(dp), intent(in)         :: time           !< Current time
    integer, intent(in)          :: s_deriv        !< State to compute derivatives from
    integer, intent(in)          :: n_prev         !< Number of previous states
    integer, intent(in)          :: s_prev(n_prev) !< Previous states
    real(dp), intent(in)         :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)          :: s_out          !< Output state
    integer, intent(in)          :: i_step         !< Step of the integrator
    integer, intent(in)          :: n_steps        !< Total number of steps
    integer                      :: n, ${IJK}$, m, level
    real(dp)                     :: inv_dr(NDIM)
    real(dp)                     :: fx(2), fy(2)
    real(dp)                     :: tmp(5)
    real(dp)                     :: dvar(@{DINDEX(bx)}@)
#:if NDIM == 3
    real(dp)                     :: fz(2)
#:endif

    call f4_update_ghostcells(f4, 1, [i_rho+s_deriv])

    !$acc parallel loop private(level, inv_dr, dvar)
    do n = 1, f4%n_blocks

       level = f4%block_level(n)
       inv_dr = 1/f4%dr_level(:, level)

#:if NDIM == 2
       !$acc loop collapse(NDIM) private(fx, fy, tmp)
#:elif NDIM == 3
       !$acc loop collapse(NDIM) private(fx, fy, fz, tmp)
#:endif
       do @{KJI_LOOP_1_to_array(f4%bx)}@

#:if NDIM == 2
          ! Compute fluxes
          tmp = f4%uu(i-2:i+2, j, i_rho+s_deriv, n)
          call muscl_flux(velocity(1), tmp, fx)

          tmp = f4%uu(i, j-2:j+2, i_rho+s_deriv, n)
          call muscl_flux(velocity(2), tmp, fy)

          ! Keep track of changes in variables
          dvar(i, j) = dt * ( &
               (fx(1) - fx(2)) * inv_dr(1) + &
               (fy(1) - fy(2)) * inv_dr(2))

          ! Store boundary fluxes
          if (i == 1) f4%bflux(j, 0, i_rho, n) = dt * fx(1)
          if (i == bx(1)) f4%bflux(j, 1, i_rho, n) = dt * fx(2)
          if (j == 1) f4%bflux(i, 2, i_rho, n) = dt * fy(1)
          if (j == bx(2)) f4%bflux(i, 3, i_rho, n) = dt * fy(2)
#:elif NDIM == 3
          ! Compute fluxes
          tmp = f4%uu(i-2:i+2, j, k, i_rho+s_deriv, n)
          call muscl_flux(velocity(1), tmp, fx)

          tmp = f4%uu(i, j-2:j+2, k, i_rho+s_deriv, n)
          call muscl_flux(velocity(2), tmp, fy)

          tmp = f4%uu(i, j, k-2:k+2, i_rho+s_deriv, n)
          call muscl_flux(velocity(3), tmp, fz)

          ! Keep track of changes in variables
          dvar(${IJK}$) = dt * ( &
               (fx(1) - fx(2)) * inv_dr(1) + &
               (fy(1) - fy(2)) * inv_dr(2) + &
               (fz(1) - fz(2)) * inv_dr(3))

          ! Store boundary fluxes
          if (i == 1)     f4%bflux(j, k, 0, i_rho, n) = dt * fx(1)
          if (i == bx(1)) f4%bflux(j, k, 1, i_rho, n) = dt * fx(2)
          if (j == 1)     f4%bflux(i, k, 2, i_rho, n) = dt * fy(1)
          if (j == bx(2)) f4%bflux(i, k, 3, i_rho, n) = dt * fy(2)
          if (k == 1)     f4%bflux(i, j, 4, i_rho, n) = dt * fz(1)
          if (k == bx(3)) f4%bflux(i, j, 5, i_rho, n) = dt * fz(2)
#:endif
       end do; ${KJI_CLOSE_LOOP}$

       ! Set output state after computations are done, since s_out can be
       ! equal to s_deriv and s_prev
       !$acc loop collapse(NDIM) private(m)
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          do m = 1, n_prev
             ! Add weighted previous states
             dvar(${IJK}$) = dvar(${IJK}$) + &
                  f4%uu(${IJK}$, i_rho+s_prev(m), n) * w_prev(m)
          end do
          f4%uu(${IJK}$, i_rho+s_out, n) = dvar(${IJK}$)
       end do; ${KJI_CLOSE_LOOP}$
    end do

    call f4_fix_c2f_flux(f4, 1, [i_rho], s_out)

  end subroutine forward_euler

  subroutine muscl_flux(v, u, flux)
    !$acc routine seq
    real(dp), intent(in)  :: v
    real(dp), intent(in)  :: u(5)
    real(dp), intent(out) :: flux(2)
    real(dp)              :: u_diff(4), uL, uR

    u_diff = u(2:5) - u(1:4)

    uL = u(2) + 0.5_dp * vanleer(u_diff(1), u_diff(2))
    uR = u(3) - 0.5_dp * vanleer(u_diff(2), u_diff(3))
    flux(1) = 0.5 * (v*uL + v*uR - abs(v) * (uR-uL))

    uL = u(3) + 0.5_dp * vanleer(u_diff(2), u_diff(3))
    uR = u(4) - 0.5_dp * vanleer(u_diff(3), u_diff(4))
    flux(2) = 0.5 * (v*uL + v*uR - abs(v) * (uR-uL))

  end subroutine muscl_flux

  elemental pure real(dp) function minmod(a, b)
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
