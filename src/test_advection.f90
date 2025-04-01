program test_adv

  use m_foap4

  implicit none
  integer, parameter :: dp = kind(0.0d0)

  integer, parameter  :: n_vars            = 2
  integer, parameter  :: i_rho             = 1
  character(len=20)   :: var_names(n_vars) = ['rho', 'tmp']
  real(dp), parameter :: velocity(2)       = [1.0_dp, 1.0_dp]

  type(foap4_t), target :: f4

  call f4_initialize(f4, "error")

  call test_advection(f4, 2, "output/test_adv")

  call f4_finalize(f4)

contains

  subroutine test_advection(f4, min_level, base_name)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: min_level
    character(len=*), intent(in) :: base_name
    integer, parameter           :: n_blocks_per_dim(2) = [1, 1]
    real(dp), parameter          :: block_length(2)     = [1.0_dp, 1.0_dp]
    integer, parameter           :: bx(2)               = [16, 16]
    integer, parameter           :: n_gc                = 2
    logical, parameter           :: periodic(2)         = [.true., .true.]
    integer, parameter           :: max_blocks          = 1000
    real(dp), parameter          :: cfl_number          = 0.5_dp
    integer                      :: n, prev_mesh_revision, n_output
    integer                      :: highest_level
    logical                      :: write_output, do_refinement
    real(dp)                     :: dt, dt_output, min_dr(2)
    real(dp)                     :: time, end_time

    time = 0.0_dp
    end_time = 1.0_dp
    dt_output = end_time / 40
    n_output = 0

    call f4_set_grid(f4, n_blocks_per_dim, block_length, bx, n_gc, &
         n_vars, var_names, periodic, min_level, max_blocks)

    call set_init_cond(f4)

    do_refinement = .true.

    if (do_refinement) then
       do n = 1, 10
          prev_mesh_revision = f4_get_mesh_revision(f4)
          call f4_update_ghostcells(f4, 1, [i_rho])
          call set_refinement_flag(f4)
          call f4_adjust_refinement(f4, .true.)
          call set_init_cond(f4)

          if (f4_get_mesh_revision(f4) == prev_mesh_revision) exit
       end do
    end if

    call f4_write_grid(f4, base_name, n_output, time)
    n_output = n_output + 1
    n = 0

    call f4_get_global_highest_level(f4, highest_level)
    min_dr = f4%dr_level(:, highest_level)

    do while (time < end_time)
       n = n + 1
       dt = cfl_number / (sum(abs(velocity)/min_dr) + epsilon(1.0_dp))
       write_output = (time + dt > n_output * dt_output)
       if (write_output) dt = n_output * dt_output - time

       call advance_heuns_method(f4, dt)
       time = time + dt

       if (write_output) then
          call f4_write_grid(f4, base_name, n_output)
          n_output = n_output + 1
       end if

       if (do_refinement) then
          call f4_update_ghostcells(f4, 1, [i_rho])
          call set_refinement_flag(f4)
          call f4_adjust_refinement(f4, .true.)

          call f4_get_global_highest_level(f4, highest_level)
          min_dr = f4%dr_level(:, highest_level)
       end if
    end do

    if (f4%mpirank == 0) then
       print *, "n_iterations: ", n
    end if

    call f4_destroy(f4)
  end subroutine test_advection

  subroutine advance_heuns_method(f4, dt)
    type(foap4_t), intent(inout) :: f4
    real(dp), intent(in)         :: dt

    call forward_euler(f4, f4%bx, f4%ilo, f4%ihi, f4%n_vars, f4%n_blocks, &
         dt, f4%uu, 0, 1, [0], [1.0_dp], 1)
    call forward_euler(f4, f4%bx, f4%ilo, f4%ihi, f4%n_vars, f4%n_blocks, &
         0.5_dp*dt, f4%uu, 1, 2, [0, 1], [0.5_dp, 0.5_dp], 0)
  end subroutine advance_heuns_method

  subroutine set_init_cond(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, i, j
    real(dp)                     :: rr(2)

    do n = 1, f4%n_blocks
       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)
             rr = f4_cell_coord(f4, n, i, j)
             f4%uu(i, j, i_rho, n) = rho_init(rr(1), rr(2))
          end do
       end do
    end do
  end subroutine set_init_cond

  pure real(dp) function rho_init(x, y)
    real(dp), intent(in) :: x, y

    if (sqrt((x-0.5_dp)**2 + (y-0.5_dp)**2) < 0.1_dp) then
       rho_init = 1.0_dp
    else
       rho_init = 0.0_dp
    end if
  end function rho_init

  subroutine set_refinement_flag(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, i, j, ref_flag, level
    real(dp)                     :: dr(2), diff
    integer, parameter           :: max_level = 4
    integer, parameter           :: min_level = 2

    do n = 1, f4%n_blocks
       level = f4%block_level(n)
       dr = f4%dr_level(:, level)

       if (f4%block_level(n) > min_level) then
          ref_flag = -1
       else
          ref_flag = 0
       end if

       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)

             diff = abs(dr(1) * (f4%uu(i+1, j, i_rho, n) + &
                  f4%uu(i-1, j, i_rho, n) - 2 * f4%uu(i, j, i_rho, n)) + &
                  dr(2) * (f4%uu(i, j+1, i_rho, n) + &
                  f4%uu(i, j-1, i_rho, n) - 2 * f4%uu(i, j, i_rho, n)))

             if (diff > 1.0e-10_dp) then
                if (f4%block_level(n) < max_level) then
                   ref_flag = 1
                   exit
                else
                   ref_flag = 0
                end if
             end if
          end do
       end do

       f4%refinement_flags(n) = ref_flag
    end do

  end subroutine set_refinement_flag

  subroutine forward_euler(f4, bx, lo, hi, n_vars, n_blocks, dt, uu, &
       s_deriv, n_prev, s_prev, w_prev, s_out)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)               :: n_blocks, bx(2), lo(2), hi(2), n_vars
    real(dp), intent(in)              :: dt
    real(dp), intent(inout)           :: uu(lo(1):hi(1), lo(2):hi(2), n_vars, n_blocks)
    integer, intent(in)               :: s_deriv        !< State to compute derivatives from
    integer, intent(in)               :: n_prev         !< Number of previous states
    integer, intent(in)               :: s_prev(n_prev) !< Previous states
    real(dp), intent(in)              :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)               :: s_out          !< Output state
    integer                           :: n, i, j, m, level
    real(dp)                          :: inv_dr(2)
    real(dp)                          :: fx(2), fy(2)
    real(dp)                          :: tmp(5)
    real(dp)                          :: dvar(bx(1), bx(2))

    call f4_update_ghostcells(f4, 1, [i_rho+s_deriv])

    !$acc parallel loop private(fx, fy, tmp, u, uprim, dvar, m, iv)
    do n = 1, n_blocks

       level = f4%block_level(n)
       inv_dr = 1/f4%dr_level(:, level)

       !$acc loop collapse(2)
       do j = 1, bx(2)
          do i = 1, bx(1)
             ! Compute x and y fluxes
             tmp = uu(i-2:i+2, j, i_rho+s_deriv, n)
             call muscl_flux(velocity(1), tmp, fx)

             tmp = uu(i, j-2:j+2, i_rho+s_deriv, n)
             call muscl_flux(velocity(2), tmp, fy)

             ! Keep track of changes in variables
             dvar(i, j) = dt * ((fx(1) - fx(2)) * inv_dr(1) + &
                  (fy(1) - fy(2)) * inv_dr(2))
          end do
       end do

       ! Set output state after computations are done, since s_out can be
       ! equal to s_deriv and s_prev
       !$acc loop collapse(3)
       do j = 1, bx(2)
          do i = 1, bx(1)
             do m = 1, n_prev
                ! Add weighted previous states
                dvar(i, j) = dvar(i, j) + &
                     uu(i, j, i_rho+s_prev(m), n) * w_prev(m)
             end do
             uu(i, j, i_rho+s_out, n) = dvar(i, j)
          end do
       end do
    end do
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
