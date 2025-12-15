subroutine feuler_finite_volume(f4, dt, dt_lim, time, s_deriv, &
     n_prev, s_prev, w_prev, s_out, i_step, n_steps)
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
  integer                      :: n, ${IJK}$, m, level, iv, ierr
  integer                      :: i_vars_deriv(n_vars)
  logical                      :: ghost_dim(NDIM), valid_cell
  real(dp)                     :: inv_dr(NDIM), cmax(NDIM), max_cfl
  real(dp)                     :: flux(n_vars, 2, NDIM)
  real(dp)                     :: tmp(1+2*n_gc, n_vars)
  real(dp)                     :: dvar(n_vars), u(n_vars)
#:if NDIM == 2
  real(dp)                     :: uprim(f4%ilo(1):f4%ihi(1), f4%ilo(2):f4%ihi(2), n_vars)
#:elif NDIM == 3
  real(dp)                     :: uprim(f4%ilo(1):f4%ihi(1), f4%ilo(2):f4%ihi(2), &
       f4%ilo(3):f4%ihi(3), n_vars)
#:endif

  i_vars_deriv = i_vars + s_deriv
  call f4_update_ghostcells(f4, n_vars, i_vars_deriv)

  max_cfl = 0.0_dp

  !$acc parallel loop private(level, inv_dr, uprim) reduction(max:max_cfl) &
  !$acc &present(f4%uu, f4%ilo, f4%ihi, f4%bx, f4%bflux, f4%bflux_ix, &
  !$acc &f4%block_level, f4%dr_level)
  do n = 1, f4%n_blocks

     level = f4%block_level(n)
     inv_dr = 1/f4%dr_level(:, level)

     !$acc loop collapse(NDIM) private(u, ghost_dim)
     do @{KJI_LOOP_array_to_array(f4%ilo, f4%ihi)}@
        ghost_dim(1) = i < 1 .or. i > f4%bx(1)
        ghost_dim(2) = j < 1 .or. j > f4%bx(2)
#:if NDIM == 3
        ghost_dim(3) = k < 1 .or. k > f4%bx(3)
#:endif

#:if NDIM == 2
        valid_cell = .not. (ghost_dim(1) .and. ghost_dim(2))
#:elif NDIM == 3
        valid_cell = .not. ( &
             (ghost_dim(1) .and. ghost_dim(2)) .or. &
             (ghost_dim(1) .and. ghost_dim(3)) .or. &
             (ghost_dim(2) .and. ghost_dim(3)))
#:endif

        if (valid_cell) then
           ! Convert to primitive, but not in corners
           u = f4%uu(${IJK}$, i_vars0+1+s_deriv:i_vars0+n_vars+s_deriv, n)
           call to_primitive(u)
           uprim(${IJK}$, :) = u
        end if
     end do; ${KJI_CLOSE_LOOP}$

     !$acc loop collapse(NDIM) private(tmp, flux, dvar, cmax, iv, m, u) &
     !$acc &reduction(max:max_cfl)
     do @{KJI_LOOP_1_to_array(f4%bx)}@

        u = uprim(${IJK}$, :)
        dvar = 0.0_dp
        call source_term(u, dvar)
        dvar = dvar * dt

#:if NDIM == 2
        ! Compute fluxes
        tmp = uprim(i-n_gc:i+n_gc, j, :)
        call flux_cell_faces(1, tmp, flux(:, :, 1), cmax(1))

        tmp = uprim(i, j-n_gc:j+n_gc, :)
        call flux_cell_faces(2, tmp, flux(:, :, 2), cmax(2))

        ! Keep track of changes in variables
        dvar = dvar + dt * ( &
             (flux(:, 1, 1) - flux(:, 2, 1)) * inv_dr(1) + &
             (flux(:, 1, 2) - flux(:, 2, 2)) * inv_dr(2))

        ! Store boundary fluxes
        if (f4%bflux_ix(0, n) > 0 .and. i == 1) &
             f4%bflux(j, i_vars0+1:i_vars0+n_vars, f4%bflux_ix(0, n)) = dt * flux(:, 1, 1)
        if (f4%bflux_ix(1, n) > 0 .and. i == f4%bx(1)) &
             f4%bflux(j, i_vars0+1:i_vars0+n_vars, f4%bflux_ix(1, n)) = dt * flux(:, 2, 1)
        if (f4%bflux_ix(2, n) > 0 .and. j == 1) &
             f4%bflux(i, i_vars0+1:i_vars0+n_vars, f4%bflux_ix(2, n)) = dt * flux(:, 1, 2)
        if (f4%bflux_ix(3, n) > 0 .and. j == f4%bx(2)) &
             f4%bflux(i, i_vars0+1:i_vars0+n_vars, f4%bflux_ix(3, n)) = dt * flux(:, 2, 2)
#:elif NDIM == 3
        ! Compute fluxes
        tmp = uprim(i-n_gc:i+n_gc, j, k, :)
        call flux_cell_faces(1, tmp, flux(:, :, 1), cmax(1))

        tmp = uprim(i, j-n_gc:j+n_gc, k, :)
        call flux_cell_faces(2, tmp, flux(:, :, 2), cmax(2))

        tmp = uprim(i, j, k-n_gc:k+n_gc, :)
        call flux_cell_faces(3, tmp, flux(:, :, 3), cmax(3))

        ! Keep track of changes in variables
        dvar = dvar + dt * ( &
             (flux(:, 1, 1) - flux(:, 2, 1)) * inv_dr(1) + &
             (flux(:, 1, 2) - flux(:, 2, 2)) * inv_dr(2) + &
             (flux(:, 1, 3) - flux(:, 2, 3)) * inv_dr(3))

        ! Store boundary fluxes
        if (f4%bflux_ix(0, n) > 0 .and. i == 1) &
             f4%bflux(j, k, i_vars0+1:i_vars0+n_vars, f4%bflux_ix(0, n)) = dt * flux(:, 1, 1)
        if (f4%bflux_ix(1, n) > 0 .and. i == f4%bx(1)) &
             f4%bflux(j, k, i_vars0+1:i_vars0+n_vars, f4%bflux_ix(1, n)) = dt * flux(:, 2, 1)
        if (f4%bflux_ix(2, n) > 0 .and. j == 1) &
             f4%bflux(i, k, i_vars0+1:i_vars0+n_vars, f4%bflux_ix(2, n)) = dt * flux(:, 1, 2)
        if (f4%bflux_ix(3, n) > 0 .and. j == f4%bx(2)) &
             f4%bflux(i, k, i_vars0+1:i_vars0+n_vars, f4%bflux_ix(3, n)) = dt * flux(:, 2, 2)
        if (f4%bflux_ix(4, n) > 0 .and. k == 1) &
             f4%bflux(i, j, i_vars0+1:i_vars0+n_vars, f4%bflux_ix(4, n)) = dt * flux(:, 1, 3)
        if (f4%bflux_ix(5, n) > 0 .and. k == f4%bx(3)) &
             f4%bflux(i, j, i_vars0+1:i_vars0+n_vars, f4%bflux_ix(5, n)) = dt * flux(:, 2, 3)
#:endif

        ! Use cmax(1) to temporarily store sum(cmax * inv_dr)
#:if NDIM == 2
        cmax(1) = cmax(1)*inv_dr(1) + cmax(2)*inv_dr(2)
#:elif NDIM == 3
        cmax(1) = cmax(1)*inv_dr(1) + cmax(2)*inv_dr(2) + cmax(3)*inv_dr(3)
#:endif
        max_cfl = max(max_cfl, cmax(1))

        ! Set output state
        do iv = 1, n_vars
           do m = 1, n_prev
              ! Add weighted previous states
              dvar(iv) = dvar(iv) + &
                   f4%uu(${IJK}$, i_vars0+iv+s_prev(m), n) * w_prev(m)
           end do
           f4%uu(${IJK}$, i_vars0+iv+s_out, n) = dvar(iv)
        end do
     end do; ${KJI_CLOSE_LOOP}$
  end do

  call f4_fix_c2f_flux(f4, n_vars, i_vars, s_out)

  dt_lim = 1/max_cfl
  call MPI_Allreduce(MPI_IN_PLACE, dt_lim, 1, MPI_DOUBLE_PRECISION, &
       MPI_MIN, f4%mpicomm, ierr)

end subroutine feuler_finite_volume
