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
  integer                      :: n, ${IJK}$, m, level, iv, ix(NDIM), ierr
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

  call f4_update_ghostcells(f4, n_vars, i_vars+s_deriv)

  max_cfl = 0.0_dp

  !$acc parallel loop private(level, inv_dr, max_cfl, uprim) &
  !$acc &reduction(max:max_cfl)
  do n = 1, f4%n_blocks

     level = f4%block_level(n)
     inv_dr = 1/f4%dr_level(:, level)

     !$acc loop collapse(NDIM) private(ix, u)
     do @{KJI_LOOP_array_to_array(f4%ilo, f4%ihi)}@
        ix = [${IJK}$]
        if (count(ix < 1 .or. ix > f4%bx) <= 1) then
           ! Convert to primitive, but not in corners
           u = f4%uu(${IJK}$, i_vars+s_deriv, n)
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
        if (i == 1) f4%bflux(j, 0, i_vars, n) = dt * flux(:, 1, 1)
        if (i == bx(1)) f4%bflux(j, 1, i_vars, n) = dt * flux(:, 2, 1)
        if (j == 1) f4%bflux(i, 2, i_vars, n) = dt * flux(:, 1, 2)
        if (j == bx(2)) f4%bflux(i, 3, i_vars, n) = dt * flux(:, 2, 2)
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
        if (i == 1)     f4%bflux(j, k, 0, i_vars, n) = dt * flux(:, 1, 1)
        if (i == bx(1)) f4%bflux(j, k, 1, i_vars, n) = dt * flux(:, 2, 1)
        if (j == 1)     f4%bflux(i, k, 2, i_vars, n) = dt * flux(:, 1, 2)
        if (j == bx(2)) f4%bflux(i, k, 3, i_vars, n) = dt * flux(:, 2, 2)
        if (k == 1)     f4%bflux(i, j, 4, i_vars, n) = dt * flux(:, 1, 3)
        if (k == bx(3)) f4%bflux(i, j, 5, i_vars, n) = dt * flux(:, 2, 3)
#:endif

        max_cfl = max(max_cfl, sum(cmax * inv_dr))

        ! Set output state
        do iv = 1, n_vars
           do m = 1, n_prev
              ! Add weighted previous states
              dvar(iv) = dvar(iv) + &
                   f4%uu(${IJK}$, i_vars(iv)+s_prev(m), n) * w_prev(m)
           end do
           f4%uu(${IJK}$, i_vars(iv)+s_out, n) = dvar(iv)
        end do
     end do; ${KJI_CLOSE_LOOP}$
  end do

  call f4_fix_c2f_flux(f4, n_vars, i_vars, s_out)

  dt_lim = 1/max_cfl
  call MPI_Allreduce(MPI_IN_PLACE, dt_lim, 1, MPI_DOUBLE_PRECISION, &
       MPI_MIN, f4%mpicomm, ierr)

end subroutine feuler_finite_volume
