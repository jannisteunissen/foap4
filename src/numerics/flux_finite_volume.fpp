#:include 'definitions_ndim.fpp'
#:include 'definitions_parallel.fpp'
subroutine feuler_finite_volume(f4, dt_in, dt_lim, s_deriv, &
     n_prev, s_prev, w_prev, s_out, i_step, n_steps)
  type(foap4_t), intent(inout) :: f4
  real(dp), intent(in)         :: dt_in
  real(dp), intent(inout)      :: dt_lim         !< Time step limit
  integer, intent(in)          :: s_deriv        !< State to compute derivatives from
  integer, intent(in)          :: n_prev         !< Number of previous states
  integer, intent(in)          :: s_prev(n_prev) !< Previous states
  real(dp), intent(in)         :: w_prev(n_prev) !< Weights of previous states
  integer, intent(in)          :: s_out          !< Output state
  integer, intent(in)          :: i_step         !< Step of the integrator
  integer, intent(in)          :: n_steps        !< Total number of steps
  integer                      :: n, ${IJK}$, m, level, iv, ierr
  logical                      :: ghost_dim(NDIM), valid_cell
  real(dp)                     :: t0, t1
  real(fp)                     :: dt, inv_dr(NDIM), cmax, cfl_sum, max_cfl
  real(fp)                     :: flux(n_tvars, 2), tmp(1+2*n_gc, n_tvars)
  real(fp)                     :: dvar(n_tvars), u(n_tvars)

  call f4_update_ghostcells(f4, n_tvars, i_tvars, s_deriv)

  t0 = MPI_Wtime()
  max_cfl = 0.0_dp
  dt = real(dt_in, fp)

  ${PARALLEL_LOOP_FLAT('collapse(NDIM+1) private(u, ghost_dim, valid_cell)')}$ ${DEFAULT_PRESENT()}$
  do n = 1, f4%n_blocks
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
           u = f4%uu(${IJK}$, i_tvars0+1:i_tvars0+n_tvars, n+s_deriv)
           call to_primitive(u)
           f4%uu_prim(${IJK}$, i_tvars0+1:i_tvars0+n_tvars, n) = u
        end if
     end do; ${KJI_CLOSE_LOOP}$
  end do

  ${PARALLEL_LOOP_FLAT('collapse(NDIM+1) private(level, inv_dr, tmp, flux, dvar, cmax, cfl_sum, iv, m, u) reduction(max:max_cfl)')}$ ${COPYIN('w_prev, s_prev')}$ ${DEFAULT_PRESENT()}$
  do n = 1, f4%n_blocks
     do @{KJI_LOOP_1_to_array(f4%bx)}@
        level = f4%block_level(n)
        inv_dr = real(1/f4%dr_level(:, level), fp)

        u = f4%uu_prim(${IJK}$, i_tvars0+1:i_tvars0+n_tvars, n)
        dvar = 0.0_fp
        call source_term(u, dvar)
        dvar = dvar * real(dt, fp)

#:if NDIM == 2
        ! Compute x-flux
        tmp = f4%uu_prim(i-n_gc:i+n_gc, j, i_tvars0+1:i_tvars0+n_tvars, n)
        call flux_cell_faces(1, tmp, flux, cmax)
        dvar = dvar + dt * inv_dr(1) * (flux(:, 1) - flux(:, 2))
        cfl_sum = cmax * inv_dr(1)

        ! Store refinement boundary fluxes
        if (f4%bflux_ix(0, n) > 0 .and. i == 1) &
             f4%bflux(j, i_tvars0+1:i_tvars0+n_tvars, f4%bflux_ix(0, n)) = dt * flux(:, 1)
        if (f4%bflux_ix(1, n) > 0 .and. i == f4%bx(1)) &
             f4%bflux(j, i_tvars0+1:i_tvars0+n_tvars, f4%bflux_ix(1, n)) = dt * flux(:, 2)

        ! Compute y-flux
        tmp = f4%uu_prim(i, j-n_gc:j+n_gc, i_tvars0+1:i_tvars0+n_tvars, n)
        call flux_cell_faces(2, tmp, flux, cmax)
        dvar = dvar + dt * inv_dr(2) * (flux(:, 1) - flux(:, 2))
        cfl_sum = cfl_sum + cmax * inv_dr(2)

        ! Store refinement boundary fluxes
        if (f4%bflux_ix(2, n) > 0 .and. j == 1) &
             f4%bflux(i, i_tvars0+1:i_tvars0+n_tvars, f4%bflux_ix(2, n)) = dt * flux(:, 1)
        if (f4%bflux_ix(3, n) > 0 .and. j == f4%bx(2)) &
             f4%bflux(i, i_tvars0+1:i_tvars0+n_tvars, f4%bflux_ix(3, n)) = dt * flux(:, 2)
#:elif NDIM == 3
        ! Compute fluxes
        tmp = f4%uu_prim(i-n_gc:i+n_gc, j, k, i_tvars0+1:i_tvars0+n_tvars, n)
        call flux_cell_faces(1, tmp, flux, cmax)
        dvar = dvar + dt * inv_dr(1) * (flux(:, 1) - flux(:, 2))
        cfl_sum = cmax * inv_dr(1)

        ! Store refinement boundary fluxes
        if (f4%bflux_ix(0, n) > 0 .and. i == 1) &
             f4%bflux(j, k, i_tvars0+1:i_tvars0+n_tvars, f4%bflux_ix(0, n)) = dt * flux(:, 1)
        if (f4%bflux_ix(1, n) > 0 .and. i == f4%bx(1)) &
             f4%bflux(j, k, i_tvars0+1:i_tvars0+n_tvars, f4%bflux_ix(1, n)) = dt * flux(:, 2)

        tmp = f4%uu_prim(i, j-n_gc:j+n_gc, k, i_tvars0+1:i_tvars0+n_tvars, n)
        call flux_cell_faces(2, tmp, flux, cmax)
        dvar = dvar + dt * inv_dr(2) * (flux(:, 1) - flux(:, 2))
        cfl_sum = cfl_sum + cmax * inv_dr(2)

        ! Store refinement boundary fluxes
        if (f4%bflux_ix(2, n) > 0 .and. j == 1) &
             f4%bflux(i, k, i_tvars0+1:i_tvars0+n_tvars, f4%bflux_ix(2, n)) = dt * flux(:, 1)
        if (f4%bflux_ix(3, n) > 0 .and. j == f4%bx(2)) &
             f4%bflux(i, k, i_tvars0+1:i_tvars0+n_tvars, f4%bflux_ix(3, n)) = dt * flux(:, 2)

        tmp = f4%uu_prim(i, j, k-n_gc:k+n_gc, i_tvars0+1:i_tvars0+n_tvars, n)
        call flux_cell_faces(3, tmp, flux, cmax)
        dvar = dvar + dt * inv_dr(3) * (flux(:, 1) - flux(:, 2))
        cfl_sum = cfl_sum + cmax * inv_dr(3)

        ! Store refinement boundary fluxes
        if (f4%bflux_ix(4, n) > 0 .and. k == 1) &
             f4%bflux(i, j, i_tvars0+1:i_tvars0+n_tvars, f4%bflux_ix(4, n)) = dt * flux(:, 1)
        if (f4%bflux_ix(5, n) > 0 .and. k == f4%bx(3)) &
             f4%bflux(i, j, i_tvars0+1:i_tvars0+n_tvars, f4%bflux_ix(5, n)) = dt * flux(:, 2)
#:endif

        max_cfl = max(max_cfl, cfl_sum)

        ! Set output state
        do iv = 1, n_tvars
           do m = 1, n_prev
              ! Add weighted previous states
              dvar(iv) = dvar(iv) + real(w_prev(m), fp) * &
                   f4%uu(${IJK}$, i_tvars0+iv, n+s_prev(m))
           end do
           f4%uu(${IJK}$, i_tvars0+iv, n+s_out) = dvar(iv)
        end do
     end do; ${KJI_CLOSE_LOOP}$
  end do

  dt_lim = 1/max_cfl
  call MPI_Allreduce(MPI_IN_PLACE, dt_lim, 1, MPI_DOUBLE_PRECISION, &
       MPI_MIN, f4%mpicomm, ierr)

  t1 = MPI_Wtime()
  f4%wtime_finite_volume = f4%wtime_finite_volume + t1 - t0

  call f4_fix_c2f_flux(f4, n_tvars, i_tvars, s_out)

end subroutine feuler_finite_volume
