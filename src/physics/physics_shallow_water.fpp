!> Convert conservative variables in-place to primitive ones
pure subroutine to_primitive(u)
  ${ROUTINE_SEQ()}$
  real(fp), intent(inout) :: u(n_tvars)
  real(fp)                :: inv_h
  integer                 :: idim

  inv_h = 1/u(i_h)
  do idim = 1, ${NDIM}$
     u(i_mom0+idim) = u(i_mom0+idim) * inv_h
  end do
end subroutine to_primitive

!> Convert primitive variables in-place to conservative ones
pure subroutine to_conservative(u)
  ${ROUTINE_SEQ()}$
  real(fp), intent(inout) :: u(n_tvars)
  integer                 :: idim

  do idim = 1, ${NDIM}$
     u(i_mom0+idim) = u(i_h) * u(i_mom0+idim)
  end do
end subroutine to_conservative

!> Compute flux (in conservative variables) from primitive variables
subroutine get_flux(flux_dim, u, flux)
  ${ROUTINE_SEQ()}$
  integer, intent(in)   :: flux_dim
  real(fp), intent(in)  :: u(n_tvars)  !< Primitive variables (h, u, v)
  real(fp), intent(out) :: flux(n_tvars)
  integer               :: idim

  ! Mass flux: h * u_flux_dim
  flux(i_h) = u(i_h) * u(i_mom0+flux_dim)

  ! Momentum flux
  do idim = 1, ${NDIM}$
     flux(i_mom0+idim) = u(i_h) * u(i_mom0+idim) * u(i_mom0+flux_dim)
  end do

  ! Add pressure term (0.5 * g * h^2) to momentum in flux direction
  flux(i_mom0+flux_dim) = flux(i_mom0+flux_dim) + &
       0.5_fp * gravity_constant * u(i_h)**2
end subroutine get_flux

!> Estimate for minimum and maximum wavespeeds for HLL-type solvers
!> Wavespeeds are u ± sqrt(g*h)
pure subroutine get_min_max_wavespeed(flux_dim, u_LR, cmin, cmax)
  ${ROUTINE_SEQ()}$
  integer, intent(in)   :: flux_dim
  real(fp), intent(in)  :: u_LR(n_tvars, 2) !< Primitive variables (L and R states)
  real(fp), intent(out) :: cmin
  real(fp), intent(out) :: cmax
  real(fp)              :: h_sqrt(2), fac, umean, cmean, cwave(2)

  ! Gravity wave speed for each state
  cwave = sqrt(gravity_constant * u_LR(i_h, :))

  ! Roe-averaged quantities
  h_sqrt = sqrt(u_LR(i_h, :))
  fac = 1/(h_sqrt(1) + h_sqrt(2))

  umean = fac * (&
       h_sqrt(1) * u_LR(i_mom0+flux_dim, 1) + &
       h_sqrt(2) * u_LR(i_mom0+flux_dim, 2))

  cmean = sqrt(0.5_fp * gravity_constant * (u_LR(i_h, 1) + u_LR(i_h, 2)))

  ! Minimum and maximum wave speeds using Davis estimate, see
  ! S.F. Davis (1988), "Simplified Second-Order Godunov-Type Methods"
  cmin = min(u_LR(i_mom0+flux_dim, 1) - cwave(1), umean - cmean)
  cmax = max(u_LR(i_mom0+flux_dim, 2) + cwave(2), umean + cmean)
end subroutine get_min_max_wavespeed
