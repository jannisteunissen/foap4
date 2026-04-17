!> Convert conservative variables in-place to primitive ones
pure subroutine to_primitive(u)
  !$acc routine seq
  real(fp), intent(inout) :: u(n_tvars)
  real(fp)                :: inv_rho, sum_v2
  integer                 :: idim

  ! Apply density floor
  u(i_rho) = max(u(i_rho), euler_rho_floor)

  inv_rho = 1/u(i_rho)
  sum_v2 = 0.0_fp
  do idim = 1, ${NDIM}$
     u(i_mom0+idim) = u(i_mom0+idim) * inv_rho
     sum_v2 = sum_v2 + u(i_mom0+idim)**2
  end do

  u(i_e) = (euler_gamma - 1.0_fp) * &
       (u(i_e) - 0.5_fp * u(i_rho) * sum_v2)

  ! Apply pressure floor
  u(i_e) = max(u(i_e), euler_p_floor)
end subroutine to_primitive

!> Convert primitive variables in-place to conservative ones
pure subroutine to_conservative(u)
  !$acc routine seq
  real(fp), intent(inout) :: u(n_tvars)
  integer                 :: idim
  real(fp)                :: sum_v2

  sum_v2 = 0.0_fp

  ! Compute momentum from density and velocity components
  do idim = 1, ${NDIM}$
     sum_v2 = sum_v2 + u(i_mom0+idim)**2
     u(i_mom0+idim) = u(i_rho) * u(i_mom0+idim)
  end do

  ! Compute energy from pressure and kinetic energy
  u(i_e) = u(i_e) * euler_inv_gamma_m1 + 0.5_fp * u(i_rho) * sum_v2
end subroutine to_conservative

!> Compute flux (in conservative variables) from primitive variables
subroutine get_flux(flux_dim, u, flux)
  !$acc routine seq
  integer, intent(in)   :: flux_dim
  real(fp), intent(in)  :: u(n_tvars)
  real(fp), intent(out) :: flux(n_tvars)
  integer               :: idim
  real(fp)              :: sum_v2

  sum_v2 = 0.0_fp

  ! Density flux
  flux(i_rho) = u(i_rho) * u(i_mom0+flux_dim)

  ! Momentum flux with pressure term
  do idim = 1, ${NDIM}$
     flux(i_mom0+idim) = u(i_rho) * u(i_mom0+idim) * u(i_mom0+flux_dim)
     sum_v2 = sum_v2 + u(i_mom0+idim)**2
  end do

  flux(i_mom0+flux_dim) = flux(i_mom0+flux_dim) + u(i_e)

  ! Energy flux
  flux(i_e) = u(i_mom0+flux_dim) * (u(i_e) * euler_inv_gamma_m1 + &
       0.5_fp * u(i_rho) * sum_v2 + u(i_e))
end subroutine get_flux

!> This implements formula (10.52) from "Riemann Solvers and Numerical Methods
!> for Fluid Dynamics" by Toro.
pure subroutine get_min_max_wavespeed(flux_dim, u_LR, cmin, cmax)
  !$acc routine seq
  integer, intent(in)   :: flux_dim
  real(fp), intent(in)  :: u_LR(n_tvars, 2) !< Primitive variables
  real(fp), intent(out) :: cmin
  real(fp), intent(out) :: cmax
  real(fp)              :: rho_sqrt(2), fac, eta2, umean, csound2(2), dmean

  rho_sqrt = sqrt(max(u_LR(i_rho, :), euler_rho_floor))
  fac = 1/(rho_sqrt(1) + rho_sqrt(2))

  umean = fac * (&
       rho_sqrt(1) * u_LR(i_mom0+flux_dim, 1) + &
       rho_sqrt(2) * u_LR(i_mom0+flux_dim, 2))

  ! Use robust sound speed calculation
  csound2(1) = get_csound2_from_prim(u_LR(:, 1))
  csound2(2) = get_csound2_from_prim(u_LR(:, 2))

  eta2 = 0.5_fp * fac**2 * rho_sqrt(1) * rho_sqrt(2)

  dmean = sqrt(fac * (rho_sqrt(1)*csound2(1) + rho_sqrt(2)*csound2(2)) + &
       eta2 * (u_LR(i_mom0+flux_dim, 2) - u_LR(i_mom0+flux_dim, 1))**2)

  cmin = umean - dmean
  cmax = umean + dmean
end subroutine get_min_max_wavespeed

!> Compute sound speed squared from primitive variable array
pure function get_csound2_from_prim(u) result(csound2)
  !$acc routine seq
  real(fp), intent(in) :: u(n_tvars)  !< Primitive variables
  real(fp)             :: csound2
  real(fp)             :: rho_safe, p_safe

  ! Apply floors for robustness
  rho_safe = max(u(i_rho), euler_rho_floor)
  p_safe = max(u(i_e), euler_p_floor)

  csound2 = euler_gamma * p_safe / rho_safe
end function get_csound2_from_prim
