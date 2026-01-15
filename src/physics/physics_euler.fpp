pure subroutine to_primitive(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(n_tvars)
  real(dp)                :: inv_rho, sum_v2
  integer                 :: idim

  inv_rho = 1/u(i_rho)
  sum_v2 = 0.0_dp
  do idim = 1, ${NDIM}$
     u(i_mom0+idim) = u(i_mom0+idim) * inv_rho
     sum_v2 = sum_v2 + u(i_mom0+idim)**2
  end do

  u(i_e) = (euler_gamma - 1.0_dp) * &
       (u(i_e) - 0.5_dp * u(i_rho) * sum_v2)
end subroutine to_primitive

pure subroutine to_conservative(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(n_tvars)
  integer                 :: idim
  real(dp)                :: sum_v2

  sum_v2 = 0.0_dp

  ! Compute momentum from density and velocity components
  do idim = 1, ${NDIM}$
     sum_v2 = sum_v2 + u(i_mom0+idim)**2
     u(i_mom0+idim) = u(i_rho) * u(i_mom0+idim)
  end do

  ! Compute energy from pressure and kinetic energy
  u(i_e) = u(i_e) * inv_gamma_m1 + 0.5_dp * u(i_rho) * sum_v2
end subroutine to_conservative

subroutine get_flux(flux_dim, u, flux)
  !$acc routine seq
  integer, intent(in)   :: flux_dim
  real(dp), intent(in)  :: u(n_tvars)
  real(dp), intent(out) :: flux(n_tvars)
  integer               :: idim
  real(dp)              :: sum_v2

  sum_v2 = 0.0_dp

  ! Density flux
  flux(i_rho) = u(i_rho) * u(i_mom0+flux_dim)

  ! Momentum flux with pressure term
  do idim = 1, ${NDIM}$
     flux(i_mom0+idim) = u(i_rho) * u(i_mom0+idim) * u(i_mom0+flux_dim)
     sum_v2 = sum_v2 + u(i_mom0+idim)**2
  end do

  flux(i_mom0+flux_dim) = flux(i_mom0+flux_dim) + u(i_e)

  ! Energy flux
  flux(i_e) = u(i_mom0+flux_dim) * (u(i_e) * inv_gamma_m1 + &
       0.5_dp * u(i_rho) * sum_v2 + u(i_e))
end subroutine get_flux

subroutine get_max_wavespeed(flux_dim, u_LR, cmax)
  !$acc routine seq
  integer, intent(in)   :: flux_dim
  real(dp), intent(in)  :: u_LR(n_tvars, 2)
  real(dp), intent(out) :: cmax
  real(dp)              :: cmin

  call get_min_max_wavespeed(flux_dim, u_LR, cmin, cmax)
  cmax = max(abs(cmin), cmax)
end subroutine get_max_wavespeed

!> This implements formula (10.52) from "Riemann Solvers and Numerical Methods
!> for Fluid Dynamics" by Toro.
pure subroutine get_min_max_wavespeed(flux_dim, u_LR, cmin, cmax)
  !$acc routine seq
  integer, intent(in)   :: flux_dim
  real(dp), intent(in)  :: u_LR(n_tvars, 2)
  real(dp), intent(out) :: cmin
  real(dp), intent(out) :: cmax
  real(dp)              :: rho_sqrt(2), fac, eta2, umean, csound2(2), dmean

  rho_sqrt = sqrt(u_LR(i_rho, :))
  fac = 1/(rho_sqrt(1) + rho_sqrt(2))

  umean = fac * (&
       rho_sqrt(1) * u_LR(i_mom0+flux_dim, 1) + &
       rho_sqrt(1) * u_LR(i_mom0+flux_dim, 1))

  ! Square of sound speed
  csound2 = euler_gamma * u_LR(i_e, :) / u_LR(i_rho, :)

  eta2 = 0.5_dp * fac**2 * rho_sqrt(1) * rho_sqrt(2)

  dmean = sqrt(fac * (rho_sqrt(1)*csound2(1) + rho_sqrt(2)*csound2(2)) + &
       eta2 * (u_LR(i_mom0+flux_dim, 2) - u_LR(i_mom0+flux_dim, 1))**2)

  cmin = umean - dmean
  cmax = umean + dmean
end subroutine get_min_max_wavespeed
