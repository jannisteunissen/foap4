pure subroutine to_primitive(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(n_vars)
  real(dp)                :: inv_rho
  integer                 :: idim

  inv_rho = 1/u(i_rho)
  do idim = 1, ${NDIM}$
     u(i_mom(idim)) = u(i_mom(idim)) * inv_rho
  end do

  u(i_e) = (euler_gamma - 1.0_dp) * &
       (u(i_e) - 0.5_dp * u(i_rho) * sum(u(i_mom)**2))
end subroutine to_primitive

pure subroutine to_conservative(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(n_vars)
  integer                 :: idim

  ! Compute energy from pressure and kinetic energy
  u(i_e) = u(i_e) * inv_gamma_m1 + 0.5_dp * u(i_rho) * sum(u(i_mom)**2)

  ! Compute momentum from density and velocity components
  do idim = 1, ${NDIM}$
     u(i_mom(idim)) = u(i_rho) * u(i_mom(idim))
  end do
end subroutine to_conservative

subroutine get_flux(flux_dim, u, flux)
  !$acc routine seq
  integer, intent(in)   :: flux_dim
  real(dp), intent(in)  :: u(n_vars)
  real(dp), intent(out) :: flux(n_vars)
  integer               :: idim

  ! Density flux
  flux(i_rho) = u(i_rho) * u(i_mom(flux_dim))

  ! Momentum flux with pressure term
  do idim = 1, ${NDIM}$
     flux(i_mom(idim)) = u(i_rho) * u(i_mom(idim)) * u(i_mom(flux_dim))
  end do

  flux(i_mom(flux_dim)) = flux(i_mom(flux_dim)) + u(i_e)

  ! Energy flux
  flux(i_e) = u(i_mom(flux_dim)) * (u(i_e) * inv_gamma_m1 + &
       0.5_dp * u(i_rho) * sum(u(i_mom)**2) + u(i_e))
end subroutine get_flux

pure subroutine get_max_wavespeed(flux_dim, u_LR, cmax)
  !$acc routine seq
  integer, intent(in)   :: flux_dim
  real(dp), intent(in)  :: u_LR(n_vars, 2)
  real(dp), intent(out) :: cmax
  real(dp)              :: cmin

  call get_min_max_wavespeed(flux_dim, u_LR, cmin, cmax)
  cmax = max(abs(cmin), cmax)
end subroutine get_max_wavespeed

pure subroutine get_min_max_wavespeed(flux_dim, u_LR, cmin, cmax)
  !$acc routine seq
  integer, intent(in)   :: flux_dim
  real(dp), intent(in)  :: u_LR(n_vars, 2)
  real(dp), intent(out) :: cmin
  real(dp), intent(out) :: cmax
  real(dp)              :: rho_sqrt(2), fac, umean, csound2(2), dmean

  rho_sqrt = sqrt(u_LR(i_rho, :))
  fac      = 1/sum(rho_sqrt)

  associate (i_v => i_mom(flux_dim))
    umean   = fac * sum(rho_sqrt * u_LR(i_v, :))
    csound2 = euler_gamma * u_LR(i_e, :) / u_LR(i_rho, :)
    dmean   = sqrt(fac * sum(rho_sqrt * csound2) + 0.5_dp * fac**2 * &
         product(rho_sqrt) * (u_LR(i_v, 2) - u_LR(i_v, 1)))
  end associate

  cmin = umean - dmean
  cmax = umean + dmean
end subroutine get_min_max_wavespeed

