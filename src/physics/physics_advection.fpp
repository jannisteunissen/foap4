pure subroutine get_flux(flux_dim, u, flux)
  !$acc routine seq
  integer, intent(in)   :: flux_dim
  real(dp), intent(in)  :: u(n_tvars)
  real(dp), intent(out) :: flux(n_tvars)
  flux(1) = advection_velocity(flux_dim) * u(1)
end subroutine get_flux

pure subroutine to_primitive(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(n_tvars)
end subroutine to_primitive

pure subroutine to_conservative(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(n_tvars)
end subroutine to_conservative

subroutine source_term(u_prim, source)
  real(dp), intent(in) :: u_prim(n_tvars)
  real(dp), intent(out) :: source(n_tvars)
  source = 0.0_dp
end subroutine source_term

pure subroutine get_min_max_wavespeed(flux_dim, u_LR, cmin, cmax)
  !$acc routine seq
  integer, intent(in)   :: flux_dim
  real(dp), intent(in)  :: u_LR(n_tvars, 2)
  real(dp), intent(out) :: cmin
  real(dp), intent(out) :: cmax
  cmin = min(advection_velocity(flux_dim), 0.0_dp)
  cmax = max(advection_velocity(flux_dim), 0.0_dp)
end subroutine get_min_max_wavespeed
