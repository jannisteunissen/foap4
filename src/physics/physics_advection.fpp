pure subroutine get_flux(flux_dim, u, flux)
  !$acc routine seq
  integer, intent(in)   :: flux_dim
  real(fp), intent(in)  :: u(n_tvars)
  real(fp), intent(out) :: flux(n_tvars)
  flux(1) = advection_velocity(flux_dim) * u(1)
end subroutine get_flux

pure subroutine to_primitive(u)
  !$acc routine seq
  real(fp), intent(inout) :: u(n_tvars)
end subroutine to_primitive

pure subroutine to_conservative(u)
  !$acc routine seq
  real(fp), intent(inout) :: u(n_tvars)
end subroutine to_conservative

subroutine source_term(u_prim, source)
  !$acc routine seq
  real(fp), intent(in) :: u_prim(n_tvars)
  real(fp), intent(out) :: source(n_tvars)
  source = 0.0_fp
end subroutine source_term

pure subroutine get_min_max_wavespeed(flux_dim, u_LR, cmin, cmax)
  !$acc routine seq
  integer, intent(in)   :: flux_dim
  real(fp), intent(in)  :: u_LR(n_tvars, 2)
  real(fp), intent(out) :: cmin
  real(fp), intent(out) :: cmax
  cmin = min(advection_velocity(flux_dim), 0.0_fp)
  cmax = max(advection_velocity(flux_dim), 0.0_fp)
end subroutine get_min_max_wavespeed
