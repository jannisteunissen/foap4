pure subroutine get_flux(flux_dim, u, flux, i0, n, ${IJK}$, f4)
  ${ROUTINE_SEQ()}$
  integer, intent(in)       :: flux_dim
  real(fp), intent(in)      :: u(n_tvars)
  real(fp), intent(out)     :: flux(n_tvars)
  integer, intent(in)       :: i0      ! 0 if lower face, 1 if upper face
  integer, intent(in)       :: n       ! block index
  integer, intent(in)       :: ${IJK}$ ! i, j, k
  type(foap4_t), intent(in) :: f4
  real(fp)                  :: v

  call get_velocity(flux_dim, v, i0, n, ${IJK}$, f4)
  flux(1) = v * u(1)
end subroutine get_flux

pure subroutine to_primitive(u)
  ${ROUTINE_SEQ()}$
  real(fp), intent(inout) :: u(n_tvars)
end subroutine to_primitive

pure subroutine to_conservative(u)
  ${ROUTINE_SEQ()}$
  real(fp), intent(inout) :: u(n_tvars)
end subroutine to_conservative

pure subroutine source_term(u_prim, source)
  ${ROUTINE_SEQ()}$
  real(fp), intent(in) :: u_prim(n_tvars)
  real(fp), intent(out) :: source(n_tvars)
  source = 0.0_fp
end subroutine source_term

pure subroutine get_min_max_wavespeed(flux_dim, u_LR, cmin, cmax, i0, n, ${IJK}$, f4)
  ${ROUTINE_SEQ()}$
  integer, intent(in)       :: flux_dim
  real(fp), intent(in)      :: u_LR(n_tvars, 2)
  real(fp), intent(out)     :: cmin
  real(fp), intent(out)     :: cmax
  integer, intent(in)       :: i0      ! 0 if lower face, 1 if upper face
  integer, intent(in)       :: n       ! block index
  integer, intent(in)       :: ${IJK}$ ! i, j, k
  type(foap4_t), intent(in) :: f4
  real(fp)                  :: v

  call get_velocity(flux_dim, v, i0, n, ${IJK}$, f4)
  cmin = min(v, 0.0_fp)
  cmax = max(v, 0.0_fp)
end subroutine get_min_max_wavespeed
