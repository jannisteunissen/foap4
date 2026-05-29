#:include 'definitions_parallel.fpp'
subroutine flux_cell_faces(flux_dim, u, flux, max_wavespeed)
  ${ROUTINE_SEQ()}$
  integer, intent(in)   :: flux_dim
  real(dp), intent(in)  :: u(1+2*n_gc, n_tvars)
  real(dp), intent(out) :: flux(n_tvars, 2)
  real(dp), intent(out) :: max_wavespeed
  real(dp)              :: cmax(2)

  call flux_tvdlf_one_side(flux_dim, 0, u, flux(:, 1), cmax(1))
  call flux_tvdlf_one_side(flux_dim, 1, u, flux(:, 2), cmax(2))
  max_wavespeed = max(cmax(1), cmax(2))
end subroutine flux_cell_faces

subroutine flux_tvdlf_one_side(flux_dim, i0, u, flux, max_wavespeed)
  ${ROUTINE_SEQ()}$
  integer, intent(in)   :: flux_dim
  integer, intent(in)   :: i0
  real(dp), intent(in)  :: u(1+2*n_gc, n_tvars)
  real(dp), intent(out) :: flux(n_tvars)
  real(dp), intent(out) :: max_wavespeed
  real(dp)              :: u_LR(n_tvars, 2), S_L, S_R
  real(dp)              :: flux_LR(n_tvars, 2)

  call reconstruct(u, i0, u_LR)

  call get_flux(flux_dim, u_LR(:, 1), flux_LR(:, 1))
  call get_flux(flux_dim, u_LR(:, 2), flux_LR(:, 2))
  call get_min_max_wavespeed(flux_dim, u_LR, S_L, S_R)
  max_wavespeed = max(abs(S_L), abs(S_R))

  ! Convert to conservative
  call to_conservative(u_LR(:, 1))
  call to_conservative(u_LR(:, 2))

  flux = 0.5_dp * (flux_LR(:, 1) + flux_LR(:, 2) - &
       max_wavespeed * (u_LR(:, 2) - u_LR(:, 1)))
end subroutine flux_tvdlf_one_side
