!> For testing purposes, simple upwind scheme with one ghost cell
pure subroutine reconstruct(u, i0, u_LR)
  !$acc routine seq
  real(fp), intent(in)  :: u(3, n_tvars)
  integer, intent(in)   :: i0
  real(fp), intent(out) :: u_LR(n_tvars, 2)
  integer               :: n

  do n = 1, n_tvars
     u_LR(n, 1) = u(i0+1, n)
     u_LR(n, 2) = u(i0+2, n)
  end do
end subroutine reconstruct
