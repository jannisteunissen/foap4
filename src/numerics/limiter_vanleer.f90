pure subroutine reconstruct(u, i0, u_LR)
  !$acc routine seq
  real(dp), intent(in)  :: u(5, n_tvars)
  integer, intent(in)   :: i0
  real(dp), intent(out) :: u_LR(n_tvars, 2)
  real(dp)              :: u_diff(3)
  integer               :: n

  do n = 1, n_tvars
     u_diff = u(i0+2:i0+4, n) - u(i0+1:i0+3, n)
     u_LR(n, 1) = u(i0+2, n) + 0.5_dp * vanleer(u_diff(1), u_diff(2))
     u_LR(n, 2) = u(i0+3, n) - 0.5_dp * vanleer(u_diff(3), u_diff(2))
  end do
end subroutine reconstruct

elemental pure real(dp) function vanleer(a, b) result(phi)
  !$acc routine seq
  real(dp), intent(in) :: a, b
  real(dp)             :: ab

  ab = a * b
  if (ab > 0) then
     phi = 2 * ab / (a + b)
  else
     phi = 0
  end if
end function vanleer

