pure subroutine reconstruct(u, i0, u_LR)
  !$acc routine seq
  real(dp), intent(in)  :: u(5, n_vars)
  integer, intent(in)   :: i0
  real(dp), intent(out) :: u_LR(n_vars, 2)
  real(dp)              :: u_diff(3)
  integer               :: n

  do n = 1, n_vars
     u_diff = u(i0+2:i0+4, n) - u(i0+1:i0+3, n)
     u_LR(n, 1) = u(i0+2, n) + 0.5_dp * minmod(u_diff(1), u_diff(2))
     u_LR(n, 2) = u(i0+3, n) - 0.5_dp * minmod(u_diff(2), u_diff(3))
  end do
end subroutine reconstruct

elemental pure real(dp) function minmod(a, b)
  !$acc routine seq
  real(dp), intent(in) :: a, b

  if (a * b <= 0) then
     minmod = 0.0_dp
  else if (abs(a) < abs(b)) then
     minmod = a
  else
     minmod = b
  end if
end function minmod
