pure subroutine reconstruct(u, i0, u_LR)
  !$acc routine seq
  real(dp), intent(in)  :: u(5, n_tvars)
  integer, intent(in)   :: i0
  real(dp), intent(out) :: u_LR(n_tvars, 2)
  real(dp)              :: u_diff(3)
  integer               :: n

  do n = 1, n_tvars
     u_diff = u(i0+2:i0+4, n) - u(i0+1:i0+3, n)
     u_LR(n, 1) = u(i0+2, n) + 0.5_dp * mc_limiter(u_diff(1), u_diff(2))
     u_LR(n, 2) = u(i0+3, n) - 0.5_dp * mc_limiter(u_diff(3), u_diff(2))
  end do
end subroutine reconstruct

elemental pure real(dp) function mc_limiter(a, b) result(phi)
  real(dp), intent(in) :: a  !< Density gradient (numerator)
  real(dp), intent(in) :: b  !< Density gradient (denominator)
  phi = limiter_gminmod(a, b, 2.0_dp)
end function mc_limiter

!> Generalized minmod limiter. The parameter theta controls how dissipative
!> the limiter is, with 1 corresponding to the minmod limiter and 2 to the MC
!> limiter.
elemental function limiter_gminmod(a, b, theta) result(phi)
  !$acc routine seq
  real(dp), intent(in) :: a, b, theta
  real(dp)             :: phi

  if (a * b > 0) then
     phi = sign(min(theta*abs(a), theta*abs(b), 0.5_dp*abs(a+b)), a)
  else
     phi = 0.0_dp
  end if
end function limiter_gminmod
