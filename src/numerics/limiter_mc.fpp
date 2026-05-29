#:include 'definitions_parallel.fpp'
pure subroutine reconstruct(u, i0, u_LR)
  ${ROUTINE_SEQ()}$
  real(fp), intent(in)  :: u(5, n_tvars)
  integer, intent(in)   :: i0
  real(fp), intent(out) :: u_LR(n_tvars, 2)
  real(fp)              :: u_diff(3)
  integer               :: n

  do n = 1, n_tvars
     u_diff = u(i0+2:i0+4, n) - u(i0+1:i0+3, n)
     u_LR(n, 1) = u(i0+2, n) + 0.5_fp * mc_limiter(u_diff(1), u_diff(2))
     u_LR(n, 2) = u(i0+3, n) - 0.5_fp * mc_limiter(u_diff(3), u_diff(2))
  end do
end subroutine reconstruct

elemental pure real(fp) function mc_limiter(a, b) result(phi)
  real(fp), intent(in) :: a  !< Density gradient (numerator)
  real(fp), intent(in) :: b  !< Density gradient (denominator)
  phi = limiter_gminmod(a, b, 2.0_fp)
end function mc_limiter

!> Generalized minmod limiter. The parameter theta controls how dissipative
!> the limiter is, with 1 corresponding to the minmod limiter and 2 to the MC
!> limiter.
elemental function limiter_gminmod(a, b, theta) result(phi)
  ${ROUTINE_SEQ()}$
  real(fp), intent(in) :: a, b, theta
  real(fp)             :: phi

  if (a * b > 0) then
     phi = sign(min(theta*abs(a), theta*abs(b), 0.5_fp*abs(a+b)), a)
  else
     phi = 0.0_fp
  end if
end function limiter_gminmod
