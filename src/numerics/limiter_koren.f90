pure subroutine reconstruct(u, i0, u_LR)
  !$acc routine seq
  real(dp), intent(in)  :: u(5, n_tvars)
  integer, intent(in)   :: i0
  real(dp), intent(out) :: u_LR(n_tvars, 2)
  real(dp)              :: u_diff(3)
  integer               :: n

  do n = 1, n_tvars
     u_diff = u(i0+2:i0+4, n) - u(i0+1:i0+3, n)
     u_LR(n, 1) = u(i0+2, n) + 0.5_dp * koren(u_diff(1), u_diff(2))
     u_LR(n, 2) = u(i0+3, n) - 0.5_dp * koren(u_diff(3), u_diff(2))
  end do
end subroutine reconstruct

!> Modified implementation of Koren limiter, to avoid division and the min/max
!> functions, which can be problematic / expensive. In most literature, you
!> have r = a / b (ratio of gradients). Then the limiter phi(r) is multiplied
!> with b. With this implementation, you get phi(r) * b
elemental pure real(dp) function koren(a, b) result(phi)
  real(dp), intent(in) :: a  !< Density gradient (numerator)
  real(dp), intent(in) :: b  !< Density gradient (denominator)
  real(dp), parameter  :: third = 1/3.0_dp
  real(dp)             :: aa, ab

  aa = a * a
  ab = a * b

  if (ab <= 0) then
     ! a and b have different sign or one of them is zero, so r is either 0,
     ! inf or negative (special case a == b == 0 is ignored)
     phi = 0
  else if (aa <= 0.25_dp * ab) then
     ! 0 < a/b <= 1/4, limiter has value 2*a/b
     phi = 2*a
  else if (aa <= 2.5_dp * ab) then
     ! 1/4 < a/b <= 2.5, limiter has value (1+2*a/b)/3
     phi = third * (b + 2*a)
  else
     ! (1+2*a/b)/6 >= 1, limiter has value 2
     phi = 2*b
  end if
end function koren
