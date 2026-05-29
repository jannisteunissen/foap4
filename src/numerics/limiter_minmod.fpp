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
     u_LR(n, 1) = u(i0+2, n) + 0.5_fp * minmod(u_diff(1), u_diff(2))
     u_LR(n, 2) = u(i0+3, n) - 0.5_fp * minmod(u_diff(3), u_diff(2))
  end do
end subroutine reconstruct

elemental pure real(fp) function minmod(a, b)
  ${ROUTINE_SEQ()}$
  real(fp), intent(in) :: a, b

  if (a * b <= 0) then
     minmod = 0.0_fp
  else if (abs(a) < abs(b)) then
     minmod = a
  else
     minmod = b
  end if
end function minmod
