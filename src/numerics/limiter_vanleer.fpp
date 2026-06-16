#:include 'definitions_parallel.fpp'
pure subroutine reconstruct(u, i0, u_LR)
  ${ROUTINE_SEQ()}$
  real(fp), intent(in)  :: u(5, n_tvars)
  integer, intent(in)   :: i0
  real(fp), intent(out) :: u_LR(n_tvars, 2)
  integer               :: n

  do n = 1, n_tvars
     u_LR(n, 1) = u(i0+2, n) + 0.5_fp * vanleer( &
          u(i0+3, n) - u(i0+2, n), &
          u(i0+2, n) - u(i0+1, n) )

     u_LR(n, 2) = u(i0+3, n) - 0.5_fp * vanleer( &
          u(i0+4, n) - u(i0+3, n), &
          u(i0+3, n) - u(i0+2, n) )
  end do
end subroutine reconstruct

elemental pure real(fp) function vanleer(a, b) result(phi)
  ${ROUTINE_SEQ()}$
  real(fp), intent(in) :: a, b

  if (a * b > 0) then
     phi = 2 * a * b / (a + b)
  else
     phi = 0
  end if
end function vanleer


