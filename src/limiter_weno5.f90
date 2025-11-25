pure subroutine reconstruct(u, i0, u_LR)
  !$acc routine seq
  real(dp), intent(in)  :: u(7, n_vars)
  integer, intent(in)   :: i0
  real(dp), intent(out) :: u_LR(n_vars, 2)
  integer               :: n

  do n = 1, n_vars
     u_LR(n, 1) = weno5_reconstruct(u(i0+1:i0+5, n), -1)
     u_LR(n, 2) = weno5_reconstruct(u(i0+2:i0+6, n), 1)
  end do
end subroutine reconstruct

!> Perform weno5 reconstruction to a cell face
!> Reference: Jiang & Shu "Efficient Implementation of Weighted ENO Schemes", JCP 1996
pure real(dp) function weno5_reconstruct(u, di)
  !$acc routine seq

  !> Cell-centered values. The face value is reconstructed at one side of
  !> the central cell at u(3), depending on the argument di.
  real(dp), intent(in) :: u(5)
  !> Direction, which can be 1 or +1. If -1, the 'left' side of the face
  !> between u(3) and u(4) is reconstructed. If +1, the 'right' side of the
  !> face between u(2) and u(3) is reconstructed.
  integer, intent(in)  :: di
  real(dp), parameter  :: c1(3)      = [2, -7, 11]/6.0_dp
  real(dp), parameter  :: c2(3)      = [-1, 5, 2]/6.0_dp
  real(dp), parameter  :: c3(3)      = [2, 5, -1]/6.0_dp
  real(dp), parameter  :: w_IS(2)    = [13/12.0_dp, 0.25_dp]
  real(dp), parameter  :: weights(3) = [0.1_dp, 0.6_dp, 0.3_dp]
  real(dp), parameter  :: weno_eps   = 1e-18_dp
  real(dp)             :: ENO(3), alpha(3), IS(3), inv_alpha_sum

  ! ENO schemes
  ENO(1) = u(3+2*di)*c1(1) + u(3+di)*c1(2) + u(3)*c1(3)
  ENO(2) = u(3+di)*c2(1)   + u(3)*c2(2)    + u(3-di)*c2(3)
  ENO(3) = u(3)*c3(1)      + u(3-di)*c3(2) + u(3-2*di)*c3(3)

  ! Smoothness indicators
  IS(1) = w_IS(1) * (u(3+2*di) - 2*u(3+di) + u(3))**2 + &
       w_IS(2) * (u(3+2*di) - 4*u(3+di) + 3*u(3))**2
  IS(2) = w_IS(1) * (u(3+di) - 2*u(3) + u(3-di))**2 + &
       w_IS(2) * (u(3+di) - u(3-di))**2
  IS(3) = w_IS(1) * (u(3) - 2*u(3-di) + u(3-2*di))**2 + &
       w_IS(2) * (3*u(3) - 4*u(3-di) + u(3-2*di))**2

  ! Compute weighted average of ENO schemes
  alpha = weights/(weno_eps + IS)**2
  inv_alpha_sum = 1/sum(alpha)
  weno5_reconstruct = inv_alpha_sum * sum(alpha * ENO)
end function weno5_reconstruct
