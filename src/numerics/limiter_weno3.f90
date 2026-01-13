pure subroutine reconstruct(u, i0, u_LR)
  !$acc routine seq
  real(dp), intent(in)  :: u(5, n_vars)
  integer, intent(in)   :: i0
  real(dp), intent(out) :: u_LR(n_vars, 2)
  integer               :: n

  do n = 1, n_vars
     u_LR(n, 1) = weno3_reconstruct(u(i0+1:i0+3, n), -1)
     u_LR(n, 2) = weno3_reconstruct(u(i0+2:i0+4, n), 1)
  end do
end subroutine reconstruct

!> Perform weno3 reconstruction to a cell face
!> Reference: Jiang & Shu "Efficient Implementation of Weighted ENO Schemes", JCP 1996
pure real(dp) function weno3_reconstruct(u, di)
  !$acc routine seq

  !> Cell-centered values. The face value is reconstructed at one side of
  !> the central cell at u(2), depending on the argument di.
  real(dp), intent(in) :: u(3)
  !> Direction, which can be 1 or +1. If -1, the 'left' side of the face
  !> between u(2) and u(3) is reconstructed. If +1, the 'right' side of the
  !> face between u(1) and u(2) is reconstructed.
  integer, intent(in)  :: di
  real(dp), parameter  :: c1(2)      = [-0.5_dp, 1.5_dp]
  real(dp), parameter  :: c2(2)      = [0.5_dp, 0.5_dp]
  real(dp), parameter  :: weights(2) = [1, 2]/3.0_dp
  real(dp), parameter  :: weno_eps   = 1e-18_dp
  real(dp)             :: ENO(2), alpha(2), IS(2), inv_alpha_sum

  ! ENO schemes
  ENO(1) = u(2+di)*c1(1) + u(2)*c1(2)
  ENO(2) = u(2)*c2(1) + u(2-di)*c2(2)

  ! Smoothness indicators
  IS(1) = (u(2+di) - u(2))**2
  IS(2) = (u(2) - u(2-di))**2

  ! Compute weighted average of ENO schemes
  alpha = weights/(weno_eps + IS)**2
  inv_alpha_sum = 1/(alpha(1) + alpha(2))
  weno3_reconstruct = inv_alpha_sum * (alpha(1) * ENO(1) + alpha(2) * ENO(2))
end function weno3_reconstruct
