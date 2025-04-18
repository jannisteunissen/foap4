module m_euler

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  real(dp), parameter :: euler_gamma  = 5.0_dp/3.0_dp
  real(dp), parameter :: inv_gamma_m1 = 1/(euler_gamma-1.0_dp)

  integer, parameter  :: n_vars_euler = 4
  integer, parameter  :: ix_rho       = 1
  integer, parameter  :: ix_momx      = 2
  integer, parameter  :: ix_momy      = 3
  integer, parameter  :: ix_mom(2)    = [2, 3]
  integer, parameter  :: ix_e         = 4

  integer, parameter :: n_variables               = 2 * n_vars_euler
  integer, parameter :: i_rho                     = 1
  integer, parameter :: i_momx                    = 2
  integer, parameter :: i_momy                    = 3
  integer, parameter :: i_e                       = 4
  integer, parameter :: i_vars_grid(n_vars_euler) = [i_rho, i_momx, i_momy, i_e]

  character(len=10), parameter :: var_names(n_variables) = [character(len=10) :: &
       "rho", "momx", "momy", "e", "cpy1", "cpy2", "cpy3", "cpy4"]

end module m_euler
