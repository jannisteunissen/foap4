module m_advection

  implicit none
  public

  integer, parameter, private  :: dp                = kind(0.0d0)
  integer, parameter           :: n_vars            = 1
  integer, parameter           :: i_rho             = 1
  integer, parameter           :: i_vars(n_vars)    = [i_rho]
  character(len=10), parameter :: var_names(n_vars) = &
       [character(len=10)      :: "rho"]
  real(dp)                     :: velocity(3)       = 1.0_dp
  !$acc declare create(velocity)

end module m_advection
