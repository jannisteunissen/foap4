module m_euler

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  real(dp), parameter :: euler_gamma  = 5.0_dp/3.0_dp
  real(dp), parameter :: inv_gamma_m1 = 1/(euler_gamma-1.0_dp)

  integer, parameter  :: n_vars = 2 + ${NDIM}$
  integer, parameter  :: i_rho       = 1
#:if NDIM == 2
  integer, parameter  :: i_mom(${NDIM}$)    = [2, 3]
  integer, parameter  :: i_e         = 4
  character(len=10), parameter :: var_names(n_vars) = [character(len=10) :: &
       "rho", "momx", "momy", "e"]
#:elif NDIM == 3
  integer, parameter  :: i_mom(${NDIM}$) = [2, 3, 4]
  integer, parameter  :: i_e         = 5
  character(len=10), parameter :: var_names(n_vars) = [character(len=10) :: &
       "rho", "momx", "momy", "momz", "e"]
#:endif

  integer, parameter :: i_vars(n_vars) = [i_rho, i_mom, i_e]

end module m_euler
