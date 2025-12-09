module m_physics_euler_${NDIM}$d

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  real(dp), parameter :: euler_gamma  = 5.0_dp/3.0_dp
  !$acc declare create(euler_gamma)
  real(dp), parameter :: inv_gamma_m1 = 1/(euler_gamma-1.0_dp)
  !$acc declare create(inv_gamma_m1)

  integer, parameter  :: n_vars = 2 + ${NDIM}$
  integer, parameter  :: i_rho = 1
  integer, parameter  :: i_mom0 = 1
  integer, parameter  :: i_e = n_vars
#:if NDIM == 2
  integer, parameter :: i_vars(n_vars) = [i_rho, i_mom0+1, i_mom0+2, i_e]
  character(len=10), parameter :: var_names(n_vars) = [character(len=10) :: &
       "rho", "momx", "momy", "e"]
#:elif NDIM == 3
  integer, parameter :: i_vars(n_vars) = [i_rho, i_mom0+1, i_mom0+2, i_mom0+3, i_e]
  character(len=10), parameter :: var_names(n_vars) = [character(len=10) :: &
       "rho", "momx", "momy", "momz", "e"]
#:endif

  real(dp) :: gravity_constant = 0.0_dp
  !$acc declare create(gravity_constant)

end module m_physics_euler_${NDIM}$d
