module m_physics_euler_${NDIM}$d

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  real(dp), parameter :: euler_gamma  = 5.0_dp/3.0_dp
  !$acc declare create(euler_gamma)
  real(dp), parameter :: inv_gamma_m1 = 1/(euler_gamma-1.0_dp)
  !$acc declare create(inv_gamma_m1)

  ! Number of temporal variables
  integer, parameter :: n_tvars = 2 + ${NDIM}$

  ! Total number of variables
  integer, parameter :: n_vars_all = n_tvars

  ! Density variable
  integer, parameter :: i_rho = 1

  ! Index offset for momentum
  integer, parameter :: i_mom0 = 1

  ! Energy variable
  integer, parameter :: i_e = n_tvars

  ! Index offset for temporal variables
  integer, parameter :: i_tvars0 = i_rho - 1

#:if NDIM == 2
  ! Indices of temporal variables
  integer, parameter :: i_tvars(n_tvars) = [i_rho, i_mom0+1, i_mom0+2, i_e]

  ! Names of variables
  character(len=10), parameter :: var_names(n_vars_all) = [character(len=10) :: &
       "rho", "momx", "momy", "e"]
#:elif NDIM == 3
  ! Indices of temporal variables
  integer, parameter :: i_tvars(n_tvars) = [i_rho, i_mom0+1, i_mom0+2, i_mom0+3, i_e]

  ! Names of variables
  character(len=10), parameter :: var_names(n_vars_all) = [character(len=10) :: &
       "rho", "momx", "momy", "momz", "e"]
#:endif

  ! Which variables are temporal
  logical, parameter :: var_temporal(n_vars_all) = .true.

  real(dp) :: gravity_constant = 0.0_dp
  !$acc declare create(gravity_constant)

end module m_physics_euler_${NDIM}$d
