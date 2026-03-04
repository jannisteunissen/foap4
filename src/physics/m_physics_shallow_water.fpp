!> Definitions for shallow water equations (2D only)
#:assert NDIM == 2
module m_physics_shallow_water_${NDIM}$d

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  ! Gravitational acceleration
  real(dp), parameter :: gravity_constant = 9.81_dp
  !$acc declare create(gravity_constant)

  ! Number of temporal variables (h, momentum components)
  integer, parameter :: n_tvars = 1 + ${NDIM}$

  ! Total number of variables
  integer, parameter :: n_vars_all = n_tvars

  ! Water height variable
  integer, parameter :: i_h = 1

  ! Index offset for momentum (h*u, h*v, ...)
  integer, parameter :: i_mom0 = 1

  ! Index offset for temporal variables
  integer, parameter :: i_tvars0 = i_h - 1

  ! Indices of temporal variables
  integer, parameter :: i_tvars(n_tvars) = [i_h, i_mom0+1, i_mom0+2]

  ! Names of variables
  character(len=10), parameter :: var_names(n_vars_all) = [character(len=10) :: &
       "h", "hu", "hv"]

  ! Which variables are temporal
  logical, parameter :: var_temporal(n_vars_all) = .true.

end module m_physics_shallow_water_${NDIM}$d
