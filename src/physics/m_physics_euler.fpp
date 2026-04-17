!> Definitions for Euler's equations of gas dynamics
module m_physics_euler_${NDIM}$d
  use m_foap4_types_${NDIM}$d

  implicit none
  public

  real(fp), protected :: euler_gamma = 5/3.0_fp
  !$acc declare create(euler_gamma)

  real(fp), protected :: euler_inv_gamma_m1 = 1/(5/3.0_fp - 1)
  !$acc declare create(euler_inv_gamma_m1)

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

  real(fp) :: euler_gravity = 0.0_dp
  !$acc declare create(euler_gravity)

  real(fp) :: euler_rho_floor = 0.0_dp
  !$acc declare create(euler_rho_floor)

  real(fp) :: euler_p_floor = 0.0_dp
  !$acc declare create(euler_p_floor)

contains

  subroutine euler_initialize(gamma, gravity, rho_floor, p_floor)
    real(dp), intent(in) :: gamma
    real(dp), intent(in) :: gravity
    real(dp), intent(in) :: rho_floor
    real(dp), intent(in) :: p_floor

    euler_gamma = real(gamma, fp)
    euler_inv_gamma_m1 = 1/(euler_gamma-1)

    euler_gravity = real(gravity, fp)

    euler_rho_floor = real(rho_floor, fp)
    euler_p_floor = real(p_floor, fp)

    !$acc update device(euler_gamma, euler_inv_gamma_m1)
    !$acc update device(euler_gravity)
    !$acc update device(euler_rho_floor, euler_p_floor)
  end subroutine euler_initialize

end module m_physics_euler_${NDIM}$d
