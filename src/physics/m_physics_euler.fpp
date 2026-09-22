!> Definitions for Euler's equations of gas dynamics
#:include 'definitions_ndim.fpp'
#:include 'definitions_parallel.fpp'
module m_physics_euler_${NDIM}$d
  use m_foap4_types_${NDIM}$d

  implicit none
  public

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
  integer, protected :: i_tvars(n_tvars) = [i_rho, i_mom0+1, i_mom0+2, i_e]
  ${DECLARE_DEVICE('i_tvars')}$

  ! Names of variables
  character(len=10), parameter :: var_names(n_vars_all) = [character(len=10) :: &
       "rho", "momx", "momy", "e"]
#:elif NDIM == 3
  ! Indices of temporal variables
  integer, protected :: i_tvars(n_tvars) = [i_rho, i_mom0+1, i_mom0+2, i_mom0+3, i_e]
  ${DECLARE_DEVICE('i_tvars')}$

  ! Names of variables
  character(len=10), parameter :: var_names(n_vars_all) = [character(len=10) :: &
       "rho", "momx", "momy", "momz", "e"]
#:endif

  ! Which variables are temporal
  logical, parameter :: var_temporal(n_vars_all) = .true.

  type euler_par_t
     real(fp) :: gamma        = 5/3.0_fp
     real(fp) :: inv_gamma_m1 = 1/(5/3.0_fp - 1)
     real(fp) :: gravity      = 0.0_dp
     real(fp) :: rho_floor    = 0.0_dp
     real(fp) :: p_floor      = 0.0_dp
  end type euler_par_t

  type(euler_par_t) :: euler_par
  ${DECLARE_DEVICE('euler_par')}$

contains

  subroutine euler_initialize(gamma, gravity, rho_floor, p_floor)
    real(dp), intent(in) :: gamma
    real(dp), intent(in) :: gravity
    real(dp), intent(in) :: rho_floor
    real(dp), intent(in) :: p_floor

    euler_par%gamma = real(gamma, fp)
    euler_par%inv_gamma_m1 = real(1/(gamma-1), fp)

    euler_par%gravity = real(gravity, fp)

    euler_par%rho_floor = real(rho_floor, fp)
    euler_par%p_floor = real(p_floor, fp)

    ${UPDATE_DEVICE('euler_par%gamma, euler_par%inv_gamma_m1')}$
    ${UPDATE_DEVICE('euler_par%gravity')}$
    ${UPDATE_DEVICE('euler_par%rho_floor, euler_par%p_floor')}$
  end subroutine euler_initialize

end module m_physics_euler_${NDIM}$d
