#:include 'definitions_ndim.fpp'
#:include 'definitions_parallel.fpp'
module m_physics_advection_${NDIM}$d
  use m_foap4_types_${NDIM}$d

  implicit none
  public

  ! Number of temporal variables
  integer, parameter           :: n_tvars = 1

  ! Total number of variables
  integer, parameter           :: n_vars_all = n_tvars + 1

  ! Variable for storing error
  integer, parameter           :: i_error = 1

  ! Density variable
  integer, parameter           :: i_rho = 2

  ! Indices of temporal variables
  integer, parameter           :: i_tvars(n_tvars) = [i_rho]

  ! Index offset for temporal variables
  integer, parameter           :: i_tvars0 = i_rho - 1

  ! Names of variables
  character(len=10), parameter :: var_names(n_vars_all) = &
       [character(len=10)      :: "error", "rho"]

  ! Velocity (in case of uniform velocity)
  real(fp)                     :: advection_velocity(${NDIM}$) = 1.0_fp
  ${DECLARE_DEVICE('advection_velocity')}$

  ! Initial location of solution
  real(dp)                     :: advection_r0(NDIM) = 0.5_dp
  ${DECLARE_DEVICE('advection_r0')}$

  ! Whether to use a Gaussian solution
  logical                      :: advection_use_gaussian = .false.
  ${DECLARE_DEVICE('advection_use_gaussian')}$

  ! Type of velocity - 1: constant, 2: rotation, 3: swirl
  integer :: advection_velocity_type = 1
  ${DECLARE_DEVICE('advection_velocity_type')}$

  ! Which variables are temporal
  logical, parameter           :: var_temporal(n_vars_all) = [.false., .true.]

contains

  subroutine advection_initialize(velocity, use_gaussian, vtype, r0)
    real(dp), intent(in) :: velocity(${NDIM}$)
    logical, intent(in)  :: use_gaussian
    integer, intent(in)  :: vtype
    real(dp), intent(in) :: r0(${NDIM}$)

    advection_velocity = real(velocity, fp)
    ${UPDATE_DEVICE('advection_velocity')}$

    advection_use_gaussian = use_gaussian
    ${UPDATE_DEVICE('advection_use_gaussian')}$

    advection_velocity_type = vtype
    ${UPDATE_DEVICE('advection_velocity_type')}$

    advection_r0 = r0
    ${UPDATE_DEVICE('advection_r0')}$
  end subroutine advection_initialize

end module m_physics_advection_${NDIM}$d
