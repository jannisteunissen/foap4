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

  type adv_params_t
     ! Velocity (in case of uniform velocity)
     real(fp) :: velocity(NDIM)  = 1.0_fp
     ! Initial location of solution
     real(dp) :: r0(NDIM)        = 0.5_dp
     ! Type of velocity - 1: constant, 2: rotation, 3: swirl
     integer  :: vtype           = 1
     ! Whether to use a Gaussian solution
     logical  :: use_gaussian    = .false.
  end type adv_params_t

  type(adv_params_t) :: adv_par
  ${DECLARE_DEVICE('adv_par')}$

  ! Which variables are temporal
  logical, parameter           :: var_temporal(n_vars_all) = [.false., .true.]

contains

  subroutine advection_initialize(velocity, use_gaussian, vtype, r0)
    real(dp), intent(in) :: velocity(${NDIM}$)
    logical, intent(in)  :: use_gaussian
    integer, intent(in)  :: vtype
    real(dp), intent(in) :: r0(${NDIM}$)

    adv_par%velocity = real(velocity, fp)
    adv_par%use_gaussian = use_gaussian
    adv_par%vtype = vtype
    adv_par%r0 = r0
    ${UPDATE_DEVICE('adv_par%velocity, adv_par%r0, adv_par%vtype, adv_par%use_gaussian')}$
  end subroutine advection_initialize

end module m_physics_advection_${NDIM}$d
