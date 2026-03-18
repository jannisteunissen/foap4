module m_physics_advection_${NDIM}$d

  implicit none
  public

  integer, parameter, private  :: dp                = kind(0.0d0)

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

  ! Velocity
  real(dp)                     :: advection_velocity(${NDIM}$) = 1.0_dp
  !$acc declare create(advection_velocity)

  ! Whether to use a Gaussian solution
  logical                      :: advection_use_gaussian = .false.
  !$acc declare create(advection_use_gaussian)

  ! Which variables are temporal
  logical, parameter           :: var_temporal(n_vars_all) = [.false., .true.]

contains

  subroutine advection_initialize(velocity, use_gaussian)
    real(dp), intent(in) :: velocity(${NDIM}$)
    logical, intent(in)  :: use_gaussian

    advection_velocity = velocity
    !$acc update device(advection_velocity)

    advection_use_gaussian = use_gaussian
    !$acc update device(advection_use_gaussian)
  end subroutine advection_initialize

end module m_physics_advection_${NDIM}$d
