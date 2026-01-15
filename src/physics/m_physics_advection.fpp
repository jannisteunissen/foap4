module m_physics_advection_${NDIM}$d

  implicit none
  public

  integer, parameter, private  :: dp                = kind(0.0d0)

  ! Number of temporal variables
  integer, parameter           :: n_tvars = 1

  ! Total number of variables
  integer, parameter           :: n_vars_all = n_tvars + 1

  ! Variable for storing AMR flags
  integer, parameter           :: i_amr_flag = 1

  ! Density variable
  integer, parameter           :: i_rho = 2

  ! Indices of temporal variables
  integer, parameter           :: i_tvars(n_tvars) = [i_rho]

  ! Index offset for temporal variables
  integer, parameter           :: i_tvars0 = i_rho - 1

  ! Names of variables
  character(len=10), parameter :: var_names(n_vars_all) = &
       [character(len=10)      :: "amr_flag", "rho"]

  ! Velocity
  real(dp)                     :: velocity(${NDIM}$) = 1.0_dp
  !$acc declare create(velocity)

  ! Which variables are temporal
  logical, parameter           :: var_temporal(n_vars_all) = [.false., .true.]

end module m_physics_advection_${NDIM}$d
