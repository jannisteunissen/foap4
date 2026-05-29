!> Module for Runge-Kutta time integration
#:include 'definitions.fpp'
#:include 'definitions_parallel.fpp'
module m_rk_${NDIM}$d

  use m_foap4_types_${NDIM}$d

  implicit none
  private

  !> Number of time integration schemes
  integer, parameter, public :: rk_num_integrators  = 6

  !> Forward Euler method
  integer, parameter, public :: rk_forward_euler    = 1
  !> Heun's method (AKA modified Euler's method, explicit trapezoidal rule), CFL
  !> coefficient of 1. See e.g. https://en.wikipedia.org/wiki/Heun%27s_method
  integer, parameter, public :: rk_heuns_method     = 2
  !> Midpoint method, see e.g. https://en.wikipedia.org/wiki/Midpoint_method
  integer, parameter, public :: rk_midpoint_method  = 3
  !> Optimal 3-stage third-order SSPRK method (Shu & Osher), CFL coefficient of
  !> 1, see e.g. https://doi.org/10.1137/S0036142902419284
  integer, parameter, public :: rk_ssprk33_method   = 4
  !> Optimal 4-stage third-order SSPRK method (Ruuth & Spiteri), CFL coefficient
  !> of 2, see e.g. https://doi.org/10.1137/S0036142902419284
  integer, parameter, public :: rk_ssprk43_method   = 5
  !> Classic 4th order Runge Kutta method, see e.g.
  !> https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
  integer, parameter, public :: rk_rk4_method       = 6

  character(len=20), public :: rk_integrator_names(rk_num_integrators) = &
       [character(len=20) :: "forward_euler", "heuns_method", &
       "midpoint_method", "ssprk33", "ssprk43", "rk4"]

  !> How many steps the time integrators take
  integer, parameter, public :: &
       rk_advance_num_steps(rk_num_integrators) = [1, 2, 2, 3, 4, 4]

  !> How many variable copies are required for the time integrators
  integer, parameter, public :: rk_advance_num_copies(rk_num_integrators) = &
       rk_advance_num_steps

  abstract interface

     !> Interface for a generic forward Euler scheme for time integration
     !>
     !> This method should advance the solution over a time dt. The method
     !> assumes that several copies are stored for the variables to be
     !> integrated. It should then operate on these different copies, which
     !> correspond to temporal states. In this way, higher-order time
     !> integration schemes can be constructed.
     !>
     !> The meaning of the temporal states is as follows. For an equation y' =
     !> f(y), the method should perform:
     !> y_out = sum(w_prev * y_prev) + dt * f(y_deriv).
     !>
     !> If the index of the variable `y` is `i`, then the index of `y_out` is
     !> `i+s_out`, etc.
     subroutine subr_feuler(f4, dt, dt_lim, s_deriv, n_prev, &
          s_prev, w_prev, s_out, i_step, n_steps)
       import
       type(foap4_t), intent(inout) :: f4
       real(dp), intent(in)         :: dt             !< Time step
       real(dp), intent(inout)      :: dt_lim         !< Computed time step limit
       integer, intent(in)          :: s_deriv        !< State to compute derivatives from
       integer, intent(in)          :: n_prev         !< Number of previous states
       integer, intent(in)          :: s_prev(n_prev) !< Previous states
       real(dp), intent(in)         :: w_prev(n_prev) !< Weights of previous states
       integer, intent(in)          :: s_out          !< Output state
       integer, intent(in)          :: i_step         !< Step of the integrator
       integer, intent(in)          :: n_steps        !< Total number of steps
     end subroutine subr_feuler
  end interface

  ! Public methods
  public :: rk_get_integrator_by_name
  public :: rk_advance

contains

  !> Get the time integrator with a given name
  integer function rk_get_integrator_by_name(name) result(n)
    character(len=*), intent(in) :: name

    do n = 1, rk_num_integrators
       if (name == rk_integrator_names(n)) exit
    end do

    if (n == rk_num_integrators + 1) then
       print *, name, " not in list of available time integrators:"
       do n = 1, rk_num_integrators
          print *, n, rk_integrator_names(n)
       end do
       error stop "Invalid name for time integrator"
    end if
  end function rk_get_integrator_by_name

  !> Advance solution over dt using time_integrator scheme, which can call
  !> forward_euler multiple times
  subroutine rk_advance(f4, dt, dt_lim, time_integrator, forward_euler)
    type(foap4_t), intent(inout) :: f4
    real(dp), intent(in)         :: dt     !< Current time step
    real(dp), intent(out)        :: dt_lim !< Time step limit
    !> One of the pre-defined time integrators (e.g. af_heuns_method)
    integer, intent(in)          :: time_integrator
    !> Forward Euler method provided by the user
    procedure(subr_feuler)       :: forward_euler
    integer                      :: n_steps, offset

    real(dp)            :: time_in
    real(dp), parameter :: third = 1/3.0_dp
    real(dp), parameter :: sixth = 1/6.0_dp

    if (time_integrator < 1 .or. time_integrator > rk_num_integrators) &
         error stop "Invalid time integrator"

    n_steps = rk_advance_num_steps(time_integrator)
    offset = f4%max_blocks
    dt_lim = 1e100_dp
    time_in = f4%time

    select case (time_integrator)
    case (rk_forward_euler)
       call forward_euler(f4, dt, dt_lim, 0, &
            1, [0], [1.0_dp], 0, 1, n_steps)
    case (rk_midpoint_method)
       call forward_euler(f4, 0.5_dp*dt, dt_lim, 0, &
            1, [0], [1.0_dp], offset, 1, n_steps)
       f4%time = time_in + 0.5_dp*dt
       call forward_euler(f4, dt, dt_lim, offset, &
            1, [0], [1.0_dp], 0, 2, n_steps)
    case (rk_heuns_method)
       call forward_euler(f4, dt, dt_lim, 0, &
            1, [0], [1.0_dp], offset, 1, n_steps)
       f4%time = time_in + dt
       call forward_euler(f4, 0.5_dp*dt, dt_lim, offset, &
            2, [0, offset], [0.5_dp, 0.5_dp], 0, 2, n_steps)
    case (rk_ssprk33_method)
       call forward_euler(f4, dt, dt_lim, 0, &
            1, [0], [1.0_dp], offset, 1, n_steps)
       f4%time = time_in + dt
       call forward_euler(f4, 0.25_dp*dt, dt_lim, offset, &
            2, [0, offset], [0.75_dp, 0.25_dp], 2*offset, 2, n_steps)
       f4%time = time_in + 0.5_dp*dt
       call forward_euler(f4, 2*third*dt, dt_lim, 2*offset, &
            2, [0, 2*offset], [third, 2*third], 0, 3, n_steps)
    case (rk_ssprk43_method)
       call forward_euler(f4, 0.5_dp*dt, dt_lim, 0, &
            1, [0], [1.0_dp], offset, 1, n_steps)
       f4%time = time_in + 0.5_dp * dt
       call forward_euler(f4, 0.5_dp*dt, dt_lim, offset, &
            1, [offset], [1.0_dp], 2*offset, 2, n_steps)
       f4%time = time_in + dt
       call forward_euler(f4, sixth*dt, dt_lim, 2*offset, &
            2, [0, 2*offset], [2*third, third], 3*offset, 3, n_steps)
       f4%time = time_in + 0.5_dp * dt
       call forward_euler(f4, 0.5_dp*dt, dt_lim, 3*offset, &
            1, [3*offset], [1.0_dp], 0, 4, n_steps)
    case (rk_rk4_method)
       ! This looks different than the standard formulation in most textbooks.
       ! The idea is to construct the states needed for the derivatives, and
       ! then take a linear combination. Note the negative coefficient used in
       ! the last step.
       call forward_euler(f4, 0.5_dp*dt, dt_lim, 0, &
            1, [0], [1.0_dp], offset, 1, n_steps)
       f4%time = time_in + 0.5_dp * dt
       call forward_euler(f4, 0.5_dp*dt, dt_lim, offset, &
            1, [0], [1.0_dp], 2*offset, 2, n_steps)
       f4%time = time_in + 0.5_dp * dt
       call forward_euler(f4, dt, dt_lim, 2*offset, &
            1, [0], [1.0_dp], 3*offset, 3, n_steps)
       f4%time = time_in + dt
       call forward_euler(f4, sixth*dt, dt_lim, 3*offset, &
            4, [0, offset, 2*offset, 3*offset], &
            [-third, third, 2*third, third], 0, 4, n_steps)
    case default
       error stop "Unknown time integrator"
    end select

    f4%time = time_in + dt
    ${UPDATE_DEVICE('f4%time')}$
  end subroutine rk_advance

end module m_rk_${NDIM}$d

