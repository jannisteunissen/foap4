!> Module for setting flags for adaptive mesh refinement
#:include 'definitions.fpp'
module m_amr_flags_${NDIM}$d

  use m_foap4_${NDIM}$d

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! Public methods
  public :: amr_flags_diff2
  public :: amr_flags_diff2_smooth

contains

  !> Set refinement flags based on an estimate of the second derivative. This
  !> is similar to the FLASH code but not equivalent, since we do not compute
  !> cross derivatives.
  subroutine amr_flags_diff2(f4, min_level, max_level, iv, &
       c_refine, c_derefine, c_eps, c_abs)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: min_level
    integer, intent(in)          :: max_level
    integer, intent(in)          :: iv
    real(dp), intent(in)         :: c_refine
    real(dp), intent(in)         :: c_derefine
    real(dp), intent(in)         :: c_eps
    real(dp), intent(in)         :: c_abs

    integer             :: n, ${IJK}$, level
    real(dp)            :: diff(NDIM), diff_norm

    !$acc parallel loop private(level, diff_norm) &
    !$acc &present(f4, f4%uu, f4%bx, f4%block_level, f4%refinement_flags)
    do n = 1, f4%n_blocks
       level = f4%block_level(n)
       diff_norm = 0.0_dp

#:if NDIM == 2
       !$acc loop collapse(2) private(diff) reduction(max: diff_norm)
       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)
             diff(1) = abs(f4%uu(i+1, j, iv, n) - 2 * f4%uu(i, j, iv, n) + &
                  f4%uu(i-1, j, iv, n)) / (c_abs + &
                  abs(f4%uu(i+1, j, iv, n) - f4%uu(i, j, iv, n)) + &
                  abs(f4%uu(i, j, iv, n) - f4%uu(i-1, j, iv, n)) + &
                  c_eps * (abs(f4%uu(i+1, j, iv, n)) + &
                  2 * abs(f4%uu(i, j, iv, n)) + abs(f4%uu(i-1, j, iv, n))))

             diff(2) = abs(f4%uu(i, j+1, iv, n) - 2 * f4%uu(i, j, iv, n) + &
                  f4%uu(i, j-1, iv, n)) / (c_abs + &
                  abs(f4%uu(i, j+1, iv, n) - f4%uu(i, j, iv, n)) + &
                  abs(f4%uu(i, j, iv, n) - f4%uu(i, j-1, iv, n)) + &
                  c_eps * (abs(f4%uu(i, j+1, iv, n)) + &
                  2 * abs(f4%uu(i, j, iv, n)) + abs(f4%uu(i, j-1, iv, n))))

             diff_norm = max(diff_norm, sqrt(diff(1)**2 + diff(2)**2))
          end do
       end do
#:elif NDIM == 3
       !$acc loop collapse(3) private(diff) reduction(max: diff_norm)
       do k = 1, f4%bx(3)
          do j = 1, f4%bx(2)
             do i = 1, f4%bx(1)
                diff(1) = abs(f4%uu(i+1, j, k, iv, n) - 2 * f4%uu(i, j, k, iv, n) + &
                     f4%uu(i-1, j, k, iv, n)) / (c_abs + &
                     abs(f4%uu(i+1, j, k, iv, n) - f4%uu(i, j, k, iv, n)) + &
                     abs(f4%uu(i, j, k, iv, n) - f4%uu(i-1, j, k, iv, n)) + &
                     c_eps * (abs(f4%uu(i+1, j, k, iv, n)) + &
                     2 * abs(f4%uu(i, j, k, iv, n)) + abs(f4%uu(i-1, j, k, iv, n))))

                diff(2) = abs(f4%uu(i, j+1, k, iv, n) - 2 * f4%uu(i, j, k, iv, n) + &
                     f4%uu(i, j-1, k, iv, n)) / (c_abs + &
                     abs(f4%uu(i, j+1, k, iv, n) - f4%uu(i, j, k, iv, n)) + &
                     abs(f4%uu(i, j, k, iv, n) - f4%uu(i, j-1, k, iv, n)) + &
                     c_eps * (abs(f4%uu(i, j+1, k, iv, n)) + &
                     2 * abs(f4%uu(i, j, k, iv, n)) + abs(f4%uu(i, j-1, k, iv, n))))

                diff(3) = abs(f4%uu(i, j, k+1, iv, n) - 2 * f4%uu(i, j, k, iv, n) + &
                     f4%uu(i, j, k-1, iv, n)) / (c_abs + &
                     abs(f4%uu(i, j, k+1, iv, n) - f4%uu(i, j, k, iv, n)) + &
                     abs(f4%uu(i, j, k, iv, n) - f4%uu(i, j, k-1, iv, n)) + &
                     c_eps * (abs(f4%uu(i, j, k+1, iv, n)) + &
                     2 * abs(f4%uu(i, j, k, iv, n)) + abs(f4%uu(i, j, k-1, iv, n))))

                diff_norm = max(diff_norm, &
                     sqrt(diff(1)**2 + diff(2)**2 + diff(3)**2))
             end do
          end do
       end do
#:endif

       if ((diff_norm > c_refine .and. level < max_level) .or. &
            level < min_level) then
          f4%refinement_flags(n) = 1
       else if ((diff_norm < c_derefine .and. level > min_level) .or. &
            level > max_level) then
          f4%refinement_flags(n) = -1
       else
          f4%refinement_flags(n) = 0
       end if
    end do

    !$acc update host(f4%refinement_flags(1:f4%n_blocks))
  end subroutine amr_flags_diff2

  !> Set refinement flags as in amr_flags_diff2, but with smoothing
  subroutine amr_flags_diff2_smooth(f4, iv_flag, n_smoothing_steps, min_level, &
       max_level, iv, c_refine, c_derefine, c_eps, c_abs)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: iv_flag !< Variable for refinement flag
    integer, intent(in)          :: n_smoothing_steps
    integer, intent(in)          :: min_level
    integer, intent(in)          :: max_level
    integer, intent(in)          :: iv !< Variable to base refinement on
    real(dp), intent(in)         :: c_refine
    real(dp), intent(in)         :: c_derefine
    real(dp), intent(in)         :: c_eps
    real(dp), intent(in)         :: c_abs

    integer             :: n, ${IJK}$, level
    real(dp)            :: diff(NDIM), diff_norm, max_flag
    real(dp), parameter :: flag_coarsen = 0.0_dp
    real(dp), parameter :: flag_stay = 1.0_dp
    real(dp), parameter :: flag_refine = 2.0_dp

    !$acc parallel loop private(level, diff_norm) &
    !$acc &present(f4, f4%uu, f4%bx)
    do n = 1, f4%n_blocks

#:if NDIM == 2
       !$acc loop collapse(2) private(diff)
       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)
             diff(1) = abs(f4%uu(i+1, j, iv, n) - 2 * f4%uu(i, j, iv, n) + &
                  f4%uu(i-1, j, iv, n)) / (c_abs + &
                  abs(f4%uu(i+1, j, iv, n) - f4%uu(i, j, iv, n)) + &
                  abs(f4%uu(i, j, iv, n) - f4%uu(i-1, j, iv, n)) + &
                  c_eps * (abs(f4%uu(i+1, j, iv, n)) + &
                  2 * abs(f4%uu(i, j, iv, n)) + abs(f4%uu(i-1, j, iv, n))))

             diff(2) = abs(f4%uu(i, j+1, iv, n) - 2 * f4%uu(i, j, iv, n) + &
                  f4%uu(i, j-1, iv, n)) / (c_abs + &
                  abs(f4%uu(i, j+1, iv, n) - f4%uu(i, j, iv, n)) + &
                  abs(f4%uu(i, j, iv, n) - f4%uu(i, j-1, iv, n)) + &
                  c_eps * (abs(f4%uu(i, j+1, iv, n)) + &
                  2 * abs(f4%uu(i, j, iv, n)) + abs(f4%uu(i, j-1, iv, n))))

             diff_norm = sqrt(diff(1)**2 + diff(2)**2)

             if (diff_norm > c_refine) then
                f4%uu(i, j, iv_flag, n) = flag_refine
             else if (diff_norm < c_derefine) then
                f4%uu(i, j, iv_flag, n) = flag_coarsen
             else
                f4%uu(i, j, iv_flag, n) = flag_stay
             end if
          end do
       end do
#:elif NDIM == 3
       !$acc loop collapse(3) private(diff)
       do k = 1, f4%bx(3)
          do j = 1, f4%bx(2)
             do i = 1, f4%bx(1)
                diff(1) = abs(f4%uu(i+1, j, k, iv, n) - 2 * f4%uu(i, j, k, iv, n) + &
                     f4%uu(i-1, j, k, iv, n)) / (c_abs + &
                     abs(f4%uu(i+1, j, k, iv, n) - f4%uu(i, j, k, iv, n)) + &
                     abs(f4%uu(i, j, k, iv, n) - f4%uu(i-1, j, k, iv, n)) + &
                     c_eps * (abs(f4%uu(i+1, j, k, iv, n)) + &
                     2 * abs(f4%uu(i, j, k, iv, n)) + abs(f4%uu(i-1, j, k, iv, n))))

                diff(2) = abs(f4%uu(i, j+1, k, iv, n) - 2 * f4%uu(i, j, k, iv, n) + &
                     f4%uu(i, j-1, k, iv, n)) / (c_abs + &
                     abs(f4%uu(i, j+1, k, iv, n) - f4%uu(i, j, k, iv, n)) + &
                     abs(f4%uu(i, j, k, iv, n) - f4%uu(i, j-1, k, iv, n)) + &
                     c_eps * (abs(f4%uu(i, j+1, k, iv, n)) + &
                     2 * abs(f4%uu(i, j, k, iv, n)) + abs(f4%uu(i, j-1, k, iv, n))))

                diff(3) = abs(f4%uu(i, j, k+1, iv, n) - 2 * f4%uu(i, j, k, iv, n) + &
                     f4%uu(i, j, k-1, iv, n)) / (c_abs + &
                     abs(f4%uu(i, j, k+1, iv, n) - f4%uu(i, j, k, iv, n)) + &
                     abs(f4%uu(i, j, k, iv, n) - f4%uu(i, j, k-1, iv, n)) + &
                     c_eps * (abs(f4%uu(i, j, k+1, iv, n)) + &
                     2 * abs(f4%uu(i, j, k, iv, n)) + abs(f4%uu(i, j, k-1, iv, n))))

                diff_norm = sqrt(diff(1)**2 + diff(2)**2 + diff(3)**2)

                if (diff_norm > c_refine) then
                   f4%uu(i, j, k, iv_flag, n) = flag_refine
                else if (diff_norm < c_derefine) then
                   f4%uu(i, j, k, iv_flag, n) = flag_coarsen
                else
                   f4%uu(i, j, k, iv_flag, n) = flag_stay
                end if

             end do
          end do
       end do
#:endif

    end do

    call smooth_flags(f4, iv_flag, n_smoothing_steps)

    !$acc parallel loop private(level, max_flag) &
    !$acc &present(f4, f4%uu, f4%bx, f4%block_level, f4%refinement_flags)
    do n = 1, f4%n_blocks
       level = f4%block_level(n)
       max_flag = 0.0_dp

       !$acc loop collapse(ndim) reduction(max:max_flag)
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          max_flag = max(max_flag, f4%uu(${IJK}$, iv_flag, n))
       end do; ${KJI_CLOSE_LOOP}$

       if ((max_flag >= flag_stay .and. level < max_level) .or. &
            level < min_level) then
          f4%refinement_flags(n) = 1
       else if ((max_flag <= flag_coarsen .and. level > min_level) .or. &
            level > max_level) then
          f4%refinement_flags(n) = -1
       else
          f4%refinement_flags(n) = 0
       end if
    end do

    !$acc update host(f4%refinement_flags(1:f4%n_blocks))
  end subroutine amr_flags_diff2_smooth

  !> Smooth flags by simple averaging, using a 'cross' stencil
  subroutine smooth_flags(f4, iv_flag, n_smoothing_steps)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: iv_flag !< Variable for refinement flag
    integer, intent(in)          :: n_smoothing_steps
    integer                      :: n, step, ${IJK}$
    real(dp)                     :: fac

#:if NDIM   == 2
    real(dp) :: fcopy(f4%ilo(1):f4%ihi(1), f4%ilo(2):f4%ihi(2))
#:elif NDIM == 3
    real(dp) :: fcopy(f4%ilo(1):f4%ihi(1), f4%ilo(2):f4%ihi(2), f4%ilo(3):f4%ihi(3))
#:endif

    ! Factor to compute average of stencil
    fac = 1.0_dp/(ndim * (2*f4%n_gc + 1))

    do step = 1, n_smoothing_steps
       call f4_update_ghostcells(f4, 1, [iv_flag])

       !$acc parallel loop
       do n = 1, f4%n_blocks
          !$acc loop collapse(ndim)
          do @{KJI_LOOP_array_to_array(f4%ilo, f4%ihi)}@
             fcopy(${IJK}$) = f4%uu(${IJK}$, iv_flag, n)
          end do; ${KJI_CLOSE_LOOP}$

          !$acc loop collapse(ndim)
          do @{KJI_LOOP_1_to_array(f4%bx)}@
#:if NDIM == 2
             f4%uu(${IJK}$, iv_flag, n) = fac * (&
                  sum(fcopy(i-f4%n_gc:i+f4%n_gc, j)) + &
                  sum(fcopy(i, j-f4%n_gc:j+f4%n_gc)))
#:elif NDIM == 3
             f4%uu(${IJK}$, iv_flag, n) = fac * (&
                  sum(fcopy(i-f4%n_gc:i+f4%n_gc, j, k)) + &
                  sum(fcopy(i, j-f4%n_gc:j+f4%n_gc, k)) + &
                  sum(fcopy(i, j, k-f4%n_gc:k+f4%n_gc)))
#:endif
          end do; ${KJI_CLOSE_LOOP}$

       end do

    end do
  end subroutine smooth_flags

end module m_amr_flags_${NDIM}$d
