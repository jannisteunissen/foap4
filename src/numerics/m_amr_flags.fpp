!> Module for setting flags for adaptive mesh refinement
#:include 'definitions.fpp'
module m_amr_flags_${NDIM}$d

  use m_foap4_types_${NDIM}$d

  implicit none
  private

  ! Public methods
  public :: amr_flags_diff2

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

    !$acc parallel loop private(level, diff_norm) default(present)
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

end module m_amr_flags_${NDIM}$d
