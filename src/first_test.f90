program first_test

  use m_foap4

  implicit none
  integer, parameter :: dp = kind(0.0d0)

  type(foap4_t) :: f4
  integer       :: n

  call f4_initialize(f4)

  n = 0
  call test_refinement(f4, [0.5_dp, 0.5_dp], "test", n)
  call test_refinement(f4, [0.99_dp, 1e-2_dp], "test", n)
  call test_refinement(f4, [1e-2_dp, 0.99_dp], "test", n)
  call test_refinement(f4, [0.99_dp, 0.99_dp], "test", n)
  call test_refinement(f4, [1e-2_dp, 1e-2_dp], "test", n)

  call f4_finalize(f4)

contains

  subroutine test_refinement(f4, refine_location, base_name, n_output)
    type(foap4_t), intent(inout) :: f4
    real(dp), intent(in)         :: refine_location(2)
    character(len=*), intent(in) :: base_name
    integer, intent(inout)       :: n_output
    integer, parameter           :: n_blocks_per_dim(2) = [1, 1]
    real(dp), parameter          :: block_length(2)     = [1.0_dp, 1.0_dp]
    integer, parameter           :: bx(2)               = [16, 16]
    integer, parameter           :: n_gc                = 1
    integer, parameter           :: n_vars              = 2
    character(len=20)            :: var_names(n_vars)   = ['rho', 'phi']
    logical, parameter           :: periodic(2)         = [.false., .false.]
    integer, parameter           :: min_level           = 1
    integer, parameter           :: max_blocks          = 10000
    integer                      :: n_refine_steps

    call f4_set_grid(f4, n_blocks_per_dim, block_length, &
         bx, n_gc, n_vars, var_names, &
         periodic, min_level, max_blocks)

    call set_init_cond(f4)
    call f4_update_ghostcells(f4, 2, [1, 2])

    n_output = n_output + 1
    call f4_write_grid(f4, base_name, n_output)

    do n_refine_steps = 1, 3
       call set_refinement_flag(f4, refine_location)
       call f4_adjust_refinement(f4)
       call f4_update_ghostcells(f4, 2, [1, 2])
       call local_average(f4)
       n_output = n_output + 1
       call f4_write_grid(f4, base_name, n_output)
    end do

    call f4_destroy(f4)
  end subroutine test_refinement

  pure real(dp) function rho_init(x, y)
    real(dp), intent(in) :: x, y
    rho_init = x + y
  end function rho_init

  pure real(dp) function phi_init(x, y)
    real(dp), intent(in) :: x, y
    phi_init = x - y
  end function phi_init

  subroutine set_init_cond(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, i, j
    real(dp)                     :: rr(2)

    do n = 1, f4%n_blocks
       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)
             rr = f4_cell_coord(f4, n, i, j)
             f4%uu(i, j, 1, n) = rho_init(rr(1), rr(2))
             f4%uu(i, j, 2, n) = phi_init(rr(1), rr(2))
          end do
       end do
    end do
  end subroutine set_init_cond

  subroutine set_refinement_flag(f4, r0)
    type(foap4_t), intent(inout) :: f4
    real(dp), intent(in)         :: r0(2)
    integer                      :: n, lvl
    real(dp)                     :: rmin(2), rmax(2)

    do n = 1, f4%n_blocks
       lvl = f4%block_level(n)
       rmin = f4%block_origin(:, n)
       rmax = rmin + f4%bx * f4%dr_level(:, lvl)

       if (all(r0 >= rmin .and. r0 <= rmax)) then
          f4%refinement_flags(n) = 1
       else
          f4%refinement_flags(n) = 0
       end if
    end do
  end subroutine set_refinement_flag

  subroutine local_average(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, i, j, iv
    real(dp)                     :: rr(2), err, sol
    real(dp), allocatable        :: tmp(:, :)

    allocate(tmp(f4%bx(1), f4%bx(2)))
    iv = 1

    do n = 1, f4%n_blocks
       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)
             rr = f4_cell_coord(f4, n, i, j)
             tmp(i, j) = 0.25_dp * ( &
                  f4%uu(i-1, j, iv, n) + &
                  f4%uu(i+1, j, iv, n) + &
                  f4%uu(i, j-1, iv, n) + &
                  f4%uu(i, j+1, iv, n))
             sol = rho_init(rr(1), rr(2))
             err = tmp(i, j) - sol
             if (abs(err) > 1e-15_dp) then
                print *, n, f4%block_level(n), i, j, err
             end if
          end do
       end do

       f4%uu(1:f4%bx(1), 1:f4%bx(2), iv, n) = tmp
    end do
  end subroutine local_average

end program first_test
