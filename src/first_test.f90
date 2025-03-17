program first_test

  use m_foap4

  implicit none
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: n_blocks_per_dim(2) = [1, 1]
  real(dp), parameter :: block_length(2) = [1.0_dp, 1.0_dp]
  integer, parameter :: bx(2) = [16, 16]
  integer, parameter :: n_gc = 1
  integer, parameter :: n_vars = 2
  character(len=20) :: var_names(n_vars) = ['rho', 'phi']
  logical, parameter :: periodic(2) = [.false., .false.]
  integer, parameter :: min_level = 1
  integer, parameter :: max_blocks = 10000

  type(foap4_t) :: f4

  call f4_initialize_grid(n_blocks_per_dim, block_length, &
       bx, n_gc, n_vars, var_names, &
       periodic, min_level, max_blocks, f4)

  call set_init_cond(f4)

  call f4_write_grid(f4, "test_0")
  call f4_update_ghostcells(f4, 2, [1, 2])
  call local_average(f4)
  call f4_write_grid(f4, "test_1")
  call f4_update_ghostcells(f4, 2, [1, 2])

  call f4_adjust_refinement(f4)
  call f4_update_ghostcells(f4, 2, [1, 2])
  call f4_write_grid(f4, "test_2")
  call local_average(f4)

  call f4_finalize_grid(f4)

contains

  pure real(dp) function rho_init(x, y)
    real(dp), intent(in) :: x, y
    rho_init = x + y
  end function rho_init

  pure real(dp) function phi_init(x, y)
    real(dp), intent(in) :: x, y
    real(dp), parameter :: pi = acos(-1.0_dp)
    phi_init = x - y !sin(x * pi) * sin(y * pi)
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

  subroutine local_average(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, i, j, iv
    real(dp)                     :: rr(2), err, sol
    real(dp), allocatable        :: tmp(:, :)

    allocate(tmp(bx(1), bx(2)))
    iv = 2

    do n = 1, f4%n_blocks
       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)
             rr = f4_cell_coord(f4, n, i, j)
             tmp(i, j) = 0.25_dp * ( &
                  f4%uu(i-1, j, iv, n) + &
                  f4%uu(i+1, j, iv, n) + &
                  f4%uu(i, j-1, iv, n) + &
                  f4%uu(i, j+1, iv, n))
             sol = f4%uu(i, j, iv, n)
             err = tmp(i, j) - sol
             if (abs(err) > 1e-15_dp) then
                print *, n, i, j, err
             end if
          end do
       end do

       f4%uu(1:bx(1), 1:bx(2), iv, n) = tmp
    end do
  end subroutine local_average

end program first_test
