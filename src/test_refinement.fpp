#:mute
#:include "fypp_definitions.fpp"
#:endmute
program test_ref

  use m_foap4

  implicit none
  integer, parameter :: dp = kind(0.0d0)

  type(foap4_t) :: f4
  integer       :: n, n_gc
  logical       :: test_coarsening
  logical       :: write_output

  call f4_initialize(f4, "error")

  write_output = .false.
  test_coarsening = .true.
  n = 0
  do n_gc = 1, 4
     if (f4%mpirank == 0) print *, "Testing uniform grid with n_gc =", n_gc
     call test_refinement(f4, n_gc, 4, 0, [1e-2_dp, 1e-2_dp], test_coarsening, &
          write_output, "output/test_unif", n)
  end do

  test_coarsening = .false.
  n = 0
  do n_gc = 1, 4
     if (f4%mpirank == 0) print *, "Testing with n_gc = ", n_gc
     call test_refinement(f4, n_gc, 3, 7, [1e-2_dp, 1e-2_dp], test_coarsening, &
          write_output, "output/test_ref", n)
     call test_refinement(f4, n_gc, 3, 7, [0.99_dp, 1e-2_dp], test_coarsening, &
          write_output, "output/test_ref", n)
     call test_refinement(f4, n_gc, 3, 7, [0.5_dp, 0.5_dp], test_coarsening, &
          write_output, "output/test_ref", n)
     call test_refinement(f4, n_gc, 3, 7, [1e-2_dp, 0.99_dp], test_coarsening, &
          write_output, "output/test_ref", n)
     call test_refinement(f4, n_gc, 3, 7, [0.99_dp, 0.99_dp], test_coarsening, &
          write_output, "output/test_ref", n)
  end do

  if (f4%mpirank == 0) call f4_print_wtime(f4)
  call f4_finalize(f4)

contains

  subroutine test_refinement(f4, n_gc, min_level, n_refine_steps, &
       refine_location, test_coarsening, write_output, base_name, n_output)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_gc
    integer, intent(in)          :: min_level
    integer, intent(in)          :: n_refine_steps
    real(dp), intent(in)         :: refine_location(2)
    logical, intent(in)          :: test_coarsening
    logical, intent(in)          :: write_output
    character(len=*), intent(in) :: base_name
    integer, intent(inout)       :: n_output
    integer, parameter           :: n_blocks_per_dim(2) = [1, 1]
    real(dp), parameter          :: block_length(2)     = [1.0_dp, 1.0_dp]
    integer, parameter           :: bx(2)               = [8, 8]
    integer, parameter           :: n_vars              = 2
    character(len=20)            :: var_names(n_vars)   = ['rho', 'phi']
    logical, parameter           :: periodic(2)         = [.false., .false.]
    integer, parameter           :: max_blocks          = 1000
    logical, parameter           :: partition           = .true.
    integer                      :: n, prev_mesh_revision

    call f4_construct_brick(f4, n_blocks_per_dim, block_length, bx, n_gc, &
         n_vars, var_names, periodic, min_level, max_blocks, &
         f4_bc_linear_extrap, 0.0_dp)

    call set_init_cond(f4)
    call f4_update_ghostcells(f4, 2, [1, 2])
    call local_average(f4)

    if (write_output) then
       n_output = n_output + 1
       call f4_write_grid(f4, base_name, n_output)
    end if

    do n = 1, n_refine_steps
       call set_refinement_flag(f4, refine_location)
       call f4_adjust_refinement(f4, partition)
       call f4_update_ghostcells(f4, 2, [1, 2])
       call local_average(f4)

       if (write_output) then
          n_output = n_output + 1
          call f4_write_grid(f4, base_name, n_output)
       end if
    end do

    if (test_coarsening) then
       do n = 1, 10
          prev_mesh_revision = f4_get_mesh_revision(f4)
          call set_coarsening_flag(f4)
          call f4_adjust_refinement(f4, partition)

          call f4_update_ghostcells(f4, 2, [1, 2])
          call local_average(f4)

          if (write_output) then
             n_output = n_output + 1
             call f4_write_grid(f4, base_name, n_output)
          end if

          if (f4_get_mesh_revision(f4) == prev_mesh_revision) exit
       end do
    end if

    call f4_destroy(f4)
  end subroutine test_refinement

  pure real(dp) function rho_init(x, y)
    @{GPU_ROUTINE_SEQ()}@
    real(dp), intent(in) :: x, y
    rho_init = x + y
  end function rho_init

  pure real(dp) function phi_init(x, y)
    @{GPU_ROUTINE_SEQ()}@
    real(dp), intent(in) :: x, y
    phi_init = x - y
  end function phi_init

  subroutine set_init_cond(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, i, j
    real(dp)                     :: rr(2)

    @{GPU_PARALLEL_LOOP()}@
    do n = 1, f4%n_blocks
       @{GPU_LOOP()}@ collapse(2) private(rr)
       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)
             rr = f4_cell_coord(f4, n, i, j)
             f4%uu(i, j, 1, n) = sum(rr) !rho_init(rr(1), rr(2))
             f4%uu(i, j, 2, n) = rr(1) - rr(2) !phi_init(rr(1), rr(2))
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

  subroutine set_coarsening_flag(f4)
    type(foap4_t), intent(inout) :: f4
    f4%refinement_flags(1:f4%n_blocks) = -1
  end subroutine set_coarsening_flag

  subroutine local_average(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, i, j, iv
    real(dp)                     :: rr(2), err, max_err, sol
    real(dp), allocatable        :: tmp(:, :)
    real(dp), parameter          :: max_difference = 1e-15_dp

    allocate(tmp(f4%bx(1), f4%bx(2)))
    iv = 1
    max_err = 0.0_dp

    @{GPU_PARALLEL_LOOP()}@ private(tmp) reduction(max:max_err)
    do n = 1, f4%n_blocks
       @{GPU_LOOP()}@ collapse(2) private(rr, sol, err) reduction(max:max_err)
       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)
             rr = f4_cell_coord(f4, n, i, j)
             tmp(i, j) = 0.25_dp * ( &
                  f4%uu(i-1, j, iv, n) + &
                  f4%uu(i+1, j, iv, n) + &
                  f4%uu(i, j-1, iv, n) + &
                  f4%uu(i, j+1, iv, n))
             sol = sum(rr) !rho_init(rr(1), rr(2))
             err = tmp(i, j) - sol
             max_err = max(abs(err), max_err)
          end do
       end do

       @{GPU_LOOP()}@ collapse(2)
       do j = 1, f4%bx(2)
          do i = 1, f4%bx(1)
             f4%uu(i, j, iv, n) = tmp(i, j)
          end do
       end do
    end do

    if (max_err > max_difference) then
       print *, f4%mpirank, max_err
       error stop "Too large error"
    end if

  end subroutine local_average

end program test_ref
