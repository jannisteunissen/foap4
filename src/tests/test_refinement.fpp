#:include '../core/definitions.fpp'
program test_ref
  use m_config
  use m_foap4_${NDIM}$d
  use m_io_${NDIM}$d

  implicit none
  integer            :: i_err = 1
  integer            :: i_phi = 2
  type(CFG_t)        :: cfg

  integer           :: min_level
  integer           :: num_refine_steps
  type(foap4_t)     :: f4
  integer           :: n, n_gc, ${IJK}$, bx(NDIM)
  integer           :: max_blocks
  logical           :: test_coarsening
  logical           :: write_output
  logical           :: abort_on_error
  real(dp)          :: r_ref(NDIM)
  real(dp)          :: load_imbalance_threshold = 1.1_dp
  character(len=80) :: output_name
  character(len=40) :: viewer = "visit"

#:if NDIM == 2
  min_level = 3
  num_refine_steps = 7
  max_blocks = 3000
#:elif NDIM == 3
  min_level = 2
  num_refine_steps = 5
  max_blocks = 10000
#:endif

  write_output    = .false.
  test_coarsening = .true.
  bx(:)           = 8
  test_coarsening = .true.
  abort_on_error  = .true.

  call CFG_update_from_arguments(cfg)
  call CFG_add_get(cfg, 'min_refinement_level', min_level, &
       'Minimum refinement level in the domain')
  call CFG_add_get(cfg, 'num_refine_steps', num_refine_steps, &
       'Number of refinement steps to take')
  call CFG_add_get(cfg, 'bx', bx, 'Size of grid blocks')
  call CFG_add_get(cfg, 'max_blocks', max_blocks, 'Max. number of blocks')
  call CFG_add_get(cfg, 'test_coarsening', test_coarsening, &
       'Test coarsening')
  call CFG_add_get(cfg, 'write_output', write_output, &
       'Write output files')
  call CFG_add_get(cfg, 'viewer', viewer, &
       'Write XDMF output for this viewer (visit or paraview)')
  call CFG_add_get(cfg, 'abort_on_error', abort_on_error, &
       'Abort when a numerical error is detected')
  call CFG_check(cfg)

  call f4_initialize(f4, "error")


  output_name = "output/test_ref_${NDIM}$D"
  n = 0
  do n_gc = 1, 4
     if (f4%mpirank == 0) print *, "Testing uniform grid with n_gc =", n_gc
     r_ref = 0.5_dp
     call test_refinement(f4, n_gc, min_level, 0, r_ref, .true., &
          test_coarsening, write_output, trim(output_name), n)
  end do

  do n_gc = 1, 4
     if (f4%mpirank == 0) print *, "Testing uniform refinement with n_gc =", n_gc
     r_ref = 0.5_dp

     ! Use only two steps of additional refinement
     call test_refinement(f4, n_gc, min_level, 2, r_ref, .true., &
          test_coarsening, write_output, trim(output_name), n)
  end do

  do n_gc = 1, 4
     if (f4%mpirank == 0) print *, "Using n_gc = ", n_gc

#:if NDIM == 2
     do j = -1, 1
        do i = -1, 1
           r_ref = [0.5_dp + i * 0.49_dp, 0.5_dp + j * 0.49_dp]
           if (f4%mpirank == 0) print *, "Refine around", r_ref
           call test_refinement(f4, n_gc, min_level, num_refine_steps, &
                r_ref, .false., test_coarsening, &
                write_output, trim(output_name), n)
        end do
     end do
#:elif NDIM == 3
     do k = -1, 1
        do j = -1, 1
           do i = -1, 1
              r_ref = [0.5_dp + i * 0.49_dp, 0.5_dp + j * 0.49_dp, &
                   0.5_dp + k * 0.49_dp]
              if (f4%mpirank == 0) print *, "Refine around", r_ref
              call test_refinement(f4, n_gc, min_level, num_refine_steps, &
                   r_ref, .false., test_coarsening, &
                   write_output, trim(output_name), n)
           end do
        end do
     end do
#:endif
  end do

  call f4_print_wtime(f4)
  call f4_finalize(f4)

contains

  subroutine test_refinement(f4, n_gc, min_level, n_refine_steps, &
       refine_location, refine_everywhere, test_coarsening, &
       write_output, base_name, n_output)
    type(foap4_t), intent(inout) :: f4
    integer, intent(in)          :: n_gc
    integer, intent(in)          :: min_level
    integer, intent(in)          :: n_refine_steps
    real(dp), intent(in)         :: refine_location(NDIM)
    logical, intent(in)          :: refine_everywhere
    logical, intent(in)          :: test_coarsening
    logical, intent(in)          :: write_output
    character(len=*), intent(in) :: base_name
    integer, intent(inout)       :: n_output
    integer, parameter           :: n_blocks_per_dim(NDIM) = 1
    real(dp), parameter          :: block_length(NDIM)     = 1.0_dp
    integer, parameter           :: n_vars                 = 2
    character(len=20)            :: var_names(n_vars)      = ['err', 'phi']
    logical, parameter           :: temporal(n_vars)       = [.false., .false.]
    integer, parameter           :: n_temporal_states      = 1
    logical, parameter           :: periodic(NDIM)         = .false.
    integer                      :: n, prev_mesh_revision

    f4%bc_callback => bc_callback

    call f4_construct_brick(f4, n_blocks_per_dim, block_length, bx, n_gc, &
         n_vars, var_names, temporal, n_temporal_states, periodic, &
         min_level, max_blocks, f4_bc_dirichlet, 0.0_dp)

    call set_init_cond(f4)
    call f4_update_ghostcells(f4, 1, [i_phi], 0)
    call compute_error(f4)

    if (write_output) then
       n_output = n_output + 1
       call io_write_grid(f4, base_name, n_output, viewer=viewer)
    end if

    do n = 1, n_refine_steps
       call set_refinement_flag(f4, refine_location, refine_everywhere)
       call f4_adjust_refinement(f4, load_imbalance_threshold)
       call f4_update_ghostcells(f4, 1, [i_phi], 0)
       call compute_error(f4)

       if (write_output) then
          n_output = n_output + 1
          call io_write_grid(f4, base_name, n_output, viewer=viewer)
       end if
    end do

    if (test_coarsening) then
       if (f4%mpirank == 0) print *, "Test coarsening"
       do n = 1, 10
          prev_mesh_revision = f4_get_mesh_revision(f4)
          call set_coarsening_flag(f4)
          call f4_adjust_refinement(f4, load_imbalance_threshold)

          call f4_update_ghostcells(f4, 1, [i_phi], 0)
          call compute_error(f4)

          if (write_output) then
             n_output = n_output + 1
             call io_write_grid(f4, base_name, n_output, viewer=viewer)
          end if

          if (f4_get_mesh_revision(f4) == prev_mesh_revision) exit
       end do
    end if

    call f4_destroy(f4)
  end subroutine test_refinement

#:if NDIM == 2
  pure real(dp) function phi_init(x, y)
    ${ROUTINE_SEQ()}$
    real(dp), intent(in) :: x, y
    phi_init = x + 2*y
  end function phi_init
#:elif NDIM == 3
  pure real(dp) function phi_init(x, y, z)
    ${ROUTINE_SEQ()}$
    real(dp), intent(in) :: x, y, z
    phi_init = x + 2*y + 3*z
  end function phi_init
#:endif

  subroutine bc_callback(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: face, i_block, ix, n, i
    real(dp)                     :: rr(NDIM)
#:if NDIM == 3
    integer                      :: j
#:endif

    ${PARALLEL()}$
    do face = 0, 2*NDIM-1
       ${LOOP('collapse(NDIM) private(i_block, ix, rr)')}$
       do n = f4%gc_phys_iface(face), f4%gc_phys_iface(face+1)-1
#:if NDIM == 2
          do i = 1, f4%bx(1)
             i_block = f4%gc_phys(n) + 1
             f4%bc_data_ix(face, i_block) = abs(f4%bc_data_ix(face, i_block))
             ix = f4%bc_data_ix(face, i_block)
             rr = f4_block_face_coord(f4, i_block, face, i)
             f4%bc_data(i, i_phi, ix) = real(phi_init(rr(1), rr(2)), fp)
             f4%bc_data_type(i, i_phi, ix) = f4_bc_dirichlet
          end do
#:elif NDIM == 3
          do j = 1, f4%bx(1)
             do i = 1, f4%bx(1)
                i_block = f4%gc_phys(n) + 1
                f4%bc_data_ix(face, i_block) = abs(f4%bc_data_ix(face, i_block))
                ix = f4%bc_data_ix(face, i_block)
                rr = f4_block_face_coord(f4, i_block, face, i, j)
                f4%bc_data(i, j, i_phi, ix) = real(phi_init(rr(1), rr(2), rr(3)), fp)
                f4%bc_data_type(i, j, i_phi, ix) = f4_bc_dirichlet
             end do
          end do
#:endif
       end do
    end do
    ${END_PARALLEL()}$
  end subroutine bc_callback

  subroutine set_init_cond(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, ${IJK}$
    real(dp)                     :: rr(NDIM)

    ${PARALLEL_LOOP('collapse(NDIM+1) private(rr)')}$ ${DEFAULT_PRESENT()}$
    do n = 1, f4%n_blocks
       do @{KJI_LOOP_1_to_array(f4%bx)}@
          rr = f4_cell_coord(f4, n, ${IJK}$)
          f4%uu(${IJK}$, i_phi, n) = real(phi_init(@{DINDEX(rr)}@), fp)
          f4%uu(${IJK}$, i_err, n) = 0.0_dp
       end do; ${KJI_CLOSE_LOOP}$
    end do
  end subroutine set_init_cond

  subroutine set_refinement_flag(f4, r0, refine_everywhere)
    type(foap4_t), intent(inout) :: f4
    real(dp), intent(in)         :: r0(NDIM)
    logical, intent(in)          :: refine_everywhere
    integer                      :: n, lvl
    real(dp)                     :: rmin(NDIM), rmax(NDIM)

    do n = 1, f4%n_blocks
       lvl = f4%block_level(n)
       rmin = f4%block_origin(:, n)
       rmax = rmin + f4%bx * f4%dr_level(:, lvl)

       if (all(r0 >= rmin .and. r0 <= rmax) .or. refine_everywhere) then
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

  subroutine compute_error(f4)
    type(foap4_t), intent(inout) :: f4
    integer                      :: n, ${IJK}$
    real(dp)                     :: rr(NDIM), sol
    logical                      :: ghost_dim(NDIM), valid_cell
    real(dp), parameter          :: max_difference = 1e-15_dp
    real(dp)                     :: max_err

    max_err = 0.0_dp

    ${PARALLEL_LOOP('collapse(NDIM+1) reduction(max:max_err) private(rr, sol, ghost_dim, valid_cell)')}$ ${DEFAULT_PRESENT()}$
    do n = 1, f4%n_blocks
       do @{KJI_LOOP_array_to_array(f4%ilo, f4%ihi)}@

          ghost_dim(1) = i < 1 .or. i > f4%bx(1)
          ghost_dim(2) = j < 1 .or. j > f4%bx(2)
#:if NDIM == 3
          ghost_dim(3) = k < 1 .or. k > f4%bx(3)
#:endif

#:if NDIM == 2
          valid_cell = .not. (ghost_dim(1) .and. ghost_dim(2))
#:elif NDIM == 3
          valid_cell = .not. ( &
               (ghost_dim(1) .and. ghost_dim(2)) .or. &
               (ghost_dim(1) .and. ghost_dim(3)) .or. &
               (ghost_dim(2) .and. ghost_dim(3)))
#:endif

          if (valid_cell) then
             rr = f4_cell_coord(f4, n, ${IJK}$)
             sol = phi_init(@{DINDEX(rr)}@)
             f4%uu(${IJK}$, i_err, n) = abs(f4%uu(${IJK}$, i_phi, n) - real(sol, fp))
          else
             f4%uu(${IJK}$, i_err, n) = 0.0_dp
          end if

          max_err = max(max_err, f4%uu(${IJK}$, i_err, n))

       end do; ${KJI_CLOSE_LOOP}$
    end do

    if (max_err > max_difference) then
       print *, f4%mpirank, "max error: ", max_err
       if (abort_on_error) then
          error stop "Too large error"
       end if
    end if

  end subroutine compute_error

end program test_ref
