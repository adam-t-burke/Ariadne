! SPDX-License-Identifier: BSD-3-Clause
! Fixture harness derived from the official L-BFGS-B 3.0 sample drivers.
! The upstream routines are downloaded by regenerate.ps1 and are not linked
! into the production Rust crate. See PROVENANCE.md and THIRD_PARTY_NOTICES.md.
program trace_generator
  implicit none
  character(len=1024) :: output_dir

  call get_command_argument(1, output_dir)
  if (len_trim(output_dir) == 0) stop 'usage: trace_generator OUTPUT_DIRECTORY'

  call run_case('driver1', trim(output_dir), 25, 5, 1.0d7, 1.0d-5, 1, .true.)
  call run_case('driver2', trim(output_dir), 25, 5, 0.0d0, 0.0d0, 2, .true.)
  call run_case('driver3-large-n', trim(output_dir), 1000, 10, 0.0d0, 0.0d0, 3, .false.)
  call run_case('edge-mixed-bounds', trim(output_dir), 4, 5, 1.0d7, 1.0d-12, 4, .true.)
  call run_case('audit-mixed-rollover', trim(output_dir), 8, 2, 1.0d7, 1.0d-10, 5, .true.)
  call run_dcsrch_cases(trim(output_dir))

contains

  subroutine run_dcsrch_cases(directory)
    character(len=*), intent(in) :: directory
    character(len=60) :: task
    character(len=2048) :: path
    integer :: unit, isave(2)
    double precision :: f, g, stp, dsave(13)

    path = trim(directory) // '/dcsrch-cases.csv'
    open(newunit=unit, file=trim(path), status='replace', action='write')
    write(unit, '(A)') 'case,step,task'

    f = 0.0d0
    g = -1.0d0
    stp = 1.0d0
    task = 'START'
    call dcsrch(f, g, stp, 1.0d-3, 0.9d0, 0.1d0, 0.0d0, 1.0d0, task, isave, dsave)
    f = -2.0d0
    g = -1.0d0
    call dcsrch(f, g, stp, 1.0d-3, 0.9d0, 0.1d0, 0.0d0, 1.0d0, task, isave, dsave)
    call write_dcsrch_record(unit, 'step-at-maximum', stp, task)

    f = 0.0d0
    g = -1.0d0
    stp = 0.0d0
    task = 'START'
    call dcsrch(f, g, stp, 1.0d-3, 0.9d0, 0.1d0, 0.0d0, 1.0d0, task, isave, dsave)
    f = 1.0d0
    g = 0.0d0
    call dcsrch(f, g, stp, 1.0d-3, 0.9d0, 0.1d0, 0.0d0, 1.0d0, task, isave, dsave)
    call write_dcsrch_record(unit, 'step-at-minimum', stp, task)

    f = 0.0d0
    g = -1.0d0
    stp = 1.0d0
    task = 'START'
    call dcsrch(f, g, stp, 1.0d-3, 0.9d0, 0.1d0, 0.0d0, 1.0d0, task, isave, dsave)
    f = -0.5d0
    g = 0.0d0
    call dcsrch(f, g, stp, 1.0d-3, 0.9d0, 0.1d0, 0.0d0, 1.0d0, task, isave, dsave)
    call write_dcsrch_record(unit, 'converged', stp, task)

    close(unit)
  end subroutine run_dcsrch_cases

  subroutine write_dcsrch_record(unit, name, step, task)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: name, task
    double precision, intent(in) :: step
    write(unit, '(A,",",ES25.17E3,",",A)') trim(name), step, trim(task)
  end subroutine write_dcsrch_record

  subroutine run_case(name, directory, n, m, factr, pgtol, mode, full_states)
    character(len=*), intent(in) :: name, directory
    integer, intent(in) :: n, m, mode
    double precision, intent(in) :: factr, pgtol
    logical, intent(in) :: full_states
    character(len=60) :: task, csave
    character(len=2048) :: path
    logical :: lsave(4), wrote_initial
    integer :: i, unit
    integer, allocatable :: nbd(:), iwa(:), isave(:)
    double precision :: f
    double precision, allocatable :: x(:), l(:), u(:), g(:), dsave(:), wa(:)

    allocate(x(n), l(n), u(n), g(n), nbd(n), iwa(3*n), isave(44), dsave(29))
    allocate(wa(2*m*n + 5*n + 11*m*m + 8*m))

    call initialize_problem(mode, n, x, l, u, nbd)
    task = 'START'
    wrote_initial = .false.
    path = trim(directory) // '/' // trim(name) // '.csv'
    open(newunit=unit, file=trim(path), status='replace', action='write')
    write(unit, '(A)') 'record,iteration,evaluations,f,projected_gradient,x,g,task'

    do
      call setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, &
                  task, -1, csave, lsave, isave, dsave)

      if (task(1:2) == 'FG') then
        call objective_gradient(mode, n, x, f, g)
        if (.not. wrote_initial) then
          call write_record(unit, 'initial', 0, 1, f, &
                            projected_gradient(n, x, g, l, u, nbd), x, g, task, .true.)
          wrote_initial = .true.
        end if
      else if (task(1:5) == 'NEW_X') then
        call write_record(unit, 'accepted', isave(30), isave(34), f, dsave(13), &
                          x, g, task, full_states)
        if (mode == 2) then
          if (isave(34) >= 99) task = 'STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
          if (dsave(13) <= 1.0d-10 * (1.0d0 + abs(f))) &
            task = 'STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'
        else if (mode == 3) then
          if (isave(34) >= 900) task = 'STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
          if (dsave(13) <= 1.0d-10 * (1.0d0 + abs(f))) &
            task = 'STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'
        end if
      else
        call write_record(unit, 'terminal', isave(30), isave(34), f, dsave(13), &
                          x, g, task, .true.)
        exit
      end if
    end do

    close(unit)
    deallocate(x, l, u, g, nbd, iwa, isave, dsave, wa)
  end subroutine run_case

  subroutine initialize_problem(mode, n, x, l, u, nbd)
    integer, intent(in) :: mode, n
    double precision, intent(out) :: x(n), l(n), u(n)
    integer, intent(out) :: nbd(n)
    integer :: i

    if (mode <= 3) then
      do i = 1, n
        x(i) = 3.0d0
        nbd(i) = 2
        if (mod(i, 2) == 1) then
          l(i) = 1.0d0
          u(i) = 1.0d2
        else
          l(i) = -1.0d2
          u(i) = 1.0d2
        end if
      end do
    else if (mode == 4) then
      ! Independent projection/mixed-bound fixture, not an upstream driver.
      x = (/ -3.0d0, 9.0d0, 7.0d0, -8.0d0 /)
      l = (/ 0.0d0, -1.0d0, 0.0d0, 0.0d0 /)
      u = (/ 2.0d0, 4.0d0, 0.0d0, 0.0d0 /)
      nbd = (/ 2, 2, 0, 1 /)
    else
      ! Audit fixture: infeasible start plus every nbd class and one fixed
      ! variable. Variables 6:8 have equal first-iteration breakpoints t=1.
      x = (/ -3.0d0, 8.0d0, 4.0d0, -6.0d0, 9.0d0, 0.0d0, 0.0d0, 0.0d0 /)
      l = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.5d0, -1.0d0, -0.5d0, 0.0d0 /)
      u = (/ 2.0d0, 0.0d0, 2.0d0, 0.0d0, 1.5d0, 1.0d0, 0.0d0, 0.5d0 /)
      nbd = (/ 2, 1, 3, 0, 2, 2, 1, 3 /)
    end if
  end subroutine initialize_problem

  subroutine objective_gradient(mode, n, x, f, g)
    integer, intent(in) :: mode, n
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: f, g(n)
    double precision :: t1, t2
    double precision, parameter :: target(4) = (/ 1.0d0, 2.0d0, -2.0d0, 3.0d0 /)
    double precision, parameter :: audit_target(8) = &
      (/ 1.0d0, -2.0d0, 4.0d0, -1.0d0, 9.0d0, 2.0d0, -1.0d0, 1.0d0 /)
    double precision, parameter :: audit_scale(8) = &
      (/ 1.0d0, 50.0d0, 0.1d0, 20.0d0, 1.0d0, 0.5d0, 0.5d0, 0.5d0 /)
    integer :: i

    if (mode <= 3) then
      ! Extended Rosenbrock objective and gradient from official drivers 1-3.
      f = 0.25d0 * (x(1) - 1.0d0)**2
      do i = 2, n
        f = f + (x(i) - x(i-1)**2)**2
      end do
      f = 4.0d0 * f

      t1 = x(2) - x(1)**2
      g(1) = 2.0d0 * (x(1) - 1.0d0) - 16.0d0 * x(1) * t1
      do i = 2, n - 1
        t2 = t1
        t1 = x(i+1) - x(i)**2
        g(i) = 8.0d0 * t2 - 16.0d0 * x(i) * t1
      end do
      g(n) = 8.0d0 * t1
    else if (mode == 4) then
      f = 0.0d0
      do i = 1, n
        f = f + (x(i) - target(i))**2
        g(i) = 2.0d0 * (x(i) - target(i))
      end do
    else
      f = 0.0d0
      do i = 1, n
        f = f + 0.5d0 * audit_scale(i) * (x(i) - audit_target(i))**2
        g(i) = audit_scale(i) * (x(i) - audit_target(i))
      end do
    end if
  end subroutine objective_gradient

  double precision function projected_gradient(n, x, g, l, u, nbd)
    integer, intent(in) :: n, nbd(n)
    double precision, intent(in) :: x(n), g(n), l(n), u(n)
    double precision :: component
    integer :: i

    projected_gradient = 0.0d0
    do i = 1, n
      component = g(i)
      if ((nbd(i) == 1 .or. nbd(i) == 2) .and. x(i) <= l(i)) &
        component = min(0.0d0, component)
      if ((nbd(i) == 2 .or. nbd(i) == 3) .and. x(i) >= u(i)) &
        component = max(0.0d0, component)
      projected_gradient = max(projected_gradient, abs(component))
    end do
  end function projected_gradient

  subroutine write_record(unit, kind, iteration, evaluations, f, pg, x, g, task, include_vectors)
    integer, intent(in) :: unit, iteration, evaluations
    character(len=*), intent(in) :: kind, task
    double precision, intent(in) :: f, pg, x(:), g(:)
    logical, intent(in) :: include_vectors

    write(unit, '(A,",",I0,",",I0,",",ES25.17E3,",",ES25.17E3,",")', advance='no') &
      trim(kind), iteration, evaluations, f, pg
    if (include_vectors) then
      call write_vector(unit, x)
      write(unit, '(A)', advance='no') ','
      call write_vector(unit, g)
    else
      write(unit, '(A)', advance='no') ','
    end if
    write(unit, '(A,A)') ',', trim(task)
  end subroutine write_record

  subroutine write_vector(unit, values)
    integer, intent(in) :: unit
    double precision, intent(in) :: values(:)
    integer :: i

    write(unit, '(A)', advance='no') '"'
    do i = 1, size(values)
      if (i > 1) write(unit, '(A)', advance='no') '|'
      write(unit, '(ES25.17E3)', advance='no') values(i)
    end do
    write(unit, '(A)', advance='no') '"'
  end subroutine write_vector

end program trace_generator
