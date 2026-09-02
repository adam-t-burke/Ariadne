! SPDX-License-Identifier: BSD-3-Clause
! Timing harness for the unmodified official L-BFGS-B 3.0 routines.
! The implementation is downloaded and hash-verified by benchmark.ps1.
program benchmark_driver
  use iso_fortran_env, only: int64
  implicit none
  integer :: n, m, repeats, run, evaluations, iterations
  integer(int64) :: started, finished, clock_rate
  integer, allocatable :: nbd(:), iwa(:), isave(:)
  double precision :: f
  double precision, allocatable :: x(:), l(:), u(:), g(:), dsave(:), wa(:)
  logical :: lsave(4)
  character(len=60) :: task, csave
  character(len=32) :: argument

  call get_command_argument(1, argument)
  read(argument, *) n
  call get_command_argument(2, argument)
  read(argument, *) m
  call get_command_argument(3, argument)
  read(argument, *) repeats

  allocate(x(n), l(n), u(n), g(n), nbd(n), iwa(3*n), isave(44), dsave(29))
  allocate(wa(2*m*n + 5*n + 11*m*m + 8*m))

  call system_clock(started, clock_rate)
  do run = 1, repeats
    call initialize()
    task = 'START'
    do
      call setulb(n, m, x, l, u, nbd, f, g, 1.0d7, 1.0d-5, wa, iwa, &
                  task, -1, csave, lsave, isave, dsave)
      if (task(1:2) == 'FG') then
        call objective_gradient()
      else if (task(1:5) /= 'NEW_X') then
        exit
      end if
    end do
    iterations = isave(30)
    evaluations = isave(34)
  end do
  call system_clock(finished)

  write(*, '(A,I0,A,I0,A,I0,A,F0.3,A,I0,A,I0)') &
      'official-fortran,n=', n, ',m=', m, ',runs=', repeats, &
      ',microseconds=', real(finished-started, kind=8)*1.0d6/real(clock_rate, kind=8)/repeats, &
      ',iterations=', iterations, ',evaluations=', evaluations

contains
  subroutine initialize()
    integer :: i
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
  end subroutine initialize

  subroutine objective_gradient()
    integer :: i
    double precision :: t1, t2
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
  end subroutine objective_gradient
end program benchmark_driver
