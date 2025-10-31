! driver program.  Reads in the integral parameters, calls lnln_pole_int, and prints the
! result.
program main
  use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
  implicit none
  complex :: lnln_pole_int
  external lnln_pole_int
  real :: klb, kub, ktemp
  
  complex :: a,c,f

  integer :: num_args, ix

  character(len=42) :: arg
  real, dimension(8) :: args

  num_args = command_argument_count()
  if(num_args.ne.8) then
    write(stderr,*) 'error; wrong number of arguments'
    return
  endif

  do ix=1,num_args
    call get_command_argument(ix,arg)
    read(arg,*) args(ix)
  enddo

  a = cmplx(args(1), args(2))
  c = cmplx(args(3), args(4))
  f = cmplx(args(5), args(6))
  klb = args(7)
  kub = args(8)

  write(*,*) 'a = ',a
  write(*,*) 'c = ',c
  write(*,*) 'f = ',f
  write(*,*) 'klb = ',klb
  write(*,*) 'kub = ',kub
  write(*,*) 'lnln_pole_int = ',lnln_pole_int(klb,kub,a,c,f)
end program main
