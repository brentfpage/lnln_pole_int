module k_crossings_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  contains

  ! Compute the k points listed in eqs. (15-18) of lnln_pole_int_doc.pdf if they are
  ! relevant. Such points that equal -a,-c, or -f (which can only be the case for pure real
  ! a,c, or f) are irrelevant because integrals with bounds spanning them are undefined. To
  ! streamline the program, this function usually returns NaN for a k_crossing that is actually
  ! -a,-c, or -f. 
  ! Also, determine if g(k) lies on (-∞,0] or if g(k) lies on [1,∞).  Here, g(k) is defined in
  ! Section 3 of lnln_pole_int_doc.pdf
  ! k_crossings : the k points listed in eqs. (15-18) of lnln_pole_int_doc.pdf
  ! discont_flag_01: whether g(k) lies on (-∞,0] AND klb<k<kub
  ! discont_flag_02: whether g(k) lies on [1,∞) AND klb<k<kub
subroutine get_k_crossings(klb, kub, a, c, f, k_crossings, discont_flag_01, discont_flag_02, acf_colinear)
  use :: utils_mod
  implicit none
  real :: klb, kub
  complex :: a, c, f
  complex :: g_k_

  real, dimension(5) :: k_crossings
  logical, dimension(5) :: discont_flag_01
  logical, dimension(5) :: discont_flag_02

  logical :: acf_colinear

  acf_colinear = acf_are_colinear(a,c,f)

  if(&
    ((aimag(a).eq.0).and.(aimag(c).eq.0)).or.&
    ((aimag(a).eq.0).and.(aimag(f).eq.0)).or.&
    ((aimag(c).eq.0).and.(aimag(f).eq.0))&
    ) then
    ! all k_crossings either nonexistent or irrelevant
    k_crossings = ieee_value(k_crossings, ieee_quiet_nan)
    discont_flag_01 = .false.
    discont_flag_02 = .false.
  else
    ! g(k) = 1-y
    call get1_k_crossing(klb, kub, a, f, k_crossings(1), discont_flag_01(1), discont_flag_02(1))

    ! g(k) = 1-xiy
    call get1_k_crossing(klb, kub, c, f, k_crossings(2), discont_flag_01(2), discont_flag_02(2))

    if((aimag(a-c).eq.0).or.(aimag(a).eq.0).or.(aimag(c).eq.0)) then
        k_crossings(3) = ieee_value(k_crossings(3), ieee_quiet_nan)
        discont_flag_01(3) = .false.
        discont_flag_02(3) = .false.
    else
      ! g(k) = 1-u
      k_crossings(3) = aimag(a*conjg(c))/aimag(c-a)
      g_k_ = (c+k_crossings(3))/(a+k_crossings(3))
      discont_flag_01(3) = (klb.lt.k_crossings(3)).and.(kub.gt.k_crossings(3)).and.(real(g_k_).lt.0.0)
      discont_flag_02(3) = (klb.lt.k_crossings(3)).and.(kub.gt.k_crossings(3)).and.(real(g_k_).gt.1.0)
    endif

     ! g(k) = z
    call quadr_sol_k_crossings(klb,kub,a,c,f,k_crossings(4:5),discont_flag_01(4:5),discont_flag_02(4:5), acf_colinear)

  endif

!   call print_k_crossing_info(k_crossings,discont_flag_01,discont_flag_02,acf_colinear)
end subroutine get_k_crossings


! for g(k) = 1-(f+k)/(f-a_or_c), determine the k_c at which Im(g(k_c)) = 0
! Case 1
!   a_or_c = a
!   g(k) = 1-y
! Case 2
!   a_or_c = c
!   g(k) = 1-ξy
! discont_flag_01 : whether klb<k_c<kub AND g(k_c)∈(-∞,0]
! discont_flag_02 : whether klb<k_c<kub AND g(k_c)∈[1,∞)
subroutine get1_k_crossing(klb, kub, a_or_c, f, k_crossing, discont_flag_01, discont_flag_02)
  implicit none
  real :: klb, kub
  complex :: a_or_c, f
  complex :: g_k_

  real :: k_crossing
  logical :: discont_flag_01
  logical :: discont_flag_02

  if((aimag(f-a_or_c).eq.0).or.(aimag(f).eq.0).or.(aimag(a_or_c).eq.0)) then
    k_crossing = ieee_value(k_crossing, ieee_quiet_nan)
    discont_flag_01 = .false.
    discont_flag_02 = .false.
  else 
    k_crossing = aimag(a_or_c*conjg(f))/aimag(f-a_or_c)
    g_k_ = -(a_or_c+k_crossing)/(f-a_or_c) 
    discont_flag_01 = (klb.lt.k_crossing).and.(kub.gt.k_crossing).and.(real(g_k_).lt.0.0)
    discont_flag_02 = (klb.lt.k_crossing).and.(kub.gt.k_crossing).and.(real(g_k_).gt.1.0)
  endif

end subroutine get1_k_crossing

! compute the k points determined by eq. (18) of lnln_pole_int_doc.pdf
! Edge cases are also described in that document.
! Also, determine if z(k) lies on (-∞,0] or if z(k) lies on [1,∞).  Here, z = (1-ξy)/(1-y), where y=(f+k)/(f-a) and ξ=(f-a)/(f-c)
! k_crossings : the k points listed in eq. (18) of lnln_pole_int_doc.pdf
! discont_flag_01: whether z(k) lies on (-∞,0] AND klb<k<kub
! discont_flag_02: whether z(k) lies on [1,∞) AND klb<k<kub
! acf_colinear : whether a, c, and f lie in a straight line in the complex plane

  subroutine quadr_sol_k_crossings(klb, kub, a, c, f, k_crossings, discont_flag_01, discont_flag_02, acf_colinear)
    implicit none
    real :: klb, kub
    complex :: a, c, f
    real, dimension(3) :: coeffs
    real :: discriminant_squared, discriminant
    complex, dimension(2) :: z
    real, dimension(2) :: k_crossings
    logical, dimension(2) :: discont_flag_01
    logical, dimension(2) :: discont_flag_02
    integer :: ic
    logical :: acf_colinear

    coeffs(1) = aimag(conjg(c-f)*(a-f)) ! given exact precision, coeffs(1).eq.0 is the same as acf_colinear

    if((aimag(a).ne.0).and.(aimag(c).ne.0).and.(aimag(f).ne.0)) then
      coeffs(2) = aimag(conjg(c-f)*(a-f)*(conjg(a)+c))
      coeffs(3) = aimag(conjg(c-f)*(a-f)*c*conjg(a))
      discriminant_squared = (coeffs(2)**2 - 4*coeffs(1)*coeffs(3))
      
      if( (discriminant_squared.gt.0).and.(.not.acf_colinear) ) then
        discriminant  = sqrt(discriminant_squared)
        k_crossings(1) = (-coeffs(2) + discriminant)/(2*coeffs(1))
        k_crossings(2) = (-coeffs(2) - discriminant)/(2*coeffs(1))
      else if ((aimag(c-a).ne.0.0).and.(aimag(c).ne.0).and.(aimag(a).ne.0).and.acf_colinear) then
        k_crossings(1) = aimag(a*conjg(c))/(aimag(c-a));
        k_crossings(2) = ieee_value(k_crossings(2), ieee_quiet_nan)
      else
        k_crossings = ieee_value(k_crossings, ieee_quiet_nan)
      endif
    else if(acf_colinear) then
      k_crossings = ieee_value(k_crossings, ieee_quiet_nan)
    else
      ! ignore the k_crossing given by -f,-c, or -a
      k_crossings(2) = ieee_value(k_crossings(2), ieee_quiet_nan)

      if ( (aimag(a).ne.0).and.(aimag(c).ne.0) ) then
        k_crossings(1) = -aimag(a*(c-f)*conjg(c-a))/aimag((c-f)*conjg(c-a))
      else if ( (aimag(a).ne.0).and.(aimag(f).ne.0) ) then
        k_crossings(1) = aimag(a*(f-c)*conjg(f-a))/coeffs(1)
      else if ( (aimag(c).ne.0).and.(aimag(f).ne.0) ) then
        k_crossings(1) = -aimag(c*(f-a)*conjg(f-c))/coeffs(1)
      else 
      ! two or more of a, c and f are pure real, so z either crosses the real axis nowhere or
      ! exclusively at two of -a, -c, and -f.  In practice, this case has already been handled
      ! by get_k_crossings
        k_crossings(1) = ieee_value(k_crossings(1), ieee_quiet_nan)
      endif
    endif

    ! it's guaranteed that discont_flag_0(1/2)=.false. for the k_crossings=-a case
    ! because integrals with bounds that span -a (with 'a' real) are undefined
    where((.not.isnan(k_crossings)).and.(k_crossings.ne.-a))
      z = (c+k_crossings)*(a-f)/(c-f)/(a+k_crossings)
      discont_flag_01 = (klb.lt.k_crossings).and.(kub.gt.k_crossings).and.(real(z).lt.0.0)
      discont_flag_02 = (klb.lt.k_crossings).and.(kub.gt.k_crossings).and.(real(z).gt.1.0)
    elsewhere
      discont_flag_01 = .false.
      discont_flag_02 = .false.
    end where

  end subroutine quadr_sol_k_crossings

  subroutine print_k_crossing_info(k_crossings, discont_flag_01, discont_flag_02, acf_colinear)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    real, dimension(5) :: k_crossings
    logical, dimension(5) :: discont_flag_01
    logical, dimension(5) :: discont_flag_02

    logical :: acf_colinear

    write(stderr,*) '---------------------'
    write(stderr,*) '---------------------'
    write(stderr,*) 'a, c, and f colinear?',acf_colinear
    write(stderr,*) '---------------------'
    write(stderr,*) 'Im(y) = 0'
    write(stderr,*) 'kc, 1-val|_kc < 0, val|_kc < 0', k_crossings(1), discont_flag_01(1), discont_flag_02(1)
    write(stderr,*) '---------------------'
    write(stderr,*) 'Im(xi*y) = 0'
    write(stderr,*) 'kc, 1-val|_kc < 0, val|_kc < 0', k_crossings(2), discont_flag_01(2), discont_flag_02(2)
    write(stderr,*) '---------------------'
    write(stderr,*) 'Im(u) = 0'
    write(stderr,*) 'k_c, 1-val|_kc < 0, val|_kc < 0', k_crossings(3), discont_flag_01(3), discont_flag_02(3)
    write(stderr,*) '---------------------'
    write(stderr,*) 'Im(z) = 0'
    write(stderr,*) 'kc, 1-val|_kc < 0, val|_kc < 0', k_crossings(4), discont_flag_02(4), discont_flag_01(4)
    write(stderr,*) 'kc, 1-val|_kc < 0, val|_kc < 0', k_crossings(5), discont_flag_02(5), discont_flag_01(5)
    write(stderr,*) '---------------------'
    write(stderr,*) '---------------------'
  end subroutine

end module
