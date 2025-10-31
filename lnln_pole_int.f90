! compute the following integral:
!  kub
! ⌠    ln(a+k)ln(c+k)
! ⎮    ────────────── dk
! ⌡         f+k
!  klb
!  using the methods outlined in lnln_pole_int_doc.pdf

function lnln_pole_int(klb,kub,a,c,f)
  use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use :: di_tri_log_mod, only : clog
  use utils_mod
  use ln_terms_mod
  use dilog_terms_mod
  use lewin_f_terms_mod
  use k_crossings_mod
  use N_type_branch_cut_corrs_mod
  implicit none
  real :: klb, kub
  complex :: a, c, f
  complex :: lnln_pole_int
  complex :: lnln_pole_int_a_eq_c

  external lnln_pole_int_a_eq_c

  real, dimension(5) :: k_crossings
  logical, dimension(5) :: discont_flag_01
  logical, dimension(5) :: discont_flag_02

  complex, dimension(2) :: y  
  complex, dimension(2) :: xiy
  complex, dimension(2) :: u
  complex, dimension(2) :: z

  logical :: acf_colinear

  integer :: int_order_sign
  real :: ktemp

  complex :: xi, val1
  logical, dimension(2) :: in_interval
  real, dimension(2) :: z_on_real_axis
  real :: u_on_real_axis
  
  int_order_sign=1
  if(klb.gt.kub) then
    ktemp = klb
    klb=kub
    kub=ktemp
    int_order_sign = -1
  endif

  if((klb.eq.kub)) then
    lnln_pole_int = (0.0, 0.0)
  else if(&
    ((aimag(a).eq.0).and.(klb.le.real(-a)).and.(kub.ge.real(-a))).or.&
    ((aimag(c).eq.0).and.(klb.le.real(-c)).and.(kub.ge.real(-c)))&
    ) then
    write(stderr,*) 'integral is undefined'
    lnln_pole_int = cmplx(ieee_value(ktemp, ieee_quiet_nan) , ieee_value(ktemp, ieee_quiet_nan))
  else if((aimag(f).eq.0).and.(klb.le.real(-f)).and.(kub.ge.real(-f))) then
    if((a.eq.(f+1)).or.(c.eq.(f+1))) then
      write(stderr,*) 'integral is unsupported - see Readme.txt'
    else
      write(stderr,*) 'integral is undefined'
    endif
    lnln_pole_int = cmplx(ieee_value(ktemp, ieee_quiet_nan) , ieee_value(ktemp, ieee_quiet_nan))
  else if ((a.eq.c).and.(c.eq.f)) then
    lnln_pole_int = (clog(kub+a)**3 - clog(klb+a)**3)/3.
  else if (a.eq.c) then
    lnln_pole_int = lnln_pole_int_a_eq_c(klb,kub,a,f)
  else if (a.eq.f) then
    ! from integration by parts
    lnln_pole_int = (clog(kub+a)**2*clog(kub+c) - clog(klb+a)**2*clog(klb+c))/2. - &
      lnln_pole_int_a_eq_c(klb,kub,a,c)/2.
  else if (c.eq.f) then
    ! from integration by parts
    lnln_pole_int = (clog(kub+c)**2*clog(kub+a) - clog(klb+c)**2*clog(klb+a))/2. - &
      lnln_pole_int_a_eq_c(klb,kub,c,a)/2.
  else
    call get_k_crossings(klb, kub, a,c,f, k_crossings, discont_flag_01, discont_flag_02, acf_colinear)
    if(acf_colinear) then
      xi = real((a-f)/(c-f))
      val1 = real((c-a)/(c-f))
    else
      xi = (a-f)/(c-f)
      val1 = (c-a)/(c-f)
    endif

    y = (f + [klb, kub])/(f-a)
    xiy = (f + [klb, kub])/(f-c)
    u = -(c-a)/([klb, kub] + a)
    z = xi*(c+[klb, kub])/(a+[klb, kub])


    lnln_pole_int = -clog(a-f) * dilog_term_corr_A(xiy, discont_flag_01(2), real((f+k_crossings(2))/(f-c)))

    lnln_pole_int = lnln_pole_int - clog(c-f) * &
      dilog_term_corr_A(y, discont_flag_01(1) , real((f+k_crossings(1))/(f-a)))

    lnln_pole_int = lnln_pole_int + clog(c-f) * clog(a - f) * ln_term_corr_A(y, discont_flag_02(1))
    lnln_pole_int = lnln_pole_int + clog(xi)**2 * ln_term_corr_A(u, discont_flag_02(3))/2.

    if( (.not.isnan(k_crossings(3)))) then
      u_on_real_axis = real((a-c)/(k_crossings(3)+a))
    else
      u_on_real_axis = ieee_value(u_on_real_axis, ieee_quiet_nan)
    endif

    lnln_pole_int = lnln_pole_int - clog(xi) * &
        dilog_term_corr_A(u, discont_flag_01(3),u_on_real_axis)

    lnln_pole_int = lnln_pole_int + &
      lewin_f_term_corr_A(y, &
        real((f + k_crossings(1))/(f-a)), &
        discont_flag_01(1))/2.

    lnln_pole_int = lnln_pole_int + &
      lewin_f_term_corr_A(xiy(:2), &
        real((f + k_crossings(2))/(f-c)), &
        discont_flag_01(2))/2.

    lnln_pole_int = lnln_pole_int + &
      lewin_f_term_corr_A(u, &
        u_on_real_axis, &
        discont_flag_01(3))/2.

    where(.not.isnan(k_crossings(4:5)))
      in_interval=((k_crossings(4:5).gt.klb).and.(k_crossings(4:5).lt.kub))
      z_on_real_axis = real(xi*(c+k_crossings(4:5))/(a+k_crossings(4:5)))
    elsewhere
      in_interval=.false.
      z_on_real_axis = ieee_value(z_on_real_axis,ieee_quiet_nan)
    endwhere

    lnln_pole_int = lnln_pole_int - &
      lewin_f_term_corr_A(1-z, &
        1 - z_on_real_axis, &
        argsort(k_crossings(4:5)),&
        in_interval, &
        discont_flag_01(4:5))/2.

    lnln_pole_int = lnln_pole_int + &
      N2_branch_cut_corr(klb,kub,a,c,f,k_crossings(3:5),discont_flag_01(3:5),discont_flag_02(3:5),acf_colinear)/2.

  ! if acf_colinear, N1 is zero for all k except for possibly a single value of k at which (1-ξy) lies on R+ and (1-y) lies on R- .
  ! in an integral over k, this single point can be ignored, so all correction terms proportional to N1 are 0.
    if(.not.acf_colinear) then
      lnln_pole_int = lnln_pole_int + &
        N1_branch_cut_corr(klb, kub, a,c,f, k_crossings([1,2,4,5]), discont_flag_01([1,2,4,5]))/2.
    endif

    lnln_pole_int = lnln_pole_int + &
      N4_and_N5_branch_cut_corr(klb, kub, a,c,f, k_crossings(:2), discont_flag_01(:2), discont_flag_02(:2), acf_colinear)

  endif

  lnln_pole_int = int_order_sign * lnln_pole_int


end function

