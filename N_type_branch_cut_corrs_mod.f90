module N_type_branch_cut_corrs_mod
  use :: di_tri_log_mod, only : clog
  contains


  ! with the help of dilog_term_corr_B, evaluate the terms in
  !  yub
  ! ⌠   [log(1-ξy) - log(1-y) + 2πi N_1]^2
  ! ⎮   ────────────────────────────────── dy ,
  ! ⌡               y
  !  ylb
  ! that involve N_1.

  ! k_crossings : k points where where 1-ξy, 1-y, and (1-ξy)/(1-y) cross the real axis
  ! y = (f+k)/(f-a)
  ! ξy = (f+k)/(f-c)
  ! discont_flag_01 : true for k crossing points k_c where 1-ξy, 1-y, and (1-ξy)/(1-y) lie
  ! specifically on the negative real axis ℝ- AND that satisfy klb<k_c<kub
  !   At k points at which 1-ξy or 1-y lies on ℝ-, both N_1 is discontinuous and, respectively, dilog(ξy) or
  !   dilog(y) is discontinuous
  !   At k points at which (1-ξy)/(1-y) lies on ℝ- , N_1 is discontinuous
  function N1_branch_cut_corr(klb, kub, a, c, f, k_crossings, discont_flag_01)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    use :: utils_mod
    use :: dilog_terms_mod
    implicit none

    real :: klb, kub
    complex :: a, c, f
    complex :: N1_branch_cut_corr
    real :: pi
    complex :: i

    real, dimension(4) :: k_crossings
    logical, dimension(4) :: discont_flag_01

    real, dimension(6) :: k_points
    real, dimension(5) :: k_midpoints
    integer, dimension(6) :: sort_idxs
    integer, dimension(6) :: orig_idxs

    integer, dimension(5) :: N1
    complex, dimension(5) :: N1_arg1
    complex, dimension(5) :: N1_arg2
    complex, dimension(6) :: xiy
    complex, dimension(6) :: y

    integer :: dilog_is_discont_point_xiy
    integer :: dilog_is_discont_point_y

    integer :: n_k_points
    
    integer :: ic

    pi = 4*atan(1.0)
    i = (0.0,1.0)

    call assemble_k_points(klb,kub,k_crossings,discont_flag_01,k_points,orig_idxs,sort_idxs,n_k_points)

    xiy(:n_k_points) = (f+k_points(:n_k_points))/(f-c)
    y(:n_k_points) = (f+k_points(:n_k_points))/(f-a)

    dilog_is_discont_point_xiy = findloc(orig_idxs(sort_idxs(:n_k_points)),2, 1)
    dilog_is_discont_point_y = findloc(orig_idxs(sort_idxs(:n_k_points)),1, 1)

    N1_branch_cut_corr = 0.0
    k_midpoints(:n_k_points-1) = (k_points(2:n_k_points) + k_points(1:n_k_points-1))/2.
    N1_arg1(:n_k_points-1) = -(c+k_midpoints(:n_k_points-1))/(f-c)
    N1_arg2(:n_k_points-1) = -(a+k_midpoints(:n_k_points-1))/(f-a)
    do ic=1,n_k_points-1
  ! handling edge cases separately to avoid bug from rounding
      if((aimag(a).eq.0).and.(aimag(c).eq.0).and.(aimag(f).eq.0)) then
        if((real(N1_arg1(ic)).gt.0).and.(real(N1_arg2(ic)).lt.0)) then
          N1(ic) = 1
        else
          N1(ic) = 0
        endif
      else if((aimag(f).eq.0).and.(aimag(c).eq.0)) then
        if ((aimag(N1_arg2(ic)).lt.0).and.(real(N1_arg1(ic)).lt.0.0)) then
          N1(ic) =-1
        else
          N1(ic) =0
        endif
      else if((aimag(f).eq.0).and.(aimag(a).eq.0)) then
        if ((aimag(N1_arg1(ic)).lt.0).and.(real(N1_arg2(ic)).lt.0.0)) then
          N1(ic) =1
        else
          N1(ic) =0
        endif
      else
        N1(ic) = split_ln_corr(N1_arg1(ic), N1_arg2(ic), .true.)
      endif
    enddo

    N1_branch_cut_corr = N1_branch_cut_corr - 4*pi*i*sum(N1(:n_k_points-1)*&
      dilog_term_corr_B(&
        xiy(:n_k_points),&
        dilog_is_discont_point_xiy)&
        )

    N1_branch_cut_corr = N1_branch_cut_corr + 4*pi*i*sum(N1(:n_k_points-1)*&
      dilog_term_corr_B(&
        y(:n_k_points),&
        dilog_is_discont_point_y)&
        )

    do ic=1,n_k_points-1
      ! N1 is guaranteed to be zero wherever y crosses the negative real axis, so no special
      ! treatment of the ∫y⁻¹ dy integral is needed
      N1_branch_cut_corr = N1_branch_cut_corr + (N1(ic)*2*pi*i)**2 * (clog(y(ic+1)) - clog(y(ic)))
    enddo
  end function


  ! with the help of dilog_term_corr_B and ln_term_corr_(A/B) , evaluate the terms in
  !  uub
  ! ⌠   [log(ξ) + log(1-u) + 2πi N_2]^2
  ! ⎮   ────────────────────────────────── du ,
  ! ⌡               u
  !  ulb
  !  that involve N_2.
  ! ξ = (f-a)/(f-c)
  ! u = -(c-a)/(k+a)

  ! k_crossings : k points where where 1-u and ξ(1-u)=(1-ξy)/(1-y) cross the real axis
  ! y = (f+k)/(f-a) 
  ! ξy = (f+k)/(f-c)
  ! for convenience, define v that can be 
  !     ⎧ 1-u
  ! v = ⎨ or
  !     ⎩ ξ(1-u)
  !   depending on the context
  ! discont_flag_01 : true for crossing points k that satisfy klb<k<kub AND for which v(k) lies on the negative real axis
  !     At k points at which 1-u lies on ℝ-, both N_2 is discontinuous and dilog(u) is discontinuous
  !     At k points at which ξ(1-u) lies on ℝ-, N_2 is discontinuous.  If acf_colinear,
  !     dilog(u) may also be discontinuous
  ! discont_flag_02 : true for crossing points k that satisfy klb<k<kub AND for which 1-v(k) lies on the negative real axis
  ! acf_colinear : whether a, c, and f lie in a straight line in the complex plane.  When this is true, ξ=xi is pure real

  function N2_branch_cut_corr(klb, kub, a,c,f,k_crossings,discont_flag_01, discont_flag_02, acf_colinear)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    use :: utils_mod
    use :: ln_terms_mod
    use :: dilog_terms_mod
    implicit none
    real :: klb, kub
    complex :: a, c, f
    complex :: N2_branch_cut_corr
    real :: pi
    complex :: i

    logical, dimension(3), intent(in) :: discont_flag_01, discont_flag_02
    real, dimension(3), intent(in) :: k_crossings

    complex :: xi

    real, dimension(5) :: k_points
    real, dimension(4) :: k_midpoints
    integer, dimension(5) :: orig_idxs
    integer, dimension(5) :: sort_idxs
    integer :: n_k_points

    integer :: ic
    logical :: acf_colinear

    complex, dimension(5) :: u
    integer :: dilog_is_discont_point
    integer, dimension(4) :: N_2
    complex, dimension(4) :: ln_u_terms

    integer :: ln_is_discont_interval, ln_is_discont_point
    
    pi = 4*atan(1.0)
    i = (0.0,1.0)
    if(acf_colinear) then
      xi = real((a-f)/(c-f))
    else
      xi = (a-f)/(c-f)
    endif

    dilog_is_discont_point=0
    if ( acf_colinear ) then ! in this case, there is maximally one k_crossing, given by k_crossings(1)=k_crossings(2) in this function
      if(discont_flag_01(1).or.discont_flag_01(2)) then
        k_points(1) = klb
        k_points(2) = k_crossings(1)
        k_points(3) = kub
        n_k_points = 3
      else
        k_points(1)=klb
        k_points(2)=kub
        n_k_points = 2
      endif
      if(discont_flag_01(1)) then
        dilog_is_discont_point = 2
      else
        dilog_is_discont_point = 0
      endif
    else
      call assemble_k_points(klb,kub,k_crossings,discont_flag_01,k_points,orig_idxs,sort_idxs,n_k_points)
      dilog_is_discont_point = findloc(orig_idxs(sort_idxs(:n_k_points)),1,1)
    endif
    
    u(:n_k_points) = (a - c)/(k_points(:n_k_points) + a)

    k_midpoints(:n_k_points-1) = (k_points(2:n_k_points) + k_points(:n_k_points-1))/2.
    do ic=1,n_k_points-1
  ! handling edge cases separately to avoid bug from rounding
      if((aimag(a).eq.0).and.(aimag(c).eq.0).and.(aimag(f).eq.0)) then
        if((real(xi).lt.0).and.(real((c+k_midpoints(ic))/(a+k_midpoints(ic))).lt.0)) then
          N_2(ic) =-1
        else
          N_2(ic) =0
        endif
      else if ((aimag(c).eq.0).and.(aimag(a).eq.0)) then
        if((real((c+k_midpoints(ic))/(a+k_midpoints(ic))).lt.0).and.(aimag(xi).gt.0)) then
          N_2(ic) =-1
        else
          N_2(ic) =0
        endif
      else
        N_2(ic) = split_ln_corr(xi,(c+k_midpoints(ic))/(a+k_midpoints(ic)),.false.)
      endif
    enddo

    N2_branch_cut_corr = -4*pi*i*sum(N_2(:n_k_points-1)*&
      dilog_term_corr_B(&
        u(:n_k_points),&
        dilog_is_discont_point)&
        )

    ln_u_terms(:n_k_points-1) = N_type_ln_terms_corr_A_or_B(u(:n_k_points),k_points(:n_k_points),&
      k_crossings(1), discont_flag_02(1), discont_flag_01(2), acf_colinear)

    N2_branch_cut_corr = N2_branch_cut_corr + &
      sum((clog(xi)*N_2(:n_k_points-1)*4*pi*i + (N_2(:n_k_points-1)*2*pi*i)**2) * &
            ln_u_terms(:n_k_points-1))
      
  end function


  ! drive computation of terms involving N_4 and N_5 in
  ! ⌠    [ln(a-f) + ln(1-y) + 2πiN_4][ln(c-f) + ln(1-ξy) + 2πiN_5]
  ! ⎮ dy ───────────────────────────────────────────────────────── .
  ! ⌡                            y
  ! Also, compute the most simple term directly in this fortran function.
  ! Here,
  ! N_4 = N_+(a-f,1-y),
  ! N_5 = N_+(c-f,1-ξy),
  ! y = (f+k)/(f-a)
  ! ξ = (f-a)/(f-c)
  ! k_crossings(1) : k point where y crosses the real axis
  ! k_crossings(2) : k point where ξy crosses the real axis

  ! discont_flag_01(1) : whether 1-y<0 at the point where it crosses the real axis AND klb<k_crossings(1)<kub
  ! In this case, N_4 and dilog(y) are discontinuous at this point

  ! discont_flag_01(2) : whether 1-ξy<0 at the point where it crosses the real axis AND klb<k_crossings(2)<kub.
  ! In this case, N_5 and dilog(ξy) are discontinuous at this point

  ! discont_flag_02(1) : whether y<0 at the point where it crosses the real axis AND klb<k_crossings(1)<kub.  In this case, ln(y) is discontinuous at this point

  ! acf_colinear : whether a, c, and f lie in a straight line in the complex plane
  function N4_and_N5_branch_cut_corr(klb, kub, a, c, f, k_crossings, discont_flag_01, discont_flag_02, acf_colinear)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    use :: utils_mod
    implicit none
    real :: klb, kub
    complex :: a, c, f
    complex :: N4_and_N5_branch_cut_corr
    integer :: N_4, N_5
    real :: pi
    complex :: i
    complex :: clog

    real, dimension(2), intent(in) :: k_crossings
    logical, dimension(2), intent(in) :: discont_flag_01
    logical, dimension(2), intent(in) :: discont_flag_02

    real :: k_midpoint

    real, dimension(4) :: k_points
    integer, dimension(4) :: sort_idxs, orig_idxs
    integer :: n_k_points
    
    integer :: ic
    
    logical :: acf_colinear
    
    i=cmplx(0.0,1.0)
    
    pi = 4*atan(1.0)

    N4_and_N5_branch_cut_corr = 0.0

    ! N_4 corrections
    N4_and_N5_branch_cut_corr = N4_and_N5_branch_cut_corr + N4_or_N5_branch_cut_corr(&
        klb, kub, a, c, a, f, k_crossings([1, 2]), discont_flag_01([1, 2]),&
        k_crossings(1), discont_flag_02(1), acf_colinear)

    ! N_5 corrections
    N4_and_N5_branch_cut_corr = N4_and_N5_branch_cut_corr + N4_or_N5_branch_cut_corr(&
        klb, kub, a, a, c, f, k_crossings([2, 1]), discont_flag_01([2, 1]),&
        k_crossings(1), discont_flag_02(1), acf_colinear)


    if(acf_colinear) then ! k_crossings(2)=k_crossings(1) in this case
      if(discont_flag_01(1).or.discont_flag_01(2)) then
        k_points(1) = klb
        k_points(2) = k_crossings(1)
        k_points(3) = kub
        n_k_points = 3
      else
        k_points(1) = klb
        k_points(2) = kub
        n_k_points = 2
      endif
    else
      call assemble_k_points(klb,kub,k_crossings,discont_flag_01,k_points,orig_idxs,sort_idxs,n_k_points)
    endif

    do ic=1,n_k_points-1
      k_midpoint = (k_points(ic) + k_points(ic+1))/2.

      N_5 = N4_or_N5_split_ln_corr(c, f, k_midpoint)
      N_4 = N4_or_N5_split_ln_corr(a, f, k_midpoint)

      ! N_4 is guaranteed to be zero wherever u crosses the negative real axis, so no special
      ! treatment of the ∫u⁻¹ du integral is needed
      N4_and_N5_branch_cut_corr = N4_and_N5_branch_cut_corr + N_4*N_5*(2*pi*i)**2* &
        (clog((f+k_points(ic+1))/(f-a)) - clog((f+k_points(ic))/(f-a)))
    enddo

  end function

  ! with the help of ln_term(s)_corr_(A/B) and dilog_term_corr_(A/B), evaluate
  !
  ! N_type_ln_terms_corr_A_or_B = 
  !  yub
  ! ⌠     ln(ζ)
  ! ⎮   N ───── dy
  ! ⌡       y
  !  ylb
  !
  ! dilog_t_terms=
  !   tub
  !  ⌠     ln(1-t)
  !  ⎮   N ─────── dt
  !  ⌡        t
  !   tlb
  !
  ! Here, 
  ! N=N_+(λ,1-s)
  ! y=(f+k)/(f-a)

  ! Case 1:
  ! a_or_c = a
  ! c_or_a = c
  ! ζ = (a-f)
  ! λ = (c-f)
  ! t = y
  ! s = ξy   , where ξ=(a-f)/(c-f)
  ! N = N_5

  ! Case 2:
  ! a_or_c = c
  ! c_or_a = a
  ! ζ = (c-f)
  ! λ = (a-f)
  ! t = ξy, where ξ=(a-f)/(c-f)
  ! s = y   
  ! N = N_4

  ! k_crossings(2) : k point where t (a fn. of k) crosses the real axis
  ! discont_flag_01(2) : whether 1-t < 0 at the point where t crosses the real axis AND whether klb<k_crossings(2)<kub

  ! k_crossings(1) : k point where s (a fn. of k) crosses the real axis
  ! discont_flag_01(1) : whether 1-s < 0 at the point where s crosses the real axis AND whether klb<k_crossings(1)<kub

  ! y_crosses_negR : whether y < 0 at the point where it crosses the real axis AND whether the k_point of this crossing point satisfies klb<k_point<kub

  ! acf_colinear : whether a, c, and f lie in a straight line in the complex plane
  function N4_or_N5_branch_cut_corr(klb, kub, a, a_or_c, c_or_a, f, k_crossings, &
      discont_flag_01, k_where_im_y_eq_0, y_crosses_negR, acf_colinear)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    use :: ln_terms_mod
    use :: dilog_terms_mod
    implicit none
    real :: klb, kub
    complex :: a, a_or_c, c_or_a, f
    complex :: N4_or_N5_branch_cut_corr
    real :: k_where_im_y_eq_0
    real :: pi
    complex :: i
    logical :: acf_colinear

    real, dimension(2), intent(in) :: k_crossings
    logical, dimension(2), intent(in) :: discont_flag_01
    logical :: y_crosses_negR

    real, dimension(2) :: k_midpoints
    complex, dimension(3) :: y_points

    real, dimension(3) :: k_points
    integer :: n_k_points

    complex, dimension(3) :: dilog_args
    integer :: dilog_is_discont_point
    integer, dimension(2) :: N4_or_N5
    complex, dimension(2) :: dilog_t_terms

    integer :: ic
    
    integer :: dilog_is_discont_interval, ln_is_discont_interval, ln_is_discont_point

    i=cmplx(0.0,1.0)
    pi = 4*atan(1.0)
    
    if(discont_flag_01(1)) then
      k_points(1) = klb
      k_points(2) = k_crossings(1)
      k_points(3) = kub
      n_k_points = 3
    else
      k_points(1) = klb
      k_points(2) = kub
      n_k_points = 2
    endif

    dilog_args(:n_k_points) = (f + k_points(:n_k_points))/(f - a_or_c)

    N4_or_N5_branch_cut_corr = 0.0

    ! if acf_colinear, then (1-ξy) and (1-y) cross the real axis at the same value of k.  As a
    ! result, referring to the definition of dilog_t_terms in the function header, both N
    ! and dilog(t) can be discontinuous at the same k.  In such cases, it is necessary to
    ! ensure that the correct value of dilog(t) is used.
    if(acf_colinear) then
      if(discont_flag_01(1)) then
        dilog_is_discont_point = 2
      else
        dilog_is_discont_point = 0
      endif
      dilog_t_terms(:n_k_points-1) = dilog_term_corr_B(&
        dilog_args(:n_k_points),&
        dilog_is_discont_point)
    else
      dilog_is_discont_interval = 0
      if(discont_flag_01(2)) then
        dilog_is_discont_interval = findloc((k_points(2:n_k_points).gt.k_crossings(2)).and.&
          (k_points(:n_k_points-1).lt.k_crossings(2)),.true.,1)
      endif

  ! for the presently treated terms, if .not.acf_colinear, then N and dilog(t) are
  ! discontinuous at separate points, so it is appropriate to use dilog...corr_A.
      dilog_t_terms(:n_k_points-1) = dilog_term_corr_A(&
        dilog_args(:n_k_points),dilog_is_discont_interval, real((f + k_crossings(2))/(f - a_or_c)))
    endif

    k_midpoints(:n_k_points-1) = (k_points(2:n_k_points) + k_points(:n_k_points-1))/2.
    do ic=1,n_k_points-1
      N4_or_N5(ic) = N4_or_N5_split_ln_corr(c_or_a, f, k_midpoints(ic))
    enddo

    N4_or_N5_branch_cut_corr = -2*pi*i*sum(N4_or_N5(:n_k_points-1)*dilog_t_terms(:n_k_points-1))

    y_points(:n_k_points) = (f+k_points(:n_k_points))/(f-a)

    N4_or_N5_branch_cut_corr = N4_or_N5_branch_cut_corr + &
      sum(N4_or_N5(:n_k_points-1)*2*pi*i*clog(a_or_c-f)*&
        N_type_ln_terms_corr_A_or_B(y_points(:n_k_points),k_points(:n_k_points),k_where_im_y_eq_0,&
        y_crosses_negR,discont_flag_01(1),acf_colinear))



  end function


  ! with the help of dilog_term_corr_B, evaluate the terms in
  !  yub                              2 
  ! ⌠    [ln(a-f) + ln(1-y) + 2πi N_4]
  ! ⎮    ───────────────────────────── dy
  ! ⌡                  y
  !  ylb
  ! that involve N_4.  The expression above is only relevant when a=c.
  ! y = (f+k)/(f-a)
  ! k_crossing : k point where y crosses the real axis
  ! discont_flag_01 : true if 1-y is on (-∞,0] at the k_crossing AND klb<k_crossing<kub.  In
  ! this case, both N_4 and dilog(y) are discontinuous at k_crossing

  function N4_branch_cut_corr_a_eq_c(klb, kub, a,f,k_crossing,discont_flag_01)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    use :: dilog_terms_mod
    implicit none
    real :: klb, kub
    complex :: a, f
    complex :: N4_branch_cut_corr_a_eq_c
    real :: pi
    complex :: i
    integer :: ic

    real :: k_crossing
    logical :: discont_flag_01

    real :: k_midpoint

    real, dimension(3) :: k_points
    integer :: n_k_points

    complex, dimension(3) :: y_points
    integer :: dilog_is_discont_point
    integer, dimension(2) :: N_4
    pi = 4*atan(1.0)
    i=(0.0,1.0)

    if ( discont_flag_01 ) then
      k_points(1) = klb
      k_points(2) = k_crossing
      k_points(3) = kub
      n_k_points = 3
      dilog_is_discont_point = 2
    else
      k_points(1) = klb
      k_points(2) = kub
      n_k_points = 2
      dilog_is_discont_point = 0
    endif

    y_points(:n_k_points) = (f+k_points(:n_k_points))/(f-a)

    do ic=1,n_k_points-1
      k_midpoint = (k_points(ic) + k_points(ic+1))/2.
      N_4(ic) = N4_or_N5_split_ln_corr(a, f, k_midpoint)
    enddo

    N4_branch_cut_corr_a_eq_c = -4*pi*i*sum(N_4(:n_k_points-1)*&
      dilog_term_corr_B(&
        y_points(:n_k_points),&
        dilog_is_discont_point)&
        )

  ! no extra correction needed b/c, in an interval where y crosses the negative real axis, N_4 is
  ! guaranteed to be 0
    do ic=1,n_k_points-1
      N4_branch_cut_corr_a_eq_c = N4_branch_cut_corr_a_eq_c + (clog(a-f)*N_4(ic)*4*pi*i + (N_4(ic)*2*pi*i)**2) * &
        (clog((f + k_points(ic+1))/(f-a)) - clog((f + k_points(ic))/(f-a)))
    enddo
      
  end function


  ! split_ln_corr with some edge case handling
  !
  ! a_or_c = a :
  !    compute N_4 = N_+(a-f,1-y)

  ! a_or_c = c :
  !    compute N_5 = N_+(c-f,1-ξy)
  ! 
  ! Here, 
  ! y=(f+k_point)/(f-a)
  ! ξ=(a-f)/(c-f)
  !
  ! N_+ is defined in utils_mod.f90
  function N4_or_N5_split_ln_corr(a_or_c, f, k_point)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    use :: utils_mod
    implicit none
    real :: k_point
    complex :: a_or_c, f
    integer :: N4_or_N5_split_ln_corr

  ! handling edge case separately to avoid bug from rounding
    if((aimag(a_or_c).eq.0).and.(aimag(f).eq.0)) then
      if((real(a_or_c-f).lt.0).and.(real(-(a_or_c+k_point)/(f-a_or_c)).lt.0)) then
        N4_or_N5_split_ln_corr=-1
      else
        N4_or_N5_split_ln_corr=0
      endif
    else if(aimag(a_or_c).eq.0) then
      if((real(a_or_c+k_point).gt.0).or.(aimag(a_or_c-f).gt.0)) then
        N4_or_N5_split_ln_corr = 0
      else
        N4_or_N5_split_ln_corr = 1
      endif
    else
      N4_or_N5_split_ln_corr = split_ln_corr(a_or_c-f, -(a_or_c+k_point)/(f-a_or_c), .false.)
    endif
  end function N4_or_N5_split_ln_corr

! compute ∫t⁻¹ dt with a possible type A or type B correction (as defined in log_terms_mod).
! ln_args : t bounds of the integral
!   Case 1 : t = y = (f+k_points)/(f-a)
!   Case 2 : t = u = (a-c)/(a+k_points)
! k_points : value of k for each t value in ln_args
! k_where_im_t_eq_0 : value of k at which the t contour crosses the real axis
! corr_flag_1 : whether t crosses specifically the negative real axis at k=k_where_im_t_eq_0
!   AND klb<k_where_im_t_eq_0<kub 
! corr_flag_2.and.acf_colinear : whether the place where t crosses the real axis is one of the
!   ln_args, in particular the second of a total of three in the relevant cases.
! results : ln_t_terms(ic) = ∫t⁻¹ dt between ln_args(ic) and ln_args(ic+1)
!   if(.not.corr_flag_1), then ln_t_terms(ic) = ln(ln_args(ic+1)) - ln(ln_args(ic))
! acf_colinear : whether a, c, and f lie in a straight line in the complex plane
  function N_type_ln_terms_corr_A_or_B(ln_args, k_points, k_where_im_t_eq_0, corr_flag_1, corr_flag_2, acf_colinear) result(ln_t_terms)
    use :: ln_terms_mod
    implicit none
    complex, dimension(:), intent(in) :: ln_args
    real, dimension(:), intent(in) :: k_points
    complex, dimension(size(ln_args)-1) :: ln_t_terms
    integer :: ln_is_discont_point
    integer :: ln_is_discont_interval
    logical :: corr_flag_1, corr_flag_2, acf_colinear
    real :: k_where_im_t_eq_0

    if(acf_colinear) then
      if(corr_flag_1.and.corr_flag_2) then
        ln_is_discont_point = 2
      else
        ln_is_discont_point = 0
      endif
      ln_t_terms = ln_terms_corr_B(ln_args,ln_is_discont_point)
    else
      if(corr_flag_1) then
        ln_is_discont_interval = findloc(k_points > k_where_im_t_eq_0,.true.,1)-1
      else
        ln_is_discont_interval = 0
      endif
      ln_t_terms = ln_term_corr_A(ln_args,ln_is_discont_interval)
    endif
  end function

end module

