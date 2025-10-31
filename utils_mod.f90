module utils_mod
  contains
! assemble into a sorted list all of:
!     klb,
!     kub,
!     each point k in k_crossings(ic) that satisfies is_discont(ic)=.true. (klb<k<kub is
!       necessary but not sufficient)

  ! out:
  ! k_points : the sorted list
  ! n_k_points : the number of points in k_points
  ! orig_idxs and sort_idxs : track the indexes of k_crossings so that it's possible to determine this index for a point in the sorted list
  subroutine assemble_k_points(klb,kub,k_crossings,is_discont, k_points, orig_idxs, sort_idxs, n_k_points)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    implicit none
    real :: klb, kub

    real, dimension(:), intent(in) :: k_crossings
    logical, dimension(:), intent(in) :: is_discont
    real, dimension(:), intent(out) :: k_points
    integer, dimension(:), intent(out) :: orig_idxs
    integer, dimension(:), intent(out) :: sort_idxs

    integer :: ic, jc, n_k_points

    orig_idxs = 0
    sort_idxs = 0

    k_points = 0.0
    k_points(1) = klb
    k_points(2) = kub
    jc = 2

    orig_idxs(1)=-1
    orig_idxs(2)=-1

    do ic=1,size(k_crossings)
      if(is_discont(ic)) then
        jc=jc+1
        k_points(jc) = k_crossings(ic)
        orig_idxs(jc) = ic
      endif
    enddo
    n_k_points = jc
    call exchange_sort(k_points(:n_k_points), sort_idxs(:n_k_points), n_k_points)
  end subroutine assemble_k_points

  subroutine exchange_sort(x, sort_idxs, x_len)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    implicit none
    real, dimension(x_len), intent(inout) :: x
    integer, dimension(size(x)), intent(out) :: sort_idxs
    integer :: x_len

    sort_idxs = argsort(x)
    x = x(sort_idxs)

  end subroutine exchange_sort

! return the indices sort_idxs such that x_in(sort_idxs) is in non-descending order with all nan
!   values at the end.  uses exchange sort, so should only be applied to short lists, size(x_in) < 10
  function argsort(x_in) result(sort_idxs)
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    implicit none

    real, dimension(:) :: x_in
    real, dimension(size(x_in)) :: x
    integer, dimension(size(x_in)) :: sort_idxs
    integer :: temp2
    real :: temp

    integer :: ic, jc

    x = x_in

    sort_idxs = [ (ic, ic=1, size(x)) ]

    do ic=1,size(x)
      if(isnan(x(ic))) then
        x(ic) = ieee_value(x(ic), ieee_positive_inf)
      endif
    enddo

    do ic=1,size(x)-1
      do jc=ic+1, size(x)
        if(x(ic).gt.x(jc)) then
          temp = x(jc)
          x(jc) = x(ic)
          x(ic) = temp

          temp2 = sort_idxs(jc)
          sort_idxs(jc) = sort_idxs(ic)
          sort_idxs(ic) = temp2
        endif
      enddo
    enddo
  end function 

  ! div=.false. case :
  ! compute
  !                  ⎧ +1  Arg(val1) + Arg(val2) <= -π
  ! N_+(val1,val2) = ⎨ -1   Arg(val2) + Arg(val2 ) > π
  !                  ⎩ 0            else

  ! div=.true. case :
  ! compute
  !                  ⎧ +1  Arg(val1) - Arg(val2) <= -π
  ! N_-(val1,val2) = ⎨ -1   Arg(val2) - Arg(val2 ) > π
  !                  ⎩ 0            else
  ! these quantities obey
  ! ln(val1 * val2) = ln(val1) + ln(val2) + 2πi * N₊(val1,val2)
  ! and
  ! ln(val1 / val2) = ln(val1) - ln(val2) + 2πi * N₋(val1,val2)

  function split_ln_corr(val1, val2, div)
      use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
      implicit none
      real :: pi
      complex :: val1, val2
      real :: ang1, ang2, angle
      integer :: split_ln_corr
      logical :: div

      if(aimag(val1).eq.0.0) then ! handling signed zero
        ang1 = atan2(0.0,real(val1))
      else
        ang1 = atan2(aimag(val1),real(val1))
      endif
      if(aimag(val2).eq.0.0) then ! handling signed zero
        ang2 = atan2(0.0,real(val2))
      else
        ang2 = atan2(aimag(val2),real(val2))
      endif

      if(div) then
        angle = ang1 - ang2
      else
        angle = ang1 + ang2
      endif

      pi = atan(1.0)*4.0

      if( angle.gt.pi) then
          split_ln_corr = -1
      else if ( angle.le.-pi) then
          split_ln_corr = 1
      else 
          split_ln_corr = 0
      endif

  end function

  ! determine whether a, c, and f lie in a straight line in the complex plane.  Assume that the real and imaginary parts of a, c, and f are accurate to within a relative tolerence given by the machine epsilon
  function acf_are_colinear(a,c,f)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    implicit none
    complex :: a, c, f 
    complex :: cf_diff, af_diff
    complex :: cf_diff_error, af_diff_error
    real :: term1_min, term1_max, term2_min, term2_max
    logical :: acf_are_colinear

    call diff_with_error(c,f,cf_diff,cf_diff_error)
    call diff_with_error(a,f,af_diff,af_diff_error)

    call mult_max_min(aimag(cf_diff),real(af_diff),aimag(cf_diff_error),real(af_diff_error), term1_min, term1_max)
    call mult_max_min(real(cf_diff),aimag(af_diff),real(cf_diff_error),aimag(af_diff_error), term2_min, term2_max)

    if ((term1_min.lt.term2_max).and.(term2_min.lt.term1_max)) then
      acf_are_colinear=.true.
    else
      acf_are_colinear=.false.
    endif

  end function acf_are_colinear

  ! get the difference of x and y and also the error on the result assuming x and y each have a relative tolerance given by the machine epsilon
  subroutine diff_with_error(x, y, xy_diff, xy_diff_error)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    implicit none

    complex :: x, y ! accurate to a relative error of TOL
    complex :: xy_diff, xy_diff_error

    complex :: i
    real :: TOL
    TOL = epsilon(real(x))
    i = (0.0,1.0)

    xy_diff = x-y
    xy_diff_error = 2*TOL*(abs(real(xy_diff))+i*abs(aimag(xy_diff))) + &
                      TOL*(abs(real(x))+i*abs(aimag(x))) + &
                      TOL*(abs(real(y))+i*abs(aimag(y)))

  end subroutine diff_with_error

  ! get the maximum and minimum of the product of x and y given that they are accurate to within ±x_err and ±y_err
  subroutine mult_max_min(x, y, x_err, y_err, min_result, max_result)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    implicit none
    real :: x, y, x_err, y_err, min_result, max_result

    real, dimension(4) :: results
    results(1) = (x-x_err)*(y-y_err)
    results(2) = (x+x_err)*(y-y_err)
    results(3) = (x-x_err)*(y+y_err)
    results(4) = (x+x_err)*(y+y_err)

    max_result = maxval(results)
    min_result = minval(results)

  end subroutine mult_max_min



  ! for integrals such as
  !  tub
  ! ⌠    ln(1-t)
  ! ⎮    ─────── N dt ,
  ! ⌡       t
  !  tlb
  ! which arise from expansion of integrands such as t^(-1) [ln(ζ) + ln(1-t) ± 2πi N]^2,
  ! it is sometimes necessary to evaluate dilog at a point, t_c, that is on the branch cut of
  ! dilog, i.e., [1,∞). Here, N is a split ln correction coefficient such as  N_+(ζ, 1-t). The
  ! code below effectively places the evaluation point t_c on the side of the real axis such
  ! that, in the relevant integration interval (a subinterval of [tlb,tub]), the t contour does
  ! not cross the real axis.

  ! in the acf_colinear case, the need for such corrections also arises for integrals with the form
  !  tub
  ! ⌠     1
  ! ⎮    ─── N dt ,
  ! ⌡     t
  !  tlb
  ! where, e.g., N = N_+(ζ, 1-t).  Typically, ∫ t^{-1} dt can be safely evaluated as ln(t) at all
  ! points of interest, i.e., where N is discontinuous or at t=tlb or tub , because t is in
  ! general not on ℝ- at these points.  In the acf_colinear case, at step discontinuities in N
  ! produced by ζ(1-t) crossing ℝ-, t is often also on ℝ-, so ∫ t^{-1} dt cannot be safely
  ! evaluated as ln(t).  In particular, if t is on ℝ-, it must be ensured that the correct branch 
  ! of ln(t) is chosen.

  ! In all cases, t_c is assumed to be the only place where the contour crosses a branch cut

  ! contour_endpoint : value of t at the upper bound of the contour, used to determine whether
  ! the contour goes from the upper half to the lower half plane or vice versa
  ! fn_eval : dilog(t_c) or ln(t_c)
  subroutine adjust_complex_part(contour_endpoint, fn_eval)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    implicit none

    complex :: contour_endpoint
    complex :: fn_eval

    integer :: sign

    if(aimag(contour_endpoint).lt.0.0) then
      sign = 1
    else
      sign = -1
    endif

    fn_eval = cmplx(real(fn_eval),sign*abs(aimag(fn_eval)))

  end subroutine adjust_complex_part
end module
