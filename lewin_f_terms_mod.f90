! evaluate 
!  t₂
! ⌠   log²(1-t′)
! ⎮   ────────── dt′
! ⌡       t′
!  t₁
! along a contour that is specified in terms of the fortran function arguments below.
! Do this using the function 𝕗 defined in eq. 8.105 of Polylogarithms and Associated Functions
! (1981, North Holland)
!         t
!        ⌠ log²(1-t′)
! 𝕗(t) = ⎮ ────────── dt′ .
!        ⌡     t′
!         0      
! 𝕗(t) is dicontinuous across t ∈ [1,∞), and these discontinuities need to be subtracted
! out if the contour crosses [1,∞).
module lewin_f_terms_mod
    interface lewin_f_term_corr_A
      module procedure lewin_f_term_corr_A_v1, lewin_f_term_corr_A_v2
    end interface
  contains

! evaluate the integral in the module header along a contour defined by
! endpoints
!     lewin_f_args(1) = t₁ 
!     lewin_f_args(2) = t₂
! lewin_f_arg_crossing : an exclusive point where the contour possibly crosses [1,∞)
! discont_flag : whether the contour does actually cross [1,∞) at lewin_f_arg_crossing
  function lewin_f_term_corr_A_v1(lewin_f_args, lewin_f_arg_crossing, discont_flag)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    use :: utils_mod
    use :: di_tri_log_mod
    implicit none

    complex, dimension(2), intent(in) :: lewin_f_args
    real :: lewin_f_arg_crossing
    logical :: discont_flag
    integer :: sign

    complex :: lewin_f_term_corr_A_v1

    lewin_f_term_corr_A_v1 = lewin_f(lewin_f_args(2)) - lewin_f(lewin_f_args(1))
    if(discont_flag) then
      if(aimag(lewin_f_args(1)).gt.0) then
        sign=-1
      else
        sign=1
      endif
      lewin_f_term_corr_A_v1 = lewin_f_term_corr_A_v1 + sign*lewin_f_discont(lewin_f_arg_crossing)
    endif
  end function

! evaluate the integral in the module header along a contour defined by
!   endpoints
!     lewin_f_args(1) = t₁ 
!     lewin_f_args(2) = t₂
!   lewin_f_arg_crossings : points where the contour crosses the real axis (and possibly other
!      points as well that are irrelevant to the considered interval)
!   sort_idxs : order in which the points lewin_f_arg_crossings are reached
!   in_interval : whether a point in lewin_f_arg_crossings is in the contour interval (t₁,t₂)
!   discont_flag : whether a point in lewin_f_arg_crossings lies specifically on [1,∞)
! The sign of the discontinuity in 𝕗(t) as the contour traverses [1,∞) depends on whether the
! contour is going from the upper half plane to the lower half plane or vice-versa,
! which is why all places where the contour crosses ℝ need to be tracked, not only those on [1,∞).
  function lewin_f_term_corr_A_v2(lewin_f_args, lewin_f_arg_crossings, sort_idxs, in_interval, discont_flag)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    use :: di_tri_log_mod
    implicit none

    complex, dimension(2), intent(in) :: lewin_f_args
    real, dimension(:), intent(in) :: lewin_f_arg_crossings
    logical, dimension(:), intent(in) :: discont_flag
    logical, dimension(:), intent(in) :: in_interval
    integer, dimension(:) :: sort_idxs

    integer :: sign
    complex :: lewin_f_term_corr_A_v2
    integer :: ic

    lewin_f_term_corr_A_v2 = lewin_f(lewin_f_args(2)) - lewin_f(lewin_f_args(1))

    if(aimag(lewin_f_args(1)).gt.0) then
      sign=1
    else
      sign=-1
    endif

    do ic=1,size(lewin_f_arg_crossings)
      if(in_interval(sort_idxs(ic))) then
        sign = -1*sign
      endif
      if(discont_flag(sort_idxs(ic))) then
        lewin_f_term_corr_A_v2 = lewin_f_term_corr_A_v2 + &
          sign*lewin_f_discont(lewin_f_arg_crossings(sort_idxs(ic)))
      endif
    enddo
  end function

end module
