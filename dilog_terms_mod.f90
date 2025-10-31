! evaluate
!  t₂
! ⌠    log(1-t′)
! ⎮    ─────────dt′
! ⌡        t′
!  t₁
! along a t' contour that is specified in terms of the fortran function arguments below.
! Do this using the dilogarithm function
!              t
!             ⌠  log(1-t′)
! dilog(t) = -⎮  ───────── dt′ .
!             ⌡      t′
!              0
! dilog(t) is discontinuous across t ∈ [1,∞), and these discontinuities need to be
! compensated for if the contour crosses [1,∞) or if a contour endpoint lies on [1,∞).

module dilog_terms_mod
  use :: di_tri_log_mod
    interface dilog_term_corr_A
      module procedure dilog_term_corr_A_single, dilog_term_corr_A_multiple
    end interface

  contains

! evaluate the integral in the module header along a t' contour that possibly crosses [1,∞) at a
! single point and that is additionally defined by
! endpoints
!     dilog_args(1) = t₁ 
!     dilog_args(2) = t₂
! is_discont : whether the contour crosses t' ∈ [1,∞)
! dilog_arg_where_im_t_eq_0 : when(is_discont), the value of t' at the point where the
!   contour crosses [1,∞)
! This function expects that if either of the endpoints lie on [1,∞), the contour coincides
! with the real axis.
  function dilog_term_corr_A_single(dilog_args, is_discont,dilog_arg_where_im_t_eq_0) result(dilog_term)
    implicit none

    complex, dimension(2) :: dilog_args
    complex, dimension(1) :: dilog_terms
    complex :: dilog_term
    
    logical :: is_discont
    real :: dilog_arg_where_im_t_eq_0

    if(is_discont) then
      dilog_terms = dilog_term_corr_A_multiple(dilog_args, 1, dilog_arg_where_im_t_eq_0)
      dilog_term = dilog_terms(1)
    else
      dilog_term = dilog(dilog_args(2)) - dilog(dilog_args(1))
    endif
  end function

  ! evaluate the same type of integral as described above dilog_term_corr_A_single but do so
  ! for all t' intervals specified by
  !     dilog_args(ic) = t₁ 
  !     dilog_args(ic+1) = t₂ ,
  ! Also, t' is assumed to cross [1,∞) in one interval at most.
  ! dilog_is_discont_interval : the interval in which the t' contour crosses [1,∞), or 0 if
  ! there is no such interval

  function dilog_term_corr_A_multiple(dilog_args, dilog_is_discont_interval,dilog_arg_where_im_t_eq_0) result(dilog_terms)
    implicit none

    complex, dimension(:) :: dilog_args
    complex, dimension(size(dilog_args)-1) :: dilog_terms

    complex :: corr_term
    integer :: sign
    
    integer :: dilog_is_discont_interval
    real :: dilog_arg_where_im_t_eq_0


    complex, dimension(size(dilog_args)) :: indef_dilog_terms
    integer :: ic

    do ic=1,size(dilog_args)
      indef_dilog_terms(ic) = dilog(dilog_args(ic))
    enddo

    do ic=1,size(dilog_args)-1
      if(ic.eq.dilog_is_discont_interval) then
        ! check whether the upper endpoint of the contour lies above or below the real axis
        if(aimag(dilog_args(size(dilog_args))).lt.0) then
          sign=-1
        else
          sign = 1
        endif
        corr_term = -sign*dilog_discont(dilog_arg_where_im_t_eq_0)
      else
        corr_term = 0.0
      endif

      dilog_terms(ic) = (indef_dilog_terms(ic+1) - indef_dilog_terms(ic) + corr_term)

    enddo

  end function

! evaluate the integral in the module header in a series of intervals specified by
!     dilog_args(ic) = t₁ 
!     dilog_args(ic+1) = t₂ ,
! where one of these interval endpoints lies on [1,∞).
! dilog_is_discont_point : specifies the endpoint that lies on [1,∞)
  function dilog_term_corr_B(dilog_args, dilog_is_discont_point) result(dilog_terms)
    use :: utils_mod
    implicit none

    complex, dimension(:) :: dilog_args
    integer :: dilog_is_discont_point
    complex, dimension(size(dilog_args)-1) :: dilog_terms

    complex, dimension(size(dilog_args)) :: indef_dilog_terms
    integer :: ic

    do ic=1,size(dilog_args)
      indef_dilog_terms(ic) = dilog(dilog_args(ic))
    enddo

    do ic=1,size(dilog_args)-1
      if(ic.eq.dilog_is_discont_point) then
        ! for dilog evaluated on the branch cut, choose the branch such that dilog is
        ! continuous in the given integration interval
        call adjust_complex_part(dilog_args(1),indef_dilog_terms(ic))
      endif
      if((ic+1).eq.dilog_is_discont_point) then
        ! for dilog evaluated on the branch cut, choose the branch such that dilog is
        ! continuous in the given integration interval
        call adjust_complex_part(dilog_args(size(dilog_args)),indef_dilog_terms(ic+1))
      endif

      dilog_terms(ic) = (indef_dilog_terms(ic+1) - indef_dilog_terms(ic))

    enddo

  end function
end module
