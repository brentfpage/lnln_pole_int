! evaluate
!  t₂
! ⌠   1
! ⎮   ── dt′
! ⌡   t′
!  t₁
! along a t' contour that is specified in terms of the fortran function arguments below.
! Do this using the logarithm function, defined by
! d ln(t)   1
! ─────── = ─
!   d t     t
! for all t∉(-∞,0].  ln(t) is discontinuous by 2πi across t∈(-∞,0], and these discontinuities
! need to be compensated for if the contour crosses (-∞,0] or if a contour endpoint lies on (-∞,0]

module ln_terms_mod
  interface ln_term_corr_A
    module procedure ln_term_corr_A_single, ln_term_corr_A_multiple
  end interface

  contains

! evaluate the integral in the module header along a t' contour that possibly crosses (-∞,0] at a
! single point and that is additionally defined by
! endpoints
!     ln_args(1) = t₁ 
!     ln_args(2) = t₂
! is_discont : whether the contour crosses t' ∈ (-∞,0]
! This function expects that if either of the endpoints lie on t∈(-∞,0], the contour coincides
! with the real axis.

function ln_term_corr_A_single(ln_args, is_discont) result(ln_term)
  use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
  use :: di_tri_log_mod, only : clog
  implicit none

  complex, dimension(2), intent(in) :: ln_args
  complex, dimension(1) :: ln_terms
  complex :: ln_term
  logical :: is_discont

  if(is_discont) then
    ln_terms = ln_term_corr_A_multiple(ln_args, 1)
    ln_term = ln_terms(1)
  else
    ln_term = clog(ln_args(2))-clog(ln_args(1))
  endif

end function

  ! evaluate the same type of integral as described above ln_term_corr_A_single but do so
  ! for all t' intervals specified by
  !     ln_args(ic) = t₁ 
  !     ln_args(ic+1) = t₂ ,
  ! Also, t' is assumed to cross t∈(-∞,0] in one interval at most.
  ! ln_is_discont_interval : the interval in which the t' contour crosses t∈(-∞,0], or 0 if
  ! there is no such interval

function ln_term_corr_A_multiple(ln_args, ln_is_discont_interval) result(ln_terms)
  implicit none

  complex, dimension(:), intent(in) :: ln_args
  complex, dimension(size(ln_args)-1) :: ln_terms
  integer :: sign
  complex :: ln_discont
  complex :: ln_t_discont_corr
  integer :: ln_is_discont_interval

  complex, dimension(size(ln_args)) :: indef_ln_terms
  integer :: ic

  ln_discont = cmplx(0,-2*4*atan(1.0)) ! -2*pi*i

  do ic=1,size(ln_args)
    indef_ln_terms(ic) = clog(ln_args(ic))
  enddo

  do ic=1,size(ln_args)-1
    ln_t_discont_corr = 0 
    if(ic.eq.ln_is_discont_interval) then
      if(aimag(ln_args(size(ln_args))).gt.0) then
        sign=-1
      else
        sign=1
      endif
      ln_t_discont_corr = -sign*ln_discont
    endif
    ln_terms(ic) = (indef_ln_terms(ic+1) - indef_ln_terms(ic) + ln_t_discont_corr)
  enddo
end function

! evaluate the integral in the module header in a series of intervals specified by
!     ln_args(ic) = t₁ 
!     ln_args(ic+1) = t₂ ,
! where one of these interval endpoints lies on (-∞,0]
! dilog_is_discont_point : specifies the endpoint that lies on (-∞,0]
function ln_terms_corr_B(ln_args, ln_is_discont_point) result(ln_terms)
  use :: utils_mod
  implicit none

  complex, dimension(:), intent(in) :: ln_args
  integer :: ln_is_discont_point

  complex, dimension(size(ln_args)) :: indef_ln_terms
  complex, dimension(size(ln_args)-1) :: ln_terms
  integer :: ic

  do ic=1,size(ln_args)
    indef_ln_terms(ic) = clog(ln_args(ic))
  enddo

  do ic=1,size(ln_args)-1
      if(ic.eq.ln_is_discont_point) then
        call adjust_complex_part(ln_args(1),indef_ln_terms(ic))
      else if ((ic+1).eq.ln_is_discont_point) then
        call adjust_complex_part(ln_args(size(ln_args)),indef_ln_terms(ic+1))
      endif
    ln_terms(ic) = (indef_ln_terms(ic+1) - indef_ln_terms(ic) )
  enddo
end function

end module
