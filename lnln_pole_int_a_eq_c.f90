! evaluate the following integral:
!  kub
! ⌠    ln²(a+k)
! ⎮    ──────── dk
! ⌡      f+k
!  klb
!  by changing variables to y=(f+k)/(f-a),
!   yub
!   ⌠   1 ⎛                     ⎞²
! = ⎮   ─ ⎝ln(1-y)+ln(a-f)+2πiN₄⎠  dy
!   ⌡   y
!   ylb
! and using lewin_f and the dilogarithm. Here, N₄ = N₊(1-y,a-f), where N₊ is defined in utils_mod.

function lnln_pole_int_a_eq_c(klb,kub,a,f)
  use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
  use :: di_tri_log_mod, only : clog
  use :: ln_terms_mod
  use :: dilog_terms_mod
  use :: lewin_f_terms_mod
  use :: k_crossings_mod
  use :: N_type_branch_cut_corrs_mod
  implicit none
  real :: klb, kub
  complex :: a, f
  complex :: lnln_pole_int_a_eq_c
  real :: k_crossing
  logical :: discont_flag_01
  logical :: discont_flag_02

  complex, dimension(2) :: y

  call get1_k_crossing(klb, kub, a, f, k_crossing, discont_flag_01, discont_flag_02)

  y = (f + [ klb, kub ])/(f-a)

  lnln_pole_int_a_eq_c = lewin_f_term_corr_A(y,real((f+k_crossing)/(f-a)),&
    discont_flag_01)

  lnln_pole_int_a_eq_c = lnln_pole_int_a_eq_c + N4_branch_cut_corr_a_eq_c(klb, kub, a,f,k_crossing,discont_flag_01)

  lnln_pole_int_a_eq_c = lnln_pole_int_a_eq_c + clog(a-f)**2 * ln_term_corr_A(y, discont_flag_02)

  lnln_pole_int_a_eq_c = lnln_pole_int_a_eq_c - 2*clog(a-f) * &
    dilog_term_corr_A(y, discont_flag_01, real((f+k_crossing)/(f-a)))
end function

