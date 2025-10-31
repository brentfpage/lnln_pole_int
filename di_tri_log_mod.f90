! math functions involving log, dilog, and trilog
module di_tri_log_mod

  interface clog
    module procedure clog_no_signed_0
  end interface

  contains
  ! if z=(real(z),-0.0), return log of z=(real(z),0.0)
  function clog_no_signed_0(z)
    implicit none
    complex :: z
    complex :: clog_no_signed_0

    if(aimag(z).ne.0.0) then
      clog_no_signed_0 = log(z)
    else
      clog_no_signed_0 = log(cmplx(real(z),0.0))
    endif
  end function clog_no_signed_0

!   return Li₂(z), i.e., the dilogarithm evaluated at z
  function dilog(z)
    implicit none
    complex :: z
    complex :: dilog
    dilog = spence(1-z)
  end function

  ! evaluate the discontinuity in dilog at dilog_arg ∈ [1,∞) as dilog_arg moves from upper to
  ! the lower half planes
  function dilog_discont(dilog_arg)
    implicit none
    complex :: dilog_discont
    real :: dilog_arg
    complex :: i
    real :: pi

    i = (0.0,1.0)
    pi = 4*atan(1.0)
    
    dilog_discont = 2*pi*i*log(dilog_arg)
  end function dilog_discont

! return Li₃(z), i.e., the trilogarithm evaluated at z
function trilog(z)
    implicit none

    complex :: z
    complex :: trilog
    ! Riemann zeta function ζ(n) at n=3
    complex, parameter :: zeta3 = (1.20205690315959428539973816151144999076498,0.0)
    real :: pi

    pi = atan(1.0) * 4.0

    if(z.eq.0.0) then
      trilog=(0.0,0.0)
    else if(z.eq.1.0) then
      trilog=zeta3
    else if (abs(z)<0.5) then ! avoid calling series1 if z^2-8z+1≈0 
      trilog = trilog_series0(z)
    else if (abs(z)<1.0) then
      trilog = trilog_series1(z)
    else if (abs(z)<2.0) then
      ! using 6.6 of lewin
      trilog = trilog_series1(z**(-1)) - 1.0/6 * clog(-z)**3 - 1.0/6 * pi**2 * clog(-z)
    else
      ! using 6.6 of lewin
      trilog = trilog_series0(z**(-1)) - 1.0/6 * clog(-z)**3 - 1.0/6 * pi**2 * clog(-z)
    endif
  end function trilog

! evaluate the trilogarithm of z by direct summation of the defining zeries for |z| < 1,
!           ∞
!           ⎲  zᵏ
! Li₃(z) =  ⎳  ── .
!          ᵏ⁼¹ kⁿ
! this function expects |z| < 0.5
  function trilog_series0(z)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    implicit none
    complex :: z
    complex :: trilog_series0
    complex :: zfac, term, res
    integer :: n, lub
    real :: TOL

    TOL = epsilon(real(z))
    lub = int(-1.1 * log(TOL)/log(2.0)) 

    res = (0.0,0.0)
    zfac = (1.0,0.0)

    do n=1,lub
      zfac = zfac * z
      term = zfac / n**3
      res = res + term
      if(abs(term).le.(TOL*abs(res))) then
        exit
      endif
    enddo
    trilog_series0 = res

    if(n.eq.(lub+1)) then
      write(stderr,*) 'trilog series 0 end reached'
    endif
  end function trilog_series0


  ! evaluate the trilogrithm of z using a similar approach as used in https://dl.acm.org/doi/10.1145/360715.360722 for the dilogarithm.  this function expects |z| < 1. also, this function likely performs badly for z^2-8z+1≈0.  judging by https://www.ams.org/journals/mcom/1979-33-146/S0025-5718-1979-0521291-X/home.html , this function also possibly performs badly for z ~ 0
  function trilog_series1(z)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    implicit none
    complex :: z
    complex :: trilog_series1
    complex :: term, res, zfac
    integer :: n, lub
    
    real :: TOL

    TOL = epsilon(real(z))
    lub=int(2.0 * (1./TOL)**(1./9.))

    zfac = (1.0,0.0)
    res = (0.0,0.0)

    do n=1,lub
      zfac = zfac * z
      term = ((zfac/n**3)/(n+1)**3)/(n+2)**3

      res = res + term
      if(abs(term).le.(TOL*abs(res))) then
        exit
      endif
    enddo

    if(n.eq.(lub+1)) then
      write(stderr,*) 'trilog series 1 lub reached'
    endif

    res = res * 8 * z**2
    res = res - 99./4*z**2 + 35./2*z &
      + 9./2*(z**2 - 1.)*dilog(z) &
      + 12*(z**2 - 2*z + 1)*clog(1-z)

    trilog_series1 = res / (z**2 - 8*z + 1)

  end function trilog_series1

  ! f(x) from 8.105 of Polylogarithms and Associated Functions - Leonard Lewin (1981, North Hooland)
  function lewin_f(x)
    implicit none

    complex :: x
    complex :: lewin_f

    lewin_f = clog(x)*clog(1-x)**2 + 2*clog(1-x)*dilog(1-x) - 2*trilog(1-x) + 2*trilog((1.0,0.0))

  end function lewin_f

  ! evaluate the discontinuity in
  !         t   2
  !        ⌠  ln (1-t′)
  ! f(t) ≡ ⎮  ───────── dt′
  !        ⌡      t′
  !         0
  ! across t=[1,∞)

  function lewin_f_discont(t_val)
    implicit none
    complex :: i
    real :: pi
    real :: t_val
    complex :: lewin_f_discont

    i = (0.0,1.0)
    pi=4*atan(1.0)

    lewin_f_discont = 4*pi*i*(log(abs(1-t_val))*log(t_val)+dilog(1-cmplx(t_val,0.0)))
  end function


! Copyright notice for functions spence, cseries_spence0.f90, and cseries_spence1.f90 below
! "
! Copyright © 2001, 2002 Enthought, Inc.
! All rights reserved.
! 
! Copyright © 2003-2019 SciPy Developers.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 
!     Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 
!     Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 
!     Neither the name of Enthought nor the names of the SciPy Developers may be used to endorse or promote products derived from this software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! "

  ! source https://github.com/scipy/scipy/blob/main/scipy/special/_spence.pxd
  !
  ! # Implement Spence's function, a.k.a. the dilogarithm, for complex
  ! # arguments. Note that our definition differs from that in the sources
  ! # by the mapping z -> 1 - z.
  ! #
  ! # Sources
  ! # [1] Zagier, "The Dilogarithm Function"
  ! # [2] functions.wolfram.com
  ! # [3] Ginsberg, Zaborowski, "The Dilogarithm Function of a Real Argument"
  ! #
  ! # Author: Josh Wilson
  ! #
  ! # Released under the same license as Scipy.
  !
  !           ⌠z  log(t)
  ! Spence(z)=⎮   ────── dt
  !           ⌡1   1-t
  !

  function spence(zin)
    implicit none
  complex :: zin
  complex :: z
  complex :: spence
  real :: pi, PISQ_6

  ! type(mp_complex) :: clog

  pi = atan(1.0)*4
  PISQ_6 = 1.64493406684822643647241516664602519

  !     """
  !     Compute Spence's function for complex arguments. The strategy is:
  !     - If z is close to 0, use a series centered at 0.
  !     - If z is far away from 1, use the reflection formula
  ! 
  !     spence(z) = -spence(z/(z - 1)) - pi**2/6 - ln(z - 1)**2/2
  ! 
  !     to move close to 1. See [1].
  !     - If z is close to 1, use a series centered at 1.
  ! 
  !     """

     z = zin 

      if (abs(z).lt.0.5) then
  !         # This step isn't necessary, but this series converges faster.
          spence = cspence_series0(z)
      else if (abs(1 - z).gt.1) then
          spence = -cspence_series1(z/(z - 1)) - PISQ_6 - 0.5*clog(z - 1.0)**2
      else
          spence = cspence_series1(z)
      endif

  end function spence

  ! source https://github.com/scipy/scipy/blob/main/scipy/special/_spence.pxd
  !
  !     """
  !     A series centered at z = 0; see
  ! 
  !     http://functions.wolfram.com/10.07.06.0005.02
  ! 
  !     """
  ! # Author: Josh Wilson
  ! #
  ! # Released under the same license as Scipy.

  ! 7/21/2025 : changed the number of terms in the series to get more precision than given by
  ! 64 bit floats; stopped using zlog1
  function cspence_series0(z)
    use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    implicit none
    real :: PISQ_6 
    real :: TOL 
    integer :: n, lub
    complex :: zfac
    complex :: sum1 
    complex :: sum2 
    complex :: term1, term2
    complex :: z
    complex :: cspence_series0
    real :: pi

    pi = atan(1.0)*4
    PISQ_6 = 1.64493406684822643647241516664602519
    TOL = epsilon(real(z))

    zfac = (1.0,0.0)
    sum1 = (0.0,0.0)
    sum2 = (0.0,0.0)

      ! expect convergence to -log10(TOL) digits before this loop upper bound is reached
    lub = int(-1.1 * log(TOL)/log(2.0) ) 
    if (z == 0.0) then
        cspence_series0 = PISQ_6
    else
        do n=1,lub
            zfac = z * zfac
            term1 = zfac/n**2
            sum1 = term1 + sum1
            term2 = zfac/n
            sum2 = term2 + sum2
            if ((abs(term1).le.(TOL*abs(sum1))).and.(abs(term2).le.(TOL*abs(sum2)))) then
                exit
            endif
        enddo
        if(n.eq.(lub+1)) then
          write(stderr,*) 'spence series 0 end reached'
        endif
        cspence_series0 = PISQ_6 - sum1 + clog(z)*sum2
    endif

  end function cspence_series0

  ! source https://github.com/scipy/scipy/blob/main/scipy/special/_spence.pxd
  ! """
  ! A series centered at z = 1 which enjoys faster convergence than
  ! the Taylor series. See [3]. The number of terms used comes from
  ! bounding the absolute tolerance at the edge of the radius of
  ! convergence where the sum is O(1).
  ! 
  ! """
  ! # Author: Josh Wilson
  ! #
  ! # Released under the same license as Scipy.

  ! 7/21/2025 : changed the number of terms in the series to get more precision than given by
  ! 64 bit floats; stopped using zlog1

  function cspence_series1(z)
      use, intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
      implicit none
      integer :: lub
      integer :: n
      complex :: zfac 
      complex :: res 
      complex :: term, z, zz, mzm1
      complex :: cspence_series1
      real :: TOL 


      TOL = epsilon(real(z))
  !     lub=int(1.1 * 10**(-log(TOL)/log(2.0)/6.0))
      lub=int(2 * (1./TOL)**(1./6.))

      zfac = (1.0,0.0)
      res = (0.0,0.0)

      if (z.eq.1) then
          cspence_series1 = (0.0,0.0)
      else
          mzm1 = 1.0 - z
          zz = mzm1**2
          do n=1,lub
              zfac = mzm1 * zfac
              term = ((zfac/n**2)/(n+1)**2)/(n+2)**2
              res = term + res
              if(abs(term).le.(TOL*abs(res))) then
                  exit
              endif
          enddo
          if(n.eq.(lub+1)) then
            write(stderr,*) 'spence_series1 lub reached, z = ', z
          endif
          res = 4.0*zz * res
          res = res + 4.0*mzm1 + 5.75*zz + 3.0*(1.0 - zz)*clog(z)
          res = res/(1.0 + 4.0*mzm1 + zz)
      endif
      cspence_series1 = res
  end function cspence_series1

end module
