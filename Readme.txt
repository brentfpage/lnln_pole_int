This program, lnln-pole-int, computes
    xub 
    ⌠    ln(a+x)ln(c+x) 
    ⎮    ────────────── dx
    ⌡         f+x 
    xlb 
for complex-valued parameters a, c, and f,
which arises in the course of evaluating
 xub         Re(a)                Re(c)
⌠        1  ⌠             1      ⌠             1
⎮    dx ─── ⎮      dy─────────── ⎮      dz─────────── .
⌡       f+x ⌡        x+y+i*Im(a) ⌡        x+z+i*Im(c)
 xlb

Details of the computation are given in lnln_pole_int_doc.pdf

Known limitation: 
    xub 
    ⌠    ln(a+x)ln(f+1+x) 
    ⎮    ──────────────── dx
    ⌡         f+x 
    xlb 
with f pure real and xlb<-f<xub is well-defined, but lnln-pole-int at present is not able to compute 
such integrals.

_________________________________________________________________________________

Installation instructions:
    Requirements: A compiler such as gfortran
    (https://fortran-lang.org/learn/os_setup/install_gfortran/)

    Download: click the green button labeled "<> Code v" on Github and then click "Download
    ZIP".  Once the zip file is downloaded, double-click it to unzip it, which will produce a
    directory lnln_pole_int-main.

    Compile: On a computer with gfortran, navigate to the directory lnln_pole_int-main on the
    command line, type 'make', and hit enter.  On a computer with a different fortran compiler, it
    will be necessary to first edit the file 'makefile'.

Usage: navigate to the directory lnln_pole_int-main on the command line and call the program as
    ./eval_int Re(a) Im(a) Re(c) Im(c) Re(f) Im(f) klb kub
The parameters that have been read in will be printed back out, followed by the result for the
integral as "lnln_pole_int = (Re(result),Im(result))".

e.g., 
    ./eval_int 3.17 -8.3234 9.209 2.3588 9.233 -4.814 1.543 9.777
         a =               (3.1699999999999999,-8.3233999999999995) 
         c =               (9.2089999999999996,2.3588000000000000) 
         f =               (9.2330000000000005,-4.8140000000000001) 
         klb =    1.5429999999999999 
         kub =    9.7769999999999992 
         lnln_pole_int = (3.7474535654559897,0.28703781175992304)




