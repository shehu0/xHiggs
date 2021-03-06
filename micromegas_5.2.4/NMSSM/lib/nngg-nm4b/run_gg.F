* run.F
* defines parameter settings for one run, then includes num.F
* this file is part of FormCalc
* last modified 14 Apr 04 th


* Whether to run in debugging mode

#define DEBUG


* Uncomment the following to check UV- and IR-finiteness

c#define DELTA 1D10
c#define MUDIM 1D100
c#define LAMBDA 1D10

* The following sets all ren. constants to zero for debugging

c#define CalcRenConst ZeroRenConst


* The model initialization code is included here.  Options for the
* XXX_ini.F file (e.g. NO_SQUARK_MIXING for mssm_ini.F) must be
* defined before the #include line.

* options for mssm_ini.F:
c#define NO_SQUARK_MIXING
c#define COMPLEX_PARAMETERS
c#define SM_ONLY
c#define NO_EXCLUSION_LIMITS
c#define TREE_LEVEL_MH

* options for thdm_ini.F:
c#define MODEL_TYPE_I
c#define MODEL_TYPE_II
c#define NO_EXCLUSION_LIMITS


* Cuts and increments.  May be defined in terms of sqrtS, the CMS energy.

#define THETACUT (0*degree)
#ifdef DEBUG
#define THETASTEP (125*degree)
#else
#define THETASTEP (1*degree)
#endif

#define K50CUT ESOFTMAX
#define K50STEP .05D0*sqrtS

#define SQRTSSTEP 10


* The integration parameters:
* - the desired relative and absolute accuracy

#define RELACCURACY 1D-3
#define ABSACCURACY 1D-10

* The multidimensional integration routines (2 -> 3) also need

#define METHOD DIVONNE
#define VERBOSE 1
#define MINEVAL 0
#define MAXEVAL 50000

* for Vegas:
#define NSTART 1000
#define NINCREASE 500

* for Suave:
#define NNEW 1000
#define FLATNESS 50

* for Divonne:
#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 1D-6
#define MAXCHISQ 10
#define MINDEVIATION .25D0

* for Cuhre:
#define KEY 0

* Now go for it:


* main.F
* the main program for the Fortran code generated by FormCalc
* this file is part of FormCalc
* last modified 15 Jul 04 th


#define SCALAR 0
#define FERMION 1
#define PHOTON 2
#define GLUON PHOTON
#define VECTOR 3

* `real' data lines are tagged with this character:
#define TAG "|"

* process.h
* defines all process-dependent parameters
* this file is part of FormCalc
* last modified 15 Jan 03 th

* Definition of the external particles.
* The TYPEn may be one of SCALAR, FERMION, PHOTON, or VECTOR.
* (PHOTON is equivalent to VECTOR, except that longitudinal
* modes are not allowed)

* Note: The initial definitions for particles 2...5 are of course
* sample entries for demonstration purposes.

#define TYPE1 FERMION
#define MASS1 MNE1
#define CHARGE1 0

#define TYPE2 FERMION
#define MASS2 MNE1
#define CHARGE2 0

#define TYPE3 PHOTON
#define MASS3 0
#define CHARGE3 0

#define TYPE4 PHOTON
#define MASS4 0
#define CHARGE4 0


#define TYPE5 PHOTON
#define MASS5 0
#define CHARGE5 0

* The combinatorial factor for identical particles in the final state:
* 1/n! for n identical particles, 1 otherwise

#define IDENTICALFACTOR 5D-1

* Possibly a colour factor, e.g.
* - an additional averaging factor if any of the incoming particles
*   carry colour,
* - the overall colour factor resulting from the external particles
*   if that cannot computed by FormCalc (e.g. if the model has no
*   colour indices, as SM.mod).

#define COLOURFACTOR 1

* Whether to include soft-photon bremsstrahlung.
* ESOFTMAX is the maximum energy a soft photon may have and may be
* defined in terms of sqrtS, the CMS energy.

#define ESOFTMAX .1D0*sqrtS

* Possibly some wave-function renormalization
* (e.g. if calculating in the background-field method)

c#define WF_RENORMALIZATION (nW*dWFW1 + nZ*dWFZ1)

#define NCOMP 2

* Include the kinematics-dependent part of the code


* 2to2.F
* the kinematics for a 2 -> 2 process
* this file is part of FormCalc
* last modified 8 Jul 04 th


* some definitions for main.F
#define LEGS 4
#define LEGS_IN 2
#define TYPES TYPE1, TYPE2, TYPE3, TYPE4

#define HEADER "2to2.h"

************************************************************************
** SetProcess sets up process-dependent constants.

	subroutine SetProcess(hel, avgfac)
	implicit none
	integer hel
	double precision avgfac

#include HEADER

	helicities = hel
	preflux = hbar_c2/(32*pi)*avgfac

	end


************************************************************************
** SetEnergy sets up energy-dependent constants.

	double precision function ThreeMom(sqrtS, ma, mb)
	implicit none
	double precision sqrtS, ma, mb

	ThreeMom = sqrt(.25D0*(sqrtS - (ma**2 - mb**2)/sqrtS)**2 -
     &    mb**2)
	end


	subroutine SetEnergy(v, *)
	implicit none
	double precision v, newsqrtS

#include HEADER
#include "util.h"

	double precision threshold, Pin, ThreeMom
	external ThreeMom

        newsqrtS=MNE1*(2+v*v/4.)
	sqrtS = newsqrtS
	Pin = ThreeMom(sqrtS, dble(MASS1), dble(MASS2))
	Pout = ThreeMom(sqrtS, dble(MASS3), dble(MASS4))
        
        if(v.gt.0.0001) then
	flux = preflux*Pout/(Pin*sqrtS**2)*v
	else
	flux = 4.0*preflux*Pout/(sqrtS**3)
	endif
c        write(*,*) Pout/Pin, 2.0/v
	call VecSet(1, dble(MASS1), Pin, 0D0, 0D0, 1D0)
	call VecSet(2, dble(MASS2), Pin, 0D0, 0D0, -1D0)

	call CalcRenConst
	reset = .TRUE.
	return

999	return 1
	end


************************************************************************
** DifferentialCS actually invokes the generated function SquaredME.
** If necessary, soft-photon bremsstrahlung corrections are added.

	subroutine DifferentialCS(result, costh)
	implicit none
	double precision result(*), costh

#include HEADER
#include "util.h"
#include "types.h"
#include "renconst.h"


	double precision sinth
	integer c

	sinth = sqrt(1 - costh**2)

	call VecSet(3, dble(MASS3), Pout, sinth, 0D0, costh)
	call VecSet(4, dble(MASS4), Pout, -sinth, 0D0, -costh)

	call SquaredME(result, helicities, reset)

	do c = 1, NCOMP
	  result(c) = result(c)*flux
	enddo


#ifdef WF_RENORMALIZATION
	result(2) = result(2) + (WF_RENORMALIZATION)*result(1)
#endif
	end



       integer function addhel(h)
       implicit none

       integer h
       
       if(h.eq.-1) addhel=2
       if(h.eq. 1) addhel=8
       if(h.eq. 0) addhel=4

       end

c**********************************************

       subroutine nnaaini(h1,h2,h3,h4,v)
       implicit none

#include "looptools.h"
#include "model.h"
       integer hel, h1,h2,h3,h4,addhel
       double precision avgfac, sqrtSfrom, sqrtSto,v
c       hel=2925       
       hel=addhel(h4)+32*(addhel(h3)+32*(addhel(h2)+32*addhel(h1)))
       avgfac=0.5D0
       if(h1.eq.0) avgfac=avgfac/2
       if(h2.eq.0) avgfac=avgfac/2

       call SetProcess(hel,avgfac)
       call ltini
C	  call ModelConstIni(*3)

#ifdef MUDIM
	  call setmudim(dble(MUDIM))
#endif
#ifdef DELTA
	  call setdelta(dble(DELTA))
#endif
#ifdef LAMBDA
	  call setlambda(dble(LAMBDA))
#endif
	  Divergence = getdelta()
           call SetEnergy(v, *3)
c	  call ModelDigest
3      return
       end
       

