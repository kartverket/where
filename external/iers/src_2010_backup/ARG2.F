      SUBROUTINE ARG2 ( IYEAR, DAY, ANGLE ) 
*+
*  - - - - - - - - -
*   A R G 2 
*  - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  The purpose of the subroutine is to compute the angular astronomical
*  argument, which depends on time, for 11 tidal argument calculations. 
*  The order of the 11 angular quantities in vector angle are given below:
*
*  01-M2, 02-S2, 03-N2, 04-K2, 05-K1, 06-O1, 07-P1, 08-Q1, 09-Mf,
*  10-Mm, 11-Ssa (See Reference 1) 
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Canonical model	
* 
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     IYEAR          i      Four digit year (Note 1)
*     DAY            d      Day of Year Greenwich Time (Note 2)
*
*  Returned:
*     ANGLE(K)       d      Angular argument for Schwiderski
*                           computation, in radians (Notes 3 and 4)
*
*  Notes:
*
*  1) This subroutine is valid only after 1973 CE.  A validation
*     test has been added to stop the subroutine if an invalid
*     year is used.
* 
*  2) Example: 32.5 for February 1 12 Noon
*     Example: 1.25 for January 1 6 AM
* 
*  3) Ocean loading phases computed from Schwiderski's models
*     refer to the phase of the associated solid Earth tide 
*     generating potential at the zero meridian according to 
*
*      OL_DR = OL_AMP ' COS (SE_PHASE" - OL_PHASE)
*
*     where OL = OCEAN LOADING TIDE,
*           SE = SOLID EARTH TIDE GENERATING POTENTIAL.
*
*     If the harmonic tide development of Cartwright, et al. 
*     (CTE) (1971, 1973) is used, make sure that SE_PHASE"
*     take into account: 
*
*     (1) the sign of SE_AMP in the tables of Cartwright et al.
*     (2) that CTE'S SE_PHASE refers to a sine rather than a 
*     cosine function if (N+M) = (DEGREE + ORDER) of the tide
*     spherical harmonic is odd.
*
*     i.e. SE_PHASE" = TAU(T) ' N1 + S(T) ' N2 + H(T) ' N3
*                 + P(T) ' N4 + N'(T) ' N5 + PS(T) ' N6
*                 + PI   If CTE'S amplitude coefficient < 0
*                 + PI/2 If (DEGREE + N1) is odd 
*
*     where TAU ... PS = astronomical arguments,
*           N1 ... N6 = CTE'S argument numbers.
*
*     Most tide generating software compute SE_PHASE" (for use
*     with cosines).
*
*  4) The double precision change from the original routine ARG.f
*     to ARG2.F yields output differences on the order of 10^-9 radians.
*
*  Called:
*     None
*
*  Test case:
*     given input: IYEAR = 2008 
*                  DAY = 311.5 (November 6 Noon)
*     expected output: ANGLE(1)  = 2.849663065753787805D0  rad
*                      ANGLE(2)  = 6.28318080000000023D0   rad
*                      ANGLE(3)  = 4.926040134021299366D0  rad
*                      ANGLE(4)  = 1.608450491115348768D0  rad
*                      ANGLE(5)  = 2.375021572352622456D0  rad
*                      ANGLE(6)  = 0.4746414933980958040D0 rad
*                      ANGLE(7)  = 3.908159227647345801D0  rad
*                      ANGLE(8)  = 2.551018561669245344D0  rad
*                      ANGLE(9)  = 5.041990012540757959D0  rad 
*                      ANGLE(10) = 4.206816878908014701D0  rad 
*                      ANGLE(11) = 1.608463638294885811D0  rad 
*
*  References:
*
*     Schwiderski, E., 1983, "Atlas of Ocean Tidal Charts and Maps, Part I:
*     The Semidiurnal Principal Lunar Tide M2," Marine Geodesy, 6, pp. 219-256. 
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2008 November  05 B.E. Stetzler    Added header and copyright
*  2008 November  07 B.E. Stetzler    Provided test case
*  2008 November  13 B.E. Stetzler    Re-defined variables 
*  2008 November  14 B.E. Stetzler    Corrected test case results 
*  2010 April     20 B.E. Stetzler    New version of routine without the
*                                     use of EQUIVALENCE
*  2011 September 20 B.E. Stetzler    Changed IYEAR variable from two
*                                     digit value to four digit value,
*                                     added year validation and IYMIN
*                                     variable, and updated test case
*  2011 October   07 B.E. Stetzler    Removed IYEAR-1900 from code and
*                                     modified ICAPD calculation 
*-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER I, ID, ICAPD, K, IYEAR, IYMIN
      PARAMETER ( K = 11 )
      PARAMETER ( IYMIN = 1974 )
      DOUBLE PRECISION ANGFAC(4,11), ANGLE(11), SPEED(11), DAY
      DOUBLE PRECISION FDAY, CAPT, H0, S0, P0, DTR
      DOUBLE PRECISION SIGM2, SIGS2, SIGN2, SIGK2, SIGK1, SIGO1,
     .                 SIGP1, SIGQ1, SIGMF, SIGMM, SIGSSA, TWOPI

*  ----------------------------------------------
*  Speed of all terms given in radians per second
*  ----------------------------------------------
      DATA SPEED(1),SIGM2/1.40519D-4,1.40519D-4/
      DATA SPEED(2),SIGS2/1.45444D-4,1.45444D-4/
      DATA SPEED(3),SIGN2/1.37880D-4,1.37880D-4/
      DATA SPEED(4),SIGK2/1.45842D-4,1.45842D-4/
      DATA SPEED(5),SIGK1/.72921D-4,.72921D-4/
      DATA SPEED(6),SIGO1/.67598D-4,.67598D-4/
      DATA SPEED(7),SIGP1/0.72523D-4,0.72523D-4/
      DATA SPEED(8),SIGQ1/0.64959D-4,0.64959D-4/
      DATA SPEED(9),SIGMF/0.053234D-4,0.053234D-4/
      DATA SPEED(10),SIGMM/0.026392D-4,0.026392D-4/
      DATA SPEED(11),SIGSSA/0.003982D-4,0.003982D-4/

      DATA ANGFAC/2.D0,-2.D0,0.D0,0.D0,4*0.D0,
     .            2.D0,-3.D0,1.D0,0.D0,2.D0,3*0.D0,
     .            1.D0,2*0.D0,.25D0,1.D0,-2.D0,0.D0,-.25D0,
     .            -1.D0,2*0.D0,-.25D0,1.D0,-3.D0,1.D0,-.25D0,
     .            0.D0,2.D0,2*0.D0,0.D0,1.D0,-1.D0,0.D0,
     .            2.D0,3*0.D0/
      PARAMETER (TWOPI = 6.283185307179586476925287D0)
      DATA DTR/0.174532925199D-1/

*  Validate year

      IF ( IYEAR.LT.IYMIN ) STOP 'Year must be after 1973 CE.'

*  Initialize day of year

      ID = DAY
* ------------------------------------------
*  Compute fractional part of day in seconds 
* ------------------------------------------
      FDAY = (DAY-ID)*86400D0
* Revision 07 October 2011: ICAPD modified 
      ICAPD = ID+365*(IYEAR-1975)+((IYEAR-1973)/4)
      CAPT = (27392.500528D0+1.000000035D0*ICAPD)/36525D0
* --------------------------------------------------
*  Compute mean longitude of Sun at beginning of day
* --------------------------------------------------
      H0=(279.69668D0+(36000.768930485D0+3.03D-4*CAPT)*CAPT)*DTR
* ---------------------------------------------------
*  Compute mean longitude of Moon at beginning of day 
* ---------------------------------------------------
      S0=(((1.9D-6*CAPT-.001133D0)*CAPT+481267.88314137D0)*CAPT
     .  +270.434358D0)*DTR
* ------------------------------------------------------------
*  Compute mean longitude of lunar perigee at beginning of day 
* ------------------------------------------------------------
      P0=(((-1.2D-5*CAPT-.010325D0)*CAPT+4069.0340329577D0)*CAPT
     .  +334.329653D0)*DTR

* Compute the tidal angle arguments

      DO 500 I = 1, K 
         ANGLE(I) = SPEED(I)*FDAY + ANGFAC(1,I)*H0 + ANGFAC(2,I)*S0
     .   + ANGFAC(3,I)*P0 + ANGFAC(4,I)*TWOPI
         ANGLE(I) = DMOD(ANGLE(I),TWOPI)
         IF (ANGLE(I).LT.0D0) ANGLE(I) = ANGLE(I) + TWOPI
  500 CONTINUE
      RETURN

* Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END
