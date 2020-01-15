      SUBROUTINE GCONV2 (A,F,X,Y,Z,PHI,LAMBDA,H)
*+
*  - - - - - - - - - - -
*   G C O N V 2
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine transforms from geocentric rectangular to 
*  geodetic coordinates.
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status:  Class 3 model
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as
*     a Class 1, 2, or 3 model.
*
*  Given:
*     a           d         Equatorial Radius of the Earth (Note 1)
*     f           d         Flattening form factor (Note 2)
*     x           d         Rectangular X coordinate (Note 3)
*     y           d         Rectangular Y coordinate (Note 3)
*     z           d         Rectangular Z coordinate (Note 3)
*
*  Returned:
*     phi         d         Latitude coordinate on the ellipsoid (Note 4)
*     lambda      d         Longitude coordinate on the ellipsoid (Note 4)
*     h           d         Height coordinate on the ellipsoid (Note 3)
*
*  Notes:
*
*  1) The parameter given is from the 1980 Geodetic Reference System, which
*     was adopted at the XVII General Assembly of the International Union
*     of Geodesy and Geophysics (IUGG).  It is expressed in meters.
*
*  2) The parameter given is from the 1980 Geodetic Reference System, which
*     was adopted at the XVII General Assembly of the International Union
*     of Geodesy and Geophysics (IUGG).  It is a dimensionless quantity.
*
*  3) Coordinates are expressed in meters.
*
*  4) Coordinates are expressed in radians.
*
*  5) This routine is closely based on the GCONV2H subroutine by
*     Toshio Fukushima (See reference 1).
*
*  6) This version of the routine uses the GRS 1980 ellipsoid parameters
*     as the given default. The user may choose to use other ellipsoids as input parameters.
*
*  Test case:
*     given input: x = 4075579.496D0 meters  Wettzell (TIGO) station
*                  y =  931853.192D0 meters
*                  z = 4801569.002D0 meters
*
*     expected output: phi    =   0.857728298603D0 radians
*                      lambda =   0.224779294628D0 radians
*                      h      = 665.9207D0 meters
*
*  References:
*
*     Fukushima, T., "Transformation from Cartesian to geodetic
*     coordinates accelerated by Halley's method", J. Geodesy (2006),
*     79(12): 689-693
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*     2006              T.   Fukushima             Original code
*     2010 March 19     B.E. Stetzler              Added header and copyright
*     2010 March 19     B.E. Stetzler              Initial standardization
*                                                  of code, capitalized
*                                                  variables for FORTRAN
*                                                  77 backwards compatibility
*     2010 March 22      B.E. Stetzler             Provided test case
*     2010 September 2   B.E. Stetzler             Corrected F to match Table
*                                                  1.2 of IERS Conventions
*                                                  (2010) and updated test case
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION A,F,X,Y,Z,PHI,LAMBDA,H,AEPS2,E2,E4T,EP,EP2,AEP,
     .                 P2,ABSZ,P,S0,PN,ZP,C0,C02,C03,S02,S03,A02,
     .                 A0,A03,D0,F0,B0,S1,CP,S12,CP2

*  Pi
      DOUBLE PRECISION DPI
      PARAMETER ( DPI = 3.141592653589793238462643D0 )

*  -------------
*  Preliminaries
*  -------------

*	Default Values: GRS1980 System (recommended in IERS Conventions)

         A = 6378137D0
         F = 0.00335281068118D0

*  Functions of ellipsoid parameters

         AEPS2 = A*A*1D-32
         E2    = (2.0D0-F)*F
         E4T   = E2*E2*1.5D0
         EP2   = 1.0D0-E2
         EP    = SQRT(EP2)
         AEP   = A*EP

* ---------------------------------------------------------
*	Compute Coefficients of (Modified) Quartic Equation
*
*	Remark: Coefficients are rescaled by dividing by 'a'
* ---------------------------------------------------------

*      Compute distance from polar axis squared

       P2 = X*X + Y*Y

*      Compute longitude lambda 

       IF ( P2.NE.0D0 ) THEN
          LAMBDA = ATAN2(Y,X)
       ELSE
          LAMBDA = 0D0
       END IF

*      Ensure that Z-coordinate is unsigned
    
       ABSZ = ABS(Z)

*      Continue unless at the poles

       IF ( P2.GT.AEPS2 ) THEN

*      Compute distance from polar axis

         P = SQRT(P2)

*      Normalize
       
         S0 = ABSZ/A
         PN = P/A
         ZP = EP*S0

*     Prepare Newton correction factors.

         C0  = EP*PN
         C02 = C0*C0
         C03 = C02*C0
         S02 = S0*S0
         S03 = S02*S0
         A02 = C02+S02
         A0  = SQRT(A02)
         A03 = A02*A0
         D0  = ZP*A03 + E2*S03
         F0  = PN*A03 - E2*C03

*     Prepare Halley correction factor.

         B0 = E4T*S02*C02*PN*(A0-EP)
         S1 = D0*F0 - B0*S0
         CP = EP*(F0*F0-B0*C0)

*     Evaluate latitude and height.

         PHI = ATAN(S1/CP)
         S12 = S1*S1
         CP2 = CP*CP
         H   = (P*CP+ABSZ*S1-A*SQRT(EP2*S12+CP2))/SQRT(S12+CP2)
      ELSE

*     Special case: pole.

         PHI = DPI/2D0
         H   = ABSZ - AEP
      END IF

*  Restore sign of latitude.
      IF ( Z.LT.0D0 ) PHI = -PHI

      RETURN

*  Finished.

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

