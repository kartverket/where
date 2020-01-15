      SUBROUTINE APG (DLAT,DLON,AZ,EL,D,GRN,GRE)
*+
*  - - - -
*   A P G
*  - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine determines the asymmetric delay d in meters caused
*  by gradients.  The north and east gradients are also provided.
*  They are based on Spherical Harmonics up to degree and order 9.
*  If the north and east gradients are used, they should be used
*  with the gradient model by Chen and Herring (1997). See Reference 1
*  and last lines of this subroutine.
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status:  Class 1 model
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
*     DLAT           d      Latitude given in radians (North Latitude)
*     DLON           d      Longitude given in radians (East Longitude)
*     AZ             d      Azimuth from north in radians
*     EL             d      Elevation angle in radians
*
*  Returned:
*     D              d      Delay in meters
*     GRN            d      North gradient in mm
*     GRE            d      East gradient in mm
*
*  Notes:
*
*  1) This a priori model cannot replace the (additional) estimation of
*     gradient parameters, if observations at elevation angles below
*     15 degrees are analyzed.  
*
*  Test case:
*     Kashima 11 Station information retrieved at:
*     ftp://ivscc.gsfc.nasa.gov/pub/config/ns/kashim11.config.txt
*
*     given input: DLAT = 0.6274877539940092D0 radians (KASHIMA 11, Japan)
*                  DLON = 2.454994088489240D0 radians
*                  AZ   = 0.2617993877991494D0 radians
*                  EL   = 0.8726646259971648D0 radians
*
*     expected output: D   = -0.9677190006296187757D-4 meters 
*                      GRN = -0.1042668498001996791D0 mm
*                      GRE = 0.4662515377110782594D-1 mm
*
*  References:
*
*     Chen, G. and Herring, T. A., 1997, ``Effects of atmospheric azimuthal
*     asymmetry on the analysis of space geodetic data,"
*     J. Geophys. Res., 102(B9), pp. 20,489--20,502, doi: 10.1029/97JB01739.
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2010 May 18       J. Boehm      Original code
*  2010 September 21 B.E. Stetzler Added header and copyright
*  2010 September 22 B.E. Stetzler Declared pi, degree, and order
*                                  as parameters
*  2010 September 22 B.E. Stetzler Provided test case
*  2010 September 29 B.E. Stetzler Chose valid elevation angle and 
*                                  corrected test case
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION     DLAT,DLON,AZ,EL,D,GRN,GRE

      DOUBLE PRECISION     V(10,10),W(10,10)
      DOUBLE PRECISION     A_N(55),B_N(55),A_E(55),B_E(55)

      INTEGER*4  I,N,M,NMAX,MMAX
      DOUBLE PRECISION     PI,X,Y,Z

      PARAMETER ( PI = 3.1415926535897932384626433D0 )    

*     degree n and order m
      PARAMETER ( NMAX = 9, MMAX = 9 )

      DATA (A_N(I),I=1,55)/ 
     .  2.8959D-02,-4.6440D-01,-8.6531D-03, 1.1836D-01,-2.4168D-02,
     . -6.9072D-05, 2.6783D-01,-1.1697D-03,-2.3396D-03,-1.6206D-03,
     . -7.4883D-02, 1.3583D-02, 1.7750D-03, 3.2496D-04, 8.8051D-05,
     .  9.6532D-02, 1.3192D-02, 5.5250D-04, 4.0507D-04,-5.4758D-06,
     .  9.4260D-06,-1.0872D-01, 5.7551D-03, 5.3986D-05,-2.3753D-04,
     . -3.8241D-05, 1.7377D-06,-4.4135D-08, 2.1863D-01, 2.0228D-02,
     . -2.0127D-04,-3.3669D-04, 8.7575D-06, 7.0461D-07,-4.0001D-08,
     . -4.5911D-08,-3.1945D-03,-5.1369D-03, 3.0684D-04, 2.4459D-05,
     .  7.6575D-06,-5.5319D-07, 3.5133D-08, 1.1074D-08, 3.4623D-09,
     . -1.5845D-01,-2.0376D-02,-4.0081D-04, 2.2062D-04,-7.9179D-06,
     . -1.6441D-07,-5.0004D-08, 8.0689D-10,-2.3813D-10,-2.4483D-10/

      DATA (B_N(I),I=1,55)/ 
     .  0.0000D+00, 0.0000D+00,-1.1930D-02, 0.0000D+00, 9.8349D-03,
     . -1.6861D-03, 0.0000D+00, 4.3338D-03, 6.1707D-03, 7.4635D-04,
     .  0.0000D+00, 3.5124D-03, 2.1967D-03, 4.2029D-04, 2.4476D-06,
     .  0.0000D+00, 4.1373D-04,-2.3281D-03, 2.7382D-04,-8.5220D-05,
     .  1.4204D-05, 0.0000D+00,-8.0076D-03, 4.5587D-05,-5.8053D-05,
     . -1.1021D-05, 7.2338D-07,-1.9827D-07, 0.0000D+00,-3.9229D-03,
     . -4.0697D-04,-1.6992D-04, 5.4705D-06,-4.4594D-06, 2.0121D-07,
     . -7.7840D-08, 0.0000D+00,-3.2916D-03,-1.2302D-03,-6.5735D-06,
     . -3.1840D-06,-8.9836D-07, 1.1870D-07,-5.8781D-09,-2.9124D-09,
     .  0.0000D+00, 1.0759D-02,-6.6074D-05,-4.0635D-05, 8.7141D-06,
     .  6.4567D-07,-4.4684D-08,-5.0293D-11, 2.7723D-10, 1.6903D-10/

      DATA (A_E(I),I=1,55)/ 
     . -2.4104D-03, 1.1408D-04,-3.4621D-04, 1.6565D-03,-4.0620D-03,
     . -6.8424D-03,-3.3718D-04, 7.3857D-03,-1.3324D-03,-1.5645D-03,
     .  4.6444D-03, 1.0296D-03, 3.6253D-03, 4.0329D-04, 3.1943D-04,
     . -7.1992D-04, 4.8706D-03, 9.4300D-04, 2.0765D-04,-5.0987D-06,
     . -7.1741D-06,-1.3131D-02, 2.9099D-04,-2.2509D-04, 2.6716D-04,
     . -8.1815D-05, 8.4297D-06,-9.2378D-07,-5.8095D-04, 2.7501D-03,
     .  4.3659D-04,-8.2990D-06,-1.4808D-05, 2.2033D-06,-3.3215D-07,
     .  2.8858D-08, 9.9968D-03, 4.9291D-04, 3.3739D-05, 2.4696D-06,
     . -8.1749D-06,-9.0052D-07, 2.0153D-07,-1.0271D-08, 1.8249D-09,
     .  3.0578D-03, 1.1229D-03,-1.9977D-04, 4.4581D-06,-7.6921D-06,
     . -2.8308D-07, 1.0305D-07,-6.9026D-09, 1.5523D-10,-1.0395D-10/

      DATA (B_E(I),I=1,55)/ 
     .  0.0000D+00, 0.0000D+00,-2.5396D-03, 0.0000D+00, 9.2146D-03,
     . -7.5836D-03, 0.0000D+00, 1.2765D-02,-1.1436D-03, 1.7909D-04,
     .  0.0000D+00, 2.9318D-03,-6.8541D-04, 9.5775D-04, 2.4596D-05,
     .  0.0000D+00, 3.5662D-03,-1.3949D-03,-3.4597D-04,-5.8236D-05,
     .  5.6956D-06, 0.0000D+00,-5.0164D-04,-6.5585D-04, 1.1134D-05,
     .  2.3315D-05,-4.0521D-06,-4.1747D-07, 0.0000D+00, 5.1650D-04,
     . -1.0483D-03, 5.8109D-06, 1.6406D-05,-1.6261D-06, 6.2992D-07,
     .  1.3134D-08, 0.0000D+00,-6.1449D-03,-3.2511D-04, 1.7646D-04,
     .  7.5326D-06,-1.1946D-06, 5.1217D-08, 2.4618D-08, 3.6290D-09,
     .  0.0000D+00, 3.6769D-03,-9.7683D-04,-3.2096D-07, 1.3860D-06,
     . -6.2832D-09, 2.6918D-09, 2.5705D-09,-2.4401D-09,-3.7917D-11/

*     unit vector
      X = DCOS(DLAT)*DCOS(DLON)
      Y = DCOS(DLAT)*DSIN(DLON)
      Z = DSIN(DLAT)
  
*     Legendre polynomials
      V(1,1) = 1.0D0
      W(1,1) = 0.0D0
      V(2,1) = Z * V(1,1)
      W(2,1) = 0.0D0

      DO N=2,NMAX
        V(N+1,1) = ((2*N-1) * Z * V(N,1) - (N-1) * V(N-1,1)) / N
        W(N+1,1) = 0.0D0
      ENDDO

      DO M=1,MMAX
        V(M+1,M+1) = (2*M-1) * (X*V(M,M) - Y*W(M,M))
        W(M+1,M+1) = (2*M-1) * (X*W(M,M) + Y*V(M,M))
        IF (M < MMAX) THEN
          V(M+2,M+1) = (2*M+1) * Z * V(M+1,M+1)
          W(M+2,M+1) = (2*M+1) * Z * W(M+1,M+1)
        ENDIF
        DO N=M+2,NMAX
          V(N+1,M+1) = ((2*N-1)*Z*V(N,M+1) - (N+M-1)*V(N-1,M+1)) / (N-M)
          W(N+1,M+1) = ((2*N-1)*Z*W(N,M+1) - (N+M-1)*W(N-1,M+1)) / (N-M)
        ENDDO
      ENDDO

*     Surface pressure on the geoid
      GRN = 0D0
      GRE = 0D0
      I = 0
      DO N=0,NMAX
        DO M=0,N
          I = I+1
          GRN = GRN + (A_N(I)*V(N+1,M+1) + B_N(I)*W(N+1,M+1))
          GRE = GRE + (A_E(I)*V(N+1,M+1) + B_E(I)*W(N+1,M+1))
        ENDDO
      ENDDO

*     calculation of the asymmetric delay in m (Chen and Herring 1997)
      D = 1.D0/(DSIN(EL)*DTAN(EL)+0.0031)*(GRN*DCOS(AZ)+GRE*DSIN(AZ)) ! mm
      D = D/1000.D0 ! m 
      
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
