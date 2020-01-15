      SUBROUTINE FCNNUT ( MJD,X,Y,DX,DY )
*+
*  - - - - - - - - - - -
*   F C N N U T
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine computes the effects of the free core nutation.  
*  Please note that the table is updated each year (see Note 4).
*  The parameter N needs to be incremented for each additional
*  year in the table.
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Class 3 model
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
*     mjd           d      Modified Julian Date, TDB (Note 1)
*
*  Returned:
*     X             d      CIP offset x component, in microas (Note 2)
*     Y             d      CIP offset y component, in microas (Note 2)
*     dX            d      Uncertainty of x component, in microas (Note 3)
*     dY            d      Uncertainty of y component, in microas (Note 3) 
*
*  Notes:
*
*  1) Though the Modified Julian Date (MJD) is strictly TDB, it is
*     usually more convenient to use TT, which makes no significant
*     difference.
*
*  1) CIP is the Celestial Intermediate Pole.  The expression
*     used is given in microarcseconds.
*  
*  2) The expression used is given in microarcseconds.
*
*  3) The updated table is maintained at the website
*     http://syrte.obspm.fr/~lambert/fcn/.
*
*  Test case: (NOT UPDATED FOR 2013 TABLE, TO BE DONE)
*     given input: MJD = 54790D0   Modified Julian Date, TDB
*                  
*     expected output:  X = -176.8012290066270680D0 microarcseconds
*                       Y = -93.51855308903756736D0 microarcseconds
*                       dX = 3.745573770491803067D0 microarcseconds
*                       dY = 3.745573770491803067D0 microarcseconds
*
*  References:
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2007 August   27 S. Lambert    Original code
*  2008 November 19 B.E.Stetzler  Added header and copyright
*  2008 November 20 B.E.Stetzler  provided a test case 
*  2009 July     10 B.E.Stetzler  Updated parameter N to 26, 
*                                 updated table to 2009, and corrected
*                                 test case results
*  2009 August   18 B.E.Stetzler  Capitalized all variables for FORTRAN
*                                 77 compatibility
*  2010 July      2 B.E.Stetzler  Updated parameter N to 27, 
*                                 updated table to 2010, and corrected
*                                 test case results
*  2011 July     13 B.E.Stetzler  Updated parameter N to 28, 
*                                 updated table to 2011, and corrected
*                                 test case results
*  2012 August   10 B.E.Stetzler  Updated parameter N to 29, 
*                                 updated table to 2012, and corrected
*                                 test case results
*  2013 December 19 G. Petit      Updated parameter N to 30, 
*                                 updated table to 2013
*                                 test case results NOT UPDATED, TO BE DONE
*-----------------------------------------------------------------------

        IMPLICIT NONE
        
        DOUBLE PRECISION MJD

        DOUBLE PRECISION X, Y    ! FCN contributions
        DOUBLE PRECISION DX, DY  ! uncertainties on X, Y
        
*       Internal variables

        INTEGER I,J,N
        PARAMETER (N=30)
        DOUBLE PRECISION PI,PER,PHI,MPE
        DOUBLE PRECISION AXC,AXS,AYC,AYS
        DOUBLE PRECISION DAXC,DAXS,DAYC,DAYS,DT,T
        DOUBLE PRECISION TABLE(4,N)
        DOUBLE PRECISION DATE(N),XC(N),XS(N),YC(N),YS(N),SX(N),SY(N)

        PARAMETER ( PI = 3.14159265358979323846D0 )
        
*       Mean prediction error

        MPE = 0.1325D0                       ! microarcseconds per day

*       FCN parameters

        PER = -430.21D0                      ! period in days
        PHI = (2D0*PI/PER)*(MJD-51544.5D0)   ! phase in radians
                
*       Block data of amplitudes for X (microas)
        
      DATA ((TABLE(I,J),I=1,4),J=1,N) /
     . 45700.D0,     4.55D0,   -36.58D0,    19.72D0, ! 1984.0
     . 46066.D0,  -141.82D0,  -105.35D0,    11.12D0, ! 1985.0
     . 46431.D0,  -246.56D0,  -170.21D0,     9.47D0, ! 1986.0
     . 46796.D0,  -281.89D0,  -159.24D0,     8.65D0, ! 1987.0
     . 47161.D0,  -255.05D0,   -43.58D0,     8.11D0, ! 1988.0
     . 47527.D0,  -210.46D0,   -88.56D0,     7.31D0, ! 1989.0
     . 47892.D0,  -187.79D0,   -57.35D0,     6.41D0, ! 1990.0
     . 48257.D0,  -163.01D0,    26.26D0,     5.52D0, ! 1991.0
     . 48622.D0,  -145.53D0,    44.65D0,     4.80D0, ! 1992.0
     . 48988.D0,  -145.12D0,    51.49D0,     5.95D0, ! 1993.0
     . 49353.D0,  -109.93D0,    16.87D0,     9.45D0, ! 1994.0
     . 49718.D0,   -87.30D0,     5.36D0,     8.25D0, ! 1995.0
     . 50083.D0,   -90.61D0,     1.52D0,     7.67D0, ! 1996.0
     . 50449.D0,   -94.73D0,    35.35D0,     4.40D0, ! 1997.0
     . 50814.D0,   -67.52D0,    27.57D0,     3.40D0, ! 1998.0
     . 51179.D0,   -44.11D0,   -14.31D0,     3.45D0, ! 1999.0
     . 51544.D0,     5.21D0,   -74.87D0,     3.26D0, ! 2000.0
     . 51910.D0,    70.37D0,  -129.66D0,     2.86D0, ! 2001.0
     . 52275.D0,    86.47D0,  -127.84D0,     2.75D0, ! 2002.0
     . 52640.D0,   110.44D0,   -42.73D0,     2.59D0, ! 2003.0
     . 53005.D0,   114.78D0,    -0.13D0,     2.53D0, ! 2004.0
     . 53371.D0,   132.96D0,    -4.78D0,     2.72D0, ! 2005.0
     . 53736.D0,   157.36D0,    28.63D0,     2.19D0, ! 2006.0
     . 54101.D0,   160.40D0,    58.87D0,     1.87D0, ! 2007.0
     . 54466.D0,   156.76D0,   101.24D0,     1.74D0, ! 2008.0
     . 54832.D0,   142.99D0,   143.01D0,     1.89D0, ! 2009.0
     . 55197.D0,    33.70D0,   184.46D0,     1.95D0, ! 2010.0
     . 55562.D0,     0.76D0,   253.70D0,     1.14D0, ! 2011.0
     . 55927.D0,    25.47D0,   271.66D0,     1.07D0, ! 2012.0
     . 56293.D0,   113.42D0,   256.50D0,     1.86D0/ ! 2013.0
          
*       Amplitudes extracted from the table

        DO I=1,N
           DATE(I)=TABLE(1,I)
           XC(I)=TABLE(2,I)
           XS(I)=TABLE(3,I)
           SX(I)=TABLE(4,I)
           YC(I)=XS(I)
           YS(I)=-XC(I)
           SY(I)=SX(I)
        END DO

*       Prediction of the amplitude at the input date

        IF (MJD.LE.DATE(1)) THEN
           AXC=XC(1)
           AXS=XS(1)
           AYC=YC(1)
           AYS=YS(1)
        ELSE IF (MJD.GE.DATE(N)) THEN
           AXC=XC(N)
           AXS=XS(N)
           AYC=YC(N)
           AYS=YS(N)
        ELSE
           DO I=1,N-1
              IF (MJD.GE.DATE(I).AND.MJD.LT.DATE(I+1)) THEN
                 T=MJD-DATE(I)
                 DT=DATE(I+1)-DATE(I)
                 DAXC=XC(I+1)-XC(I)
                 DAXS=XS(I+1)-XS(I)
                 DAYC=YC(I+1)-YC(I)
                 DAYS=YS(I+1)-YS(I)
                 AXC=XC(I)+(DAXC/DT)*T
                 AXS=XS(I)+(DAXS/DT)*T
                 AYC=YC(I)+(DAYC/DT)*T
                 AYS=YS(I)+(DAYS/DT)*T
              END IF
           END DO
        END IF
        
*       Computation of X and Y

        X=AXC*DCOS(PHI)-AXS*DSIN(PHI)  ! microas
        Y=AYC*DCOS(PHI)-AYS*DSIN(PHI)  ! microas
        
*       Prediction of the uncertainty at the input date

        IF (MJD.LE.DATE(1)) THEN
           AXC=SX(1)+MPE*(DATE(1)-MJD)
           AXS=SX(1)+MPE*(DATE(1)-MJD)
           AYC=SY(1)+MPE*(DATE(1)-MJD)
           AYS=SY(1)+MPE*(DATE(1)-MJD)
        ELSE IF (MJD.GE.DATE(N)) THEN
           AXC=SX(N)+MPE*(MJD-DATE(N))
           AXS=SX(N)+MPE*(MJD-DATE(N))
           AYC=SY(N)+MPE*(MJD-DATE(N))
           AYS=SY(N)+MPE*(MJD-DATE(N))
        ELSE
           DO I=1,N-1
              IF (MJD.GE.DATE(I).AND.MJD.LT.DATE(I+1)) THEN
                 T=MJD-DATE(I)
                 DT=DATE(I+1)-DATE(I)
                 DAXC=SX(I+1)-SX(I)
                 DAXS=SX(I+1)-SX(I)
                 DAYC=SY(I+1)-SY(I)
                 DAYS=SY(I+1)-SY(I)
                 AXC=DABS(SX(I)+(DAXC/DT)*T)
                 AXS=DABS(SX(I)+(DAXS/DT)*T)
                 AYC=DABS(SY(I)+(DAYC/DT)*T)
                 AYS=DABS(SY(I)+(DAYS/DT)*T)
              END IF
           END DO
        END IF
        
*       Computation of the uncertainties

        DX=AXC+AXS  ! microas
        DY=AYC+AYS  ! microas

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
