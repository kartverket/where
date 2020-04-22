module hfeop_xyu
! Hold tidal information.
!  Maximum tides is 200.
!  So far the model with the largest is JPL which has ~160
! History
! 2017Nov22 JMGipson. Original version
! 2019Oct09 Two bugs found and fixed by Michael Gerstl. Rotines calc_gmst and calc_nut_arg.
!           Effect is ~1.0 uas in PM, ~0.01 us in UT1
      
!      character*180  lhfeop_xyu_file                !this holds the name of the file. 
      logical khfeop_xyu
      integer , parameter :: max_tide=200
      integer num_tide             !number we find   
      integer*2 itide_arg(6,max_tide)    !tidal arguments. GST+pi followed by nutation.             
      character*8 ctide(max_tide)
      character*8 cdoodson(max_tide)             
      real*8 tide_period(max_tide)       !period
      real*8 tide_coef(2,4,max_tide)     !Coefficients of model in order: (sin,cosine)x(X,Y,UT1,LOD)x Max_tide                 

Contains
!***********************************************************************************
      subroutine import_tides_xyu(lhfeop_file)
      implicit none 
! Written by:
!    JMGipson 2017Nov22
!    Import HFEOP model in IERS format 
!
! input
      character*80 lhfeop_file     !file that contains the data. 

! local
      integer ilu 
      character*128   ldum          !temporary buffer to hold input line.  
      logical kexist
      
      inquire(exist=kexist, file=lhfeop_file)
      if(.not. kexist) then
         write(*,*) "Did not find HF_EOP_XYU file "//trim(lhfeop_file)
         stop
      endif 
  
      num_tide=0
      ilu=35
      open(ilu, file=lhfeop_file)     
! Skip any lines that  begin with "#"     
      
! get the number of components in model.
      do while(num_tide .lt. max_tide) 
        read(ilu,'(a)',end=100) ldum       !read in a line, exit when reach EOF.     
        if(ldum(1:1) .ne. "#") then     
          num_tide=num_tide+1  
          read(ldum,*,err=110) ctide(num_tide),itide_arg(1:6,num_tide),cdoodson(num_tide), &
    &     tide_period(num_tide), tide_coef(1:2,1:4,num_tide)
!         write(*,*) itide_arg(:,num_tide), tide_coef(:,:, num_tide)
         endif 
       end do 
       write(*,*) "ran out of space in import_iers_tides ", num_tide
      
100   continue      
!      write(*,*) "read in ", num_tide
      close(ilu)        
      return 
      
110   continue
      write(*,*) "Error reading in line: "
      write(*,*) trim(ldum)
      stop
      
      end subroutine
!***********************************************************************************
      SUBROUTINE calc_hf_eop_xyu(RMJD_TT,Delta_T,EOP)   
      IMPLICIT NONE                           
!     Routine to compute diurnal and semi-diurnal HF-EOP when model is in IERS format.
!     
! Passed. 
       real*8 RMJD_TT     !modified julian day TT 
       real*8 Delta_T     !TT-RMJD
 
! returned
      real*8    EOP(4,2)     !X,Y, UT1, LOD and their derivatives        
! Function      
 
! Local                        
      integer   itide     !loop counter over tides.
      integer   jeop  !loop counter over EOP
 
! arg,argd -- argument, argument_dot.
      real*8  arg,argd, fund_arg(6,2)
      logical knew_order                  !if true, then GMST+pi is first arguement. if false, then last. 
      data knew_order/.true./
      
!
!**** Get the fundamental arguments at this epoch
!
      call tide_angles(RMJD_TT, Delta_T, knew_order, fund_arg)
    
      EOP=0.d0
!     Now loop over the tidal contribution.  
      do itide = 1, num_tide
!         Get the argument and the time_derivative. 
          arg=dotarg(itide_arg(1,itide), fund_arg(1,1))
!         arg = mod(arg, 2.d0*pi)
          argd=dotarg(itide_arg(1,itide), fund_arg(1,2))      
!          write(*,*) tide_coef(1:2,1:3,itide) 
          !stop
          do jeop=1,4
             EOP(jeop,1)=eop(jeop,1)+sin(arg)*tide_coef(1,jeop,itide)+cos(arg)*tide_coef(2,jeop,itide)
             EOP(jeop,2)=eop(jeop,2)+argd*(cos(arg)*tide_coef(1,jeop,itide)-sin(arg)*tide_coef(2,jeop,itide))
          end do           
        end do      
        
!***** That is all.
        return
        end subroutine calc_hf_eop_xyu
! ****************************************************************************************8
! Calculate the tide angles. 
      subroutine tide_angles(rmjd_TT,Delta_T,knew_order,fund_arg)
      implicit none 
! Passed
      real*8 rmjd_TT       !rmjd (TT)
      real*8 Delta_T       !TT-UT1 (seconds)
      logical knew_order    !True: Put the arguments in the order GMST+pi, fund_arg. (2010 IERS conventions.)
                            !False: Order is nutation, GMST+pi
! returned
      real*8 fund_arg(6,2)  !Argument and their time derivatives. 

! local
      real*8 nut_arg(5,2)    !Fundamental arguments of nutation. Second index is time derivative.
      real*8 gmst(2)         !GMST. Second index is time derivative. 
      real*8 RMJD_UT         !Modified julian date TT.

      real*8 pi, twopi 
      parameter ( pi    = 3.1415926535897932D0 )
      parameter (twopi = 2.d0*pi)
     

      rmjd_UT = rmjd_TT-Delta_T/86400.d0
      call calc_gmst(rmjd_UT,gmst)
           
      call calc_nut_arg(rmjd_tt,nut_arg)
! Add pi to GMST 
      gmst(1)=dmod(gmst(1)+pi,twopi)

! If knew_order=.true. then put GMST+pi as the first argument
!              =.false. use 1996 IERS order where GMST+pi was last. 

!      knew_order=.false. 
      if(knew_order) then 
        fund_Arg(1,1:2)=gmst 
        fund_arg(2:6,1:2)=nut_arg 
      else
        fund_Arg(6,1:2)=gmst
        fund_arg(1:5,1:2)=nut_arg 
      endif   

      return
      end subroutine tide_angles
     
! ********************************************************************** 
! Compute GMST
      subroutine calc_gmst(RMJD_UT,gmst)

      implicit none 
      real*8 RMJD_UT          
      real*8 gmst(2)     !GMST and it's derivative. 
      
      logical kuse_mod   !Do the calculation in a way that is more numerically stable. 
   
! local
      real*8  T          !UT Julian centuries since J2000. 

      real*8 pi, twopi 
      real*8 sec2rad         
      real*8 sec_per_circ    
      parameter ( pi    = 3.1415926535897932D0 )
      parameter (twopi = 2.d0*pi)
      parameter (sec_per_circ=1296000d0)           !how many arcseconds in a circle. 
      parameter (sec2rad=twopi/sec_per_circ)       !conversion from seconds to radians.
   
      real*8 tmp_big, tmp_small

! coefficients of expansion for GMST 

! An alternate form of co(1,2) is (876600d0*3600d0 + 8640184.812866d0)
! Below this is expanded.  
      real*8 co(4)
      data co/67310.54841d0,   3164400184.812866d0, 0.093104d0, -6.2d-6/

!Rewrite the second coefficient in another way. 
!    co(2)=cot(1,j)*sec_per_circ+cot(2,j)
      real*8 cot(2)
      data cot/2442.d0, -431815.18734d0/ 

      t=(RMJD_UT-51544.5d0)/36525.d0        
  
      kuse_mod=.true.
      if(kuse_mod) then
        tmp_big=co(1)+dmod(cot(1)*T,1.d0)*sec_per_circ+cot(2)*T
        tmp_big=dmod(tmp_big,sec_per_circ)
      else 
         tmp_big=co(1)+co(2)*T       
      endif 

!%%  tmp_small = co(3)*T*2 + co(4)*T**3    ! original
     tmp_small = co(3)*T**2 + co(4)*T**3   ! correction. Michael Gerstl

    
! Convert from time-seconds to arc-seconds
      tmp_big=tmp_big*15.d0
      tmp_small=tmp_small*15.d0
      

      gmst(1)=dmod(tmp_big+tmp_small,sec_per_circ)    !If we are the boundary, may overflow
      if(gmst(1) .lt. 0.d0) gmst(1)=gmst(1)+sec_per_circ 
      gmst(1)=dmod(gmst(1),sec_per_circ)*sec2rad     !convert to radians   

      gmst(2)=co(2)+2.d0*co(3)*T+3.d0*co(4)*T**2
      gmst(2)=gmst(2)*15*sec2rad/36525.d0/86400.d0  !convert to radians/sec
      return
      end subroutine calc_gmst
!******************************************************************************************
      subroutine calc_nut_arg(RMJD_TT, arg)
      real*8 RMJD_TT               ! Modified Julian Day (TT) 
      real*8 arg(5,2)

      logical kuse_mod          !use the mod calculation 
      real*8 t
      real*8 pi, twopi
      real*8 sec2rad         
      real*8 sec_per_circ    
      parameter ( pi    = 3.1415926535897932D0 )
      parameter (twopi = 2.d0*pi)
      parameter (sec_per_circ=1296000d0)
      parameter (sec2rad=twopi/sec_per_circ)       !conversion from seconds to radians.
      real*8 tmp_big, tmp_small
      real*8 tmp

! Written by JMGipson
!   2017Oct02. 
!
! Compute the fundamental arguments and their derivatives.
! Order is: l, lp, f,d, omega 

! The coefficients for the data. 
! The order is constant, T, T^2, T^3, T^4

! The values in the table below come from the routine fundarg.f 
!  http://iers-conventions.obspm.fr/2010/2010_official/chapter8/software/FUNDARG.F

! the version below expands this out.
      real*8 co(5,5), cot(2,5)
      data co / &
     & 485868.249036d0, 1717915923.2178d0,  31.8792d0, 0.051635d0,  -0.00024470d0, &
     & 1287104.79305d0, 129596581.0481d0,   -0.5532d0, -0.000136d0, 0.00001149d0,  &
     & 335779.526232d0, 1739527262.8478d0, -12.7512d0, -0.001037d0, 0.00000417d0,  &
     & 1072260.70369d0, 1602961601.2090d0,  -6.3706, 0.006593d0, -0.00003169d0,    &
     & 450160.398036d0, -6962890.2665d0,    0.007702d0, 7.4722d0, -0.00005939d0 /

! The values in COT are related to co(2,*) by
!    co(2,j)=cot(1,j)*sec_per_circ+cot(2,j)
      data cot /  &
     & 1326.d0, -580076.722d0,  &
     &  100.d0, -3418.9519d0,   &
     & 1342.d0,  295262.8478d0, &
     & 1237.d0, -190398.791d0,  &
     &   -5.d0, -482890.2665d0/

      t=(RMJD_TT-51544.5)/36525.d0   
!      write(*,*) "NEW CENT ", T
! The general formula is:
!       arg(j,1)=co(1,j)+co(2,j)*T+Co(3,j)*T**2+co(4,j)*T**3+co(5,j)*T**4  
!       arg(j,2)=co(2,j)+co(3,j)*T+Co(4,j)*T**2+Co(5,j)*T**3


! If this is true, do the calculation in a way that is slighltly more numerically stable. 
      kuse_mod=.true.

! do it the following way for numerical stability.
      do i=1,5
        if(kuse_mod) then
          tmp_big=co(1,i)+dmod(cot(1,i)*T,1.d0)*sec_per_circ+cot(2,i)*T
          tmp_big=dmod(tmp_big,sec_per_circ)
        else 
          tmp_big=co(1,i)+co(2,i)*T       
        endif 

!%%   tmp_small = co(3,i)*T*22 + co(4,i)*T**3 + co(5,i)*T**4   ! original
      tmp_small = co(3,i)*T**2 + co(4,i)*T**3 + co(5,i)*T**4   ! Michael Gerstl

! Conversion from time-seconds to angle-seconds 
        arg(i,1)=dmod(tmp_big+tmp_small,sec_per_circ)    !If we are the boundary, may overflow
        if(arg(i,1) .lt. 0.d0) arg(i,1)=arg(i,1)+sec_per_circ 

        arg(i,1)=dmod(arg(i,1),sec_per_circ)*sec2rad     !convert to radians
 
! Now do the rates. 
        arg(i,2)=co(2,i)+2.d0*co(3,i)*T+3.d0*co(4,i)*T**2+4*co(5,i)*T**3
!         arg(i,2)=co(2,i)    
! Conversion from time-seconds to angle-seconds  
        arg(i,2)=arg(i,2)*sec2rad/36525.d0/86400.d0    !Convert to radians/per sec      
      end do 
   
      return
      end subroutine calc_nut_arg
!
!********************************************************
      FUNCTION dotarg(iarg,angles)
      IMPLICIT NONE                         
      real*8 dotarg
! Calculate dot product
      integer*2 iarg(6)   !multipliers of fundametal  arguments. 
      real*8 angles(6)    !Values of GMST+pi, 5 fundamental arguments
      integer*2 i         !local do index
!
       dotarg=0.
       do i=1,6
         dotarg=dotarg+iarg(i)*angles(i)
       end do
       return
       end function dotarg
! ****************************************************
End module 


