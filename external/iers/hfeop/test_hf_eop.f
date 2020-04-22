      program test_hf_eop
      use hfeop_xyu      !module that contains data and calcuation routines.
! Program to read in HF-EOP files in IERS format, compute HF-EOP for a two week period and write out the results.
! The models are presumed to be in the directory ../models
!    John Gipson
!    2016Oct20 
!    301-614-6876
!   john.m.gipson@nasa.gov
!

     implicit none 
!
! functions
      integer*4 julday           !julian date   
      character*21 jd_to_date_1  !convert fjday to ascii string. Assumes FJDAY is UTC. 


      character*128 lhfeop_file     !input file that contains model information
      character*128 leop_out       !output file 
      real*8 delta_T
      real*8 MJDAY_UTC         !UTC time tag
      real*8 FJDAY_UTC         !Julian version of above
      real*8 MJDAY_TT          !TT tag
      real*8 EOP(4,2)          !XP,YP,UT1, LOD  and their rates. Units are (uas,uas,uts, uas/Day) 
                       
      integer i                !counter 
      integer j                !counter over models
      integer ind

      Delta_T=  68.931118         !Difference on 2017Nov28
   
! Read in several models and compute HF-EOP for each one. 
! k1_xyu and s2_xyu just have a single tide and are for testing purposes. 
    
      do j=1,6 
      if(j .eq. 1) then
        lhfeop_file= "../models/k1_xyu.txt"
        leop_out="k1_xyu.out"
      else if(j .eq. 2) then
        lhfeop_file= "../models/s2_xyu.txt"
        leop_out="s2_xyu.out"  
      else if(j .eq. 3) then 
        lhfeop_file="../models/abn_comb_2012_xyu.txt"
        leop_out="abn_comb_2012_xyu.out"
      else if(j .eq. 4) then
        lhfeop_file="../models/2017a_astro_xyu.txt"
        leop_out="2017a_astro_xyu.out"
      else if(j .eq. 5) then
        lhfeop_file="../models/fes2012_xyu.txt"
        leop_out="fes2012_xyu.out"
      else if(j .eq. 6) then
        lhfeop_file= "../models/madzak_xyu.txt"
        leop_out="madzak_xyu.out"        
      endif 

      write(*,*) "Doing "//trim(lhfeop_file) 
     
      call import_tides_xyu(lhfeop_file)    
! 
      open(2,file=leop_out)
! write out the results for the CONT17 time period. 2017NOV28 through 2017DEC12 
      write(2,*) "HFEOP predictions from model "//trim(lhfeop_file) 
      write(2,*) "X, Y, UT1, LOD units= (uas,uas,us, us/day) "
      write(2,*) "Time derivatives multiplied by 1000" 
  
      write(2,'("    MJD          Time            |     X      Y       UT1   LOD     |    Xdot    Ydot   UT1dot  LODdot |")')  
      MJDAY_UTC=58085           ! 2017NOV28   
! write out the EOP for two weeks
      do  i=1, 15*24  
         MJDAY_TT=MJDAY_UTC+Delta_T/86400.d0 
         call calc_hf_eop_xyu(MJDAY_TT, Delta_T, eop)
         fjday_UTC=mjday_UTC+2400000.5                !convert to modified Julian day  for use by jd_to_date_1
         write(2,'(f10.3,1x,a," | ",2(4f8.2," | ") )') mjday_UTC, jd_to_date_1(fjday_UTC) , &
      &       eop(1:4,1), eop(1:4,2)*1000.d0    
    
         mjday_UTC=mjday_UTC+1.d0/24.d0
      end do
      close(2)

      end do      !end looping over models. 
      end program 

