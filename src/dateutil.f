*************************************************************************
*	This file contains two FORTRAN 77 routines for converting between
*	Julian and calendar dates. Routine JTODMY converst from 
*	Julian day number to calendar date (Gregorian); DMYTOJ
*	goes the other way round. Algorithms obtained from 
*	http://hermetic.magnet.ch/cal_stud/jdn.htm on 28th November
*	2001 (Julian day 2452242).
*************************************************************************
      SUBROUTINE JTODMY(JULIAN,D0,M0,Y0,DD,MM,YY,IFAIL)
*
*     Converts a date in JULIAN format, where day 0 is D0/M0/Y0, 
*     to DD/MM/YY in the Gregorian calendar. If Y0, M0 and D0 are all 
*     zero, the origin defaults to 24th November, 4713BC (which is when 
*     the calendar starts!). The origin is actually 1st January 4712BC 
*     according to the *Julian* calendar. This algorithm is only valid
*     for dates after 1st March, 4900BC. Algorithm due to Henry F. Fliegel
*     and Thomas C. Van Flandern. Tested against Splus version, and agrees.
*
      INTEGER JULIAN,Y0,M0,D0,YY,MM,DD,IFAIL
      INTEGER L,N,I,J,JD
      
      CHARACTER*72 MESSAGE(2)

      IFAIL = 0

      IF ((Y0.EQ.0).AND.(M0.EQ.0).AND.(D0.EQ.0)) THEN
       JD = JULIAN
      ELSE
       CALL DMYTOJ(D0,M0,Y0,JD,IFAIL)
       IF (IFAIL.NE.0) RETURN
       JD = JULIAN + JD
      ENDIF
      L = JD + 68569
      N = (4*L) / 146097
      L = L - ((146097*N) + 3) / 4
      I = (4000 * (L+1)) / 1461001
      L = L - ((1461*I) / 4) + 31
      J = (80*L)/2447
      DD = L - ((2447*J) / 80)
      L = J/11
      MM = J + 2 - (12*L)
      YY = (100 * (N-49)) + I + L

      IF (((10000*YY)+(100*MM)+DD).LT.-49000301) THEN
       IFAIL = 1
       WRITE(MESSAGE,1)
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
       RETURN
      ENDIF

 1    FORMAT('****ERROR**** in Julian-to-date conversion: ',/,
     +     'Algorithm only works for dates before 1st March 4900BC!')
      END
*************************************************************************
*************************************************************************
*************************************************************************
      SUBROUTINE DMYTOJ(DD,MM,YY,JULIAN,IFAIL)
*
*     Converts date DD/MM/YY (on the Gregorian calendar) into a 
*     Julian day number (origin: 24th November 4713BC). This algorithm
*     is valid only for dates after 1st March, 4800BC. Algorithm due 
*     to Henry F. Fliegel and Thomas C. Van Flandern. Tested against 
*     Splus version, and agrees.
*
      INTEGER YY,MM,DD,JULIAN,IFAIL
      
      CHARACTER*72 MESSAGE(2)

      IFAIL = 0
      IF (YY.EQ.0) THEN
       IFAIL = 1
       WRITE(MESSAGE,1)
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       RETURN
      ELSEIF (((10000*YY)+(100*MM)+DD).LT.-48000301) THEN
       IFAIL = 2
       WRITE(MESSAGE,2)
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
       RETURN
      ENDIF

      JULIAN = (1461 * (YY + 4800 + (MM - 14) / 12)) / 4 +
     +         (367 * (MM - 2 - 12 * ((MM - 14 ) / 12))) / 12 -
     +         (3 * ((YY + 4900 + (MM - 14) / 12) / 100)) / 4 +
     +         DD - 32075      

 1    FORMAT('****ERROR**** in date-to-Julian conversion: ',
     +     'Gregorian year 0 doesn''t exist!')
 2    FORMAT('****ERROR**** in date-to-Julian conversion: ',/,
     +     'Algorithm only works for dates before 1st March 4800BC!')
      END
******************************************************************
******************************************************************
******************************************************************
      Subroutine Timediff(Date1,Time1,Date2,Time2,Elapsed,Ifail)
******************************************************************
*     Returns the number of seconds between two times as returned
*     by a call to DATE_AND_TIME. This is to be used primarily 
*     for determining run times. Date? and Time? (input, character)  
*     are the values of DATE and TIME obtained from two separate  
*     calls to DATE_AND_TIME; Elapsed (output, double) is the 
*     difference between them in seconds. Ifail is an error flag.
******************************************************************
      Character, intent(in) :: Date1*8,Time1*10,Date2*8,Time2*10
      Double precision, intent (out) :: Elapsed
      Integer, intent(out) :: Ifail
******************************************************************
*     Extra Integers
*     ~~~~~~~~~~~~~~
*     CN?,YR?,  Centuries, years, months and days corresponding 
*     MO?,DY?   to Date1 and Date2
*     HR?,MN?,  Hours and minutes
*     J1,J2     Julian day numbers
******************************************************************
      Integer CN1,CN2,YR1,YR2,MO1,MO2,DY1,DY2,HR1,HR2,MN1,MN2,J1,J2
******************************************************************
*     Extra Doubles
*     ~~~~~~~~~~~~~
*     SC?       Seconds corresponding to Date1 and Date2
******************************************************************
      Double precision SC1,SC2
      
      Ifail = 0
      
      Read(Date1,'(4I2)') CN1,YR1,MO1,DY1 
      Read(Date2,'(4I2)') CN2,YR2,MO2,DY2
      Read(Time1,'(2I2,F6.3)') HR1,MN1,SC1 
      Read(Time2,'(2I2,F6.3)') HR2,MN2,SC2
      
      Call DMYTOJ(DY1,MO1,(100*CN1+YR1),J1,Ifail)
      If (Ifail.NE.0) Return
      Call DMYTOJ(DY2,MO2,(100*CN2+YR2),J2,Ifail)
      If (Ifail.NE.0) Return
  
      Elapsed = Dble(J2-J1)*8.64D4 + Dble(HR2-HR1)*3.6D3 + 
     +          Dble(MN2-MN1)*6.0D1 + SC2-SC1   
     
      End subroutine Timediff
      
