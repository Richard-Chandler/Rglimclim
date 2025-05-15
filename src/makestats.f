      Subroutine DailyStats(FilNo,Nsites,Nvars,Nregs,WantedSites,
     +                      WantedVars,Thresholds,RegionDefs,
     +                      TmpFlNo,MoStats,Ifail)
*
*     To compute summary statistics from a data file containing a daily 
*     simulation. Arguments:
*
*       FilNo   - Number of file from which to read the simulated data (which
*                 should already be connected for input on entry)
*       Nsites  - Number of sites for which data is presented in the data file
*       Nvars   - Number of variables with data in file
*       Nregs   - Number of regions for which summaries have been requested
*       WantedSites     - Binary vector indicating which sites to calculate
*                         summaries for
*       WantedVars      - Ditto, variables
*       Thresholds      - Thresholds for use in calculation of conditional
*                         summaries (threshold exceedances)
*       RegionDefs      - A binary matrix indicating which sites are included
*                         in each of the regions for which summaries have
*                         been requested (rows correspond to sites, columns
*                         to regions)
*       TmpFlNo         - Number of temporary file containing site codes.
*       MoStats         - array of summary statistics for each month, site and
*                         variable. The elements are:
*
*                               1 Mean
*                               2 Standard deviation
*                               3 Maximum
*                               4 Minimum
*                               5 Lag 1 autocorrelation
*                               6 Lag 2 autocorrelation
*                               7 Lag 3 autocorrelation
*                               8 Proportion of threshold exceedances
*                               9 Mean conditional on threshold exceedance
*                              10 Std dev conditional on threshold exceedance
*                      11 onwards Correlations with remaining variables
*
*       Ifail           - Error code
*
      Integer, intent(in) :: FilNo,Nsites,Nvars,Nregs,TmpFlNo
      Integer, intent(in) :: WantedSites(Nsites),WantedVars(Nvars)
      Integer, intent(in) :: RegionDefs(Nsites,Nregs)
      Double precision, intent(in) :: Thresholds(Nvars)
      Double precision, intent(out) :: 
     +                  MoStats(12,Nsites+Nregs,Nvars,10+Nvars)
      Integer, intent(out) :: Ifail
******************************************************************************
*       Additional INTEGERs
*       ^^^^^^^^^^^^^^^^^^^
*       YY,MM,DD        - Year, month and day
*       OLDYR,OLDMO,OLDDY - Year, month and day yesterday
*       NDAYS           - Vector containing sample sizes for each month
*       DAYS_IN_MONTH   - Function to return number of days in a month
*       SitesPerRegion  - Vector containing the number of sites in each region
*       WantedLocs      - Vector indicating for which locations (including
*                         both sites and regions) summaries are required
*       ACFN            - Numbers of pairs of observations contributing to
*                         each estimated autocorrelation
*       DONE            - Indicator for whether we've reached EOF
*       I,J             - counters
******************************************************************************
      Integer YY,MM,DD,OLDYR,OLDMO,OLDDY,NDAYS(12),DONE
      Integer DAYS_IN_MONTH,SitesPerRegion(Nregs)
      Integer WantedLocs(Nsites+Nregs)
      Integer ACFN(12,Nsites+Nregs,Nvars,3)
      Integer I,J
******************************************************************************
*       Additional DOUBLEs
*       ^^^^^^^^^^^^^^^^^^
*       DatArray        - array holding current day's data
*       Modata          - array of data for this month.
*       ACFMeans        - array holding means for calculations of 
*                         autocorrelations (this because ACFs are
*                         calculated as correlations between shifted
*                         series rather than using the usual formula,
*                         as an easy way to handle the replication
*                         of months that will usually be present). 
******************************************************************************
      Double precision DatArray(Nsites,Nvars)
      Double precision MoData(31,Nsites+Nregs,Nvars)
      Double precision ACFMeans(12,Nsites+Nregs,Nvars,3,2)
******************************************************************************
*       Additional CHARACTERs
*       ^^^^^^^^^^^^^^^^^^
*       Scodes  - vector of short site codes
******************************************************************************
      Character Scodes(Nsites)*4
***************************************************************************
*
*       Initialise (3rd and 4th columns of MOSTATS are for max and min)
*
      MoStats = 0.0D0
      MoStats(1:12,1:(Nsites+Nregs),1:Nvars,3) = -1.0D100
      MoStats(1:12,1:(Nsites+Nregs),1:Nvars,4) = 1.0D100
      DatArray = 0.0D0
      Ifail = 0
      OLDMO = 0
      DONE = 0
      YY = 0 
      MM = 0
      DD = 0
      ACFN = 0
      Ndays = 0
      WantedLocs = 1
      WantedLocs(1:Nsites) = WantedSites
      ACFMeans = 0.0D0
*
*       Read site codes from temporary file
*
      CALL GetScodes(TmpFlNo, NSites, Scodes, Ifail)
      IF (IFAIL.NE.0) RETURN
*
*       Increment date unless this is the first time through, in which
*       case we don't know what the first date is going to be
*
 10   IF (YY.NE.0) THEN
       DD = DD + 1
       IF (DD.GT.DAYS_IN_MONTH(MM,YY)) THEN
        DD = 1
        MM = MM + 1
        IF (MM.GT.12) THEN
         MM = 1
         YY = YY + 1
        ENDIF
       ENDIF
      ENDIF

      Call DATRD(FilNo,DatArray,-1.0D100,Scodes,Nsites,Nvars,
     +                                        YY,MM,DD,DONE,Ifail)
      IF (Ifail.NE.0) RETURN
*
*     New month - compute stats
*
      IF ((MM.NE.OLDMO).AND.(OLDMO.NE.0)) THEN
       CALL STUPDT(MoData,MoStats,ACFMeans,ACFN,Nsites+Nregs,Nvars,
     +             WantedLocs,WantedVars,Thresholds,OLDMO,OLDDY)
       NDAYS(OLDMO) = NDAYS(OLDMO) + OLDDY
      ENDIF
*
*     Update data matrix
*
      MoData(DD,1:Nsites,1:Nvars) = DatArray
*
*       Now calculate regional mean time series
*
      MoData(DD,(Nsites+1):(Nsites+Nregs),1:Nvars) = 0.0D0
      SitesPerRegion = SUM(RegionDefs,Dim=1)
      
      Do I=1,Nregs
       Do J=1,Nsites
        If (RegionDefs(J,I).EQ.1) then
         MoData(DD,Nsites+I,1:Nvars) = MoData(DD,Nsites+I,1:Nvars) + 
     +                                 DatArray(J,1:Nvars)    
        Endif
       End Do
       MoData(DD,Nsites+I,1:Nvars) = MoData(DD,Nsites+I,1:Nvars) /
     +                                      Dble(SitesPerRegion(I))
      End Do

      OLDMO = MM
      OLDYR = YY
      OLDDY = DD
      IF (DONE.EQ.0) GOTO 10

      NDAYS(OLDMO) = NDAYS(OLDMO) + OLDDY
      
      CALL STUPDT(MoData,MoStats,ACFMeans,ACFN,Nsites+Nregs,Nvars,
     +                WantedLocs,WantedVars,Thresholds,OLDMO,OLDDY)
      CALL STCALC(MoStats,ACFMeans,Nsites+Nregs,Nvars,WantedLocs,
     +                                          WantedVars,NDAYS,ACFN)

      END
************************************************************************
************************************************************************
************************************************************************
      SUBROUTINE STUPDT(MoData,MoStats,ACFMeans,ACFN,Nlocs,Nvars,
     +                  WantedLocs,WantedVars,Thresholds,MONTH,NDAYS)

*
*     Updates contributions to monthly summary statistics. NB NDAYS
*     here is the number of days in a single month, whereas in other
*     routines it's the total number of days in all instances of the
*     current month. 
*
      Integer, intent(in) :: Nlocs,Nvars,NDAYS,MONTH
      Integer, intent(in) :: WantedLocs(Nlocs),WantedVars(Nvars)
      Integer, intent(inout) :: ACFN(12,Nlocs,Nvars,3)
      Double precision, intent (in) :: MODATA(31,Nlocs,Nvars)
      Double precision, intent (in) :: Thresholds(Nvars)
      Double precision, intent(inout) :: 
     +          MoStats(12,Nlocs,Nvars,10+Nvars),
     +          ACFMeans(12,Nlocs,Nvars,3,2)
      INTEGER I,J,J2,K,LAG
*
*     Columns of MOSTATS: 1=Mean,2=SD, 3=Max, 4=Min, 5=ACF1,6=ACF2,7=ACF3,
*     8=P(exceed), 9=E(X|exceedance), 10=SD(X|exceedance), 11 onwards=
*     correlations with other variables
*
      Do I=1,Nlocs
       If (WantedLocs(I).EQ.0) goto 100 
       Do J=1,Nvars
        If (WantedVars(J).EQ.0) goto 99
        Do K=1,NDAYS
         MOSTATS(MONTH,I,J,1) = MOSTATS(MONTH,I,J,1) + MODATA(K,I,J)
         MOSTATS(MONTH,I,J,2) = MOSTATS(MONTH,I,J,2) + 
     +                                              (MODATA(K,I,J)**2)
         IF (MODATA(K,I,J).GT.MOSTATS(MONTH,I,J,3)) THEN
          MOSTATS(MONTH,I,J,3) = MODATA(K,I,J)
         ENDIF
         IF (MODATA(K,I,J).LT.MOSTATS(MONTH,I,J,4)) THEN
          MOSTATS(MONTH,I,J,4) = MODATA(K,I,J)
         ENDIF
*
*       Conditional statistics only where thresholds are defined (columns
*       8, 9 and 10)
*
         If ((Thresholds(J).GT.-1.0D100).AND.
     +       (MODATA(K,I,J).GT.Thresholds(J))) then
          MOSTATS(MONTH,I,J,8) = MOSTATS(MONTH,I,J,8) + 1.0D0
          MOSTATS(MONTH,I,J,9) = MOSTATS(MONTH,I,J,9) + MODATA(K,I,J)
          MOSTATS(MONTH,I,J,10) = MOSTATS(MONTH,I,J,10) + 
     +                                               (MODATA(K,I,J)**2)
         End if
        End Do
*
*       Autocorrelations (columns 5, 6 and 7),and updates to elements of
*       ACFMeans
*
        Do LAG=1,3
         Do K=1,NDAYS-LAG
          MOSTATS(MONTH,I,J,4+LAG) = MOSTATS(MONTH,I,J,4+LAG) + 
     +                            ( MODATA(K,I,J) * MODATA(K+LAG,I,J) )
          ACFMEANS(MONTH,I,J,LAG,1) = ACFMEANS(MONTH,I,J,LAG,1) + 
     +                                                    MODATA(K,I,J)
          ACFMEANS(MONTH,I,J,LAG,2) = ACFMEANS(MONTH,I,J,LAG,2) + 
     +                                                MODATA(K+LAG,I,J)
          ACFN(MONTH,I,J,LAG) = ACFN(MONTH,I,J,LAG) + 1
         End Do
        End Do
*
*       Correlations with other variables
*
        If (Nvars.GT.1) then
         Do K=1,Ndays
          Do J2=1,J-1
           If (WantedVars(J2).EQ.0) goto 98
           MoStats(Month,I,J,10+J2) = MoStats(Month,I,J,10+J2) + 
     +                                (MoData(K,I,J)*MoData(K,I,J2))      
  98      End Do
         End Do 
        End If
  99   End Do
 100  End Do

      END
************************************************************************
************************************************************************
************************************************************************
      SUBROUTINE STCALC(MOSTATS,ACFMeans,Nlocs,Nvars,WantedLocs,
     +                                         WantedVars,Ndays,ACFN)
*
*     Finalises monthly summary statistics
*
      INTEGER, intent(in) :: Nlocs,Nvars,Ndays(12)
      Integer, intent(in) :: WantedLocs(Nlocs),WantedVars(Nvars)
      Double precision, intent(inout) :: 
     +                         MOSTATS(12,Nlocs,Nvars,10+Nvars),
     +                         ACFMeans(12,Nlocs,Nvars,3,2)
      Integer, intent(in) :: ACFN(12,Nlocs,Nvars,3)

      INTEGER I,J,J2,LAG,M
*
*     Columns of MOSTATS: 1=Mean,2=SD, 3=Max, 4=Min, 5=ACF1,6=ACF2,7=ACF3,
*     8=P(exceed), 9=E(X|exceedance), 10=SD(X|exceedance), 11 onwards=
*     correlations with other variables
*
      Do 100 M=1,12
       If (Ndays(M).EQ.0) then
        Mostats(M,1:Nlocs,1:Nvars,1:(10+Nvars)) = -1.0d100
        Goto 100
       End if
       Do 110 I=1,Nlocs
        If (WantedLocs(I).EQ.0) goto 110 
        Do 120 J=1,Nvars 
         If (WantedVars(J).EQ.0) goto 120
         MOSTATS(M,I,J,1) = MOSTATS(M,I,J,1) / DBLE(Ndays(M))
         If (Ndays(M).LE.1) then
          MOSTATS(M,I,J,2) = -1.0D100
          MOSTATS(M,I,J,5:7) = -1.0D100
          Goto 120
         End if 
         MOSTATS(M,I,J,2) = MOSTATS(M,I,J,2) - 
     +                        (DBLE(Ndays(M))*(MOSTATS(M,I,J,1)**2))
         MOSTATS(M,I,J,2) = DSQRT(MOSTATS(M,I,J,2) / DBLE(NDAYS(M)-1))
*
*       Correlations with other variables (NB only the ones with indices
*       below J have the required quantities calculated at present)
*
         If (Nvars.GT.1) then
          Do J2=1,J-1
           If (WantedVars(J2).EQ.0) goto 98
            MoStats(M,I,J,10+J2) = MoStats(M,I,J,10+J2) -
     +                             (DBLE(Ndays(M))*MoStats(M,I,J,1)*
     +                                             MoStats(M,I,J2,1))
            MoStats(M,I,J,10+J2) = MoStats(M,I,J,10+J2) /
     +                              (DBLE(Ndays(M)-1) * 
     +                                    MoStats(M,I,J,2) *
     +                                    MoStats(M,I,J2,2))
            MoStats(M,I,J2,10+J) = MoStats(M,I,J,10+J2)
  98      End Do
          MoStats(M,I,J,10+J) = 1.0D0
         End If
*
*       Autocorrelations. NB calculation is as sample correlation between
*       pairs (with separate means for the "before" and "after" elements
*       of each pair), rather than the usual formula for an ACF - this 
*       is to ensure the results are valid correlations when combining 
*       information from multiple years. 
*
         DO 160 LAG=1,3
          ACFMeans(M,I,J,LAG,1) = ACFMeans(M,I,J,LAG,1) /
     +                                       DBLE(ACFN(M,I,J,LAG))
          ACFMeans(M,I,J,LAG,2) = ACFMeans(M,I,J,LAG,2) /
     +                                       DBLE(ACFN(M,I,J,LAG))
          MOSTATS(M,I,J,4+LAG) = MOSTATS(M,I,J,4+LAG) -
     +                           (Dble(ACFN(M,I,J,LAG)) * 
     +                              ACFMeans(M,I,J,LAG,1) * 
     +                              ACFMeans(M,I,J,LAG,2))
          IF (MOSTATS(M,I,J,2).GT.0.0D0) THEN
           MOSTATS(M,I,J,4+LAG) = MOSTATS(M,I,J,4+LAG) / 
     +                            (Dble(Ndays(M)-1) * 
     +                             (MOSTATS(M,I,J,2)**2))
          ELSE
           MOSTATS(M,I,J,4+LAG) = -1.0D100
          ENDIF
 160     CONTINUE
*
*       Threshold exceedance statistics if requested
*
         IF (MOSTATS(M,I,J,8).GT.0.0D0) THEN
          MOSTATS(M,I,J,9) = MOSTATS(M,I,J,9) / MOSTATS(M,I,J,8)
         ELSE
          MOSTATS(M,I,J,9) = -1.0D100
         ENDIF
         IF (MOSTATS(M,I,J,8).GT.1.0D0) THEN
          MOSTATS(M,I,J,10) = MOSTATS(M,I,J,10) - 
     +                        (MOSTATS(M,I,J,8)*(MOSTATS(M,I,J,9)**2))
          MOSTATS(M,I,J,10) = DSQRT(MOSTATS(M,I,J,10) / 
     +                                       (MOSTATS(M,I,J,8)-1.0D0))
         ELSE
          MOSTATS(M,I,J,10) = -1.0D100
         ENDIF
         MOSTATS(M,I,J,8) = MOSTATS(M,I,J,8) / DBLE(NDAYS(M))
 120    Continue
 110   Continue
 100  Continue

      END
