*       PROGRAM AUTHOR:         RICHARD E. CHANDLER
******************************************************************************
*                               (email: richard@stats.ucl.ac.uk)
*       AFFILIATION:            DEPARTMENT OF STATISTICAL SCIENCE
*                                UNIVERSITY COLLEGE LONDON
*                                 GOWER STREET
*                                  LONDON WC1E 6BT
*                                   UK
******************************************************************************
*  This file contains a collection of FORTRAN routines used in all of the
*  programs relating to Generalised Linear Modelling of daily climate
*  data. For more details of these programs, see the headers to the 
*  files fit_logi.f, fit_gamm.f and simrain.f
*
*  REVISION HISTORY: see manual.
*  © UCL 2002
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE TAUDEF(GLBVAL,GLBCOD,NGLB,THRESH,TRACE,THRTYP,
     +                                                  MODEL,MXP)
*
*     To define a level for `soft', `hard' or `trace' thresholding of 
*     non-negative data. Arguments:
*
*     GLBVAL    - Array of `global' quantities. The first element 
*                 contains the threshold. DOUBLE PRECISION, input
*     GLBCOD    - Indexing array, linking the quantities in the 
*                 model to the elements of GLBVAL. INTEGER, input
*     NGLB      - Number of global quantities defined in the model.
*                 INTEGER, input.
*     THRESH    - Threshold. DOUBLE PRECISION, output.
*     TRACE     - Threshold, if `trace' thresholding is to be used.
*                 DOUBLE PRECISION, output.
*     THRTYP    - Indicator for type of thresholding. 1 = `soft',
*                 2 = `hard', 0 = `none'. INTEGER, output.
*     MODEL     - Code for the type of model we're dealing with (10 is
*                 Gaussian, for which threshold should be large and
*                 negative)
*     MXP       - Used for dimensioning. INTEGER, input
******************************************************************************
      INTEGER MXP
      INTEGER GLBCOD(MXP),NGLB,THRTYP,MODEL
      DOUBLE PRECISION GLBVAL(MXP),THRESH,TRACE
******************************************************************************
*     Additional INTEGERs
*     ^^^^^^^^^^^^^^^^^^^
*     I         - counter
******************************************************************************
      INTEGER I
*
*     Basic strategy: threshold is zero for non-negative quantities and
*     something large and negative for real-valued quantities unless 
*     instructed otherwise. If threshold is defined, appropriate element 
*     of GLBCOD will be coded as 1000 + x, where x is 1 for `trace' 
*     treatment of small values, 2 for `soft' thresholding and 3 for 
*     `hard' thresholding. Since `trace' treatment does not involve 
*     thresholding, we get THRTYP = x - 1.
*
      IF ((MODEL.EQ.10).OR.(MODEL.EQ.20)) THEN
       THRESH = -1.0D100
       TRACE = -1.0D100
      ELSE      
       THRESH = 0.0D0
       TRACE = 0.0D0
      ENDIF
      THRTYP = 0
      DO 100 I=1,NGLB
       IF (GLBCOD(I)/1000.NE.1) GOTO 100
       THRTYP = MOD(GLBCOD(I),1000) - 1
       IF (THRTYP.EQ.0) THEN
        TRACE = GLBVAL(1)
       ELSE
        THRESH = GLBVAL(1)
       ENDIF
       RETURN
 100  CONTINUE

      END
******************************************************************************
******************************************************************************
******************************************************************************
      INTEGER FUNCTION DAYS_IN_MONTH(MONTH,YEAR)
      IMPLICIT NONE
C
C Author       : Kevin Black
C Date written : January 1985
C Abstract     :
C
C Integer function returns the number of days in the month MONTH for the
C year YEAR (leap years are accounted for)
C
      INTEGER MONTH,YEAR
C
C Local variables
C
      INTEGER NDAYS(12)
      DATA NDAYS/31,28,31,30,31,30,31,31,30,31,30,31/
C
C Get number of days for months other than February
C
      IF(MONTH.NE.2)THEN
         DAYS_IN_MONTH=NDAYS(MONTH)
C
C To be a leap year the year number must be divisble by 4
C but not by 100. However if the year is divisble by 4
C and by 400 then it is a leap year (this being the exception
C to the 100 divisibility rule).
C
      ELSE IF(YEAR/4*4.EQ.YEAR.AND.
     *        (YEAR/400*400.EQ.YEAR.OR.
     *         YEAR/100*100.NE.YEAR))THEN
              DAYS_IN_MONTH=29
           ELSE
              DAYS_IN_MONTH=28
           ENDIF
      RETURN
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE COVSET(COVS,ExtUnits,NP,NSITES,NVARS,RespIdx,COVCODE,
     +                  TWO,THREE,SITINF,Distance,YEAR,MONTH,DAY,
     +                  DatArray,AllowIncAvge,PWTIDX,MXLAG,BETA,THETA,
     +                  TRACE,RECALC,FORCE,MISSFL,ICHECK,IFAIL,MXP,SIM,
     +                  MISSST)
*******************************************************************************
*       Calculates covariates in model for 1 days' data.
*******************************************************************************
*       INTEGER variables
*       ^^^^^^^^^^^^^^^^^
*       ExtUnits- Contains numbers of units connected to input files 
*                 containing external data (1=annual,2=monthly,3=daily) 
*       NP      - No. of parameters in the model. Input
*       NSITES  - No. of sites. Input
*       NVARS   - No. of variables in input data file. Input
*       RespIdx - Index of response variable in input data file
*       AllowIncAvge    - Indicator for whether or not weighted averages
*                 of non-response variables will be considered as 
*                 valid in model fitting if the variable is missing at
*                 the current site 
*       COVCODE - Codes for covariates in the model. Input
*       TWO     - list of 2-way interactions. Input
*       THREE   - list of 3-way interactions. Input
*	MXLAG	- We'll chuck out any observations for which there are
*                 missing values during any of the MXLAG immediately 
*                 preceding days.
*	YEAR  	- current year. Input
*	MONTH   - current month. Input
*       DAY	- current day. Input
*       PWTIDX  - Index for locating parameters relating to nonlinear 
*                 weighting schemes for previous days' values.
*       WTSCHM  - Weighting scheme required.
*	OLDYR	- Value of YEAR last time through (do we have to recalculate?)
*	OLDMN	- Value of MONTH last time through (ditto)
*	RECALC	- Indicates whether we're calculating everything (0) or
*		  just things associated with nonlinear functions (>0). *NB
*                 if > 0, it is assumed that all required predictors are
*                 present for this case*
*	FORCE	- Forces recalculation of everything even if the year and
*		  month are the same as on previous call (for initialising
*		  COVS in programs where more than one model is involved)
*	MISSFL	- Array of flags indicating whether there's any
*		  missing data at any of the sites today. Output.
*       MXN     - Maximum no. of sites. Input
*       MXP     - Maximum no. of parameters. Input
*       SIM     - 1 if we're simulating, 0 otherwise (used to halt
*                 simulation if data on external predictors is missing, and 
*                 also to determine whether we need to calculate derivative
*                 information for nonlinear parameters).
*       MISSST  - Indicator - should we set MISSFL in this routine?
*	IDX1	- Covariate no. for a nonlinear parameter we're 
*		  currently dealing with
*	IDX2	- parameter no, relative to the covariate
*       EXTLAG  - lag, for external predictors
*       YMISS } - indicators telling us whether any required external
*       MMISS }   predictors are missing
*       DMISS
*       VARNUM  - Number of variable with which we're currently working
*       ICHECK  - Indicator to tell us whether parametrisation has been 
*                 checked (only need to carry out checks on first call
*                 to this routine)
*       IFAIL   - Error flag
*******************************************************************************
      INTEGER WTSCHM,MXP,SIM,MISSST
      INTEGER, intent(in) :: ExtUnits(3),NP(10),NSITES,NVARS,RespIdx
      INTEGER, intent(in) :: AllowIncAvge
      INTEGER COVCODE(MXP),TWO(MXP,2),THREE(MXP,3)
      INTEGER YEAR,MONTH,DAY,OLDYR,OLDMN,PWTIDX(MXP,0:3,NVARS)
      INTEGER RECALC,FORCE,MISSFL(NSITES),MXLAG(NVARS)
      INTEGER I,J,K,IDX1,IDX2,EXTLAG,YMISS,MMISS,DMISS,VARNUM
      INTEGER ICHECK,IFAIL
*******************************************************************************
*       DOUBLE PRECISION variables
*       ^^^^^^^^^^^^^^^^^^^^^^^^^^
*       COVS    - Values of covariates at each site. Input/output. The routine
*		  modifies only those elements which need to be modified NB THIS
*		  REQUIRES THAT ON ENTRY, THE ELEMENTS OF COVS 
*		  CORRESPONDING TO YEAR AND MONTHLY EFFECTS HAVE NOT BEEN
*		  UPDATED SINCE THE PREVIOUS CALL, IF WE'RE IN THE SAME YEAR/
*		  MONTH AS WE WERE FOR THE PREVIOUS CALL!
*       SITINF - Raw information regarding each site
*	DatArray - current and previous days' values
*	TRACE	- trace threshold
*	BETA	- Parameter vector
*	THETA	- array containing parameters in nonlinear transformations,
*		  arranged to line up with their covariates. Values in
*		  positions NP(6)+1 to NP(7) actually contain derivatives of 
*		  nonlinear functions wrt the parameters.
*	TRENDFN	- Value of long-term trend, + its derivatives wrt any
*		  nonlinear parameters we're trying to estimate. Element 0
*		  is the value, 1 and 2 are derivatives.
*	MONFN	- FUNCTION used for monthly contributions
*	DYPRED	- Predictor (& possible derivatives wrt nonlinear parameters)
*                 corresponding to previous days.
*       EXTDY   - value of external forcing variable for today
*       Distance- array of inter-site distances. First slice is X-separation,
*                 second slice is Y-separation and third slice is Euclidean
*                 separation
*	TMP	- temporary storage
*******************************************************************************
      DOUBLE PRECISION COVS(NSITES,MXP),SITINF(NSITES,MXP,0:3)
      DOUBLE PRECISION DatArray(NSITES,NVARS,0:10)
      DOUBLE PRECISION TRENDFN(0:2),MONFN,DYPRED(0:3),EXTDY
      DOUBLE PRECISION TRACE,BETA(0:MXP),THETA(MXP,3)
      Double precision, intent(in) :: Distance(3,NSITES,NSITES)
      DOUBLE PRECISION TMP
*
*	Remember year and month effects in case we want to use them again
*	(and start them all off as zeroes). Also remember whether we've
*       checked the model.
*
      SAVE OLDYR,OLDMN,YMISS,MMISS,DMISS
      DATA OLDYR,OLDMN,YMISS,MMISS,DMISS /5*0/
*
*	If we're checking for missing data, initialise MISSFL
*       (assume no missing data); and if we're FORCEing recalculation
*       of everything, ensure that other missing flags are reinitialised.
*       Also calculate inter-site distances, if these are needed by 
*      a weighting scheme subsequently
*

      IF (MISSST.EQ.1) THEN
       MISSFL(1:NSITES) = 0
       IF ((FORCE.EQ.1).OR.(RECALC.GT.0)) THEN 
        YMISS = 0
        MMISS = 0
        DMISS = 0
       ENDIF
      ENDIF
*
*		This routine does all sites for 1 day; therefore we only need 
*       to compute year and month effects, and `external' daily effects, 
*       once (save on flops). Don't do any parameters that we
*		don't have to! Values are only set if the year/month changes.
*       This means that these elements of COVS must NOT be changed 
*       between calls. First years (trends and `external forcings'
*       are dealt with differently. NB the need to deal with missing
*       data when external forcings are involved - effectively, we
*       just go onto the next year).
*

      DO 50 J=NP(1)+1,NP(2)
       IF ((RECALC.GT.0).AND.(THETA(J,1).GT.1.0D8)) GOTO 50
       IF ((YEAR.NE.OLDYR).OR.(FORCE.EQ.1)) THEN
        IF (MOD(IABS(COVCODE(J)),1000).LE.50) THEN
         CALL TRENDSB(TRENDFN,YEAR,J,COVCODE(J),THETA,ICHECK,MXP,IFAIL)
         IF (IFAIL.NE.0) RETURN
        ELSE
         EXTLAG = COVCODE(J)/1000
         CALL EXTSET(ExtUnits(1),1,MOD(IABS(COVCODE(J)),1000)-50,
     +                    DAY,MONTH,YEAR,EXTLAG,TRENDFN(0),YMISS,IFAIL)
         IF (IFAIL.NE.0) RETURN
        ENDIF
        COVS(1:NSITES,J) = TRENDFN(0)
*
*	There may be nonlinear functions floating round here,
*	so we'll take the opportunity to set the derivatives for their
*	parameters if we're fitting models (rather than simulating). NB 
*       we know whether or not a nonlinear parameter is to do with the 
*       current predictor because the predictor is given by COVCODE/1000 
*       (integer division). NB ALSO: for the moment, we put the value of 
*       the derivative in the row of THETA corresponding to the nonlinear 
*       parameter (this does not conflict with anything else, since 
*       actual parameter values of THETA are stored in the row of the 
*       covariate to which they correspond). This is necessary to keep 
*       track of what the derivatives actually are when interaction terms 
*       are involved, without having to recalculate every time.
*       
        IF (SIM.EQ.0) THEN
         DO 57 K=NP(6)+1,NP(7)
          IF (COVCODE(K)/1000.EQ.J) THEN
           IDX2 = MOD(COVCODE(K),1000)
           THETA(K,IDX2) = TRENDFN(IDX2)
          ENDIF
 57      CONTINUE
        ENDIF
       ENDIF
 50   CONTINUE
      OLDYR = YEAR
*
*       Monthly, then daily, effects. Skip if we're only RECALCulating 
*       nonlinear contributions. Again, `external' monthly contributions
*       are dealt with differently.
*
      IF ((RECALC.GT.0).AND.(FORCE.EQ.0)) GOTO 99

      DO 60 J=NP(2)+1,NP(3)
       IF ((MONTH.NE.OLDMN).OR.(FORCE.EQ.1)) THEN
        IF (MOD(IABS(COVCODE(J)),1000).LE.50) THEN
         CALL MONSB(MONFN,MONTH,COVCODE(J),IFAIL)
         IF (IFAIL.NE.0) RETURN
        ELSE
         EXTLAG = COVCODE(J)/1000
         CALL EXTSET(ExtUnits(2),2,MOD(IABS(COVCODE(J)),1000)-50,
     +                     DAY,MONTH,YEAR,EXTLAG,MONFN,MMISS,IFAIL)
         IF (IFAIL.NE.0) RETURN
        ENDIF
        COVS(1:NSITES,J) = MONFN
       ENDIF
 60   CONTINUE

      DO 70 J=NP(3)+1,NP(4)
       IF (.NOT.(MOD(IABS(COVCODE(J)),1000).LE.10)) THEN
*
*     Have to pass DatArray array over here, even though we're not going
*     to use it.
*
        IF (MOD(IABS(COVCODE(J)),1000).LE.50) THEN
         CALL DYSET(DAY,MONTH,DatArray,0,NSITES,NVARS,RespIdx,
     +     AllowIncAvge,TRACE,COVCODE(J),THETA,PWTIDX,
     +     Distance,ICHECK,IFAIL,DYPRED,MXP)
         IF (IFAIL.NE.0) RETURN
         EXTDY = DYPRED(0)
        ELSE
         EXTLAG = COVCODE(J)/1000
         CALL EXTSET(ExtUnits(3),3,MOD(IABS(COVCODE(J)),1000)-50,
     +                     DAY,MONTH,YEAR,EXTLAG,EXTDY,DMISS,IFAIL)
         IF (IFAIL.NE.0) RETURN
        ENDIF
        COVS(1:NSITES,J) = EXTDY
       ENDIF
 70   CONTINUE

 99   OLDMN = MONTH
*
*     Missing external predictors? If so, no point continuing - set
*     everything to missing (do the entire MISSFL array at once) and 
*     return EXCEPT if we're just recalculating - in which case the 
*     missing flags are probably the result of a SAVE on a previous 
*     call that has not been reset
*
      IF( (RECALC.EQ.0).AND.
     +    ((YMISS.EQ.1).OR.(MMISS.EQ.1).OR.(DMISS.EQ.1))) THEN 
       IF (MISSST.EQ.1) MISSFL = 1
       IF (SIM.EQ.1) IFAIL = 70
       RETURN
      ENDIF
*
*       That's it for effects that just have to be computed once. Now 
*       the remainder. Start with main effects. If we're recomputing
*       at a site which has missing predictor data according to the
*       MISSFL array, skip the whole section.
*
      DO 100 I=1,NSITES
       IF ((MISSST.EQ.0).AND.(MISSFL(I).EQ.1)) GOTO 100
*
*	Site effects. May be some nonlinear parameters to estimate here, so
*	as before we copy the derivatives into the appropriate place in
*	THETA (NB we're now doing a site at a time, so these elements of
*	THETA will get changed around, but that doesn't matter since they're
*	not used anywhere else).
*
       DO 110 J=1,NP(1)
        IF ((RECALC.GT.0).AND.(THETA(J,1).GT.1.0D8)) GOTO 110
        COVS(I,J) = SITINF(I,COVCODE(J),0)
        DO 111 K=NP(6)+1,NP(7)
         IF (COVCODE(K)/1000.EQ.J) THEN
          IDX2 = MOD(COVCODE(K),1000)
          THETA(K,IDX2) = SITINF(I,COVCODE(J),IDX2)
         ENDIF
 111     CONTINUE
 110   CONTINUE
*
*       Daily effects. At this point we're just setting previous days' 
*       values - 'external' daily effects don't vary between sites and 
*       have already been set. If any previous values (or current values
*       of other required variables) are missing, set a MISSING flag for 
*       this site and go on to the next. The effect will be not to alter 
*       anything that was previously in COVS, which should therefore be
*       initialised in the calling program if used for simulation. 
*
       Do J=1,NVARS
        Do K=0,MXLAG(J)
         IF (DatArray(I,J,K).LT.-1.0D99) THEN
          IF (MISSST.EQ.1) MISSFL(I) = 1
          GOTO 100
         ENDIF
        End Do
       End Do

       DO 130 J=NP(3)+1,NP(4)
        IF (MOD(IABS(COVCODE(J)),1000).LE.10) THEN
*
*     There may be nonlinear parameters to estimate here, associated 
*     (for example) with weighted averages of previous days' values.
*     We can tell whether or not this is the case, because there will
*     be non-zero pointers in PWTIDX. If so, we copy their derivatives 
*     into a spare row of THETA (this will get overwritten as we move 
*     through the sites, but that's OK as we're basically using it for 
*     temporary storage, a site at a time). There's an additional 
*     complication here, which is that we will generally have a single
*     parametrisation that covers various different lagged values. This
*     will give rise to several different sets of derivatives, so we need
*     more than one spare row of THETA. Creativity required: if the
*     current predictor is the one that the weighting scheme is nominally
*     attached to, use a spare row in the `nonlinear' part of THETA. 
*     Otherwise, use the row corresponding to the nonlinear part (it
*     doesn't get used for anything else).
*
         WTSCHM = MOD(COVCODE(J)/10000,10)
         VARNUM = COVCODE(J) / 1000000
         IF (RECALC.GT.0) THEN
*
*     Need to break this down into 2 bits, to avoid referencing
*     PWTIDX(0,*,*) which isn't defined
*
          IF (WTSCHM.EQ.0) GOTO 130
          IF ((PWTIDX(WTSCHM,1,VARNUM).EQ.0).AND.
     +        (PWTIDX(WTSCHM,2,VARNUM).EQ.0).AND.
     +        (PWTIDX(WTSCHM,3,VARNUM).EQ.0)) GOTO 130
         ENDIF
         CALL DYSET(DAY,MONTH,DatArray,I,NSITES,NVARS,RespIdx,
     +              AllowIncAvge,TRACE,COVCODE(J),THETA,PWTIDX,
     +              Distance,ICHECK,IFAIL,DYPRED,MXP)
         IF (IFAIL.NE.0) RETURN
*
*       DYSET can pick up missing values that haven't yet been marked -
*       need to trap these
*
         IF (DYPRED(0).LT.-1.0D99) THEN
          IF (MISSST.EQ.1) MISSFL(I) = 1
          GOTO 100
         ENDIF
         COVS(I,J) = DYPRED(0)
         IF (WTSCHM.EQ.0) GOTO 130
         DO 135 K=1,3
          IF (PWTIDX(WTSCHM,K,VARNUM).GT.0) THEN
           IF (J.EQ.PWTIDX(WTSCHM,0,VARNUM)) THEN
            THETA(PWTIDX(WTSCHM,K,VARNUM),K) = DYPRED(K)
           ELSE
            THETA(J,K) = DYPRED(K)
           ENDIF
          ENDIF
 135     CONTINUE
        ENDIF
 130   CONTINUE
*
*       Now 2-way interactions. Just do the ones for which main
*	effects have been defined (not all when we're just recalculating).
*	If we're recalculating, we only do it when there's a nonlinear 
*	parameter in an associated main effect.
*
       DO 140 J=NP(4)+1,NP(5)
        IF ((RECALC.GT.0).AND.
     +      (THETA(TWO(J,1),1).GT.1.0D8).AND.
     +      (THETA(TWO(J,2),1).GT.1.0D8) ) GOTO 140
        COVS(I,J) = COVS(I,TWO(J,1))*COVS(I,TWO(J,2))
 140   CONTINUE  
*
*       And 3-way interactions
*
       DO 150 J=NP(5)+1,NP(6)
        IF ((RECALC.GT.0).AND.
     +      (THETA(THREE(J,1),1).GT.1.0D8).AND.
     +      (THETA(THREE(J,2),1).GT.1.0D8).AND.
     +      (THETA(THREE(J,3),1).GT.1.0D8) ) GOTO 150
        COVS(I,J) = COVS(I,THREE(J,1))
     +              *COVS(I,THREE(J,2))
     +              *COVS(I,THREE(J,3))
 150   CONTINUE  
*
*	Now the 'dummy' elements corresponding to derivatives of 
*	nonlinear parameters. We need to take the derivatives of the
*	nonlinear function itself (set into spare bits of THETA earlier),
*	and scale as appropriate: they appear in the linear predictor as 
*	multiples of other parameters and (where interactions are 
*	involved) covariates. NB this only needs to be done during 
*       model fitting, not during simulation. 
*	
       IF (SIM.EQ.0) THEN
        DO 160 J=NP(6)+1,NP(7)
         IDX1 = COVCODE(J)/1000
         IDX2 = MOD(COVCODE(J),1000)
         CALL NLCNTR(IDX1,NP,BETA,TWO,THREE,COVS,I,NSITES,MXP,TMP,IFAIL)
         IF (IFAIL.NE.0) RETURN
         COVS(I,J) = THETA(J,IDX2)*TMP
*
*       Parameters in weighting schemes are fiddly, because they may 
*       relate to several main effects rather than just one. So 
*       now we have to go through all the *other* previous days' 
*       predictors looking for contributions to this parameter ...
*
         WTSCHM = 0
         IF ((IDX1.GT.NP(3)).AND.(IDX1.LE.NP(4))) THEN
          IF (MOD(IABS(COVCODE(J)),1000).LE.10) THEN
           WTSCHM = MOD(COVCODE(IDX1)/10000,10)
           VARNUM = COVCODE(IDX1) / 1000000
          ENDIF
         ENDIF
         IF (WTSCHM.NE.0) THEN
          DO 165 IDX1 = NP(3)+1,NP(4)
*
*	Only calculate contributions from the current weighting scheme ...
*
           IF (MOD(COVCODE(IDX1)/10000,10).NE.WTSCHM) GOTO 165
*
*	... and the current transformation ...
*
           IF (MOD(IABS(COVCODE(IDX1)),1000).GT.10) GOTO 165
*
*	... and the current variable number
*
           IF (COVCODE(IDX1)/1000000 .NE. VARNUM) GOTO 165
*
*	COVCODE(J)/1000 is the covariate to which the parameter
*	was originally attached, so we don't need it again
*
           IF (IDX1.EQ.COVCODE(J)/1000) GOTO 165
           CALL NLCNTR(IDX1,NP,BETA,TWO,THREE,COVS,I,NSITES,MXP,
     +                                                      TMP,IFAIL)
           IF (IFAIL.NE.0) RETURN
*
*     For these contributions, the derivatives of THETA are stored
*     in the rows corresponding to main effects
*
           COVS(I,J) = COVS(I,J) + (THETA(IDX1,IDX2)*TMP)
 165      CONTINUE
         ENDIF
 160    CONTINUE
       ENDIF
 100  CONTINUE

      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE NLCNTR(PREDNO,NP,BETA,TWO,THREE,COVS,SITE,
     +                                        MXN,MXP,XOUT,IFAIL)
*******************************************************************************
*     To calculate contributions to a `dummy' predictor representing 
*     derivative wrt a nonlinear parameter, that are associated with 
*     a given main effect in a GLM. Arguments:
*
*     PREDNO -  Number of main effect required. INTEGER, input
*     NP     -  Index of model predictors in different categories.
*               INTEGER, input.
*     BETA   -  Curremt estimate of coefficient vector. DOUBLE, input.
*     TWO -     Table of 2-way interactions. INTEGER, input.
*     THREE -   Guess.
*     COVS   -  Array of covariate values. DOUBLE, input.
*     SITE   -  Current site.
*     MXN,MXP - Dimensioning. INTEGER, input.
*     XOUT   -  Required contribution. DOUBLE, output.
*     IFAIL -   Error flag
*******************************************************************************
      INTEGER PREDNO,NP(10),MXN,MXP,SITE,IFAIL
      INTEGER TWO(MXP,2),THREE(MXP,3)
      DOUBLE PRECISION BETA(0:MXP),COVS(MXN,MXP),XOUT
*******************************************************************************
*     Additional INTEGERs
*     ^^^^^^^^^^^^^^^^^^^
*     I    - counter
*******************************************************************************
      INTEGER I
*******************************************************************************
*     Additional CHARACTERs
*     ^^^^^^^^^^^^^^^^^^^^^
*     MESSAGE   - Error messages
*******************************************************************************
      CHARACTER MESSAGE*200

      IFAIL = 0
*
*     Check that PREDNO is a main effect
*
      IF (PREDNO.GT.NP(4)) GOTO 990

      XOUT = BETA(PREDNO)
      DO 100 I=NP(4)+1,NP(5)
       IF ( (TWO(I,1).EQ.PREDNO) ) THEN
        XOUT = XOUT + (BETA(I) * COVS(SITE,TWO(I,2)))
       ELSEIF ( (TWO(I,2).EQ.PREDNO) ) THEN
        XOUT = XOUT + (BETA(I) * COVS(SITE,TWO(I,1)))
       ENDIF
 100  CONTINUE

      DO 120 I=NP(5)+1,NP(6)
       IF ( (THREE(I,1).EQ.PREDNO) ) THEN
        XOUT = XOUT + 
     +      (BETA(I) * COVS(SITE,THREE(I,2)) * 
     +                 COVS(SITE,THREE(I,3)))
       ELSEIF ( (THREE(I,2).EQ.PREDNO) ) THEN
        XOUT = XOUT + 
     +      (BETA(I) * COVS(SITE,THREE(I,1)) * 
     +                 COVS(SITE,THREE(I,3)))
       ELSEIF ( (THREE(I,3).EQ.PREDNO) ) THEN
        XOUT = XOUT + 
     +      (BETA(I) * COVS(SITE,THREE(I,1)) *
     +                 COVS(SITE,THREE(I,2)))
       ENDIF
 120  CONTINUE

      RETURN

*
*     Error trapping
*
 990  MESSAGE = '****ERROR**** in routine NLCNTR, PREDNO is not '//
     +          'a main effect.'
      CALL INTPR(TRIM(MESSAGE),-1,0,0)
      IFAIL = 999
      
      END
******************************************************************************
******************************************************************************
******************************************************************************
      Subroutine GetHandle(UnitNumber)
*******************************************************************************
*     Finds an available file handle. Search starts at 20, to avoid 
*     potential clashes with processor-dependent reserved handles
*******************************************************************************
      Integer, intent(out) :: UnitNumber
*******************************************************************************
*     LOGICAL variables
*     ^^^^^^^^^^^^^^^^^
*     IsOpen    TRUE if unit is already connected, FALSE otherwise
*******************************************************************************
      Logical IsOpen
      
      UnitNumber = 20
 10   INQUIRE(UNIT=UnitNumber,OPENED=IsOpen)
      IF (IsOpen) THEN
       UnitNumber = UnitNumber + 1
       GOTO 10
      ENDIF
      
      End Subroutine GetHandle
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE SCROPN(ACODE,FCODE,NFIELDS,FILNO)
*******************************************************************************
*       Finds an available file handle and opens a SCRATCH file
*       connected to that handle.
*******************************************************************************
*       INTEGER variables
*       ^^^^^^^^^^^^^^^^^
*	ACODE	Integer, input: access code (1 for direct, 2 for sequential)
*	FCODE   Integer, input: format code (1 for formatted, 2 for unformatted)
*       NFIELDS	Number of fields per record - mandatory if ACCESS = 
*               'DIRECT', ignored otherwise. If used, each record is 
*               long enough for NFIELDS double precision entries. Input.
*       FILNO	Number of file when it is opened. This is determined automatically
*               by subroutine GetHandle.
*       BYTELN    Record length required to store a single double precision
*                 number
*       RECL      Total record length required for NFIELDS numbers
*******************************************************************************
      INTEGER, intent(in) :: ACODE, FCODE, NFIELDS
      Integer, intent(out) :: FILNO
      Integer BYTELN, RECL
*******************************************************************************
*       CHARACTER variables
*       ^^^^^^^^^^^^^^^^^^^
*       FORM	Format status - either 'FORMATTED' or 'UNFORMATTED'
*******************************************************************************
      CHARACTER FORM*11
*******************************************************************************

*
*       Define FORM
*
      If (FCODE.EQ.1) then
       FORM = '  FORMATTED'
      Else If (FCODE.EQ.2) then
       FORM = 'UNFORMATTED'
      End If
*
*       Locate a file handle
*
      CALL GetHandle(FILNO)
*
*       Figure out record lengths
*
      INQUIRE(iolength=BYTELN) 1.0D0
      RECL = NFIELDS*BYTELN
*
*       Most compilers ignore the RECL argument if ACCESS=SEQUENTIAL,
*       but I've had problems with gfortran - hence slightly 
*       clunky IF structure
*
      IF (ACODE.EQ.1) THEN
       OPEN(FILNO,STATUS='SCRATCH',ACCESS='DIRECT',FORM=FORM,RECL=RECL)
      ELSE IF (ACODE.EQ.2) THEN 
       OPEN(FILNO, STATUS='SCRATCH', ACCESS='SEQUENTIAL', FORM=FORM)
      ENDIF

      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE BCALC(P,ALPHA,A,B)
***************************************************************************
*	To find the value of B in the conditionally independent 2 weather
*	state model for rainfall occurrence (called from routine WDALLOC),
*	such that the overall probability of a site being wet is P (input).
*	ALPHA is the probability of a wet day (input). log(A) is the
*	coefficient of x (0/1 weather state) in the model (input). B is 
*	output.
*	There is a unique value of B which makes the marginal probability
*	correct. It is obtained by taking one of the roots (we know
*	which one) of a quadratic. We use the numerically robust algorithm
*	suggested by Press et al, Numerical Recipes, Section 5.6.
*
*	This routine is called by fit_logi.f and simrain.f
***************************************************************************
      DOUBLE PRECISION P,ALPHA,A,B
***************************************************************************
*	Extra DOUBLEs
*	^^^^^^^^^^^^^
*	TMP	- Temporary storage
*	QUADA }	- Coefficients in quadratic equation to determine B
*	QUADB }	
*	QUADC }
***************************************************************************
      DOUBLE PRECISION QUADA,QUADB,QUADC,TMP
***************************************************************************
      QUADA = A*P
      QUADB = 1.0D0 + (ALPHA*(A-1.0D0)) - (P*(A+1.0D0))
      QUADC = P-1.0D0
      TMP = DSQRT( (QUADB**2)-(4.0D0*QUADA*QUADC) )
      TMP = -( QUADB + DSIGN(TMP,QUADB) )/ 2.0D0
      IF (QUADB.LT.0.0D0) THEN
       B = TMP/QUADA
      ELSE
       B = QUADC/TMP
      ENDIF
      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE ANSPRM(NU,MEAN,VAR,IFAIL)
***************************************************************************
*       Finds mean and variance of a normal distribution such that when
*       cubed, the result will behave like a Gamma random variable with
*	a mean of 1 and a variance of 1/NU. This basically involves 
*	(i) working with a linear transformation of the underlying normal 
*	and finding the mean and variance of that - there is an explicit 
*	formula for the square of the mean at each value of the variance 
*	so we then have to just solve numerically for the variance (ii) 
*	backtransforming to get the numbers we want. 
***************************************************************************
*       DOUBLE PRECISION variables
*       ^^^^^^^^^^^^^^^^^^^^^^^^^^
*       NU      - shape parameter of gamma distribution
*       MEAN    - required mean of Anscombe residuals
*       VAR     - required variance
*       LLIM }  - lower and upper limits used when solving for MEAN and
*       ULIM }    VAR
*       TOL     - convergence criterion when solving for MEAN and VAR
*       F1,F2   - values of function used to solve for MEAN and VAR
*       MID     - new estimate of VAR (using bisection - tried regula falsi
*                 but it was hugely unstable)
*       FMID    - function value at MID
***************************************************************************
      DOUBLE PRECISION NU,MEAN,VAR
      DOUBLE PRECISION LLIM,ULIM,TOL,F1,F2,MID,FMID
***************************************************************************
*       INTEGER variables
*       ^^^^^^^^^^^^^^^^^
*       IFAIL   Error flag
***************************************************************************
      Integer, intent(out) :: IFAIL      
***************************************************************************
*               CHARACTER variables
*               ^^^^^^^^^^^^^^^^^^^
*       MESSAGE - messages to screen
***************************************************************************
      CHARACTER MESSAGE*255
***************************************************************************

      IFAIL = 0
      TOL = 1.0D-8
*
*       For a solution to exist there is an upper bound on VAR. 
*
      LLIM = 1.0D-8
      ULIM = (1.0D0/(15.0D0*NU))**(1.0D0/3.0D0)
*
*       Evaluate the function at the boundaries
*
      VAR = LLIM
      CALL ANSFN(F1,MEAN,LLIM,NU,IFAIL)
      IF (IFAIL.NE.0) RETURN
      CALL ANSFN(F2,MEAN,ULIM,NU,IFAIL)
      IF (IFAIL.NE.0) RETURN
*
*       Check that the root is bracketed ...
*
      IF (F1*F2.GT.0.0D0) THEN
       WRITE(MESSAGE,1)
       CALL INTPR(TRIM(MESSAGE),-1,0,0)
       IFAIL = 999
       RETURN
      ENDIF
*
*       Compute midpoint (unsophisticated but stable, which regula 
*       falsi WASN'T)
*
 20   MID = ( LLIM+ULIM )/2.0D0
      CALL ANSFN(FMID,MEAN,MID,NU,IFAIL)
      IF (IFAIL.NE.0) RETURN
      IF (FMID*F1.LT.0.0D0) THEN
       ULIM = MID
       F2 = FMID
      ELSEIF (FMID*F2.LT.0.0D0) THEN
       LLIM = MID
       F1 = FMID
      ENDIF
      VAR = MID
      IF ( ((ULIM-LLIM).GT.TOL).OR.(DABS(F2-F1).GT.TOL) ) GOTO 20

      RETURN

 1    FORMAT('****ERROR**** Can''t solve for the mean and variance',
     +' of Anscombe residuals.')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE ANSFN(F,MEAN,VAR,NU,IFAIL)
***************************************************************************
*       Evaluates the function F, which we're trying to set to zero
*       when calculating the mean and variance of Anscombe residuals
*       for a given value of NU. F is strictly a function of VAR only;
*       we also calculate the corresponding value of MEAN in this
*       subroutine. A,B and C are coefficients in a quadratic equation
*       for MU squared, which always has positive roots providing 
*       VAR is positive. 
***************************************************************************
      DOUBLE PRECISION F,MEAN,VAR,NU,A,B,C
***************************************************************************
*       INTEGER variables
*       ^^^^^^^^^^^^^^^^^
*       IFAIL   Error flag
***************************************************************************
      Integer, intent(out) :: IFAIL      
***************************************************************************
*               CHARACTER variables
*               ^^^^^^^^^^^^^^^^^^^
*       MESSAGE - messages to screen
***************************************************************************
      CHARACTER MESSAGE(2)*255
***************************************************************************
      
      IFAIL = 0
*
*       First calculate MEAN. Context (+probably maths later on, 
*       if I thought about it) dictates MEAN > 0
*
      A = 9.0D0*VAR
      B = 36.0D0*(VAR**2)
      C = (15.0D0*(VAR**3)) - (1.0D0/NU)
*
*       Compute square of mean by solving quadratic
*
      MEAN = (-B+DSQRT( (B**2)-(4.0D0*A*C) ) )/(2.0D0*A)
*
*       If < 0, it's probably a rounding error and should be zero
*
      IF (MEAN.LT.0.0D0) THEN
       IF (DABS(MEAN).LT.1.0D-10) THEN
        MEAN = 0.0D0
       ELSE
        WRITE(MESSAGE,1)
        CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
        IFAIL = 999
        RETURN
       ENDIF
*
*       Otherwise we just need to square root whatever we've got
*
      ELSE
       MEAN = DSQRT(MEAN)
      ENDIF
*
*     Now calculate the corresponding function value
*
      F = MEAN*( (3.0D0*VAR) + (MEAN**2) ) - 1.0D0

 1    FORMAT('****ERROR**** I''ve just tried to square root a ',
     +'negative number in routine ANSFN.',/,
     +'It''s too big to be a rounding error.')
      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE PRGERR(ROUTIN,SRCFIL,ERRNO,IFAIL)
*
*	To spit out an error message when programming errors are
*	detected. Arguments (all input):
*	ROUTIN	- Name of routine where error was detected. CHARACTER*8
*	SRCFIL	- Name of source file containing that routine. CHARACTER*12
*	ERRNO	- Identifies error. Possible values:
*			1 - Incompatible dimensions of routine arguments
*			    and internal variables. Probably indicates 
*			    a PARAMETER statement which is wrong in the
*			    routine.
*			2 - Insufficient internal storage in a routine
*			    to cope with size of arguments
***************************************************************************
      INTEGER ERRNO,I,IFAIL
      CHARACTER ROUTIN*8,SRCFIL*12,MESSAGE(4)*255
***************************************************************************

      IFAIL = 997
      IF (ERRNO.EQ.1) THEN
       WRITE(MESSAGE,1) ROUTIN,SRCFIL
       Do I=1,3
        CALL INTPR(TRIM(MESSAGE(I)),-1,0,0)
       End Do
      ELSEIF (ERRNO.EQ.2) THEN
       WRITE(MESSAGE,2) ROUTIN,SRCFIL
       Do I=1,4
        CALL INTPR(TRIM(MESSAGE(I)),-1,0,0)
       End Do
      ELSE
       WRITE(MESSAGE,99) ROUTIN,SRCFIL
       Do I=1,2
        CALL INTPR(TRIM(MESSAGE(I)),-1,0,0)
       End Do
       IFAIL = 999
      ENDIF

 1    FORMAT('*****ERROR***** In routine ',A8,' (source file ',
     +A12,'),',/,'some internal array dimensions are incompatible ',
     +'with the routine''s',/,'arguments. Check the PARAMETER ',
     +'statements in this routine, and recompile.')
 2    FORMAT('*****ERROR***** In routine ',A8,' (source file ',
     +A12,'),',/,'there is insufficient internal storage to cope ',
     +'with this problem.',/,'You should increase the dimension ',
     +'of some or all of the internal arrays for',
     +/,'this routine, and recompile.')
 99   FORMAT('*****ERROR***** Unidentified programming error ',
     +'found in routine ',A8,/,'(source file ',A12,')')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE CORSET(FILNO,SPMOD,RHO,CORMAT,NSITES,
     +                           SITINF,NATTR,SCODES,MXP,IFAIL)
****************************************************************************
*	Sets up a spatial dependence matrix according to model SPMOD.
*       NB it is assumed that all required parameters etc. have been 
*       checked, so no checking is done internally here. NB also 
*       the correlation between a site and itself is ALWAYS 1. 
*       Arguments:
*
*	FILNO	(integer, input). Channel number of file from which 
*               to read dependences and distances, if necessary.
*	SPMOD	(integer, input). Identifies which model to use for the
*		spatial dependence.
*	RHO	(double precision, input). Vector of spatial model 
*		parameters.
*	CORMAT	(double precision, output). Matrix of computed inter-site
*		dependences.
*	NSITES	(integer, input) No. of sites in simulation.
*	SITINF  (double precision, input). Array of site attributes. Used
*		to calculate distance-dependent dependences.
*       NATTR   (integer, input). Number of attributes defined in SITINF
*	SCODES 	(character, input). Array of short site codes.
*	MXP	(integer,input). Dimensioning.
*       IFAIL   (integer, output). Error code.
****************************************************************************
      INTEGER FILNO,SPMOD,NSITES,NATTR,MXP
      DOUBLE PRECISION RHO(MXP),CORMAT(NSITES,NSITES)
      Double precision, intent(in) :: SITINF(NSITES,MXP,0:3)
      CHARACTER SCODES(NSITES)*4
****************************************************************************
*	Additional INTEGERs
*	-------------------
*	I,J,K	counters
*       IFAIL   error flag
****************************************************************************
      INTEGER I,J,K,IFAIL
****************************************************************************
*	Additional DOUBLEs
*	------------------
*	CORR	Correlation read from input file
*       XSEP }  Inter-site separation in terms of the first two 
*       YSEP }  attributes defined in SITINF
*       DIST    Inter-site distance
****************************************************************************
      DOUBLE PRECISION CORR,XSEP,YSEP,DIST
****************************************************************************
*	Additional CHARACTERs
*	---------------------
*	SCODE1} Site codes read from input file
*	SCODE2}
*       MESSAGE For screen output
****************************************************************************
      CHARACTER*4 SCODE1,SCODE2,MESSAGE(5)*255
      
      IFAIL = 0
*
*	'Independence' case. We treat it as a special case of the 
*	other models.
*
      IF (SPMOD.EQ.0) THEN
       DO 90 I=1,NSITES
        CORMAT(I,I) = 1.0D0
        DO 91 J=I+1,NSITES
         CORMAT(I,J) = 0.0D0
         CORMAT(J,I) = 0.0D0
 91     CONTINUE
 90    CONTINUE
*
*	Here's the `empirical correlation' case. We make sure all
*	the pairs of sites we're simulating for are represented in
*	the file we're reading from - initialise all correlations
*	to something impossible, then read in from file, and then 
*	pass through again to check that nothing impossible is left.
*       The value of -9.999 is the same as that used in the fitting
*       programs to indicate a pair of sites with no observations
*       in common. 
*
      ELSEIF (SPMOD.EQ.1) THEN
       READ(FILNO,*)
       DO 100 I=1,NSITES
        CORMAT(I,I) = 1.0D0
        DO 101 J=I+1,NSITES
         CORMAT(I,J) = -9.999D0
         CORMAT(J,I) = -9.999D0
 101    CONTINUE
 100   CONTINUE
 
 110   READ(FILNO,'(A4,4X,A4,4X,F7.4)',END=199) SCODE1,SCODE2,CORR
       DO 120 I=1,NSITES
        DO 121 J=I+1,NSITES
         IF ( ((SCODE1.EQ.SCODES(I)).AND.(SCODE2.EQ.SCODES(J))).OR.
     +      ((SCODE1.EQ.SCODES(J)).AND.(SCODE2.EQ.SCODES(I))) ) THEN
          CORMAT(I,J) = CORR
          CORMAT(J,I) = CORR
          GOTO 110
         ENDIF
 121    CONTINUE
 120   CONTINUE
       GOTO 110

 199   CLOSE(FILNO)
       DO 130 I=1,NSITES
        DO 131 J=I+1,NSITES
         IF (DABS(CORMAT(I,J)).GT.1.0D0) THEN
          WRITE(MESSAGE,10) SCODES(I),SCODES(J)
          Do K=1,5
           CALL INTPR(TRIM(MESSAGE(K)),-1,0,0)
          End Do
          IFAIL = 200
          RETURN
         ENDIF
 131    CONTINUE
 130   CONTINUE
*
*	Here's the `constant dependence' case
*
      ELSEIF (SPMOD.EQ.2) THEN
       DO 200 I=1,NSITES
        CORMAT(I,I) = 1.0D0
        DO 201 J=I+1,NSITES
         CORMAT(I,J) = RHO(1)
         CORMAT(J,I) = RHO(1)
 201    CONTINUE
 200   CONTINUE
*
*       Exponential and powered exponential decay to zero or 
*       a threshold, possibly with a nugget (but not both
*       nugget and threshold at present!)
*
      ELSEIF ((SPMOD.GE.3).AND.(SPMOD.LE.8)) THEN
       XSEP = 0.0D0
       YSEP = 0.0D0
       DO 220 I=1,NSITES
        CORMAT(I,I) = 1.0D0
        DO 221 J=I+1,NSITES
         IF (NATTR.GE.1) XSEP=SITINF(I,1,0)-SITINF(J,1,0)
         IF (NATTR.GE.2) YSEP=SITINF(I,2,0)-SITINF(J,2,0)
         DIST = DSQRT( (XSEP*XSEP) + (YSEP*YSEP) ) 
         IF (SPMOD.EQ.3) THEN
          CORMAT(I,J) = DEXP(-RHO(1)*DIST)
         ELSEIF (SPMOD.EQ.4) THEN
          CORMAT(I,J) = RHO(2) + ( (1.0D0-RHO(2))*DEXP(-RHO(1)*DIST))
         ELSEIF (SPMOD.EQ.5) THEN
          CORMAT(I,J) = DEXP(-RHO(1)*(DIST**RHO(2)))
         ELSEIF (SPMOD.EQ.6) THEN
          CORMAT(I,J) = RHO(3) + 
     +            ( (1.0D0-RHO(3)) * DEXP(-RHO(1)*(DIST**RHO(2))) )
         ELSE IF (SPMOD.EQ.7) THEN
          CORMAT(I,J) = RHO(2) * DEXP(-RHO(1)*DIST)
         ELSEIF (SPMOD.EQ.8) THEN
          CORMAT(I,J) = RHO(3) * DEXP(-RHO(1)*(DIST**RHO(2)))
         ENDIF
         CORMAT(J,I) = CORMAT(I,J)
 221    CONTINUE
 220   CONTINUE
*
*	No other valid options defined at present
*
      ELSE
       IFAIL = 999
       RETURN
      ENDIF
      
      RETURN
            
 10   FORMAT('****ERROR**** I don''t have a valid spatial ',
     +'correlation between sites',
     +/A4, ' and ',A4,'. This is either because there is no ',
     +'entry in the input file ',
     +/'for this site pair, or because the correlation in the ',
     +'file exceeds 1 in',
     +/'magnitude (which itself could be a missing value flag ',
     +'of -9.999).',/,'Please sort this out!')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE DATRD(FILNO,CURRENT,MissVal,SCODES,NSITES,NVARS,
     +                                        YY,MM,DD,DONE,IFAIL)
*
*       Reads in rows of data from file connected to unit FILNO 
*       and allocates them to array CURRENT if they correspond to 
*       date DD/MM/YY. If DD,MM and YY are all zero, read the first 
*       record
*
      INTEGER FILNO,YY,MM,DD,NSITES,NVARS,DONE,IFAIL
      DOUBLE PRECISION MissVal,CURRENT(NSITES,NVARS)
      CHARACTER*4 SCODES(NSITES)
******************************************************************************
*       Additional INTEGERS
*       ^^^^^^^^^^^^^^^^^^^
*       DF      - Day, as read from input file
*       MF      - Month, ditto
*       YF      - Year, ditto
*       IOERR   - Input file error flag - will be < 0 when EOF is attained
******************************************************************************
      INTEGER DF,MF,YF,IOERR,I,J
******************************************************************************
*       Additional DOUBLEs
*       ^^^^^^^^^^^^^^^^^^
*       FromFile   - values read from file
******************************************************************************
      DOUBLE PRECISION FromFile(NVARS)
******************************************************************************
*       Additional CHARACTERs
*       ^^^^^^^^^^^^^^^^^^^^^
*       CODEF           - site code as defined in input file
*       MESSAGE         - Error messages
*       InFormat        - Format for reading a record from input file
******************************************************************************
      CHARACTER CODEF*4,MESSAGE*255,InFormat*72

      YF=0
      MF=0
      DF=0
      IFAIL=0
*
*       Read a line of the input file
*
      WRITE(InFormat,'(''(I4,I2,I2,A4,'',I12,''F6.2)'')') NVARS
 19   READ(FILNO,InFormat,ERR=91,END=25,IOSTAT=IOERR)
     +                           YF,MF,DF,CODEF,(FromFile(I),I=1,NVARS)
*
*       If we've got the right day, update CURRENT according to the site code. 
*       We will have the right day EITHER if it's the first time of
*       reading the file OR if YF,MF and DF correspond to YY,MM and DD
*
      IF ( ((10000*(YF-YY))+(100*(MF-MM))+DF-DD.EQ.0)
     +                                  .OR.(YY.EQ.0) ) THEN
       Do I=1,NSITES
        IF (CODEF.EQ.SCODES(I)) THEN
         Do J=1,NVARS
          CURRENT(I,J) = FromFile(J)
          IF (DABS(FromFile(J)-MissVal).LT.1.0D-6) 
     +                                CURRENT(I,J)=-1.0D101
         End Do
        ENDIF
       End Do
*
*       If this is the first time through, set YY,MM and DD correctly
*
       IF (YY.EQ.0) THEN
        YY = YF
        MM = MF
        DD = DF
       ENDIF
*
*       If we overshoot, rewind the input file 1 line and go back
*
      ELSEIF ((10000*(YF-YY))+(100*(MF-MM))+DF-DD.GT.0) THEN
       BACKSPACE (UNIT=FILNO)
       RETURN
*
*       If the input file is the wrong order, stop with an error 
*       message
*
      ELSEIF ((10000*(YF-YY))+(100*(MF-MM))+DF-DD.LT.0) THEN
       GOTO 92
      ENDIF

      GOTO 19
*
*       End of file
*       
 25   IF (IOERR.LT.0) THEN
       DONE = 1
      ENDIF

      RETURN
*
*       Error trapping
*
 91   WRITE(MESSAGE,1) DF,MF,YF
      CALL INTPR(TRIM(MESSAGE),-1,0,0)
      IFAIL = 5
      RETURN
 92   WRITE(MESSAGE,2) DF,MF,YF,CODEF
      CALL INTPR(TRIM(MESSAGE),-1,0,0)
      IFAIL = 5
      RETURN

 1    FORMAT('Input file error: last successful access was for ',
     +       'data on ',I2,'/',I2,'/',I4,'.')
 2    FORMAT('Input file out of order for ',I2,'/',I2,'/',I4,
     +       ', at site ',A4,'.')

      END
******************************************************************************
******************************************************************************
******************************************************************************
      Subroutine DistCalc(SiteInfo,NSITES,NATTR,MXP,Distance)
******************************************************************************
*
*       To calculate an array of inter-site distances. Arguments:
*
*       SiteInfo        Array containing attribute information for each site.
*                       Distances will be calculated from the first two 
*                       attributes (if present).
*       NSITES          Number of sites
*       NATTR			Number of attributes defined for each site
*       MXP             For dimensioning Siteinfo
*       Distance        Array of inter-site distances. First slice is 
*                       X-separation, second slice is Y-separation and 
*                       third slice is Euclidean separation.
******************************************************************************
      Integer, intent(in) :: NSITES, NATTR, MXP
      Double precision, intent(in) :: Siteinfo(NSITES,MXP,0:3)
      Double precision, intent(out) :: Distance(3,NSITES,NSITES)
******************************************************************************
*       Additional Integers
*       -------------------
*       I,J             Counters
******************************************************************************
      Integer I,J
      
      If (NATTR.EQ.0) then
       Distance = 0.0d0
      Else if (NATTR.EQ.1) then
       Do I=1,NSITES
        Distance(1,I,I) = 0.0d0
        Do J=I+1,NSITES
         Distance(1,I,J) = Siteinfo(I,1,0)-Siteinfo(J,1,0)
         Distance(1,J,I) = -Distance(1,I,J)
        End Do
       End Do
       Distance(2,1:NSITES,1:NSITES) = 0.0d0
       Distance(3,1:NSITES,1:NSITES) = Distance(1,1:NSITES,1:NSITES)
      Else
       Do I=1,NSITES
        Distance(1:3,I,I) = 0.0d0
        Do J=I+1,NSITES
         Distance(1,I,J) = Siteinfo(I,1,0)-Siteinfo(J,1,0)
         Distance(2,I,J) = Siteinfo(I,2,0)-Siteinfo(J,2,0)
         Distance(3,I,J) = Dsqrt( (Distance(1,I,J)**2) + 
     +                            (Distance(2,I,J)**2) )
         Distance(1,J,I) = -Distance(1,I,J)
         Distance(2,J,I) = -Distance(2,I,J)
         Distance(3,J,I) = Distance(3,I,J)
        End Do
       End Do
      End if
         
      End subroutine DistCalc
