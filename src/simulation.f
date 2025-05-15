      SUBROUTINE SIMULATE(MODEL,IsPrecip,NSITES,NATTR,NVARS,
     +                    CurVar,MXP1,MXP2,SITINF1,SITINF2,
     +                    TmpFlNo,NP1,NP2,
     +                    COVCODE1,COVCODE2,TWO1,TWO2,THREE1,THREE2,
     +                    FOUID1,FOUID2,LEGID1,LEGID2,PWTID1,PWTID2,
     +                    SPMOD1,SPMOD2,SITXFM1,SITXFM2,GLBCOD,MaxLag,
     +                    BETA1,BETA2,THETA1,THETA2,PHI1,PHI2,RHO1,
     +                    RHO2,GLBVAL,InitVals,DateLims,UnitNos,
     +                    ReadFrom,IFAIL)
******************************************************************************
*       To simulate from a GLM fitted to daily climate time series.
*
*      INPUT ARGUMENTS
*      ^^^^^^^^^^^^^^^
*       MODEL   - Code(s) for model to be simulated. NB some models e.g.
*               - precipitation models and joint mean-variance models have 
*                 two components. These are the codes as defined in R
*                 routine model.num(). 
*       IsPrecip- 1 if this is a precipitation-like variable with an 
*                 "occurrence" and "amounts" component, 0 otherwise
*       NSITES  - Number of sites in study region
*       NATTR   - Number of attributes defined for each site
*       NVARS   - Number of variables in input data file
*       CurVar  - The index of the variable in input data file that we 
*                 are currently simulating.
*       MXP1,   - Maximum numbers of parameters in the models.
*       MXP2
*       MaxPars - The larger of MXP1 and MXP2. Used for dimensioning SITINF.
*       SITINF1,- Double precision arrays containing information 
*       SITINF2   regarding each site (Eastings, Northings etc.). 
*                 Column 0 holds the actual value, 1 and 2 hold 
*                 derivatives wrt any parameters of nonlinear
*                 transformations. Strictly speaking it isn't 
*                 necessary to have one of these for each model, but
*                 it makes the dimensioning easier. 
*       TmpFlNo - Number of temporary file containing 4-character site codes
*       NP1,NP2 - Integer arrays indicating numbers of parameters 
*                 in each model component. 
*                 See documentation for NP in fitting.f header for details.
*       BETA1,  - Parameter vectors in each model component (double 
*       BETA2     precision; element 0 is the constant term in each case).
*       COVCODE1 - Code numbers of covariates in the two model components
*       COVCODE2 
*       GLBCOD   - Definitions of `global' quantities
*       MaxLag  - Only keep track of last MaxLag days
*	FOUID1 }- Index arrays for Fourier representation of site effects. 
*	FOUID2 }  See file funcdefs.f for full details
*	LEGID1 }- Index arrays for Legendre polynomial representation of
*	LEGID2 }  site effects.
*       PWTID1 }- Indexing for weighted averages of previous days' values
*       PWTID2 }  ((I,0)th entry is the number of the covariate to which the
*                 parametrisation for the Ith weighting scheme is attached,
*                 (I,J)th is the model number of the Jth parameter in the
*                 scheme).
*	SPMOD1	 - Spatial correlation models for each component (SPMOD2
*	SPMOD2	   used only for precipitation models)
*	SITXFM1	} Select nonlinear transformations of site attributes in
*	SITXFM2 } each model. First column chooses the attribute, second 
*		  the transformation
*       TWO1    - indices of 2-way interactions in logistic regression model
*       TWO2    - ditto, intensity model
*       THREE1  - guess ...
*       THREE2  -    "
*	THETA1	- parameters in nonlinear transformations, for each model
*	THETA2	  component.
*	RHO1	- parameters defining spatial correlation matrix in each 
*	RHO2	  model component
*       PHI1,   - Dispersion parameters (double precision)
*       PHI2
*       GLBVAL  - Double precision array containing values of `global'
*                 quantities such as trace thresholds
*       InitVals - Values to use for initialising each variable at the
*                  start of simulation in case values are missing. 
*       DateLims - Two-element vector of integers containing the start
*                 and end months of simulation, in the form YYYYMM
*       UnitNos - Numbers of units attached to input and output files
*                 (which should be opened prior to calling this routine).
*                 For a description of the unit numbers, see the header
*                 for the CheckFiles() routine in the R code.
*       ReadFrom- Number of unit from which to read data (this will be
*                 95 for the first variable of any simulation, 96 for 
*                 the remaining variables)
*       IFAIL   - error flag
***************************************************************************
      Integer, intent(in) :: MODEL(2),IsPrecip,Nsites,NATTR,NVARS
      Integer, intent(in) :: CurVar,MXP1,MXP2,NP1(10),NP2(10),TmpFlNo
      Integer, intent(in) :: COVCODE1(MXP1),COVCODE2(MXP2)
      Integer, intent(in) :: TWO1(MXP1,2),TWO2(MXP2,2)
      Integer, intent(in) :: THREE1(MXP1,3),THREE2(MXP2,3)
      Integer, intent(in) :: FOUID1(MXP1),FOUID2(MXP2)
      Integer, intent(in) :: LEGID1(MXP1),LEGID2(MXP2)
      Integer, intent(in) :: SITXFM1(MXP1,2),SITXFM2(MXP2,2)
      Integer, intent(in) :: SPMOD1,SPMOD2
      Integer, intent(in) :: PWTID1(MXP1,0:3,NVARS),
     +                       PWTID2(MXP2,0:3,NVARS)
      Integer, intent(in) :: MaxLag,GLBCOD(MXP1),DateLims(2)
      Integer, intent(in) :: UnitNos(100),ReadFrom
      Integer, intent(out) :: IFAIL
      Double precision, intent(in) :: GLBVAL(MXP1),InitVals(NVARS)
      Double precision, intent(in) :: BETA1(0:MXP1),BETA2(0:MXP2)
      Double precision, intent(in) :: THETA1(MXP1,3),THETA2(MXP2,3)
      Double precision, intent(in) :: PHI1,PHI2
      Double precision, intent(in) :: RHO1(MXP1),RHO2(MXP2)
      Double precision, intent (inout) :: SITINF1(NSITES,MXP1,0:3),
     +                                    SITINF2(NSITES,MXP2,0:3)
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>    PARAMETERS      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*       MXPTXT          Maximum number of parameters for CHARACTER 
*                       ARRAYS IN THE DESCRIBE COMMON BLOCK. THIS CAN BE 
*                       REMOVED AFTER SHIFTING ALL CHARACTER MANIPULATIONS 
*                       TO R
*
*	NB IF THIS VALUE IS CHANGED, THE CORRESPONDING ENTRY IN
*       OTHER SOURCE FILES WILL ALSO NEED CHANGING. 
******************************************************************************
        INTEGER MXPTXT
        PARAMETER (MXPTXT=80)
******************************************************************************

*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>      VARIABLES     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
***************************************************************************
*       INTEGER variables
*       ^^^^^^^^^^^^^^^^^
*       FLAGS   - array indicating required status of each entry of current
*                 data array:
*                 0 = Value known or not needed, don't change it
*                 1 = Value unknown, simulate from the model
*       OFlags  - Value of FLAGS yesterday
*       OWDFlags- Value of FLAGS yesterday at the stage of allocating wet/dry
*                 sites
*       MXLArr  - array containing values of MaxLag for calls to COVSET
*                 (for simulation purposes these are all equal to -1)
*	MISSFL	- Flags indicating whether or not covariates are missing 
*                 at each site
*       WetDry  - array holding 0/1 indicators for precipitation 
*                 occurrence
*	OWetDry - value of WetDry yesterday
*       ChangeFlags - Indicates whether (1) or not (0) any element of
*                 the FLAGS array has changed since yesterday for a 
*                 particular day. Required as an argument to COVSET.
*       Nknown  - Number of sites whose values are known
*       NKWet   - Precip only: number of sites that are known to be wet
*       ModToPass - Code of model to pass to DailySim routine
*	FORCE	- flag to indicate whether to force recalculation of all 
*	   	  covariates (1) or just those which have changed since last
*		  calculation (0)
*       YY      - current year
*       YRLIMS  - first and last years to do
*       MM      - current month
*       MOLIMS  - first and last months to do
*       DD      - current day
*       DYLIMS  - first and last days
*       DAYS_IN_MONTH   - no. of days in current month (integer function)
*       SITE    - Current site
*       DATE    - current month, in YYYYMM format
*       JSTART  - Julian day number corresponding to start date of 
*                 simulation
*       JFINIT  - Julian day number for first record in scratch file
*       JCUR    - Julian day number for current date
*       RECNO   - Number of record from which to read data
*       YF,MF,DF- Year, month and day read from scratch file
*       SEED    - seed for random number generator
*       WTSCH1 }- Indexing for weighted average of previous days' values
*       WTSCH2 }  (Ith element is the weighting scheme used in Ith
*                 covariate).
*	PUPDT	- indicates which components of the linear predictors are
*		  to be updated. If PUPDT = 1, all components are updated;
*		  if 2, year, seasonal and daily effects, if 3, just
*		  seasonal and daily and if 4 just daily
*       THRTYP  - Type of thresholding to be carried out on small
*                 positive values. 0: no thresholding; 1: `soft'
*                 thresholding (zero everything below threshold, 
*                 subtract threshold from everything above); 2:
*                 `hard' thresholding (just zero everything below
*                 threshold)
*       ICHECK  - Indicates whether model parametrisations have been 
*                 checked
*       WDRecalc- Indicator for whether we need to recalculate covariance
*                 matrix inverses etc. in WDALLOC.
***************************************************************************
      INTEGER, dimension (NSITES) :: FLAGS,OFlags,OWDFlags,MISSFL,
     +                                               WetDry,OWetDry
      INTEGER MXLarr(NVARS),ChangeFlags,Nknown,NKWet,ModToPass
      INTEGER YY,YRLIMS(2),MM,MOLIMS(2),DD,DYLIMS(2),DAYS_IN_MONTH
      INTEGER SITE,DATE,JSTART,JFINIT,JCUR,RECNO,YF,MF,DF
      INTEGER WTSCH1(MXP1),WTSCH2(MXP2)
      INTEGER FORCE,PUPDT,THRTYP,ICHECK,WDRecalc
      INTEGER I,J,K
***************************************************************************
*       DOUBLE PRECISION variables
*       ^^^^^^^^^^^^^^^^^^^^^^^^^^
*       DatArray- array holding all current and previous values
*       ETA1, ETA2      - linear predictors for the two GLMs
*       MU1     - mean for the first GLM
*       ANSMEAN - vector of mean Anscombe residuals for gamma model
*       MeanToPass      - Vector of means to pass to DailySim
*       SDToPass        - Vector of standard deviations, similarly
*       WDCSD   - Vector of conditional standard deviations for use in 
*                 latent Gaussian correlation structures
*       COVS1   - Values of covariates in each model component
*       COVS2   
*	CORMAT1	- Inter-site latent correlations for logistic model (for
*                 passing to WDALLOC)
*	CorToPass - Inter-site residual correlations, all other models
*                 (for passing to DailySim)
*       SIGMA11 } Partitions of residual covariance matrix in DailySim 
*       SIGMA22 } routine (passed back here so that their values 
*       SIGMA12 }  can be saved between calls)
*       WDSIG11 } Ditto for residual covariance matrix of latent Gaussian
*       WDSIG22 } variables in WDALLOC
*       WDSIG12 } 
*       CVAR    - Covariance matrix of latent Gaussian variables conditioned
*                 on observed sites, in WDALLOC
*       SIGMA   - Conditional covariance matrix for DailySim
*       MOTOT   - Monthly totals at each site
*       YRTOT   - Annual totals at each site
*       PRED$   - contributions to linear predictors, potentially for both
*                 model components ($=1,2). The elements indicate which
*                  set of explanatory variables is represented:
*                       1 - all the bits we know once the site is fixed
*                       2 - the additional bits once site and year are fixed
*                       3 - site, year and month
*                       4 - site, year, month and previous day
*       Distance- array of inter-site distances. First slice is X-separation,
*                 second slice is Y-separation and third slice is Euclidean
*                 separation
******************************************************************************
      Double precision, dimension(NSITES) :: ETA1,ETA2,MU1,ANSMEAN
      DOUBLE PRECISION, dimension(NSITES) :: MeanToPass,SDToPass,WDCSD
      DOUBLE PRECISION CORMAT1(NSITES,NSITES)
      DOUBLE PRECISION CorToPass(NSITES,NSITES)
      DOUBLE PRECISION DatArray(NSITES,NVARS,0:10)
      DOUBLE PRECISION MOTOT(NSITES,12),YRTOT(NSITES)
      DOUBLE PRECISION PRED1(NSITES,4),PRED2(NSITES,4)
      DOUBLE PRECISION COVS1(NSITES,MXP1),COVS2(NSITES,MXP2)
      Double precision, dimension(NSITES,NSITES) ::
     +    SIGMA11,SIGMA22,SIGMA12,SIGMA,WDSIG11,WDSIG12,WDSIG22,WDCVAR
      Double precision Distance(3,NSITES,NSITES)
      DOUBLE PRECISION THRESH,TRACE
***************************************************************************
*       CHARACTER variables
*       ^^^^^^^^^^^^^^^^^^^
*       SCODES  - array of short site codes
*       ATTRTXT - Text for site attributes
*       MOTXT   - Text for seasonal model components
*       TRNDTXT - Text for trend components
*       DYTXT   - Text for daily components
*       PWTTXT  - Text for parameters in nonlinear weighting schemes for
*                 averaging previous days' values
*	TRPTXT	- Text for nonlinear trend transformations
*	SXPTXT	- Text for nonlinear site transformations
*	SPATXT	- Text for spatial correlation models
*	SPPTXT	- Text for parameters in spatial correlation 
*		  model (dimension to MXP1 for convenience)
*       CORTXT  - array describing quantities to which inter-site 
*                 correlations correspond
*       GLBTXT  - Text for global quantities
*       MESSAGE - For error messages
***************************************************************************
      CHARACTER SCODES(NSITES)*4
      CHARACTER*70 ATTRTXT1(MXP1),ATTRTXT2(MXP2)
      CHARACTER*70 PWTTXT(MXPTXT,3),MOTXT(MXPTXT)
      CHARACTER*70 TRNDTXT(MXPTXT),DYTXT(MXPTXT)
      CHARACTER*70 TRPTXT(MXPTXT,3),SXPTXT(MXPTXT,3)
      CHARACTER*70 SPATXT(0:MXPTXT),SPPTXT(MXPTXT,4)
      CHARACTER*70 CORTXT(MXPTXT),GLBTXT(MXPTXT)
      CHARACTER*70 MESSAGE
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>    LINE NUMBERS    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*       The following line numbering convention is used:
*       LABEL RANGE     DESCRIPTION
*       -----------     -----------
*       1 - 99          FORMAT statements
*       500 - 599       Daily simulation        
*       600-699         Multiple simulations
*       900-999         Error trapping
***************************************************************************
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>    COMMON BLOCKS    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*       DESCRIBE        Descriptive text
**************************************************************************
      COMMON /DESCRIBE/ SXPTXT,MOTXT,TRNDTXT,TRPTXT,DYTXT,
     +                  SPATXT,SPPTXT,CORTXT,PWTTXT,GLBTXT

*
*       Initializing (Fortran90 syntax). NB initialise OFlags to -1, 
*       so that on the first day ChangeFlags is guaranteed to be set
*       to 1 and everything is initialised when simulating to 
*       correspond to the first day of simulation. 
*
      OFlags = -1
      OWDFlags = -1
      OWetDry = -1
      COVS1 = 0.0D0
      COVS2 = 0.0D0
      WTSCH1 = 0
      WTSCH2 = 0
      ICHECK = 0
      MXLarr = -1      
      WDRecalc = 1
*
*       Read site codes from temporary file (they are not used directly anywhere, 
*       but they appear in some error messages)
*
      CALL GetScodes(TmpFlNo, NSites, Scodes, Ifail)
      IF (IFAIL.NE.0) RETURN
*
*       Next lines define the model code to pass to the main daily 
*       simulation routine. This is equal to MODEL(1) except in the
*       case when we're doing precipitation, in which case the main
*       routine simulates gamma variates as defined in MODEL(2). 
*
      ModToPass = MODEL(1)
      IF ( (ModToPass.EQ.1).AND.(IsPrecip.Eq.1) ) ModToPass = 11
      DO I=1,NATTR
       IF (I.LE.MXP1) THEN 
        WRITE(ATTRTXT1(I),'(''~[[site attribute#'',I4,''#]]~'')') I
       ENDIF
       IF (I.LE.MXP2) THEN
        WRITE(ATTRTXT2(I),'(''~[[site attribute#'',I4,''#]]~'')') I
       ENDIF
      END DO
*
*       Ensure that external files and pointers are positioned correctly
*
      Call FileReset(UnitNos(3:5))
*
*       Calculate any requested nonlinear transformations of site attributes,
*       along with inter-site distances
*
      CALL ATTRXFM(SITINF1,NSITES,ATTRTXT1,NATTR,COVCODE1,
     +                          SITXFM1,FOUID1,LEGID1,THETA1,NP1,
     +                          ICHECK,IFAIL,MXP1)
      IF (IFAIL.NE.0) RETURN
      IF (MODEL(2).GT.0) THEN 
       CALL ATTRXFM(SITINF2,NSITES,ATTRTXT2,NATTR,COVCODE2,
     +                           SITXFM2,FOUID2,LEGID2,THETA2,NP2,
     +                           ICHECK,IFAIL,MXP2)
       IF (IFAIL.NE.0) RETURN
      ENDIF
      
      Call DistCalc(SITINF1,NSITES,MXP1,Distance)      
*
*	Calculate spatial correlation / covariance matrices (& convert
*       to covariances for gamma & homoscedastic normal models). This is 
*       a bit messy: for logistic models we need to pick up a correlation
*       matrix to pass to WDALLOC (and possibly another to pass to 
*       DailySim for precipitation models), but for everything else we 
*       just need one for DailySim. Note also that the "independence" 
*       structure is treated as a special case, so a model still needs
*       to be defined (hence SPMOD1.GE.0 rather than SPMOD1.GT.0).
*
      IF ((SPMOD1.GE.0).AND.(SPMOD1.LT.20)) THEN
       IF (MODEL(1).EQ.1) THEN 
        CALL CORSET(UnitNos(6),SPMOD1,RHO1,CORMAT1,NSITES,SITINF1,
     +                              NATTR,SCODES,MXP1,IFAIL)
        IF (IFAIL.NE.0) RETURN
       ELSE
        CALL CORSET(UnitNos(6),SPMOD1,RHO1,CorToPass,NSITES,SITINF1,
     +                              NATTR,SCODES,MXP1,IFAIL)
        IF (IFAIL.NE.0) RETURN
        IF (MODEL(1).EQ.11) THEN
         CALL ANSET(ANSMEAN,CorToPass,PHI1,NSITES,IFAIL)
         IF (IFAIL.NE.0) RETURN
        ENDIF
       ENDIF
      ENDIF
*
*       This is the "intensity" correlation structure for precip-like
*       models
*      
      IF ((IsPrecip.EQ.1).AND.(SPMOD2.GE.0).AND.(SPMOD2.LT.20)) THEN
       CALL CORSET(UnitNos(7),SPMOD2,RHO2,CorToPass,NSITES,SITINF2,
     +                              NATTR,SCODES,MXP2,IFAIL)
       IF (IFAIL.NE.0) RETURN
       IF (MODEL(2).EQ.11) THEN
        CALL ANSET(ANSMEAN,CorToPass,PHI2,NSITES,IFAIL)
        IF (IFAIL.NE.0) RETURN
       ENDIF
      ENDIF
*
*       Now define start and end dates etc. 
*
      Do I=1,2
       YRLIMS(I) = Datelims(I) / 100
       MOLIMS(I) = MOD(Datelims(I),100)
       DYLIMS(I) = 1
      End Do
      DYLIMS(2) = DAYS_IN_MONTH(MOLIMS(2),YRLIMS(2))
*
*     Define thresholds for small positive values, if required
*
      CALL TAUDEF(GLBVAL,GLBCOD,NP1(8),THRESH,TRACE,THRTYP,
     +                                           MODEL(1),MXP1)
*
*       Initialise the data array and vector of simulated values.
*       Start by finding the Julian day number of the first record
*       in the scratch file holding the data; and use this to figure 
*       out the number of the first record that needs to be read. Then
*       read this record along with the required number of previous
*       days' values. Note that the initialisation of the data array
*       itself is for the day before the start of simulation (hence
*       loop from MaxLag-1 to 0 because we'll never need the values
*       corresponding to MaxLag); this ensures that when the date is 
*       incremented on the first day of simulation, the data all end 
*       up in the right place. Note also that the spin-up values are 
*       also written to the output scratch file, because this will be 
*       the input file for the next variable in a multivariate setting.
*
      READ(UnitNos(ReadFrom),REC=1) SITE,YY,MM,DD,
     +                                 (DatArray(1,J,0),J=1,NVARS)
      CALL DMYTOJ(DD,MM,YY,JFINIT,IFAIL)
      IF (IFAIL.NE.0) RETURN
      CALL DMYTOJ(DYLIMS(1),MOLIMS(1),YRLIMS(1),JSTART,IFAIL)
      Do I=MaxLag-1,0,-1
       Do J=1,NSITES
        RECNO = NSITES*(JSTART-JFINIT-I-1)+J
        READ(UnitNos(ReadFrom),REC=RECNO) 
     +                   SITE,YY,MM,DD,(DatArray(J,K,I),K=1,NVARS)
        Do K=1,NVARS
         If (DatArray(J,K,I).LT.-1.0D99) DatArray(J,K,I) = InitVals(K) 
        End Do
        WRITE(UnitNos(96),REC=RECNO) 
     +                   SITE,YY,MM,DD,(DatArray(J,K,I),K=1,NVARS)
       End Do
      End Do
*
*       Simulation/interpolation starts here. Loop over years, months 
*	and days of simulation. Don't need to compute site contributions
*	again, so set PUPDT = 1 here; and set the current Julian date
*       to the day before the start of simulation so that it increments
*       correctly.
*
      PUPDT = 1

      JCUR = JSTART-1
      DO 500 YY = YRLIMS(1),YRLIMS(2)
*
*       A new year - need to update trend predictors (unless we're
*	already calculating site predictors)
*
       IF (PUPDT.GE.2) PUPDT = 2
*
*       Initialise the current yearly and monthly totals to zero.
*
       YRTOT = 0.0D0
       MOTOT = 0.0D0
       DO 520 MM = 1,12
*
*       Skip if we're outside requested range
*
        DATE = (100*YY) + MM
        IF ( (DATE.LT.DATELIMS(1)).OR.(DATE.GT.DATELIMS(2)) ) GOTO 520
*
*       New month - need to update seasonal predictors (unless we're already 
*	doing trend predictors). Also FORCE recalculation of all
*       covariates
*
        IF (PUPDT.GE.3) PUPDT = 3
        FORCE = 1
        DO 540 DD = 1,DAYS_IN_MONTH(MM,YY)
         Nknown = 0
*
*       First thing is to update the data values (including 
*       thresholding if requested), FLAGS etc. 
*
         JCUR = JCUR + 1
         DatArray(1:NSITES,1:NVARS,1:MaxLag) = 
     +                       DatArray(1:NSITES,1:NVARS,0:(MaxLag-1))
         Do J=1,NSITES
          RECNO = NSITES*(JCUR-JFINIT)+J
          READ(UnitNos(ReadFrom),REC=RECNO) 
     +                     SITE,YF,MF,DF,(DatArray(J,K,0),K=1,NVARS)
          If (DatArray(J,CurVar,0).GE.-1.0D100) then
           Nknown = Nknown + 1
           FLAGS(J) = 0
           If (THRTYP.GT.0) then
            If (DatArray(J,CurVar,0).LT.0.0D0) goto 20
            If (DatArray(J,CurVar,0).LT.THRESH) then
             DatArray(J,CurVar,0) = 0.0D0
            Elseif (THRTYP.EQ.1) then
             DatArray(J,CurVar,0) = DatArray(J,CurVar,0) - THRESH
            Endif
           Endif
          Else
           FLAGS(J) = 1
          Endif
 20      End Do
*
*     Calculate linear predictors, and set FORCE to zero (initially it's 1).
*     Also set PUPDT to 4 so that next time through we'll just update
*     day-to-day quantities. 
*
         CALL COVSET(COVS1,UnitNos(3:5),NP1,NSITES,NVARS,CurVar,
     +               COVCODE1,TWO1,THREE1,SITINF1,Distance,YY,MM,DD,
     +               DatArray,1,PWTID1,MXLarr,BETA1,THETA1,TRACE,0,
     +               FORCE,MISSFL,ICHECK,IFAIL,MXP1,1,1)

         IF (IFAIL.NE.0) RETURN
         CALL SETPRED(PRED1,COVS1,NP1,NSITES,TWO1,THREE1,BETA1,
     +                                       PUPDT,NSITES,MXP1)
         IF (MODEL(2).GT.0) THEN
          CALL COVSET(COVS2,UnitNos(3:5),NP2,NSITES,NVARS,CurVar,
     +                COVCODE2,TWO2,THREE2,SITINF2,Distance,YY,MM,DD,
     +                DatArray,1,PWTID2,MXLarr,BETA2,THETA2,TRACE,0,
     +                FORCE,MISSFL,ICHECK,IFAIL,MXP2,1,1)
          IF (IFAIL.NE.0) RETURN
          CALL SETPRED(PRED2,COVS2,NP2,NSITES,TWO2,THREE2,BETA2,
     +                                        PUPDT,NSITES,MXP2)
         ENDIF

         FORCE = 0
         PUPDT = 4
*
*     If we got through that lot, model parametrisation must be OK -
*     no further need to check
*
         ICHECK = 1
         If (Nknown.eq.NSITES) goto 554
*
*       Calculate linear predictors, means and standard deviations to 
*       pass to WDALLOC and DailySim later
*
         ETA1(1:NSITES) = PRED1(1:NSITES,4)
         ETA2(1:NSITES) = PRED2(1:NSITES,4)

         DO 550 I=1,NSITES
*
*       If MODEL(1) is 1 we're doing a logistic regression; this 
*       could be stand-alone or coupled with precip. 
*
          IF (MODEL(1).EQ.1) THEN
           MU1(I) = 1.0D0 / (1.0D0 + DEXP(-ETA1(I)))
           IF (MODEL(2).EQ.11) THEN
            MeanToPass(I) = DEXP(ETA2(I))
           ENDIF
          ELSE IF ((MODEL(1).EQ.10).OR.(MODEL(1).EQ.20)) THEN
           MeanToPass(I) = ETA1(I)
           IF (MODEL(1).EQ.10) THEN
            SDToPass(I) = DSQRT(PHI1)
           ELSE
*
*       Joint mean variance-model: calculate standard deviations
*       by halving the linear predictor before applying exponent
*       (otherwise have to do SQRT as well as EXP).
*
            IF (MODEL(2).EQ.21) THEN
             SDToPass(I) = DEXP(ETA2(I)/2.0D0)
            ELSE
             IFAIL = 999
             RETURN
            ENDIF
           ENDIF
*
*       Standard gamma GLM
*
          ELSE IF (MODEL(1).EQ.11) THEN
           MeanToPass(I) = DEXP(ETA1(I))
          ELSE
           IFAIL = 999
           RETURN
          ENDIF
 550     CONTINUE
*
*       Generate today's values at all required sites. The logic is:
*
*       1. If (precip or a logistic model) {
*           Generate wet/dry pattern using imputation if necessary
*           Amend FLAGS so that it is 0 only for "observed wet" sites;
*              then the zeroes will correspond to the values we're 
*              conditioning on when generating intensities
*          } else {
*           Generate a "dummy" wet / dry pattern that is 1 everywhere,
*           to ensure that other variables are all generated at every
*           site (this is for coding convenience, because it allows
*           us to use the same basic code structure for precipitation
*           as for other variables)
*          }
*       2. If any elements of FLAGS have changed since the last time
*              around, recalculate required matrix inverses and 
*              Cholesky decompositions
*       3. Sample a residual vector from the appropriate multivariate
*              normal distribution, and back-transform to the required 
*              vector of values
*
         IF ((IsPrecip.EQ.1).OR.(ModToPass.EQ.1)) THEN
          WetDry(1:NSITES) = -1
          NKWet = 0
          WDRecalc = 0
          Do I=1,NSITES
           IF (FLAGS(I).EQ.0) THEN
            IF (DABS(DatArray(I,CurVar,0)).LT.1.0D-6) THEN
             WetDry(I) = 0
            ELSE
             WetDry(I) = 1
             NKWet = NKWet + 1
            ENDIF
           ENDIF
           If (WDRecalc.Eq.0) then
            If (FLAGS(I).NE.OWDFlags(I)) WDRecalc = 1
           Endif
          End Do
          OWDFlags = FLAGS
          CALL WDALLOC(WetDry,NSITES,Nknown,NKWet,FLAGS,SPMOD1,
     +                 Eta1,RHO1,CORMAT1,WDSIG11,WDSIG12,WDSIG22,
     +                 WDCVAR,WDCSD,WDRecalc,IFAIL)
          IF (IFAIL.NE.0) RETURN
         ELSE
          WetDry(1:NSITES) = 1
         ENDIF
*
*       Next line: result will be 0 if either FLAGS or WETDRY is zero 
*       (i.e. if the value is known or if it has just been set to a 'dry' 
*       zero), otherwise result will be 1. Result: FLAGS contains ones only
*       in those elements that need to be generated in the DailySim routine.
*	For non-precipitation variables, covariance matrices etc. will need
*       to be recalculated only if FLAGS has changed since yesterday. For
*	the intensity component of precipitation variables however, these
*	quantities will also need to be recalculated if the wet-dry 
*	pattern has changed since yesterday.
*
         FLAGS(1:NSITES) = FLAGS(1:NSITES) * WetDry(1:NSITES)
         ChangeFlags = 0
         Do I=1,NSITES
          If ((Flags(I).NE.OFlags(I)).OR.
     +        (WetDry(I).NE.OWetDry(I))) THEN
           ChangeFlags = 1
           Goto 551
          End If
         End Do

 551     CALL DailySim(ModToPass,DatArray(1:NSITES,CurVar,0),WetDry,
     +                FLAGS,NSITES,NKWET,MeanToPass,
     +                SDToPass,CorToPass,ANSMEAN,Sigma11,Sigma12,
     +                Sigma22,Sigma,ChangeFlags,IFAIL)
         IF (IFAIL.NE.0) THEN
          IF (IFAIL.EQ.400) THEN
           WRITE(MESSAGE,'(A14,I8)') 'Error on date ',100*DATE+DD
           CALL INTPR("",-1,0,0)
           CALL INTPR("",-1,0,0)
           CALL INTPR(TRIM(MESSAGE),-1,0,0)
          ENDIF
          RETURN
         ENDIF

         OFlags(1:NSITES) = FLAGS(1:NSITES)
         OWetDry(1:NSITES) = WetDry(1:NSITES)

*
*     If the models were fitted to *soft thresholded* data, we need to add
*     the threshold back in to non-zero amounts (unless we're just 
*     simulating occurrences, for some reason). Hard-thresholded data
*     are left as simulated. And write the data back to a second 
*     scratch file
*
 554     DO 555 I = 1,NSITES
          IF ((ModToPass.NE.1).AND.(THRTYP.EQ.1)) THEN
           IF (DatArray(I,CurVar,0).GT.0.0D0) 
     +           DatArray(I,CurVar,0) = DatArray(I,CurVar,0) + THRESH
          ENDIF
          RECNO = NSITES*(JCUR-JFINIT)+I
          WRITE(UnitNos(96),REC=RECNO) 
     +                    I,YY,MM,DD,(DatArray(I,J,0),J=1,NVARS)
 555     CONTINUE
*
*     Now: if models were fitted to *soft thresholded* data, we need to
*     re-apply the threshold so that, when today's values become
*     covariates for tomorrow's distributions, they mean the same as
*     they did when the models were fitted. At this point there won't
*     be any positive values less than the threshold (see line 555 above)
*     so we don't need to worry about them.
*
         IF ((ModToPass.NE.1).AND.(THRTYP.EQ.1)) THEN
          DO 565 I = 1,NSITES
           IF (DatArray(I,CurVar,0).GT.0.0D0) 
     +         DatArray(I,CurVar,0) = DatArray(I,CurVar,0) - THRESH
 565      CONTINUE
         ENDIF

 540    CONTINUE
 520   CONTINUE

 500  CONTINUE
 
      RETURN

      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE ANSET(ANSMEAN,ANSVAR,PHI,NSITES,IFAIL)
*
*       Takes correlation matrix of Anscombe residuals and converts
*        to covariance matrix, based on dispersion parameter PHI of the
*        gamma distributions. ANSMEAN is vector of mean Anscombe 
*        residuals (output), ANSVAR the correlation (input) and
*        covariance (output) matrix, PHI the dispersion parameter (input),
*        and NSITES the number of sites.
******************************************************************************
      INTEGER NSITES,IFAIL
      DOUBLE PRECISION ANSMEAN(NSITES),ANSVAR(NSITES,NSITES),PHI
******************************************************************************
*       Additional DOUBLES
*       ^^^^^^^^^^^^^^^^^^
*       MEAN    - mean Anscombe residual
*       VAR     - Variance of Anscombe residuals
*       NU      - Shape parameter of gamma dbn (=1/PHI)
******************************************************************************
      DOUBLE PRECISION MEAN,VAR,NU
******************************************************************************
*       Additional INTEGERs
*       ^^^^^^^^^^^^^^^^^^^
*       I,J        - counters
******************************************************************************
      INTEGER I,J
*
*       First find mean and variance of Anscombe residuals
*
      NU = 1.0D0 / PHI
      
      CALL ANSPRM(NU,MEAN,VAR,IFAIL)
*
*        Now calculate the covariance matrix
*
      DO 100 I=1,NSITES
       ANSMEAN(I) = MEAN
       ANSVAR(I,I) = VAR
       DO 105 J=I+1,NSITES
        ANSVAR(I,J) = ANSVAR(I,J) * VAR
        ANSVAR(J,I) = ANSVAR(I,J)
 105   CONTINUE
 100  CONTINUE

      RETURN
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE SETPRED(PRED,COVS,NP,NSITES,TWOWAY,THREEWAY,
     +                                       BETA,PUPDT,MXN,MXP)
*******************************************************************************
*       Calculates linear predictor in a model - updates only the
*        bits requested by the status of PUPDT
*******************************************************************************
*       INTEGER variables
*       ^^^^^^^^^^^^^^^^^
*       NP      - No. of parameters in the model. Input
*       NSITES  - No. of sites. Input
*       TWOWAY  - list of 2-way interactions. Input
*       THREEWAY- list of 3-way interactions. Input
*        PUPDT	- Indicates which predictors to update. Input/output
*       MONTH   - the current month. Input
*       MXN     - Maximum no. of sites. Input
*       MXP     - Maximum no. of parameters. Input
*******************************************************************************
      INTEGER NP(10),NSITES,TWOWAY(MXP,2),THREEWAY(MXP,3)
      INTEGER MXN,MXP
      INTEGER PUPDT
      INTEGER I,J
*******************************************************************************
*       DOUBLE PRECISION variables
*       ^^^^^^^^^^^^^^^^^^^^^^^^^^
*       PRED    - site contribution to linear predictor. Output
*       COVS    - Values of covariates at each site
*       BETA    - Vector of parameters
*******************************************************************************
      DOUBLE PRECISION PRED(MXN,4),COVS(MXN,MXP),BETA(0:MXP)

*
*       Site effects ...
*
      IF (PUPDT.LE.1) THEN
       DO 100 I=1,NSITES
        PRED(I,1) = BETA(0)
        DO 105 J=1,NP(1)
          PRED(I,1) = PRED(I,1) + (BETA(J)*COVS(I,J))
 105    CONTINUE
 100   CONTINUE
*
*       2-way interactions - only those between site effects
*
       DO 120 J=NP(4)+1,NP(5)
        IF ( (TWOWAY(J,1).LE.NP(1)).AND.(TWOWAY(J,2).LE.NP(1)) ) THEN
         DO 125 I=1,NSITES
           PRED(I,1) = PRED(I,1) + (BETA(J)*COVS(I,J))
 125     CONTINUE
        ENDIF
 120   CONTINUE  
*
*       And 3-way interactions
*
       DO 130 J=NP(5)+1,NP(6)
        IF ( (THREEWAY(J,1).LE.NP(1)).AND.
     +       (THREEWAY(J,2).LE.NP(1)).AND.
     +       (THREEWAY(J,3).LE.NP(1))) THEN
         DO 135 I=1,NSITES
           PRED(I,1) = PRED(I,1) + (BETA(J)*COVS(I,J))
 135     CONTINUE
        ENDIF
 130   CONTINUE  
      ENDIF
*
*       Year effects ...
*
      IF (PUPDT.LE.2) THEN
       DO 250 I=1,NSITES
        PRED(I,2) = PRED(I,1)
 250   CONTINUE
       DO 300 J=NP(1)+1,NP(2)
        DO 305 I=1,NSITES
         PRED(I,2) = PRED(I,2) + (BETA(J)*COVS(I,J))
 305    CONTINUE
 300   CONTINUE
*
*       2-way interactions (year-year and year-site) ...
*
       DO 320 J=NP(4)+1,NP(5)
        IF ( (TWOWAY(J,1).LE.NP(2)).AND.(TWOWAY(J,2).LE.NP(2)) ) THEN
         IF ( (TWOWAY(J,1).LE.NP(1)).AND.(TWOWAY(J,2).LE.NP(1)) )
     +                                                   GOTO 320
         DO 325 I=1,NSITES
          PRED(I,2) = PRED(I,2) + (BETA(J)*COVS(I,J))
 325     CONTINUE
        ENDIF
 320   CONTINUE  
*
*       And 3-way interactions
*
       DO 330 J=NP(5)+1,NP(6)
        IF ( (THREEWAY(J,1).LE.NP(2)).AND.
     +       (THREEWAY(J,2).LE.NP(2)).AND.
     +       (THREEWAY(J,3).LE.NP(2))) THEN
         IF ( (THREEWAY(J,1).LE.NP(1)).AND.
     +        (THREEWAY(J,2).LE.NP(1)).AND.
     +        (THREEWAY(J,3).LE.NP(1))) GOTO 330
         DO 335 I=1,NSITES
           PRED(I,2) = PRED(I,2) + (BETA(J)*COVS(I,J))
 335     CONTINUE
        ENDIF
 330   CONTINUE  
      ENDIF
*
*        Seasonal effects ...
*
      IF (PUPDT.LE.3) THEN
       DO 450 I=1,NSITES
        PRED(I,3) = PRED(I,2)
 450   CONTINUE
       DO 500 J=NP(2)+1,NP(3)
        DO 505 I=1,NSITES
          PRED(I,3) = PRED(I,3) + (BETA(J)*COVS(I,J))
 505    CONTINUE
 500   CONTINUE
*
*       2-way interactions (month-month, month-year and month-site)
*
       DO 520 J=NP(4)+1,NP(5)
        IF ( (TWOWAY(J,1).LE.NP(3)).AND.(TWOWAY(J,2).LE.NP(3)) ) THEN
         IF ( (TWOWAY(J,1).LE.NP(2)).AND.(TWOWAY(J,2).LE.NP(2)) )
     +                                                   GOTO 520
         DO 525 I=1,NSITES
           PRED(I,3) = PRED(I,3) + (BETA(J)*COVS(I,J))
 525     CONTINUE
        ENDIF
 520   CONTINUE  
*
*       And 3-way interactions
*
       DO 530 J=NP(5)+1,NP(6)
        IF ( (THREEWAY(J,1).LE.NP(3)).AND.
     +       (THREEWAY(J,2).LE.NP(3)).AND.
     +       (THREEWAY(J,3).LE.NP(3))) THEN
         IF ( (THREEWAY(J,1).LE.NP(2)).AND.
     +        (THREEWAY(J,2).LE.NP(2)).AND.
     +        (THREEWAY(J,3).LE.NP(2))) GOTO 530
         DO 535 I=1,NSITES
           PRED(I,3) = PRED(I,3) + (BETA(J)*COVS(I,J))
 535     CONTINUE
        ENDIF
 530   CONTINUE  
      ENDIF
*
*        Previous day - always done so don't need to condition on PUPDT.
*
      DO 650 I=1,NSITES
       PRED(I,4) = PRED(I,3)
 650  CONTINUE
      DO 700 J=NP(3)+1,NP(4)
       DO 705 I=1,NSITES
        PRED(I,4) = PRED(I,4) + (BETA(J)*COVS(I,J))
 705   CONTINUE
 700  CONTINUE
*
*       2-way interactions ...
*
      DO 720 J=NP(4)+1,NP(5)
       IF ( (TWOWAY(J,1).LE.NP(4)).AND.(TWOWAY(J,2).LE.NP(4)) ) THEN
        IF ( (TWOWAY(J,1).LE.NP(3)).AND.(TWOWAY(J,2).LE.NP(3)) )
     +                                                   GOTO 720
        DO 725 I=1,NSITES
         PRED(I,4) = PRED(I,4) + (BETA(J)*COVS(I,J))
 725    CONTINUE
       ENDIF
 720  CONTINUE  
*
*       And 3-way interactions
*
      DO 730 J=NP(5)+1,NP(6)
       IF ( (THREEWAY(J,1).LE.NP(4)).AND.
     +      (THREEWAY(J,2).LE.NP(4)).AND.
     +      (THREEWAY(J,3).LE.NP(4))) THEN
        IF ( (THREEWAY(J,1).LE.NP(3)).AND.
     +       (THREEWAY(J,2).LE.NP(3)).AND.
     +       (THREEWAY(J,3).LE.NP(3))) GOTO 730
        DO 735 I=1,NSITES
         PRED(I,4) = PRED(I,4) + (BETA(J)*COVS(I,J))
 735    CONTINUE
       ENDIF
 730  CONTINUE  

      END
******************************************************************************
******************************************************************************
******************************************************************************
      Subroutine ReadData2(UnitNos,MissVal,NSITES,NVARS,TmpFlNo,
     +                    DateRange,ImputeTill,MaxLag,IFAIL)
******************************************************************************
*       To read a dataset and write a binary scratch file (which has
*       already been opened) containing all of the data along with 
*       markers for values to be simulated. Arguments:
*
*       UnitNos Contains channel numbers of input data files. See GLMFIT
*               header for details. The scratch file here corresponds to 
*               the unit in UnitNos(95)
*       MissVal Value corresponding to missing data in the input file. 
*       NSITES  Number of sites
*       NVARS   Number of variables in data file
*       TmpFlNo Number of temporary file containing 4-character site codes
*       DateRange First and last months to simulate, in format YYYYMM.
*       ImputeTill - Last date for imputation. 
*       MaxLag  Number of previous days' observations needed to form 
*               a linear predictor - used to determine the earliest 
*               date for which data should be stored. 
*       IFAIL   Error flag
******************************************************************************
      Integer, intent(in) :: UnitNos(100),NSITES,NVARS,DateRange(2)
      Integer, intent(in) :: ImputeTill,Maxlag,TmpFlNo
      Double precision, intent(in) :: MissVal
      Integer, intent(out) :: IFAIL
******************************************************************************
*       Additional INTEGERs
*       -------------------
*	YF,MF,} Year, month and day read from file
*	DF    } 
*	YY,MM,} Year, month and day to write out
*	DD    } 
*       SITE    Current site
*       DONE    Indicator for whether we have finished reading data
*       JSTART  Julian day number at which to start storing output
*       JEND    Julian day number at which to end
*       JCUR    Number of current Julian day
*       JF      Julian day number corresponding to date read from file    
*       JIT     Julian day number corresponding to ImputeTill
*       J       Counter
*       DAYS_IN_MONTH   Function to return # of days in a given month 
******************************************************************************
      Integer YF,MF,DF,YY,MM,DD,SITE,DONE,J
      INTEGER JSTART,JEND,JCUR,JF,JIT,DAYS_IN_MONTH
******************************************************************************
*       Additional DOUBLEs
*       ------------------
*       DatFromFile Array of observations read from file
*       DatToWrite  Array of data to be written out 
******************************************************************************
      Double precision DatFromFile(NSITES,NVARS)
      Double precision DatToWrite(NSITES,NVARS)
***************************************************************************
*       Additional CHARACTERs
*       ---------------------
*       SCODES  Array of 4-character site identifiers
***************************************************************************
      Character SCODES(NSITES)*4
***************************************************************************
*
*       Initialise (NB use of Fortran 90 array capabilities)
*
      IFAIL = 0
      YF = 0
      MF = 0
      DF = 0
      JF = 0
      DONE = 0
*
*       Ensure that all files and pointers are positioned correctly
*
      REWIND UnitNos(1)
*
*       Read site codes from temporary file
*
      CALL GetScodes(TmpFlNo, NSites, Scodes, Ifail)
      IF (IFAIL.NE.0) RETURN
*
*       Find Julian day number for first day to be stored ...
*
      CALL DMYTOJ(1,MOD(DateRange(1),100),DateRange(1)/100,
     +                                              JSTART,IFAIL)
      IF (IFAIL.NE.0) RETURN
      IF (MaxLag.GE.0) JSTART = JSTART - MaxLag
*
*       ... and last day ...
*
      MM = MOD(DateRange(2),100)
      YY = DateRange(2) / 100
      CALL DMYTOJ(DAYS_IN_MONTH(MM,YY),MM,YY,JEND,IFAIL)
      IF (IFAIL.NE.0) RETURN
*
*       ... and for last day of imputation ...
*
      MM = MOD(ImputeTill,100)
      YY = ImputeTill / 100
      CALL DMYTOJ(DAYS_IN_MONTH(MM,YY),MM,YY,JIT,IFAIL)
      IF (IFAIL.NE.0) RETURN
*
*       Now loop over each day: read from input file unless (a) we've 
*       already got to the end of it (b) the last data we read from it
*       were for some date in the future (c) we've got beyond the last
*       date for which imputation is required. 
*
      Do JCUR=JSTART,JEND
       DatToWrite = -1.0D101
       IF ((DONE.EQ.0).AND.(JF.LT.JCUR).AND.(JCUR.LE.JIT)) THEN
*
*       If we're in the middle of reading the input file, need to increment
*       the date in order that DATRD will function correctly
*
 50     IF (JF.NE.0) THEN
         CALL JTODMY(JF+1,0,0,0,DF,MF,YF,IFAIL)
         IF (IFAIL.NE.0) RETURN
        ENDIF
        DatFromFile = -1.0D101
        CALL DATRD(UnitNos(1),DatFromFile,MissVal,
     +                   SCODES,NSITES,NVARS,YF,MF,DF,DONE,IFAIL)
        IF (IFAIL.NE.0) RETURN
        CALL DMYTOJ(DF,MF,YF,JF,IFAIL)
        IF ((JF.LT.JCUR).AND.(DONE.EQ.0)) GOTO 50
       ENDIF
       IF (JF.EQ.JCUR) DatToWrite = DatFromFile
*
*       Find the year, month and day corresponding to the current date
*
       CALL JTODMY(JCUR,0,0,0,DD,MM,YY,IFAIL)
       IF (IFAIL.NE.0) RETURN
       Do SITE=1,NSITES
        WRITE(UnitNos(95),REC=NSITES*(JCUR-JSTART)+SITE,ERR=920) 
     +                   SITE,YY,MM,DD,(DatToWrite(SITE,J),J=1,NVARS)
       End Do
      End Do
*
*       Reset input file
*
      Call FileReset(UnitNos(1))
      
      Return
      
 920  CALL INTPR('**** ERROR: while writing to scratch file. '//
     +            'You are probably out of disk space.',-1,0,0)

      End Subroutine ReadData2
******************************************************************************
******************************************************************************
******************************************************************************
      Subroutine WriteData(UnitNos,MissVal,NSITES,NVARS,
     +                    TmpFlNo,OUTOPT,DateLims,OutLims,
     +                    WHICHREGS,REGCODE,NREGS,IFAIL)
******************************************************************************
*       To read simulated data from a binary scratch file and
*       produce required output files / objects. Arguments:
*
*       UnitNos Contains channel numbers of input data files. See GLMFIT
*               header for details. The scratch file here corresponds to 
*               the unit in UnitNos(96)
*       MissVal Value to write for missing data.
*       NSITES  Number of sites
*       NVARS   Number of variables in data file
*       TmpFlNo Number of temporary file containing 4-character site codes
*       OutOpt  selects output options:
*               1 = Output summary information only (monthly and annual 
*                   totals, averaged over the whole region and/or various
*                   subregions)
*               2 = Write an output file with one record for every day,
*		    containing daily values for each site.
*                   Do not produce summary measures
*               3 = Do both.
*       OutLims	Two-element vector containing the start and end dates for 
*		simulation (scratch file may contain initialisation values
*               prior to start date). Format is YYYYMM
*       OutLims	Two-element vector containing the start and end dates for 
*		which daily output should be produced if it has been 
*		requested. Format is YYYYMM
*       WHICHREGS- Chooses which regions we want to generate summary statistics
*                  for. Element 0 is the whole area
*       NREGS   - Total number of regions (excepting the whole area)
*       IFAIL   Error flag
******************************************************************************
      Integer, intent(in) :: UnitNos(100),NSITES,NVARS
      Integer, intent(in) :: TmpFlNo,OutOpt
      Integer, intent(in) :: DateLims(2),OutLims(2)
      Integer, intent(in) :: REGCODE(NSITES),WHICHREGS(0:NSITES),NREGS
      Double precision, intent(in) :: MissVal
      Integer, intent(out) :: IFAIL
******************************************************************************
*       Additional INTEGERs
*       -------------------
*	YF,MF,} Year, month and day read from file
*	DF    } 
*	YY,MM,} Year, month and day to write out
*	DD    } 
*       SITE    Current site
*       RECNO   Number of record to read from scratch file
*       NRECS   Number of records in scratch file
*       CurDate Current date, in format YYYYMMDD
*       N       Numbers of observations contributing to monthly 
*               means / totals
*       J       Counter
*************************************Date*****************************************
      Integer YF,MF,DF,YY,MM,DD,SITE,RECNO,NRECS,CurDate,J
      Integer N(NSITES,NVARS,12)
******************************************************************************
*       Additional DOUBLEs
*       ------------------
*       DatFromFile Array of observations read from file
*       MoTot       Array of monthly totals
*       YrTot       Array of annual totals
******************************************************************************
      Double precision DatFromFile(NVARS)
      Double precision MoTot(NSITES,NVARS,12),YrTot(NSITES,NVARS)
***************************************************************************
*       Additional CHARACTERs
*       ---------------------
*       SCODES  Array of 4-character site identifiers
*       DLYFMT  Output format for daily data
***************************************************************************
      Character SCODES(NSITES)*4,DLYFMT*26
***************************************************************************
*
*       Initialise
*
      IFAIL = 0
      N = 0
      YY = -9999
      MoTot = 0.0D0
      YrTot = 0.0D0
      WRITE(DLYFMT,'(A13,I8,A5)') '(I4,I2,I2,A4,',NVARS,'F6.2)'
*
*       Read site codes from temporary file
*
      CALL GetScodes(TmpFlNo, NSites, Scodes, Ifail)
      IF (IFAIL.NE.0) RETURN
*
*       Now read through the scratch file a line at a time
*
      INQUIRE(UNIT=UnitNos(96),NEXTREC=NRECS)
      NRECS=NRECS-1
      
      REWIND UnitNos(20)
      Do RECNO=1,NRECS
       READ(UnitNos(96),REC=RECNO,ERR=999) 
     +                     SITE,YF,MF,DF,(DatFromFile(J),J=1,NVARS)
       If (YF.NE.YY) then
*
*       Write this year's summary stats to output file if required
*
        If ((YY.NE.-9999).AND.((OutOpt.Eq.1).OR.(OutOpt.Eq.3))) THEN
         YrTot = SUM(MoTot,DIM=3)
         CurDate = (100*YY) + MM

         If (CurDate.GE.DateLims(1).AND.CurDate.LE.DateLims(2)) then
          CALL OUTSUM(YY,YRTOT,MOTOT,N,NSITES,NVARS,REGCODE,
     +                         WHICHREGS,NREGS,Missval,UnitNos(21))
         End If
        End If
        N = 0
        MoTot = 0.0D0
        YrTot = 0.0D0
       End If
       YY = YF
       MM = MF
       DD = DF
*
*       Depending on output option selected, EITHER update 
*       summary stats for each site for this month, OR write 
*       the day's results to the output file (if within the 
*       range requested), OR both
*
       Do J=1,NVARS
        If (DatFromFile(J).LE.-1.0D100) then 
         DatFromFile(J) = MissVal
        Else if ( (OutOpt.EQ.1).OR.(OutOpt.EQ.3) ) then
         MoTot(SITE,J,MM) = MoTot(SITE,J,MM) + DatFromFile(J)
         N(SITE,J,MM) = N(SITE,J,MM) + 1
        End if
       End Do 
       If ( (OutOpt.EQ.2).OR.(OutOpt.EQ.3) ) THEN
        CurDate = (100*YY) + MM
        IF (CurDate.GE.OutLims(1).AND.CurDate.LE.OutLims(2)) THEN
         WRITE(UnitNos(20),DLYFMT) YY,MM,DD,SCODES(SITE),
     +                                   (DatFromFile(J),J=1,NVARS)
        End if
       End if
      End Do
*
*       Write final set of summary stats to output file if required
*
      If ((OutOpt.Eq.1).OR.(OutOpt.Eq.3)) THEN
       YrTot = SUM(MoTot,DIM=3)
       CurDate = (100*YY) + MM
       If (CurDate.GE.DateLims(1).AND.CurDate.LE.DateLims(2)) then
        CALL OUTSUM(YY,YRTOT,MOTOT,N,NSITES,NVARS,REGCODE,
     +                             WHICHREGS,NREGS,Missval,UnitNos(21))
       End if
      End If

      Return
      
 999  IFAIL = 999
      Return
      
      End Subroutine WriteData
*******************************************************************************
*******************************************************************************
*******************************************************************************
      SUBROUTINE OUTSUM(YEAR,YRTOT,MOTOT,N,NSITES,NVARS,REGCODE,
     +                             WHICHREGS,NREGS,Missval,FILNO)
*
*       Calculates annual and monthly summary statistics from the site stats
*       contained in YRTOT and MOTOT, and writes records to summary output
*       file FILNO for each region selected in WHICHREGS. NREGS is the 
*       total no. of regions defined; N contains numbers of observations
*       from which means are calculated; and Missval is the number to write
*       for a missing value
*
******************************************************************************
*       INTEGER variables
*       ^^^^^^^^^^^^^^^^^
*       YEAR    - current year
*       NSITES  - number of sites
*       NVARS   - number of variables
*       REGCODE- indicates which region each site is in
*       WHICHREGS - indicates which regions we want summary measures for
******************************************************************************
      Integer, intent(in) :: YEAR,NSITES,NVARS,NREGS,FILNO
      Integer, intent(in) :: REGCODE(NSITES),WHICHREGS(0:NSITES)
      Integer, intent(in) :: N(NSITES,NVARS,12)
      INTEGER I,J,K
******************************************************************************
*       DOUBLE PRECISION variables
*       ^^^^^^^^^^^^^^^^^^^^^^^^^^
*       YRTOT   - annual totals for each site
*       MOTOT   - monthly  "    "     "
*       MOMEAN  - monthly areal means } we have an element for each region.
*       YRMEAN  - annual areal means  } The number of regions is < NSITES
*       WEIGHTS - weights for averaging in each region (NB the dimensioning
*                 is for NSITES, the assumption being that there will always
*                  be fewer regions than sites)
*       YRWT    - weights for calculating yearly averages
******************************************************************************
      Double precision, intent(in) :: MOTOT(NSITES,NVARS,12)
      Double precision, intent(in) :: YRTOT(NSITES,NVARS),Missval
      DOUBLE PRECISION MOMEAN(0:NSITES,NVARS,12),YRMEAN(0:NSITES,NVARS)
      DOUBLE PRECISION WEIGHTS(0:NSITES,NVARS,12),YRWT
******************************************************************************
*       CHARACTER variables
*       ^^^^^^^^^^^^^^^^^^^
*       OutFormat       Format string for output
******************************************************************************
      Character OutFormat*40
******************************************************************************
*
*       Initialize ...
*
      YRMEAN = 0.0D0
      WEIGHTS = 0.0D0
      MOMEAN = 0.0D0
      
      WRITE(OutFormat,'(A13,I8,A19)') 
     +      '(I4,1X,I3,1X,',NVARS,'(13(F6.2,1X)))'
*
*       and calculate totals and weights for each region. This is a little bit
*       cunning ... for each site, the element of the ??MEAN array corresponding
*       to the region code is updated. The WEIGHTS for this region are only
*       incremented if this region has been chosen. If no region is defined
*       for a site (i.e. REGCODE is 0) we just update the whole area stats
*
      DO 150 I=1,NSITES
       DO 151 J=1,NVARS
        IF (REGCODE(I).NE.0) THEN
         YRMEAN(REGCODE(I),J) = YRMEAN(REGCODE(I),J) + YRTOT(I,J)
        ENDIF
        YRMEAN(0,J) = YRMEAN(0,J) + YRTOT(I,J)
        DO 155 K=1,12
         IF (REGCODE(I).NE.0) THEN
          WEIGHTS(REGCODE(I),J,K) = WEIGHTS(REGCODE(I),J,K)
     +                        + (DBLE(N(I,J,K)*WHICHREGS(REGCODE(I))))
          MOMEAN(REGCODE(I),J,K) = MOMEAN(REGCODE(I),J,K)+MOTOT(I,J,K)
         ENDIF
         WEIGHTS(0,J,K) = WEIGHTS(0,J,K) + (DBLE(N(I,J,K)*WHICHREGS(0)))
         MOMEAN(0,J,K) = MOMEAN(0,J,K) + MOTOT(I,J,K)
 155    CONTINUE
 151   CONTINUE
 150  CONTINUE
*
*       Now write a record for each region with non-zero weights. 
*
      DO 200 I=0,NREGS
       DO 201 J=1,NVARS
        YRWT = SUM(WEIGHTS(I,J,1:12))
        DO 205 K=1,12
         IF (WEIGHTS(I,J,K).GT.0.0D0) THEN
          MOMEAN(I,J,K) = MOMEAN(I,J,K) / WEIGHTS(I,J,K)
         ELSE
          MOMEAN(I,J,K) = Missval
         ENDIF
 205    CONTINUE
        IF (YRWT.GT.0.0D0) THEN
         YRMEAN(I,J) = YRMEAN(I,J) / YRWT
        ELSE
         YRMEAN(I,J) = MISSVAL
        ENDIF
 201   CONTINUE
       IF (WHICHREGS(I).GT.0) THEN
        WRITE(FILNO,OutFormat) 
     +        YEAR,I,((MOMEAN(I,J,K),K=1,12),YRMEAN(I,J),J=1,NVARS)
       ENDIF
 200  CONTINUE

      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE DailySim(Model,YSIM,WetDry,FLAGS,NSITES,NKREQ,
     +                          MEANS,StdDevs,CovMat,ANSMEAN,SIGMA11,
     +                          SIGMA12,SIGMA22,SIGMA,Recalc,IFAIL)
***************************************************************************
*       © UCL 2002-2018
*
*       This subroutine is the core of the simulation procedure. It 
*       generates a set of daily values for a subset of NSITES sites
*       (depending on the value of the FLAGS) and outputs the results 
*       in the array YSIM, which is is also used in input - if values 
*       are known for any site (indicated by a zero element of FLAGS) 
*       they will be preserved. 
*
*       Arguments:
*
*       Model   - Code defining type of distribution from which to 
*                 simulate. These are the codes as defined in R
*                 routine model.num(). 
*       YSIM    - array of daily values
*       WetDry  - array of zeroes and ones; values will be generated 
*                 only for positions containing ones. The idea is that
*                 for precipitation, the ones correspond to "wet" sites;
*                 for other variables, all elements of WetDry should be
*                 set to 1. 
*       FLAGS   - array indicating required status of each entry of YSIM:
*                 0 = Value known, don't change it (including because it
*                     may have been set to zero already if we're 
*                     dealing with a precipitation-like variable)
*                 1 = Value unknown, simulate from the model
*       NSITES  - no. of sites required (i.e. index of last entry in YSIM)
*       NKREQ   - Number of sites that are of interest (because WetDry=1)
*                 and where YSIM is recorded
*       MEANS   - vector of modelled means for each site.
*       StDevs  - vector of modelled standard deviations
*       ANSMEAN - vector of mean Anscombe residuals
*       CovMat  - Residual covariance matrix (correlation matrix for 
*                 joint mean-covariance models)
*       SIGMA11 } Partitions of residual covariance matrix 
*       SIGMA22 } (passed back to calling routine so that their 
*       SIGMA12 } values can be saved between calls)
*       SIGMA   Covariance matrix if conditioning on known values
*       Recalc  - Indicator for whether to recalculate matrix inverses
*       IFAIL   - Error flag
***************************************************************************
      Integer, intent(in) :: Model,NSITES,Recalc
      Integer, intent(in) :: FLAGS(NSITES),WETDRY(NSITES)
      Double precision, intent(in) :: MEANS(NSITES),StdDevs(NSITES),
     +                         ANSMEAN(NSITES),CovMat(NSITES,NSITES)
      Integer, intent(inout) :: NKREQ
      Double precision, intent(inout) :: YSIM(NSITES)
      Double precision, intent(inout), dimension(NSITES,NSITES) ::
     +                         SIGMA11,SIGMA22,SIGMA12,SIGMA
      Integer, intent(out) :: IFAIL
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>             INTERNAL VARIABLES        >>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
***************************************************************************
*       INTEGER variables
*       ^^^^^^^^^^^^^^^^^
*       NNeeded  - No. of sites we've got to invent something for
*       CurIdx   - Index of current site in INDEX ordering (see below)   
*       Ntries   - Number of attempts so far to generate a vector of 
*                  non-negative Anscombe residuals in gamma models
*       INDEX    - contains site numbers of elements of ZRESIDS etc. when
*                  we partition the wet sites.
***************************************************************************
      INTEGER NNeeded,CurIdx,Ntries,INDEX(NSITES)
      INTEGER I,J,K
***************************************************************************
*       DOUBLE PRECISION variables
*       ^^^^^^^^^^^^^^^^^^^^^^^^^^
*       ZRESIDS - Anscombe residuals for wet sites
*       ZBQLNOR - normal random number generator
*       MU      Mean residual vector if conditioning on known values
*       TMPARR? Temporary storage
***************************************************************************
      DOUBLE PRECISION ZRESIDS(NSITES),MU(NSITES),TMPARR1(NSITES)
      DOUBLE PRECISION ZBQLNOR
***************************************************************************
*       CHARACTER variables
*       ^^^^^^^^^^^^^^^^^^^
*       MESSAGE Error messages
***************************************************************************
      CHARACTER*255 MESSAGE(4)

***************************************************************************
*
*	If the response variable is binary then just copy the WETDRY array 
*	into YSIM and exit
*
      IF (MODEL.EQ.1) THEN
       YSIM(1:NSITES) = DBLE(WETDRY(1:NSITES))
       RETURN
      END IF
*
*       How many "unknown" sites do we not need to bother with (e.g.
*       because they are known to be dry)? Also set known dry sites to zero.
*
      NNeeded = 0
      NKREQ = 0

      DO 195 I = 1,NSITES
       IF (WETDRY(I).EQ.0) THEN
        YSIM(I) = 0.0D0
       ELSE IF (FLAGS(I).EQ.1) THEN
        NNeeded = NNeeded + 1
       ENDIF
 195  CONTINUE

      IF (NNeeded.EQ.0) RETURN
*
*       Now compute residuals for sites where values of YSIM are known. 
*       It is convenient to re-order things: the first NNeeded elements 
*       of the array correspond to the sites we're trying to infill, and 
*       the remainder are the sites which are known. We ignore sites 
*       where WetDry=0 (whether with a FLAG of zero or not), as these 
*       are excluded from the analysis when fitting precipitation 
*       intensity models and will tend to result in underestimation. 
*       For bookkeeping, we store the ordering in the array INDEX. NB if 
*       there aren't any unknown sites because WetDry=0 everywhere, bail 
*       out. We use NKREQ as a counter - which involves recalculating it, 
*       but it will end up the same as it was to start with.
*     
      CurIdx = 0
      DO 200 I = 1,NSITES
*
*       (i) Unknown sites for which simulation is required.
*
       IF ( (FLAGS(I).NE.0).AND.(WetDry(I).EQ.1) ) THEN
        CurIdx = CurIdx + 1
        INDEX(CurIdx) = I
        IF ((Model.EQ.10).OR.(Model.EQ.20)) THEN
         MU(CurIdx) = 0.0D0
        ELSE IF (MODEL.EQ.11) THEN
         MU(CurIdx) = ANSMEAN(I)
        ENDIF
*
*       (ii) Sites that must be considered and have known values
*
       ELSEIF ( (FLAGS(I).EQ.0).AND.(WetDry(I).EQ.1) ) THEN
        NKREQ = NKREQ + 1
        INDEX(NNeeded+NKREQ) = I
        IF ((Model.EQ.10).OR.(Model.EQ.20)) THEN
         ZRESIDS(Nneeded+NKREQ) = (YSIM(I)-MEANS(I)) / StdDevs(I)
         MU(Nneeded+NKREQ) = 0.0D0
        ELSE IF (Model.EQ.11) THEN
         ZRESIDS(Nneeded+NKREQ) = (YSIM(I)/MEANS(I))**(1.0D0/3.0D0)
         MU(Nneeded+NKREQ) = ANSMEAN(I)
        ELSE
         IFAIL = 999
         RETURN
        ENDIF
       ENDIF
 200  CONTINUE
*
*       Now calculate the components of the partitioned covariance matrix,
*       if we need to recalculate it (otherwise skip and go straight to
*       generation). Note: SIGMA11 has dimension Nneeded*Nneeded 
*	(corresponding to the unknown sites); SIGMA12 is Nnneeded*NKREQ; 
*	and Sigma22 is NKREQ*NKREQ.
*
      IF (Recalc.EQ.1) THEN
       DO 210 I=1,NNeeded
        DO 211 J=1,NNeeded
         SIGMA11(I,J) = CovMat(INDEX(I),INDEX(J))
 211    CONTINUE
        DO 212 J=1,NKREQ
         SIGMA12(I,J) = CovMat(INDEX(I),INDEX(J+Nneeded))
 212    CONTINUE
 210   CONTINUE
       DO 215 I=1,NKREQ
        DO 216 J=1,NKREQ
         SIGMA22(I,J) = CovMat(INDEX(I+Nneeded),INDEX(J+Nneeded))
 216    CONTINUE
 215   CONTINUE
*
*       Invert SIGMA22 if there are any knowns (this plays a part in both
*       the conditional mean and the conditional variance of the unknowns)
*
       IF (NKREQ.GT.0) CALL MATINV(SIGMA22,NKREQ,1,NSITES,1,IFAIL)
       IF (IFAIL.NE.0) RETURN
      ENDIF
*
*
*       Now we can calculate the conditional mean and covariance matrix of
*       the unknown quantities. TMPARR1(J) holds the (I,J)th element of
*       SIGMA12 * SIGMA22. Note that the ordering of everything in MU? and
*       SIGMA? is with unknowns first, whereas things like CovMat are in 
*       original site order. Note also that this will work if there are
*       no knowns - MU and SIGMA are effectively just copied then.
*

      DO 220 I=1,NNeeded
*
*       Mean ...
*
       DO 225 J=1,NKREQ
        TMPARR1(J) = 0.0D0
        DO 226 K=1,NKREQ
         TMPARR1(J) = TMPARR1(J) + (SIGMA12(I,K)*SIGMA22(K,J))
 226    CONTINUE
        MU(I) = MU(I) + 
     +          (( ZRESIDS(NNeeded+J)-MU(NNeeded+J) )*TMPARR1(J))
 225   CONTINUE
*
*       ... and covariance matrix. Formula requires SIGMA21, but this is
*       just the transpose of SIGMA12 so we use SIGMA12(J,K) in place of
*       SIGMA21(K,J). NB later SIGMA gets overwritten, so we only do this
*       bit if a RECALC is required
*
       IF (Recalc.eq.1) THEN
        DO 230 J=1,NNeeded
         SIGMA(I,J) = SIGMA11(I,J)
         DO 231 K=1,NKREQ
          SIGMA(I,J) = SIGMA(I,J) - (SIGMA12(J,K)*TMPARR1(K))
 231     CONTINUE
 230    CONTINUE
       ENDIF
 220  CONTINUE
*      
*       If the Cholesky decomposition of SIGMA needs to be recomputed 
*       because of a RECALC, do it now
*
      IF (Recalc.eq.1) THEN 
       CALL CHOLDEC(SIGMA,NNeeded,NSITES,IFAIL)
       IF (IFAIL.NE.0) THEN
        WRITE(MESSAGE,2) INDEX(IFAIL)
        CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(3)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(4)),-1,0,0)
        IFAIL = 200
        RETURN
       ENDIF
      ENDIF

*
*       Almost there! Next step is to generate a random multivariate
*       normal vector for the unknowns, if this is required. 
*       We do it by generating a vector of iid normals
*       in TMPARR1, then premultiplying by any matrix A such that A.A'
*       = SIGMA. We use the transpose of the upper triangular matrix in
*       the Cholesky decomposition, which overwrites SIGMA and gets
*       stored for future use - this also speeds things up as half its
*       elements are zero so we don't have to do the sums for them.
*       Result goes into ZRESIDS, the appropriate elements of which will
*       then be back-transformed into YSIM.
*
*       NB only update if FLAGS and WetDry are both 1 - otherwise, 
*       the value is already known
* 
      Ntries = 0
 250  DO 255 I=1,NNeeded
       ZRESIDS(I) = MU(I)
       TMPARR1(I) = ZBQLNOR(0.0D0,1.0D0)
 255  CONTINUE
*
*       ZBQLNOR generates *pairs* of deviates using the Box-Muller
*       algorithm, and saves the unused one for reuse on the next
*       call. This is done for efficiency, but it means that 
*       simulations potentially are not reproducible unless you
*       can guarantee that there's no spare deviate sitting there
*       at the end of the last run. The way to guarantee this is
*       to ensure that you always generate an even number of deviates.
*
      IF (2*(NNeeded/2).NE.NNeeded) 
     +                  TMPARR1(NNeeded) = ZBQLNOR(0.0D0,1.0D0) 
*
*       Multiply TMPARR by transpose of upper triangular matrix and 
*       add to ZRESIDS
*
      DO 256 I=1,NNeeded
       DO 257 K=1,I
        ZRESIDS(I) = ZRESIDS(I)+(SIGMA(K,I)*TMPARR1(K))
 257   CONTINUE
 256  CONTINUE

*
*       Et enfin - the back-transformation. 
*
      DO 300 I=1,NNeeded
       IF ( (FLAGS(INDEX(I)).NE.0).AND.(WetDry(INDEX(I)).EQ.1) ) THEN
        IF ((Model.EQ.10).OR.(Model.EQ.20)) THEN
         YSIM(INDEX(I)) = MEANS(INDEX(I)) + 
     +                         (StdDevs(INDEX(I)) * ZRESIDS(I))
        ELSE IF (Model.EQ.11) THEN
*
*       For gamma distributions, condition on the values being strictly positive
*       as this is what the data yield. If any negative values are obtained, 
*       the correct thing to do is to throw away the entire current simulation 
*       and start again. For small numbers of sites this is usually entirely
*       feasible; however, with large numbers of sites it can sometimes be
*       *very* difficult to generate a multivariate normal realisation with
*       no negative values (particularly when carrying out imputations when 
*       several other sites record zeroes). To avoid getting stuck in 
*       infinite loops here therefore, if we've had no luck after 10+sqrt(Nneeded) 
*       trials then just return a zero. This is slightly unsatisfactory but 
*       should not affect results noticeably (because if it keeps generating
*       zeroes then clearly it's trying to put a very small value in there
*       somewhere)
*
         IF (ZRESIDS(I).LT.0.0D0) THEN 
          Ntries = Ntries + 1
          If (Ntries.lt.int(10+dsqrt(dble(Nneeded)))) goto 250
          ZRESIDS(I) = 0.0D0
         END IF
         YSIM(INDEX(I)) = MEANS(INDEX(I))*( ZRESIDS(I)**3 )
        ELSE
         IFAIL = 999
         RETURN
        ENDIF
*
*       Trap values that won't fit in fixed-format output files - 
*       they are symptomatic of exploding simulations. The
*       final "not" in the IF clause is to trap things that
*       are not a number e.g. due to illegal arithmetic
*       operations. 
*
        IF ((YSIM(INDEX(I)).LE.-1.0D2) .OR. 
     +      (YSIM(INDEX(I)).GE.1.0D3) .OR.
     +      (.NOT.(DABS(YSIM(INDEX(I))).LT.1.0D8))) THEN
         IFAIL = 400
         RETURN
        ENDIF
       ENDIF
 300  CONTINUE
            
      RETURN      
      
 2    FORMAT(5X,'****ERROR**** input correlation matrix isn''t ',
     +'positive definite.',/5X,'Problem found at site ',I4,'. ',
     +'If you are using an empirical',/5X,'correlation matrix, you ',
     +'may need to fit a valid correlation model to',/5X,
     +'get round this problem. Otherwise, check for duplicate sites.')
      END
c***************************************************************************
c***************************************************************************
c***************************************************************************
c      SUBROUTINE MAKEMO(NSITES)
c*
c*     To make monthly time series of areal rainfall totals from GLM
c*     simulated daily data files
c*
c      INTEGER NSITES,NWANT
c      PARAMETER(NWANT=10)
c
c      INTEGER SIM,YEAR,MONTH,DAY,OLDMO,OLDYR,OLDSIM,I
c      INTEGER SITLST(NWANT)
c      REAL RAIN(NSITES),RAINMO
c
c      DATA(SITLST(I),I=1,NWANT) /5,6,13,15,16,18,19,20,27,33/
c
c      OPEN(10,FILE="daily.sim")
c      OPEN(11,FILE="monthly.sim")
c
c      OLDMO = 0
c*
c*     The extra 2 entires of RAIN are for 3- and 10-site average
c*
c 10   READ(10,'(I3,1X,I4,1X,2(I2,1X),100(F6.2,3X))',END=99)
c     +                    SIM,YEAR,MONTH,DAY,(RAIN(I),I=1,NSITES)
c      IF (MONTH.NE.OLDMO) THEN
c       IF (OLDMO.NE.0) THEN
c        RAINMO = RAINMO/REAL(NWANT)
c        WRITE(11,'(I3,1X,I4,1X,I2,1X,F6.2)') 
c     +                           OLDSIM,OLDYR,OLDMO,RAINMO
c       ENDIF
c       IF (SIM.NE.OLDSIM) WRITE(*,'(''Simulation '',I3,''...'')') SIM
c       RAINMO = 0.0
c      ENDIF
c
c      DO 20 I=1,NWANT
c       RAINMO = RAINMO + RAIN(SITLST(I))
c 20   CONTINUE
c      OLDSIM = SIM
c      OLDYR = YEAR
c      OLDMO = MONTH
c      GOTO 10
c
c 99   CLOSE(10)
c      RAINMO = RAINMO/REAL(NWANT)
c      WRITE(11,'(I3,1X,I4,1X,I2,1X,F6.2)') 
c     +                           OLDSIM,OLDYR,OLDMO,RAINMO
c
c      END

