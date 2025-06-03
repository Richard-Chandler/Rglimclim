      SUBROUTINE GLMFIT(MODEL,MAXIT,YCHECK,NSITES,NVARS,MissVal,
     +                  RespIdx,AllowIncAvge,MXP,NATTR,SITINF,
     +                  TmpFlNo,NP,BETA,COVCODE,
     +                  TWO,THREE,THETA,CovNaive,CovRobust,LogLik,
     +                  Deviance,NOBS,P,DF,SITXFM,FOUIDX,LEGIDX,MXLAG,
     +                  PWTIDX,PHI,YBAR,SY,ERRBAR,MSE,RSQ,PRBAR,SPR,
     +                  SEPR,ARBARO,SARO,ARBARE,SARE,BYPRED,PRMNTH,
     +                  PRYEAR,PRSITE,NYEARS,FIRSTYR,SPMOD,RHO,GLBCOD,
     +                  GLBVAL,DoFit,DoResid,UnitNos,NeedRead,
     +                  DateRange,Verbosity,IFAIL)
******************************************************************************
*       To fit a GLM to daily climate time series.
*
*      INPUT ARGUMENTS
*      ^^^^^^^^^^^^^^^
*       MODEL   - Code for model to be fitted. Options:
*                 1     Logistic regression model for binary responses
*                 10    Gaussian model with identity link
*                 11    Gamma GLM with log link
*	MAXIT	- Maximum no of iterations permitted. Negative values
*                 lead to unlimited iterations
*       YCHECK  - Code to determine whether to check the values of the
*                 response variable. Options:
*                 0     No checking performed
*                 1     Checking performed and response variables 
*                       calculated as follows:
*                       - If model=1, any strictly positive value is
*                         set to 1 and all other values are set to
*                         zero
*                       - If model=11, only strictly positive values 
*                         are retained
*       NSITES  - Number of sites in study region
*       NVARS   - Number of variables in input data file
*       MissVal - Value used in the data file to denote a missing observation
*       RespIdx - Index of response variable in input data file
*       AllowIncAvge    - Indicator for whether or not weighted averages
*                 of non-response variables will be considered as 
*                 valid in model fitting if the variable is missing at
*                 the current site 
*       NATTR   - Number of attributes defined for each site
*       SITINF  - Double precision array containing information 
*                 regarding each site (Eastings, Northings etc.). 
*                 Column 0 holds the actual value, 1 and 2 hold 
*                 derivatives wrt any parameters of nonlinear
*                 transformations
*       TmpFlNo - Number of temporary file containing list of 4-character
*                 site codes
*       NP      - Integer array indicating numbers of parameters 
*                 in each model component. Elements are:
*                 1 - No. of main effects associated with site
*                 2 - NP(1) + No. of main effects associated with year
*                 3 - NP(2) + No. of main effects associated with month
*                 4 - NP(3) + No. of main effects associated with previous
*                               days' response variable (not trace values)
*                 5 - NP(4) + No. of 2-way interactions
*                 6 - NP(5) + No. of 3-way interactions
*		  7 - NP(6) + No. of nonlinear parameters
*                 8 - Number of `global' parameters such as 
*                     trace thresholds
*		  9 - Number of dispersion parameters
*		 10 - No. of spatial correlation parameters
*       MXP       Maximum number of parameters in each model (integer)
*       BETA    - Current estimate of parameter vector (double 
*                 precision; element 0 is the constant term)
*       COVCODE - Integer array containing code numbers of covariates
*                 (NB this potentially gest modified in routine 
*                 ATTRXFM)
*       TWO     - Integer array giving indices of two-way interactions 
*                 in model
*       THREE   - guess ...
*	THETA	- Double precision array of parameters in nonlinear
*                 functions of basic covariates. Each function has up 
*                 to 2 parameters. The elements of THETA are arranged 
*                 to line up with the covariates to which they
*                 correspond. For parameters we're estimating rather 
*                 than assuming to be known, THETAs duplicate the 
*                 last few elements of BETA - basically because its
*                 easier to do the bookkeeping if they're in BETA as
*                 well, as they can then be updated as part of the
*		  IWLS procedure.
*       CovNaive  - Double precision array containing covariance matrix
*                 of parameter estimates
*       LogLik  - Double precision: Independence log-likelihood for 
*                 the fitted model
*       Deviance- Double precision: Deviance computed from the 
*                 independence log-likelihood
*	NOBS	- Integer: no. of observations
*       P       - Integer: number of parameters estimated (excluding
*                 constant and dispersion parameters)
*       DF      - Integer: residual degrees of freedom for model
*	SITXFM	- Integer array selecting nonlinear transformations of
*                 site attributes. First column chooses the attribute,
*                 second the xfrmation.
*	FOUIDX	- Integer array indexing Fourier representation of 
*                 site effects. See file funcdefs.f for full details
*	LEGIDX	- Integer array indexing Legendre polynomial
*                 representation of site effects.
*	MXLAG	- Integer array, indicating number of immediately previous 
*                 days' values of each variable required to be present at 
*                 a site in order for it to be included in the fitting.
*                 This allows comparison of models with different 
*                 numbers of lags.
*       PWTIDX  - Integer array indexing  weighted averages of previous
*                 days' values ((I,0,K)th entry is the number of the
*                 covariate to which the parametrisation for the Ith
*                 weighting scheme for variable K is attached; (I,J,K)th 
*                 is the model number of the Jth parameter in the scheme).
*       YBAR    - Mean of Ys
*       SY      - Std Dev of Ys
*       ERRBAR  - Mean prediction error
*       MSE     - Mean squared error
*       RSQ     - R-squared (1 - Error sum of squares/Data sum of squares)
*       PRBAR   - Mean Pearson residual
*       SPR     - Std Dev of pearson residuals
*	SEPR	- Standard Error of mean Pearson residual
*       PRMNTH  - Pearson residuals by month } Column 0 is N, 1 is mean,
*       PRSITE  - And by site                } 2 is Std Dev, 3 is S.e. (mean)
*       PRYEAR  - And by year                } if model is correct
*       ARBARO  - Observed mean Anscombe residual (gamma model)
*       SARO    - Observed Std dev of Anscombe residuals (gamma model)
*       ARBARE  - Expected mean Anscombe residual (gamma model)
*       SARE    - Expected Std dev of Anscombe residuals (gamma model)
*	BYPRED	- Table of occurrence frequencies by predicted       } Logistic
*		- probability. Row 1 is observed proportion in each  } model
*		- category, Row 2 is expected, and Row 3 is N.       } 
*       PHI     - Dispersion parameter (double precision)
*	SPMOD   - Integer, indicating choice of structure for 
*                 inter-site dependence
*	RHO	- Double precision array containing parameters defining
*                 inter-site dependence structure. Dimensioned to
*		  MXP purely for convenience
*       GLBCOD  - Integer array of codes defining `global' quantities 
*                 in GLBVAL
*       GLBVAL  - Double precision array containing values of `global'
*                 quantities such as trace thresholds
*       DoFit   - perform parameter estimation? (1=yes,0=no)
*       DoResid - Do a residual analysis? (1=yes, 2=no)
*       UnitNos - Numbers of units attached to input and output files
*                 (which should be opened prior to calling this routine).
*                 For a description of the unit numbers, see the header
*                 for the CheckFiles() routine in the R code. 
*       NeedRead- Do we need to read the data from file (1) or have the 
*                 required scratch files already been produced (0)? 
*       DateRange 2-element vector (output) containing the dates of the 
*                 last observation used in fitting, in the form YYYYMMDD
*       Verbosity Controls verbosity of output (larger values: more verbose)
*       IFAIL   - error flag
***************************************************************************
      INTEGER, INTENT(IN) :: MODEL,MAXIT,YCHECK,AllowIncAvge
      Integer, intent(in) :: NSITES,NATTR,Verbosity,TmpFlNo
      DOUBLE PRECISION SITINF(NSITES,MXP,0:3)
      INTEGER, INTENT(IN) :: MXP, NP(10), NVARS, RespIdx
      INTEGER, INTENT(INOUT) :: DateRange(2),COVCODE(MXP)
      INTEGER, INTENT(IN) :: SPMOD,MXLAG(NVARS), GLBCOD(MXP)
      INTEGER, INTENT(IN) :: SITXFM(MXP,2),TWO(MXP,2),THREE(MXP,3)
      INTEGER, INTENT(IN) :: FOUIDX(MXP),LEGIDX(MXP)
      INTEGER, INTENT(IN) :: PWTIDX(MXP,0:3,NVARS), DoFit, DoResid
      INTEGER, INTENT(IN) :: FIRSTYR,NYEARS,UnitNos(100),NeedRead
      INTEGER, INTENT(OUT) :: NOBS,P,DF,IFAIL
      DOUBLE PRECISION, INTENT(INOUT) :: BETA(0:MXP), THETA(MXP,3)
      DOUBLE PRECISION, INTENT(OUT) :: CovNaive(0:MXP,0:MXP)
      DOUBLE PRECISION, INTENT(OUT) :: CovRobust(0:MXP,0:MXP)
      DOUBLE PRECISION, INTENT(OUT) :: LogLik, Deviance
      DOUBLE PRECISION, INTENT(INOUT) :: PHI, RHO(MXP)
      Double precision, intent(out) :: YBAR,SY,ERRBAR,MSE,RSQ
      Double precision, intent(out) :: PRBAR,SPR,SEPR
      Double precision, intent(out) :: ARBARO,SARO,ARBARE,SARE
      Double precision, intent(out) :: BYPRED(3,0:9),PRMNTH(12,0:3),
     +      PRSITE(NSITES,0:3),PRYEAR(FIRSTYR:FIRSTYR+NYEARS-1,0:3)
      DOUBLE PRECISION, INTENT(IN) :: MissVal,GLBVAL(MXP)
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>    PARAMETERS      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*       MXP2            Maximum number of parameters. THIS IS NEEDED
*                       ONLY FOR DIMENSIONING CHARACTER ARRAYS IN THE
*                       DESCRIBE COMMON BLOCK, AND CAN BE REMOVED
*                       AFTER SHIFTING ALL CHARACTER MANIPULATIONS TO
*                       R
*
*	NB IF THESE VALUES ARE CHANGED, THE CORRESPONDING ENTRY IN
*       OTHER SOURCE FILES WILL ALSO NEED CHANGING. 
******************************************************************************
        INTEGER MXP2
        PARAMETER (MXP2=80)
******************************************************************************
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>      VARIABLES     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
***************************************************************************
*               DOUBLE PRECISION variables
*               ^^^^^^^^^^^^^^^^^^^^^^^^^^
*       Y       Vector of all daily responses
*       ETA     Linear predictor for the current case
*       MU      Predicted mean for current case
*       SD      Standard deviation for current case
*       PI      Pi. Used for calculating seasonal components.
*       THRESH  Threshold used for dealing with small positive values
*       TRACE   Ditto, when we're trating them as `trace' values
*       Distance- array of inter-site distances. First slice is X-separation,
*                 second slice is Y-separation and third slice is Euclidean
*                 separation
******************************************************************************
      DOUBLE PRECISION Y,THRESH,TRACE,ETA,MU,SD
      Double precision Distance(3,NSITES,NSITES)
******************************************************************************
*               INTEGER variables
*               ^^^^^^^^^^^^^^^^^
*       I       Counter
*       TRENDSEL        Selects long-term trend function.
*       TRENDINTR       Selects which components long-term
*                       trend will interact with. When they're
*                       selected it contains the numbers of the
*                       components
*       TOGGLE          To toggle elements of TRENDINTR
*       PARMIDX         To partition the covariates into main, trace &
*                       trend segments
*	SITE	- number of current site
*	YY,MM,DD - year, month and day (used for calculation of
*		   previous days' response variables)
*       DAYS_IN_MONTH   - no. of days in current month (integer function)
*       THRTYP  - Type of thresholding to be carried out on small
*                 positive values. 0: no thresholding; 1: `soft'
*                 thresholding (zero everything below threshold, 
*                 subtract threshold from everything above); 2:
*                 `hard' thresholding (just zero everything below
*                 threshold)
*       ICHECK   - Indicator for whether model parametrisation has been 
*                  checked
***************************************************************************
      INTEGER SITE,YY,MM,DD,THRTYP
      INTEGER I,ICHECK
******************************************************************************
***************************************************************************
***************************************************************************
*               CHARACTER variables
*               ^^^^^^^^^^^^^^^^^^^
*       TRENDLBL- Label for long-term trend function used
*	TRPTXT	- Describes nonlinear parameters in trend functions
*       OPTION  - Used for 'Y/N' input
*       SCODES  - array of short site codes
*       ATTRTXT - Text for site attributes
*	SXPTXT	- Text for parameters in nonlinear site transformations
*       MOTXT   - Text for monthly model components
*       TRNDTXT - Text for trend components
*	TRPTXT	- Text for nonlinear parameters in trend components
*       DYTXT   - Text for previous day components
*       PWTTXT  - Text for parameters in nonlinear weighting schemes for
*                 averaging previous days' values
*	SPATXT	- Text for spatial correlation models
*	SPPTXT	- Text for parameters in spatial correlation 
*		  model (dimension to MXP for convenience)
*       CORTXT  - array describing quantities to which inter-site 
*                 correlations correspond
*       GLBTXT  - Text for global quantities
*       MESSAGE - Used for on-screen messages
***************************************************************************
      CHARACTER SCODES(NSITES)*4
      CHARACTER*70 ATTRTXT(MXP)
      CHARACTER*70 MOTXT(MXP2),TRNDTXT(MXP2),DYTXT(0:MXP2)
      CHARACTER*70 PWTTXT(MXP2,3),TRPTXT(MXP2,3)
      CHARACTER*70 SXPTXT(MXP2,3),GLBTXT(MXP2)
      CHARACTER*70 SPATXT(0:MXP2),SPPTXT(MXP2,4)
      CHARACTER*70 CORTXT(MXP2)
      CHARACTER MESSAGE*255
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>    COMMON BLOCKS   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /DESCRIBE/ SXPTXT,MOTXT,TRNDTXT,TRPTXT,DYTXT,
     +                  SPATXT,SPPTXT,CORTXT,PWTTXT,GLBTXT
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>    FILE NUMBERS    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*       The following FORTRAN file numbers are allocated in this program:
*       UNIT    DESCRIPTION
*       ----    -----------
*	99	Temporary scratch file used for doing character/numeric
*		conversion
******************************************************************************
******************************************************************************
******************************************************************************
**********                  EXECUTABLES START HERE                  **********
******************************************************************************
*
*       Initialise
*
      ICHECK = 0
      DO I=1,NATTR
       WRITE(ATTRTXT(I),'(''~[[site attribute#'',I4,''#]]~'')') I
      END DO
*
*       Read site codes from temporary file
*
      CALL GetScodes(TmpFlNo, NSites, Scodes, Ifail)
      IF (IFAIL.NE.0) RETURN

*
*       Calculate any requested nonlinear transformations of site attributes,
*       and inter-site distance information
*
      CALL ATTRXFM(SITINF,NSITES,ATTRTXT,NATTR,COVCODE,SITXFM,
     +                      FOUIDX,LEGIDX,THETA,NP,ICHECK,IFAIL,MXP)

      IF (IFAIL.NE.0) GOTO 990

      Call DistCalc(SITINF,NSITES,NATTR,MXP,Distance)
*
*	P is the number of parameters estimated. If NP(7) is 0,
*	there are no nonlinear parameters in the model and there
*	are NP(6) parameters to be estimated.
*
      P = NP(7)
      IF (P.EQ.0) P = NP(6)
*
*     Define thresholds for small positive values, if required
*
      CALL TAUDEF(GLBVAL,GLBCOD,NP(8),THRESH,TRACE,THRTYP,MODEL,MXP)
*
*       Read data, and set up design matrix. We read the data a day at a time,
*	and then copy the appropriate bits to the X and Y arrays.
*

      If (NeedRead.Eq.1) then 
       If (Verbosity.Gt.0) CALL INTPR('Reading data ...',-1,0,0)
       Call ReadData(UnitNos,MODEL,MissVal,RespIdx,AllowIncAvge,
     +               YCHECK,THRESH,THRTYP,NSITES,MXP,MXLAG,NVARS,
     +               P,NP,COVCODE,TWO,THREE,PWTIDX,SCODES,SITINF,
     +               Distance,BETA,THETA,TRACE,ICHECK,NOBS,DateRange,
     +               IFAIL)
     
       IF (IFAIL.NE.0) GOTO 990
       IF (NOBS.EQ.0) GOTO 930
      EndIf  
*
*       And start estimating.
*
      If (DoFit.Eq.1) then 
       If (Verbosity.Gt.0) CALL INTPR('Beginning estimation ...',-1,0,0)
       If (Verbosity.Gt.1) then
        WRITE(MESSAGE,'(A15,A20,A30)') 'Iteration','Log-likelihood',
     +                         'Largest standardised score'
        CALL INTPR(TRIM(MESSAGE),-1,0,0)
        WRITE(MESSAGE,'(A15,A20,A30)') '---------','--------------',
     +                         '--------------------------'
        CALL INTPR(TRIM(MESSAGE),-1,0,0)
       End if
       CALL IWLS(MODEL,RespIdx,AllowIncAvge,NVARS,MissVal,NP,P,
     +          COVCODE,TWO,THREE,FOUIDX,LEGIDX,PWTIDX,MXP,MAXIT,
     +          MXLAG,NOBS,NSITES,NATTR,SITXFM,GLBCOD,BETA,THETA,
     +          CovNaive,CovRobust,PHI,LogLik,Deviance,GLBVAL,
     +          SITINF,Distance,ATTRTXT,SCODES,UnitNos,Verbosity,IFAIL)
       IF (IFAIL.NE.0) GOTO 990
*
*	Do estimation of spatial dependence structure if required,
*	then do final model output.
*
       IF (SPMOD.NE.0) THEN
        IF (MODEL.EQ.1) THEN
         IF ((SPMOD.EQ.21).OR.(SPMOD.EQ.22)) THEN 
          CALL SPEST(UnitNos(91),NSITES,NOBS,SPMOD,RHO,MXP,
     +                                          Verbosity,IFAIL)
          IF (IFAIL.NE.0) GOTO 990
         ELSEIF (SPMOD.GE.1) THEN
          CALL CORLAT(UnitNos(91),NSITES,NOBS,SPMOD,SCODES,RHO,
     +                   UnitNos(6),SITINF,NATTR,MXP,Verbosity,IFAIL)
          IF (IFAIL.NE.0) GOTO 990
         ENDIF
        ELSEIF ((MODEL.EQ.10).OR.(MODEL.EQ.11)) THEN
         CALL COREST(UnitNos(91),NSITES,NOBS,MODEL,SPMOD,SCODES,RHO,
     +                   UnitNos(6),SITINF,NATTR,MXP,Verbosity,IFAIL)
         IF (IFAIL.NE.0) GOTO 990
        ENDIF
       ENDIF
*
*       Reset dispersion parameter for logistic model (-ve value in 
*       calling R code means that it's fixed)
*
       IF (MODEL.EQ.1) PHI = -1.0D0
      End If
*
*       At this stage, we should be able to calculate DF regardless of
*       whether or not the estimation was done
*
      DF = NOBS - (P+1)
*
*	Return unless residual analysis was requested
*
      IF (DoResid.EQ.0) Goto 920

      CALL RESID(MODEL,DABS(PHI),P+1,NOBS,UnitNos(91),NSITES,
     +                     YBAR,SY,ERRBAR,MSE,RSQ,PRBAR,SPR,SEPR,
     +                     ARBARO,SARO,ARBARE,SARE,BYPRED,PRMNTH,
     +                     PRSITE,PRYEAR,FIRSTYR,NYEARS,IFAIL)
      IF (IFAIL.NE.0) GOTO 990

*
*	Return unless user asked for output file for further residual 
*       analysis. NB we have already ascertained that any existing file 
*       can be overwritten; the final column of scratch file 91 contains
*       the weights for each case, which are proportional to the inverse
*       of the variances.
*
      If (DoResid.EQ.1) Goto 920

      REWIND(UnitNos(2))
      WRITE(UnitNos(2),7)
      DO 500 I=1,NOBS
       READ(UnitNos(91),REC=I) SITE,YY,MM,DD,Y,MU,ETA,SD
       SD = DSQRT(PHI / SD)
       WRITE(UnitNos(2),8) SCODES(SITE),YY,MM,DD,Y,MU,SD
 500  CONTINUE
*
*       Reset pointers for external files and close any remaining 
*       scratch files relating to inter-site correlations
*
 920  Call FileReset(UnitNos)
      RETURN
*
*	Error trapping
*

 930  IFAIL = 7

 990  Call FileReset(UnitNos) 
      RETURN

 7    FORMAT('SITE YEAR MONTH DAY OBSERVED PREDICTED        SD')
 8    FORMAT(A4,1X,I4,T14,I2,T17,I3,T21,F8.4,T31,F8.4,T41,F8.4)
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE SPEST(FILNO,NSITES,NOBS,SPMOD,RHO,MXP,Verbosity,IFAIL)
*
*   To estimate spatial dependence structure so as to match the 
*   variance of the histogram of proportions of wet sites.
*
******************************************************************************
*               INTEGER variables
*               ^^^^^^^^^^^^^^^^^
*   I,J - Counters
*   FILNO - File number of scratch file with residual info in
*   MXP - Maximum number of parameters
*   MXN - Maximum number of sites
*   NOBS - no. of observations
*   NSITES - no. of sites required
*   NDAYS - No of days
*   SPMOD - spatial dependence model selected
*   CASEID - array of identifiers for each case, for purposes of
*            residual analysis. Col.1 contains site number, 2 
*            contains year, 3 contains month and 4 day
*   OLDDAT }- Used to determine when the date changes
*   NEWDAT }  
*   SCRFNO	- Channel number of scratch file to store daily stats
*   BRACK	- Indicator in iteration to find out whether we've
*		  yet bracketed the root of an equation.
*   ITER	- Iteration number
*   Verbosity   - Controls verbosity of output
*   IFAIL       - Error code
***************************************************************************
      INTEGER NSITES,NOBS,FILNO,SPMOD,MXP
      INTEGER CASEID(4),OLDDAT,NEWDAT,NDAYS
      INTEGER SCRFNO,Verbosity,IFAIL
      INTEGER BRACK,ITER
      INTEGER I,J
***************************************************************************
*               DOUBLE PRECISION variables
*               ^^^^^^^^^^^^^^^^^^^^^^^^^^
*   Y       Daily response variable
*   PROBS   Predicted probabilities
*   MSY	    Depends on the model:
*              Model 21 - Mean square of observations
*              Model 22 - Sum of squared `Pearson residuals'
*   MSP	    Model-dependent:
*              Model 21 - Mean square of probs
*              Model 22 - Sum(No. of sites - 1).Just using it because
*              it's convenient storage that happens to exist.
*   YBAR    Mean of observations each day
*   PBAR    Mean of predictions each day
*   NT	    No of observations each day
*   RHO	    Vector of parameter estimates for spatial dependence
*   ONEDAY  Today's amounts. Col 1 is observed, 2 is predicted
*   VAR0    Theoretical variance of proportion of wet days under
*           independence
*   SPVAR   Expected value of mean square proportion of wet sites
*   LB,UB   Lower and upper bounds for parameter
*   FL,FB   Objective function evaluated at LB and UB
*   FCUR    Objective function evaluated at current point.
*   TOL	    Accuracy required.
*   VPWET   Theoretical variance of proportion of wet sites under
*           current model parametrisation.
*   TMP     Temporary storage
******************************************************************************
      DOUBLE PRECISION Y,PROBS
      DOUBLE PRECISION MSY,MSP,YBAR,PBAR,NT
      DOUBLE PRECISION ONEDAY(NSITES,2)
      DOUBLE PRECISION RHO(MXP),VAR0,SPVAR
      DOUBLE PRECISION LB,UB,FL,FU,FCUR,TOL
      DOUBLE PRECISION TMP
***************************************************************************
*     CHARACTER variables
*     ^^^^^^^^^^^^^^^^^^^
*     MESSAGE   - For screen output
***************************************************************************
      CHARACTER MESSAGE*255
*
*	Find an available file handle and open a scratch file on it 
*	(start searching at number 85)
*
      IFAIL = 0
      SCRFNO = 85
      CALL SCROPN(1,2,2,SCRFNO)
*
*       Initialise arrays (easier this way than in a DATA statement)
*
      TOL = 1.0D-8
      BRACK = 0
      ITER = 0
      OLDDAT = 0
      MSY = 0.0D0
      MSP = 0.0D0
      VAR0 = 0.0D0
      FU = 0.0D0
      NDAYS = 0
      DO 90 I=1,NSITES
       ONEDAY(I,1) = -1.0D12
       ONEDAY(I,2) = -1.0D12
 90   CONTINUE
*
*   Read data & calculate the mean squares of both observations
*   and predictions. Also calculate necessary theoretical 
*   quantities, including starting points for iteration where
*   necessary (assume independence). Going to NOBS+1 guarantees we 
*   include the final day's measurements. All sums of squares
*   for model 21 are weighted according to the number of contributing
*   observations.
*
      If (Verbosity.Gt.0) then
       WRITE(MESSAGE,2)
       CALL INTPR(TRIM(MESSAGE),-1,0,0)
      End if
      DO 100 I=1,NOBS+1
       IF (I.LE.NOBS) THEN
        READ(FILNO,REC=I) (CASEID(J),J=1,4),Y,PROBS
       ENDIF
*
*	Now, if the date has changed OR if we're on the last observation, 
*	update yesterday's contributions. This scheme
*	relies on the fact that the data are in chronological order. As
*	we go, we write each day's mean prediction and number of 
*	active sites to a scratch file. NB all sums of squares are
*	weighted by number of contributing observations.
*
       NEWDAT = (10000*CASEID(2))+(100*CASEID(3))+CASEID(4)
       IF (((NEWDAT.NE.OLDDAT).AND.(I.NE.1)).OR.(I.EQ.NOBS+1)) THEN
        OLDDAT=NEWDAT
        YBAR = 0.0D0
        PBAR = 0.0D0
        NT = 0.0D0
        TMP = 0.0D0
        DO 200 J=1,NSITES
         IF (ONEDAY(J,1).GT.-0.9D12) THEN
          YBAR = YBAR + ONEDAY(J,1)
          PBAR = PBAR + ONEDAY(J,2)
          TMP = TMP + ( ONEDAY(J,2)*(1.0D0-ONEDAY(J,2)) )
          NT = NT + 1.0D0
          ONEDAY(J,1) = -1.0D12
          ONEDAY(J,2) = -1.0D12
         ENDIF
 200    CONTINUE
*
*     Some model-dependent stuff ...
*
        IF (SPMOD.EQ.21) THEN
         YBAR = YBAR/NT
         PBAR = PBAR/NT
         MSY = MSY + (NT*(YBAR**2))
         MSP = MSP + (NT*(PBAR**2))
*
*	MODEL 21: Divide by NT instead of NT**2 to allow for weighting
*
         VAR0 = VAR0 + (TMP/NT)
         NDAYS = NDAYS + 1
         WRITE(SCRFNO,REC=NDAYS) PBAR,NT
        ELSEIF (SPMOD.EQ.22) THEN
*
*     MODEL 22: not much to do.
*
         PBAR = PBAR/NT 
         TMP = PBAR * (1.0D0 - PBAR)
         MSP = MSP + NT - 1.0D0
         MSY = MSY + (((YBAR-(NT*PBAR))**2)/(NT*TMP)) - 1.0D0
        ELSE
         WRITE(MESSAGE,1)
         CALL INTPR(TRIM(MESSAGE),-1,0,0)
         CLOSE (SCRFNO)
         IFAIL = 999
         RETURN
        ENDIF
       ENDIF
*
*	And copy the latest entry into today's array
*
       ONEDAY(CASEID(1),1) = Y
       ONEDAY(CASEID(1),2) = PROBS
 100  CONTINUE

      IF (SPMOD.EQ.21) THEN
*
*	MODEL 21: choose parameters to fit to mean and variance 
*	of proportion of wet sites each day. Use the secant method to
*	initally bracket the root of the equation (Observed variance) -
*	(Theoretical variance) = 0, then regula falsi/bisection thereafter.
*	NB if user has requested RHO(1)=0 we reset it to 1.
*
       MSY = MSY/DBLE(NOBS)
       MSP = MSP/DBLE(NOBS)
       VAR0 = VAR0/DBLE(NOBS)
       LB = 0.0D0
       IF (RHO(1).LT.1.0D-12) RHO(1) = 1.0D0
       UB = RHO(1)
       FL = MSY - (VAR0 + MSP)
*
*	For each day, read mean predictions, individual predictions,
*	compute Bs and hence theoretical variances. We now use the
*	first column of ONEDAY to store the probabilities, and
*	the second column to store the Bs. NB SPVAR returns 
*       something large and negative if an error condition occurs. 
*
 299   ITER = ITER + 1
       FCUR = MSY - MSP - SPVAR(FILNO,SCRFNO,NDAYS,NSITES,RHO,SPMOD)
       IF (FCUR.GT.1.0D10) THEN 
        IFAIL = 999
        CLOSE(SCRFNO)
        RETURN
       ENDIF
*
*	That's the value of the difference in variances at
*	the current objective function value. At first we check
*	to make sure the root is actually bracketed (use the 
*	fact that we know the theoretical variance increases
*	with A so that FCUR decreases). Once the root is bracketed
*	we use regula falsi to get to the solution.
*
       IF (FCUR.GT.0.0D0) THEN
        IF (BRACK.EQ.0) THEN
         TMP = RHO(1)+( FCUR*(RHO(1)-LB)/(FL-FCUR) )
         UB = TMP
         LB = RHO(1)
         RHO(1) = UB
         FL = FCUR
        ELSE
         LB = RHO(1)
         FL = FCUR
        ENDIF
       ELSE
        BRACK = 1
        UB = RHO(1)
        FU = FCUR
       ENDIF
       IF (BRACK.EQ.1) THEN
        RHO(1) = ( (FL*UB)-(FU*LB) )/(FL-FU)
       ENDIF
       IF (DABS(FCUR).GT.TOL) GOTO 299
      ELSEIF (SPMOD.EQ.22) THEN
*
*     MODEL 22: All is very straightforward - analytical expression
*
       RHO(1) = (MSP/MSY) - 1.0D0
      ELSE
       WRITE(MESSAGE,1)
       CALL INTPR(TRIM(MESSAGE),-1,0,0)
       CLOSE (SCRFNO)
       IFAIL = 999
       RETURN
      ENDIF
      CLOSE (SCRFNO)

 1    FORMAT('****ERROR**** Invalid value of SPMOD in routine SPEST.')
 2    FORMAT('Iteration complete - now computing spatial ',
     +'dependence structure ...')

      END
******************************************************************************
******************************************************************************
******************************************************************************
      FUNCTION SPVAR(FILNO1,FILNO2,NDAYS,NSITES,RHO,SPMOD)
*
*	To evaluate the variance of daily proportion of wet sites
*	under spatial dependence models for binary data. FILNO1 
*       contains individual 
*	observations and predictions, FILNO2 contains numbers of sites 
*	for each day, and mean predictions. In fact, what is returned
*	is the expected value of 
*	[ \sum_{t=1}^{T} n_{t} * (YBAR_{t})^{2} ]/ \sum_{t=1}^{T} n_{t}
******************************************************************************
*               INTEGER variables
*               ^^^^^^^^^^^^^^^^^
*	FILNO1	- File number of random access file containing observations
*		  and predictions
*	FILNO2	- File number of random access file containing daily
*		  summary info
*	NDAYS	- Number of days in database
*	NTODAT	- Number of cases read to date
*	SPMOD	- Spatial structure chosen
*	I,J,K	- Counters
*	CASEID	- Contains date & site identifiers for each case
*       IFAIL   - Error flag
******************************************************************************
      INTEGER FILNO1,FILNO2,NDAYS,NSITES,SPMOD,CASEID(4)
      INTEGER NTODAT,I,J,K
******************************************************************************
*               DOUBLE PRECISION variables
*               ^^^^^^^^^^^^^^^^^^^^^^^^^^
*	RHO	- Array of spatial dependence parameters
*	SPVAR	- Final function value
*	A	- Parameter to estimate in dependence model 21
*	B	- Value of B corresponding to A in model 21
*	PBAR	- Probability that a day is `wet' in model 21
*	Y	- Current observation
*	PROBS	- Current prediction
*	CURVAR	- Contribution to overall variance from current day
*	NT	- Number of sites active today
*	ONEDAY	- Vector of rainfall occurrences for today
*	TMP	- temporary storage
******************************************************************************
      DOUBLE PRECISION RHO(4),SPVAR,PBAR
      DOUBLE PRECISION A,B,Y,PROBS,CURVAR,NT
      DOUBLE PRECISION ONEDAY(NSITES,2)
      DOUBLE PRECISION TMP
***************************************************************************
*     CHARACTER variables
*     ^^^^^^^^^^^^^^^^^^^
*     MESSAGE   - For screen output
***************************************************************************
      CHARACTER MESSAGE*255
******************************************************************************
*
*	Initialise:
*
      SPVAR = 0.0D0
*
*	Dependence model 21 (independence given binary weather state)
*
      IF (SPMOD.EQ.21) THEN
       NTODAT = 0
       A = DEXP(RHO(1))
       DO 300 I=1,NDAYS
        DO 302 J=1,NSITES
         ONEDAY(J,1) = -1.0D12
 302    CONTINUE
        READ(FILNO2,REC=I) PBAR,NT
        DO 305 J=1,NINT(NT)
         READ(FILNO1,REC=NTODAT+J) (CASEID(K),K=1,4),Y,PROBS
         ONEDAY(CASEID(1),1) = PROBS
         CALL BCALC(PROBS,PBAR,A,B)
         ONEDAY(CASEID(1),2) = B
 305    CONTINUE
        NTODAT = NTODAT+NINT(NT)
*
*	So now I know the Ps and the Bs, and can calculate
*	the covariances for each site pair, hence the 
*	variance. Covariances go in TMP since we don't need them
*	otherwise. Variances are calculated from marginal probabilities,
*	since covariance formula is inapplicable for i=j.
*
        CURVAR = 0.0D0
        DO 350 J=1,NSITES
         IF (ONEDAY(J,1).GT.-0.9D12) THEN
          TMP = ONEDAY(J,1)*(1.0D0-ONEDAY(J,1))
          CURVAR = CURVAR + TMP
          DO 351 K=J+1,NSITES
           IF (ONEDAY(K,1).GT.-0.9D12 ) THEN
            TMP = (1.0D0-PBAR)*ONEDAY(J,1)*ONEDAY(K,1)*
     +           (1.0D0-ONEDAY(J,1))*(1.0D0-ONEDAY(K,1))*
     +           (1.0D0-ONEDAY(J,2))*(1.0D0-ONEDAY(K,2))
            TMP = TMP/( PBAR*
     +           (1.0D0-(ONEDAY(J,1)*(1.0D0-ONEDAY(J,2))))*
     +           (1.0D0-(ONEDAY(K,1)*(1.0D0-ONEDAY(K,2)))) )
            CURVAR = CURVAR + (2.0D0*TMP)
           ENDIF
 351      CONTINUE
         ENDIF
 350    CONTINUE
*
*	Weighting by # of observations means divide by NT, not NT**2
*
        CURVAR = CURVAR/NT
        SPVAR = SPVAR + CURVAR
 300   CONTINUE
       SPVAR = SPVAR/DBLE(NTODAT)
      ELSE
       WRITE(MESSAGE,1)
       CALL INTPR(TRIM(MESSAGE),-1,0,0)
       SPVAR = -1.0D10
       CLOSE (FILNO1)
       CLOSE (FILNO2)
       RETURN
      ENDIF

 1    FORMAT('****ERROR**** Invalid value of SPMOD in routine SPVAR')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE CORLAT(FILNO,NSITES,NOBS,SPMOD,SCODES,
     +                  RHO,OUTFNO,SITINF,NATTR,MXP,Verbosity,IFAIL)
*
*       To estimate latent Gaussian inter-site correlations
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMETERS      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*       MXNOBS  - Maximum number of days for which two sites have
*                 simultaneous observations. 50000 corresponds to 
*                 roughly 150 years. 
******************************************************************************
c        INTEGER MXNOBS
c        PARAMETER (MXNOBS=50000)
******************************************************************************
*               INTEGER variables
*               ^^^^^^^^^^^^^^^^^
*   FILNO   - File number of scratch file with residual info in
*   MXP	    - Maximum number of parameters
*   NOBS    - no. of observations
*   NSITES  - no. of sites required
*   SPMOD   - spatial dependence model selected
*   OUTFNO  - File number of output file (for writing observed
*             correlation matrix)
*   NATTR   - Number of site attributes that have been defined
*   Verbosity Controls verbosity of output
***************************************************************************
      INTEGER NSITES,NOBS,FILNO,SPMOD,MXP,OUTFNO,NATTR,Verbosity
***************************************************************************
*               Additional INTEGERs
*               ^^^^^^^^^^^^^^^^^^^
*   I,J,K,L - Counters
*   CASEID  - array of identifiers for each case, for purposes of
*             residual analysis. Col.1 contains site number, 2 
*             contains year, 3 contains month and 4 day
*   NSOFAR  - Number of observations read so far for a pair of sites
*   OLDDAT  - Date of previous observation read from file
*   NEWDAT  - Date of current observation
*   ITER    - Iteration number in regula falsi search
*   IFAIL   - Error flag 
***************************************************************************
      INTEGER CASEID(4)
      INTEGER I,J,K,L,NSOFAR,OLDDAT,NEWDAT,ITER,IFAIL
***************************************************************************
*               DOUBLE PRECISION variables
*               ^^^^^^^^^^^^^^^^^^^^^^^^^^
*   SITINF  Array containing site information. See calling program
*           for details
*   RHO	    Vector of parameter estimates for spatial dependence
*   JNTPRB  Matrix of overall joint occurrence probabilities at pairs of
*           sites
*   CORMAT  Estimated spatial correlation matrix of latent Gaussian field
*   CDFBVN  Joint upper tail probabilities for a standard bivariate 
*           normal distribution
*   QNORM   Quantile function for a univariate normal distribution
*   LB, UB  Lower and upper bounds for latent correlation
*   FL, FU  Values of objective function at current lower and upper bounds
*   FCUR    Value of objective function at current estimate
*   TOL     Tolerance for judging convergence of correlation search
*   JNTN    Ns going to make up entries in JNTPRB
*   Y	    Binary rainfall occurrence indicator read from file
*   P	    Modelled probability of occurrence, again read from file
*   TAU     Array of Z-scores corresponding to marginal 
*           occurrence probabilities for each pair of sites
*   XSEP }  Separation between current site pair in the first 
*   YSEP }  two attributes contained in SITINF
******************************************************************************
      DOUBLE PRECISION SITINF(NSITES,MXP,0:3)
      DOUBLE PRECISION JNTPRB(NSITES,NSITES),CORMAT(NSITES,NSITES)
      DOUBLE PRECISION JNTN(NSITES,NSITES),Y,P
      Double precision, dimension(:,:), allocatable :: TAU
      DOUBLE PRECISION CDFBVN,QNORM,FCUR,FL,FU,LB,UB,TOL
      DOUBLE PRECISION RHO(MXP),XSEP,YSEP
***************************************************************************
*               CHARACTER variables
*               ^^^^^^^^^^^^^^^^^^^
*       SCODES  - array of short site codes
*       MESSAGE - Messages to screen
***************************************************************************
      CHARACTER SCODES(NSITES)*4,MESSAGE(2)*255
***************************************************************************
*
*	Allocate storage for TAU
*
***************************************************************************
      allocate (TAU(1:NOBS, 1:2), Stat=IFAIL)
      if (ifail.ne.0) then
       write(message(1),6)
       call intpr(trim(message(1)),-1,0,0)
       ifail = 999 
       return
      end if
*
*       Calculate matrix of observed joint probabilites
*
      If (Verbosity.Gt.0) then
       WRITE(MESSAGE,3)
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      End if
      CALL EMPCOR(FILNO,NSITES,NOBS,1,0,3,JNTPRB,JNTN,IFAIL)
      IF (IFAIL.NE.0) goto 999
*
*       Now, for each pair of sites, find the latent Gaussian 
*       correlation corresponding to the observed joint occurrence
*       probability. This is done using regula falsi / bisection. 
*       It's a bit tedious - have to read all the probabilities
*       for each pair of sites every time. 
*
      DO 100 I=1,NSITES-1
       CORMAT(I,I) = 1.0D0
       DO 101 J=I+1,NSITES
        If (Verbosity.Gt.1) then
         WRITE(MESSAGE(1),4) SCODES(I),SCODES(J)
         CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
        End if
        IF (JNTN(I,J).GT.DBLE(NOBS)) THEN
         CALL PRGERR(' CORLAT ',' fitting.f  ',2,IFAIL)
         goto 999
        ENDIF
        IF (JNTN(I,J).GT.0.0D0) THEN
*
*       Read marginal probabilities for current site pair. Going to NOBS+1
*       guarantees we include the final day's measurements
*
         NSOFAR = 0
         OLDDAT = 0
         TAU(1,1) = -1.0D12
         TAU(1,2) = -1.0D12
         DO 120 K=1,NOBS+1
          IF (K.LE.NOBS) THEN
           READ(FILNO,REC=K) (CASEID(L),L=1,4),Y,P
           NEWDAT = (10000*CASEID(2))+(100*CASEID(3))+CASEID(4)
           IF ((NEWDAT.NE.OLDDAT).OR.(I.EQ.NOBS+1)) THEN
            IF ((TAU(NSOFAR+1,1).GT.-1.0D12).AND.
     +          (TAU(NSOFAR+1,2).GT.-1.0D12)) THEN 
             NSOFAR = NSOFAR + 1
             IF (INT(JNTN(I,J)).EQ.NSOFAR) GOTO 199
            ENDIF
            OLDDAT=NEWDAT
            TAU(NSOFAR+1,1) = -1.0D12
            TAU(NSOFAR+1,2) = -1.0D12
           ENDIF
           IF (CASEID(1).EQ.I) 
     +                   TAU(NSOFAR+1,1) = -QNORM(P,0.0D0,1.0D0,IFAIL)
           IF (CASEID(1).EQ.J) 
     +                   TAU(NSOFAR+1,2) = -QNORM(P,0.0D0,1.0D0,IFAIL)
          ENDIF
 120     CONTINUE
*
*     At this point we have an array containing the thresholds corresponding
*     to marginal probabilities at each pair of sites in the first NSOFAR
*     rows. We can now find the latent Gaussian correlation corresponding 
*     to the overall probability of joint occurrence. Start by bracketing
*     the correlation in either the range (0,1) or (-1,0) depending on
*     whether the observed joint probability is larger or smaller than 
*     that expected under independence. 
*
 199     LB = -1.0D0
         UB = 1.0D0
         FU = 0.0D0
         FL = 0.0D0
         DO 200 K=1,NSOFAR
          FU = FU + CDFBVN(TAU(K,1),TAU(K,2),UB,1.0D-5)
          FL = FL + CDFBVN(TAU(K,1),TAU(K,2),LB,1.0D-5)
 200     CONTINUE
         FU = JNTPRB(I,J) - (FU / DBLE(NSOFAR))
         FL = JNTPRB(I,J) - (FL / DBLE(NSOFAR))
*
*       It is conceivably possible (and has occurred!) that the
*       observed joint probability falls outside the interval 
*       dictated by the modelled margins. In this case issue a
*       warning, use whichever extreme correlation is most 
*       appropriate, and carry on. 
*
         IF (FL*FU.GT.0.0D0) THEN
          IF (FL.GT.0.0D0) THEN
           CORMAT(I,J) = UB
          ELSE
           CORMAT(I,J) = LB
          ENDIF
          If (Verbosity.Gt.0) then
           WRITE(MESSAGE,5) CORMAT(I,J)
           CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
           CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
          End if
          GOTO 101
         ENDIF
*
*       Regula falsi starts here. NB solving for the correlations
*       is computationally expensive, and in any case there isn't 
*       much point trying to get it more accurate than the 
*       precision of the joint probability estimate itself. 
* 
         ITER = 0
         TOL = DSQRT(JNTPRB(I,J)*(1.0D0-JNTPRB(I,J))/DBLE(NSOFAR))
         TOL = DMAX1(1.0D-5,TOL / 1.0D1)
 210     ITER = ITER + 1
         CORMAT(I,J) = ( (FL*UB)-(FU*LB) )/(FL-FU)
         FCUR = 0.0D0
         DO 211 K=1,NSOFAR
          FCUR = FCUR + CDFBVN(TAU(K,1),TAU(K,2),
     +                         CORMAT(I,J),TOL/1.0D2)
 211     CONTINUE
         FCUR = JNTPRB(I,J) - (FCUR / DBLE(NSOFAR))
         IF (FL*FCUR.LE.0.0D0) THEN
          UB = CORMAT(I,J)
          FU = FCUR
         ELSEIF (FU*FCUR.LT.0.0D0) THEN
          LB = CORMAT(I,J)
          FL = FCUR
         ELSE 
          WRITE(MESSAGE(1),2)
          CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
          IFAIL = 999 
          goto 999
         ENDIF

         IF ((DABS(FCUR).GT.TOL).AND.(DABS(UB-LB).GT.TOL)) GOTO 210
        ELSE
*
*       Missing correlations - code as -9.999
*
         CORMAT(I,J) = -9.999D0
        ENDIF
 101   CONTINUE
 100  CONTINUE
      CORMAT(NSITES,NSITES) = 1.0D0
*
*	And produce estimates of correlation structure. For all
*       fitted models, empirical correlations are written to file
*       (missing values coded as -9.999) since this provides the
*       opportunity to assess the fit visually afterwards. 
*
      IF ((SPMOD.GT.0).AND.(SPMOD.LT.20)) THEN
       WRITE(OUTFNO,18)
       XSEP = 0.0D0
       YSEP = 0.0D0
       DO 500 I=1,NSITES-1
        DO 501 J=I+1,NSITES
         IF (NATTR.GE.1) XSEP=SITINF(I,1,0)-SITINF(J,1,0)
         IF (NATTR.GE.2) YSEP=SITINF(I,2,0)-SITINF(J,2,0)
         WRITE(OUTFNO,19) SCODES(I),SCODES(J),CORMAT(I,J),
     +                                    INT(JNTN(I,J)),XSEP,YSEP
 501    CONTINUE
 500   CONTINUE
       IF (SPMOD.GT.1) CALL CORFIT(SPMOD,CORMAT,JNTN,SITINF,NATTR,
     +                                   SCODES,NSITES,RHO,MXP,IFAIL)
       IF (IFAIL.NE.0) goto 999
      ELSE
       WRITE(MESSAGE(1),1)
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       IFAIL = 999 
       goto 999
      ENDIF
      
 999  deallocate (TAU)
      return

 1    FORMAT('****ERROR**** Invalid value of SPMOD in routine CORLAT')
 2    FORMAT('****ERROR**** Can''t bracket correlation in routine ',
     +'CORLAT.')
 3    FORMAT('Iteration complete - now computing spatial ',
     +'dependence structure',/,'(this could take some time) ...')
 4    FORMAT(10X,'--- Sites ',A4,' and ',A4,' ---')
 5    FORMAT('****WARNING**** Observed joint occurrence frequency ',
     +'is incompatible',/,'with modelled marginal probabilities. ',
     +'Setting latent correlation to ',F4.1)
 6    FORMAT('****ERROR**** unable to allocate storage in routine ',
     +       'CORLAT.')
 18   FORMAT('Site1  Site2  Corr      N   Xsep    Ysep ') 
 19   FORMAT(A4,4X,A4,4X,F7.4,1X,I5,2(1X,F7.3))

      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE COREST(FILNO,NSITES,NOBS,MODEL,SPMOD,SCODES,
     +                  RHO,OUTFNO,SITINF,NATTR,MXP,Verbosity,IFAIL)
*
*       To estimate residual correlation structure for normal 
*       or gamma observations: raw residuals for Gaussian obs,
*       Anscombe for gamma
******************************************************************************
*               INTEGER variables
*               ^^^^^^^^^^^^^^^^^
*   I,J	   - Counters
*   FILNO  - File number of scratch file with residual info in
*   MXP	   - Maximum number of parameters
*   NOBS   - no. of observations
*   MODEL  - code indicating which family of distributions we're using
*   NSITES - no. of sites required
*   SPMOD  - spatial dependence model selected
*   OUTFNO - File number of output file (for writing observed
*            correlation matrix)
*   NATTR  - Number of site attributes that have been defined
*   RESTYP - Type of residual for which we're calculating correlations,
*            as follows:
*	        0 - raw responses
*	        1 - Pearson residuals (=raw residuals for Gaussian
*                   observations)
*	        2 - Anscombe residuals for gamma observations
*   VerbosityControls verbosity of on-screen output
*   IFAIL  - Error flag
***************************************************************************
      INTEGER NSITES,NOBS,FILNO,MODEL,SPMOD,MXP,OUTFNO,NATTR
      INTEGER RESTYP,Verbosity,IFAIL
      INTEGER I,J
***************************************************************************
*               DOUBLE PRECISION variables
*               ^^^^^^^^^^^^^^^^^^^^^^^^^^
*   SITINF Array containing site information. See calling program
*          for details
*   RHO    Vector of parameter estimates for spatial dependence
*   CORMAT Spatial correlation matrix/matrix of cross-products
*   CORRN  Ns going to make up entries in CORMAT
*   XSEP } Separations between current site pair in the first 
*   YSEP } two attributes contained in SITINF
******************************************************************************
      DOUBLE PRECISION SITINF(NSITES,MXP,0:3)
      DOUBLE PRECISION CORMAT(NSITES,NSITES),CORRN(NSITES,NSITES)
      DOUBLE PRECISION RHO(MXP),XSEP,YSEP
***************************************************************************
*               CHARACTER variables
*               ^^^^^^^^^^^^^^^^^^^
*       SCODES  - array of short site codes
*       MESSAGE - messages to screen
***************************************************************************
      CHARACTER SCODES(NSITES)*4,MESSAGE*255
***************************************************************************
*	Check dimensioning:
*
      IFAIL = 0
*
*       Calculate estimated empirical correlation matrix
*
      If (Verbosity.Gt.0) then
       CALL INTPR('Iteration complete - now computing spatial '//
     +            'dependence structure ...',-1,0,0)
      End if
      RESTYP = 1
      IF (MODEL.EQ.11) RESTYP = 2
      CALL EMPCOR(FILNO,NSITES,NOBS,MODEL,RESTYP,2,CORMAT,CORRN,IFAIL)
      IF (IFAIL.NE.0) RETURN
*
*	And produce estimates of correlation structure. For all
*       fitted models, empirical correlations are written to file
*       (missing values coded as -9.999) since this provides the
*       opportunity to assess the fit visually afterwards. 
*
      IF ((SPMOD.GT.0).AND.(SPMOD.LT.20)) THEN
       REWIND(OUTFNO)
       WRITE(OUTFNO,18)
       XSEP = 0.0D0
       YSEP = 0.0D0
       DO 100 I=1,NSITES-1
        DO 101 J=I+1,NSITES
         IF (NATTR.GE.1) XSEP=SITINF(I,1,0)-SITINF(J,1,0)
         IF (NATTR.GE.2) YSEP=SITINF(I,2,0)-SITINF(J,2,0)
         IF (CORRN(I,J).LE.0.0D0) CORMAT(I,J) = -9.999D0
         WRITE(OUTFNO,19) SCODES(I),SCODES(J),CORMAT(I,J),
     +                                   INT(CORRN(I,J)),XSEP,YSEP
 101    CONTINUE
 100   CONTINUE
       IF (SPMOD.GT.1) CALL CORFIT(SPMOD,CORMAT,CORRN,SITINF,NATTR,
     +                                   SCODES,NSITES,RHO,MXP,IFAIL)
       IF (IFAIL.NE.0) RETURN
      ELSE
       WRITE(MESSAGE,1)
       CALL INTPR(TRIM(MESSAGE),1,0,0)
       IFAIL = 999
      ENDIF

 1    FORMAT('****ERROR**** Invalid value of SPMOD in routine COREST')
 18   FORMAT('Site1  Site2  Corr      N   Xsep    Ysep ') 
 19   FORMAT(A4,4X,A4,4X,F7.4,1X,I5,2(1X,F7.3))

      END
******************************************************************************
******************************************************************************
******************************************************************************
      Subroutine ReadData(UnitNos,MODEL,MissVal,RespIdx,AllowIncAvge,
     +                    YCHECK,THRESH,THRTYP,NSITES,MXP,MXLAG,NVARS,
     +                    P,NP,COVCODE,TWO,THREE,PWTIDX,SCODES,SITINF,
     +                    Distance,BETA,THETA,TRACE,ICHECK,NOBS,
     +                    DateRange,IFAIL)
******************************************************************************
*       To read a dataset and write a binary scratch file (which has
*       already been opened) containing response and covariate information.
*       Arguments:
*
*       UnitNos Contains channel numbers of input data files. See GLMFIT
*               header for details. The scratch file here corresponds to 
*               the unit in UnitNos(90)
*       MODEL   Code for model being fitted
*       MissVal Value corresponding to missing data in the input file. 
*       YCHECK  Code to determine whether to check the values of the
*               response variable. See GLMFIT header for details. 
*       THRESH  Threshold used for dealing with small positive values
*       THRTYP  Type of thresholding to be carried out on small
*               positive values. See GLMFIT header for more details. 
*       NSITES  Number of sites
*       MXP     Maximum number of parameters
*       MXLAG   Only use observations for which at least MXLAG previous
*               days' observations are available at the same site.
*       P       Number of parameters estimated
*       NP      Array indicating numbers of parameters in each model 
*               component. See GLMFIT header for details.
*       COVCODE Array containing code numbers of covariates.
*       TWO     Array giving indices of two-way interactions in model
*       THREE   Guess ...
*       PWTIDX  Array indexing  weighted averages of previous days' 
*               values. See GLMFIT header for details.
*       SCODES  Array of 4-character site identifiers
*       SITINF  Array containing site information. See GLMFIT header 
*               for details
*       Distance- array of inter-site distances. First slice is X-separation,
*                 second slice is Y-separation and third slice is Euclidean
*                 separation
*       BETA    Current estimate of parameter vector       } For compatibility
*	THETA	Array of parameters in nonlinear functions } with COVSET
*       TRACE   Equivalent of THRESH for 'trace value' handling
*       ICHECK  Indicator for whether model parametrisation has been checked
*       NOBS    Total number of observations read
*       DateRange Dates of first and last observations, in the format 
*               YYYYMMDD. If non-negative on entry, only observations in 
*               this window will be used for fitting.
*       IFAIL   Error flag
*         
******************************************************************************
      Integer, intent(in) :: UnitNos(100),MODEL,YCHECK,THRTYP,NSITES
      Integer, intent(in) :: MXP,NVARS,MXLAG(NVARS),P,NP(10)
      Integer, intent(in) :: RespIdx,COVCODE(MXP),PWTIDX(MXP,0:3,NVARS)
      Integer, intent(in) :: TWO(MXP,2),THREE(MXP,3),AllowIncAvge
      Double precision, intent(in) :: MissVal
      Double precision, intent(in) :: THRESH, SITINF(NSITES,MXP,0:3)
      Double precision, intent(in) :: BETA(0:MXP),THETA(MXP,3),TRACE
      Double precision, intent(in) :: Distance(3,NSITES,NSITES)
      Character, intent(in) :: SCODES(NSITES)*4
      Integer, intent(inout) :: DateRange(2),ICHECK
      Integer, intent(out) :: NOBS,IFAIL
******************************************************************************
*       Additional INTEGERs
*       -------------------
*	YY,MM,} Year, month and day (used for calculation of previous
*	DD    } days' response variables)
*       CurDate Current date, in format YYYYMMDD
*       SITE    Current site
*	MISSFL	Flags for sites which have missing covariate data
*       WARNED  Indicator for whether or not a warning has been issued
*               about strange response variable values
*       DONE    Indicator for whether we have finished reading data
*       FORCE   Indicates whether or not to force recalculation of 
*               all covariate values when reading data (this should 
*               be done for the first observation, in case there is 
*               anything left in a SAVE statement from a previous call)
*       J       Counter
******************************************************************************
      Integer YY,MM,DD,CurDate,SITE,MISSFL(NSITES),WARNED,DONE,FORCE,J
******************************************************************************
*       Additional DOUBLEs
*       ------------------
*       RSPNSE  Vector of response variables for current day
*       Y       Response for current case
*       CaseWt  Weight for current case (used for joint mean-variance
*               models)
*       X       Covariates for current case
*	COVS    Values of covariates at each site for 1 day
*       DatArray Array of current and previous days' observations: Jth column
*               corresponds to J days previously (J=0 is current day)
******************************************************************************
      Double precision DatArray(NSITES,NVARS,0:10),RSPNSE(NSITES)
      DOUBLE PRECISION Y,CaseWt,X(0:MXP),COVS(NSITES,MXP)
***************************************************************************
*       Additional CHARACTERs
*       ---------------------
*       MESSAGE - Used for on-screen messages
***************************************************************************
      Character MESSAGE*255
***************************************************************************
*
*       Initialise (NB use of Fortran 90 array capabilities)
*
      DatArray = -1.0D101
      NOBS=0
      IFAIL = 0
      WARNED = 0
      YY = 0
      MM = 0
      DD = 0
      DONE = 0
      FORCE = 1
      X(0) = 1.0D0
*
*       Ensure that all files and pointers are positioned correctly
*

      REWIND(UnitNos(1))
      Call FileReset(UnitNos(3:5))

*
*       Read through the input file, retaining only data from the 
*       requested time period if the elements of DateRange are 
*       non-negative on input 
*      
 50   CALL UPDATE(UnitNos(1),DatArray,THRESH,THRTYP,MissVal,RespIdx,
     +            NSITES,NVARS,YY,MM,DD,MXLAG,SCODES,DONE,IFAIL)
      IF (IFAIL.NE.0) RETURN
      CurDate = (10000*YY) + (100*MM) + DD
      If (DateRange(1).GE.0) then
       If (CurDate.GT.DateRange(2)) then
        Done = 1
        Goto 60
       Else If ( (DONE.EQ.0) .AND. (CurDate.LT.DateRange(1)) ) Then
        GOTO 50
       End If
      End if 

      CALL COVSET(COVS,UnitNos(3:5),NP,NSITES,NVARS,RespIdx,COVCODE,
     +     TWO,THREE,SITINF,Distance,YY,MM,DD,DatArray,AllowIncAvge,
     +     PWTIDX,MXLAG,BETA,THETA,TRACE,0,FORCE,MISSFL,ICHECK,IFAIL,
     +     MXP,0,1)

      IF (IFAIL.NE.0) RETURN
      FORCE = 0
*
*	Use this value only if all required predictors are present
*
      DO 51 SITE=1,NSITES
       RSPNSE(SITE) = DatArray(SITE,RespIdx,0)
       IF ((MISSFL(SITE).EQ.0).AND.(RSPNSE(SITE).GT.-1.0D99)) THEN
        IF ((YCHECK.EQ.0).OR.(YCHECK.EQ.2)) THEN
         IF (MODEL.EQ.1) THEN
          IF ( (DABS(RSPNSE(SITE)).GT.1.0D-6).AND.
     +         (DABS(RSPNSE(SITE)-1.0D0).GT.1.0D-6) ) THEN
           IF (YCHECK.EQ.2) GOTO 910
           IF ((YCHECK.GE.0).AND.(WARNED.EQ.0)) THEN 
            MESSAGE = '**** WARNING **** data file contains non-'//
     +                'binary values of response variable.'
            CALL INTPR(TRIM(MESSAGE),-1,0,0)
            MESSAGE = 'Maybe you used the wrong value of the '//
     +                'response.check argument to GLCfit()?'
            CALL INTPR(TRIM(MESSAGE),-1,0,0)
            WARNED = 1
           ENDIF
          ENDIF
         ELSEIF (MODEL.EQ.11) THEN
          IF (RSPNSE(SITE).LE.0.0D0) THEN
           IF (YCHECK.EQ.2) GOTO 911
           IF ((YCHECK.GE.0).AND.(WARNED.EQ.0)) THEN 
            MESSAGE = '**** WARNING **** data file contains non-'//
     +                'positive values of response variable.'
            CALL INTPR(TRIM(MESSAGE),-1,0,0)
            MESSAGE = 'Maybe you used the wrong value of the '//
     +                'response.check argument to GLCfit()?'
            CALL INTPR(TRIM(MESSAGE),-1,0,0)
            WARNED = 1
           ENDIF
          ENDIF
         ENDIF
        ENDIF
        Y = RSPNSE(SITE)
        IF ((YCHECK.EQ.1).AND.(MODEL.EQ.1)) THEN
         IF (Y.GT.0.0D0) THEN
          Y = 1.0D0
         ELSE
          Y = 0.0D0
         ENDIF 
        ELSEIF ((YCHECK.EQ.1).AND.(MODEL.EQ.11)) THEN
         IF (Y.LE.0.0D0) GOTO 51
        ENDIF         
        NOBS = NOBS + 1
        DO 52 J=1,P
         X(J) = COVS(SITE,J)
 52     CONTINUE
        CaseWt = 1.0D0
*
*	Write this case to first scratch file, and flag to indicate
*       no further parameter checking is required
*
        WRITE(UnitNos(90),REC=NOBS,ERR=920) SITE,YY,MM,DD,
     +                                      Y,CaseWt,(X(J),J=1,P)
        IF (ICHECK.EQ.0) ICHECK = 1
       ENDIF 
 51   CONTINUE
      IF (DONE.EQ.0) GOTO 50
      
*
*       Reset all files and pointers, and DateRange
*
 60   Call FileReset(UnitNos(3:5))
      If (NOBS.GT.0) then 
       READ(UnitNos(90),REC=1) SITE,YY,MM,DD
       DateRange(1) = (10000*YY) + (100*MM) + DD
       READ(UnitNos(90),REC=NOBS) SITE,YY,MM,DD
       DateRange(2) = (10000*YY) + (100*MM) + DD
      End If
      
      RETURN

 910  IFAIL = 90
      RETURN
 911  IFAIL = 91
      RETURN
 920  CALL INTPR('**** ERROR: while writing to scratch file. '//
     +            'You are probably out of disk space.',-1,0,0)
      IFAIL = 1
      RETURN

      End Subroutine ReadData
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE UPDATE(FILNO,DatArray,THRESH,THRTYP,MissVal,RespIdx,
     +                  NSITES,NVARS,YY,MM,DD,MXLAG,SCODES,DONE,IFAIL)

*
*       This updates the values held in the DatArray array. If 
*       thresholding of the response variable has been requested 
*       (THRTYP > 0) it's done here.
*
      INTEGER FILNO,YY,MM,DD,MXLAG(NVARS),NSITES,NVARS,DONE
      Integer THRTYP,RespIdx,IFAIL
      Double precision DatArray(NSITES,NVARS,0:10),THRESH,MissVal
      CHARACTER SCODES(NSITES)*4
******************************************************************************
*       Additional INTEGERS
*       ^^^^^^^^^^^^^^^^^^^
*       DAYS_IN_MONTH   - needed to calculate DP
******************************************************************************
      INTEGER DAYS_IN_MONTH
      INTEGER I,J,K
      
      IFAIL = 0
*
*       First update the DatArray values. 
*
      Do I=1,NSITES
       Do J=1,NVARS 
        Do K=MXLAG(J),1,-1
         DatArray(I,J,K) = DatArray(I,J,K-1)
        End Do
        DatArray(I,J,0) = -1.0D101
       End Do
      End Do
*
*       Increment date unless this is the first time through, in which
*       case we don't know what the first date is going to be
*
      IF (YY.NE.0) THEN
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
*
*       Read data for today
*    
      CALL DATRD(FILNO,DatArray(1:NSITES,1:NVARS,0),MissVal,
     +                  SCODES,NSITES,NVARS,YY,MM,DD,DONE,IFAIL)
      IF (IFAIL.NE.0) RETURN
*
*       Threshold if requested
*
      IF (THRTYP.GT.0) THEN
       DO 20 I = 1,NSITES
        IF (DatArray(I,RespIdx,0).LT.0.0D0) GOTO 20
        IF (DatArray(I,RespIdx,0).LT.THRESH) THEN
         DatArray(I,RespIdx,0) = 0.0D0
        ELSEIF (THRTYP.EQ.1) THEN
         DatArray(I,RespIdx,0) = DatArray(I,RespIdx,0) - THRESH
        ENDIF
 20    CONTINUE
      ENDIF
    
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE IWLS(MODEL,RespIdx,AllowIncAvge,NVARS,MissVal,NP,P,
     +                COVCODE,TWO,THREE,FOUIDX,LEGIDX,PWTIDX,MXP,
     +                MAXIT,MXLAG,NOBS,NSITES,NATTR,SITXFM,GLBCOD,
     +                BETA,THETA,CovNaive,Covrobust,PHI,LOGL,DVNCE,
     +                GLBVAL,SITINF,Distance,ATTRTXT,SCODES,UnitNos,
     +                Verbosity,IFAIL)
******************************************************************************
*	To fit a GLM, possibly involving nonlinear transformations of
*	covariates, using Iterative Weighted Least Squares. Parameter
*	vectors, standard errors and a log-likelihood are returned.
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*////////////	NB THIS ROUTINE CONTAINS QUITE A FEW SECTIONS WHICH  ///////
*////////////	ARE MODEL-SPECIFIC. THESE SECTIONS ARE ENCLOSED BY   ///////
*////////////   LINES OF CHEVRONS LIKE THIS, FOR EASY RECOGNITION IF ///////
*////////////   YOU SUBSEQUENTLY WANT TO AMEND THEM OR ADD EXTRA     ///////
*////////////   MODELS.                                              ///////
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

******************************************************************************
*	ARGUMENTS
*	~~~~~~~~~
*	MODEL	- INTEGER, input. Defines the distribution being fitted.
*		  Values:
*			1  -	Logistic regression for binary responses
*			11 -	Gamma model with log link
*       RespIdx - INTEGER, input. Index of response variable in input file.
*       NVARS   - INTEGER, input. Number of variables in input file.
*       MissVal - Double precision, input. Value indicating a missing 
*                 observation in input file. 
*	NP	- INTEGER, input. Table indicating numbers of parameters
*		  in each category of predictor.
*	P	- INTEGER, input. Number of parameters being estimated.
*	COVCODE	- INTEGER, input. Table defining covariates.
*	TWO 	- INTEGER, input. Table defining 2-way interactions.
*	THREE 	- INTEGER, input. Table defining 3-way interactions.
*	FOUIDX	- INTEGER, input. Table pointing to boundaries of
*		  orthogonal series representation of site effects
*		  (Fourier representations).
*	LEGIDX	- INTEGER, input. Table pointing to boundaries of
*		  orthogonal series representation of site effects
*		  (Legendre polynomial representations).
*       PWTIDX  - INTEGER, input. Table pointing to locations of 
*                 parameters/covariates relating to weighted averages of
*                 previous days' values.
*	MXP	- INTEGER, input. Max no of parameters in model.
*	MAXIT	- INTEGER, input. Maximum number of iterations required.
*	MXLAG	- INTEGER, input. Array of flags to indicate how many previous
*		  days' values must be available in order for a 
*		  case to be included.
*	NOBS	- INTEGER, input. Number of observations.
*	NSITES	- INTEGER, input. Number of sites used.
*	NATTR	- INTEGER, input. Number of site attributes.
*       SITXFM  - INTEGER, input. Selects nonlinear transformations 
*		  of site attributes.
*	GLBCOD  - INTEGER, input. Codes defining global quantities.
*	LOGL	- DOUBLE, output. Final log-likelihood.
*       DVNCE   - Double, output. Final deviance
*	BETA	- DOUBLE, input/output. Parameter estimates.
*	THETA	- DOUBLE, input/output. Estimates of parameters in
*		  nonlinear transformations of predictors.
*       CovNaive- DOUBLE, output. Naive covariance matrix of parameter estimates
*       CovRobust- DOUBLE, output. Robust covariance matrix of parameter
*		  estimates                 
*	PHI	- DOUBLE, output. Dispersion parameter, where needed.
*	GLBVAL  - DOUBLE, input. Values of global quantities.
*	SITINF- DOUBLE, input. Site attributes.
*       Distance- array of inter-site distances. 
*       ATTRTXT - CHARACTER, input. Text for site attributes
*       SCODES  - CHARACTER, input. 4-character site codes
*       UnitNos - INTEGER, input. Numbers of connections for input data
*                 files. See GLMFIT header for details. 
*       Verbosity Controls verbosity of output. Higher values => more verbose
*       IFAIL   - INTEGER, output. For error trapping
******************************************************************************
      INTEGER, intent(in) :: MXP,P,NOBS,NSITES,MODEL,MAXIT,UnitNos(100)
      INTEGER, intent(in) :: RespIdx,NVARS,AllowIncAvge,Verbosity
      INTEGER, intent(in) :: NP(10),COVCODE(MXP),MXLAG(NVARS)
      INTEGER TWO(MXP,2),THREE(MXP,3),FOUIDX(MXP),LEGIDX(MXP)
      INTEGER PWTIDX(MXP,0:3,NVARS),NATTR,SITXFM(MXP,2),GLBCOD(MXP)
      INTEGER, intent(out) ::  IFAIL
      Double precision, intent(in) :: MissVal,Distance(3,NSITES,NSITES)
      DOUBLE PRECISION BETA(0:MXP),THETA(MXP,3)
      DOUBLE PRECISION, intent(out) :: CovNaive(0:MXP,0:MXP)
      DOUBLE PRECISION, intent(out) :: CovRobust(0:MXP,0:MXP)
      DOUBLE PRECISION, intent(out) :: LOGL, DVNCE
      DOUBLE PRECISION GLBVAL(MXP),PHI
      DOUBLE PRECISION SITINF(NSITES,MXP,0:3)
      CHARACTER*70 ATTRTXT(MXP)
      CHARACTER SCODES(NSITES)*4
******************************************************************************
*	Extra INTEGERs - those marked # are just dummies that are required
*	~~~~~~~~~~~~~~	 by other routines such as COVSET.
*	ITER	- current iteration number
*	SITE	- number of current site
*	YY	- Year for current case
*	MM	- Month for current case
*	DD	- Day for current case
*       OY,OM,OD- Year, month and day for previous case
*	CASEID	- array of identifiers for each case.
*	RECNO	- Record number in file scratch file - used when updating 
*		  covariate attributes
* 	MISSFL	- Flags for sites which have missing covariate data
*       CNVRGE  - convergence status: 1=initial value, 2=iterating,
*                 3=converged
*       IDMAXU  - Parameter number with largest score
*       RECALC  - indicator for whether or not we have to calculate nonlinear
*                 functions at the beginning of an iteration. Codes:
*                 0 Calculate all covariates
*                 1 Calculate only covariates associated with nonlinear 
*                   parameters
*                 2 We've converged!
*       RESWRT  - indicator for whether residual information has been 
*                 written to scratch file. 
*	XRC	- indicator for whether we have to recompute X'WX (i.e. if
*		  either X matrix or weights have changed).
*       PRVRST  - Indicator for whether we need to reset the DatArray
*                 array while calculating nonlinear predictors
*	OVRSHT	- Used to keep track of number of times we correct for
*		  `overshooting' in an IWLS iteration. If this happens
*		  so that the log-likelihood decreases substantially,
*		  we backtrack up to 10 times. Rounding errors and changes
*		  to data values in EM part of algorithm may mean that we
*		  never get exactly back to where we were.
*	TRWRN	- Indicates if we've issued a warning message about trace 
*		  values on this iteration
*	FORCE	- Used to force COVSET to recalculate seasonal and yearly
*		  contributions for the first case in the dataset on each
*		  iteration - to be on the safe side.
*       THRTYP  - Type of thresholding to be used
*       Y0,M0,D0- Year, month and day read from original daily data file
*                 when recalculating the DatArray array
*       DONE    - Argument to UPDATE
*       UFILNO  - Number of scratch file for daily scores.
*       NDAYS   - Number of days of data (= number of rows in UFILNO)
*       P1      - Index of last parameter estimated - equal to P for
*                 models without a dispersion parameter, P+1 otherwise.
*	I,J,K	- Counters
******************************************************************************
      INTEGER ITER,SITE,YY,MM,DD,OY,OM,OD,CASEID(4),RECNO
      INTEGER MISSFL(NSITES),CNVRGE,RECALC,RESWRT,XRC,PRVRST,DONE
      INTEGER IDMAXU,OVRSHT,TRWRN,FORCE,THRTYP,Y0,M0,D0,UFILNO,NDAYS
      INTEGER P1,I,J,K
******************************************************************************
*	Extra DOUBLEs - those marked as # are just dummies which are
*	~~~~~~~~~~~~~   required by other routines such as COVSET
*	OLDLOGL	- Log-likelihood from previous iteration
*       LLINVL  - `Invalid' log-likelihood
*       LLINC   - Increase in Log L on last iteration
*       LLTOL   - Value of new logL which is deemed significantly
*                 worse than the old one
*	TOL	- Small number used for determining convergence
*	XWX,XWZ	- The arrays X'WX and X'Wz in IWLS (see McCullagh &
*                 Nelder, 1989, for example)
*       XWSWX   - The filling the the sandwich, for the sandwich estimator
*                 of variance in the presence of spatial dependence.
*       SCORE   - Score vector
*       SUMSCR  - Sum of scores for current day
*       SCSTD   - Standard deviations of scores
*       MAXU    - Largest component of score vector
*	BETADJ	- Adjustment to BETA from this iteration
*       BETA0   - Starting value of BETA
*       LLR0    - Log likelihood ratio corresponding to BETA0
*       LAMBDA  - Multiple of step size to use
*       LLTAB   - Table of log-likelihoods for different values of LAMBDA
*       MU	- Fitted value for current case
*       SUMSQ	- Used to calculate dispersion parameters in some models.
*       ETA	- Linear predictor for the current case
*       Z	- Adjusted dependent variate for the current case
*       U       - Contribution to score vector from current case
*       W	- Weight for the current case
*       WU,WZ   - W*U and W*Z - to save on repeated calculations
*       Y	- Observed value for the current case
*       CaseWt  - Externally-defined weight for the current case (by 
*                 contrast with W, which combines this with the weight
*                 from the IWLS algorithm)
*       X	- Covariates for the current case
*       RSPNSE	- Vector of responses for today
*#	DatArray- Observations for current and previous days (returned 
*#                from Covset)
*       WtVec   - Vector of case weights for today
*	COVS	- Covariate values for current case.
*       SUMLNY  - Sum of Ln(Y)s - used in calculating full likelihood.
*       THRESH  - Threshold for dealing with small positive values
*	TRACE   - Ditto, treating them as `trace' values
*       MLPHI   - Maximum likelihood estimate of dispersion parameter
*       DIGAMM  - Digamma function
*       TRIGAMM - Trigamma function
*	TMP	- Temporary storage
******************************************************************************
      DOUBLE PRECISION OLDLOGL,TOL,LLTOL,LLINVL,LLINC
      DOUBLE PRECISION XWX(0:MXP,0:MXP),XWZ(0:MXP)
      DOUBLE PRECISION XWSWX(0:MXP,0:MXP),SCORE(0:MXP),MAXU
      DOUBLE PRECISION MU,SUMSQ,ETA,Z,W,U,WU,WZ,SUMSCR(0:MXP)
      DOUBLE PRECISION SCSTD(0:MXP),Y,CaseWt,X(0:MXP),LLR0
      DOUBLE PRECISION BETA0(0:MXP),BETADJ(0:MXP),LAMBDA(3),LLTAB(3)
      DOUBLE PRECISION RSPNSE(NSITES),DatArray(NSITES,NVARS,0:10)
      DOUBLE PRECISION WtVec(NSITES),THRESH,TRACE,COVS(NSITES,MXP)
      DOUBLE PRECISION SUMLNY,MLPHI,DIGAMM,TRIGAMM,TMP
******************************************************************************
*	Extra CHARACTERs
*	~~~~~~~~~~~~~~~~
*       MESSAGE - Used for writing messages to screen
******************************************************************************
      CHARACTER MESSAGE(4)*200
******************************************************************************
*
*	Initialise ... (NB set X(0) to 1 as it will always represent a
*	constant term in a model).
*
      IFAIL = 0
      TOL = 1.0D-4
      LLINVL = -1.0D32
      CovNaive = 0.0D0
      CovRobust = 0.0D0
      IF (NOBS.GT.100000) TOL=TOL*DBLE(NOBS)/1.0D5
      ITER = 0
      LOGL = LLINVL
      LLTOL = LOGL - 1.0D6
      LLINC = -LOGL
      OLDLOGL = -(2.0D0*DABS(LOGL))
      X(0) = 1.0D0
      U = 0.0D0
      W = 0.0D0
      SUMSQ = 0.0D0
      CNVRGE = 2
      OVRSHT = 0
      TRWRN = 0
      PHI = 1.0D0
      P1 = P
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF (MODEL.EQ.1) PHI = -1.0D0
      IF ((MODEL.EQ.10).OR.(MODEL.EQ.11)) P1 = P+1
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      MLPHI = PHI
      CALL TAUDEF(GLBVAL,GLBCOD,NP(8),THRESH,TRACE,THRTYP,MODEL,MXP)
*
*     For robust covariance matrix calculations, we'll need to write 
*     daily sums of scores to a scratch file as we go along. Start 
*     looking for a free file handle at 85.
*
      CALL SCROPN(2,2,0,UFILNO)
      
      NDAYS = 0
*
*     Store the initial value of BETA, for e.g. Wald tests later
*     on (Fortran 90 syntax)
*
      BETA0 = BETA
*
*     We'll need to recalculate the DatArray array if there are
*     any nonlinearities relating to previous days' responses
*
      PRVRST = 0
      DO 30 J=NP(6)+1,NP(7)
       I = COVCODE(J)/1000
       IF ((I.GT.NP(3)).AND.(I.LE.NP(4))) THEN
        PRVRST = 1
        GOTO 31
       ENDIF
 30   CONTINUE
*
*	By the time we get here, all covariates are set so we only
*	need to recalculate things if there are nonlinear effects
*
 31   RECALC = 1
*
*	... and start iterating. We only update the
*	old log-likelihood when it was produced by an IWLS step
*	rather than a `corrected overshoot'. We also set
*       a threshold for judging whether we're going to 
*       overshoot.
*      
 100  IF (OVRSHT.EQ.0) THEN
       OLDLOGL = LOGL
       LLTOL = OLDLOGL
      ENDIF
      
      LOGL = 0.0D0
      SUMSQ = 0.0D0
      NDAYS = 0
      REWIND(UFILNO)
      READ(UnitNos(90),REC=1) SITE,OY,OM,OD
      IF ((ITER.EQ.0).OR.(TRACE.GT.0.0D0)) SUMLNY = 0.0D0
*
*	Initialise X'WX and X'Wz for this iteration. For some models
*	we can speed things up by only calculating X'WX once. 
*	Exceptions are when the weights depend nontrivially on the
*	fitted values (eg. in logistic regression, and when there
*	are nonlinear parameters to estimate so the design matrix 
*	itself changes from iteration to iteration. The XRC flag 
*	indicates which case we're dealing with.
*
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF ( (ITER.EQ.0).OR.(NP(7).GT.NP(6)).OR.(MODEL.EQ.1) ) THEN
       XRC = 1
      ELSE
       XRC = 0
      ENDIF
*
*       For gamma model, calculate terms involving PHI in the 
*       score for the dispersion parameter
*
      IF (MODEL.EQ.11) THEN
       TMP = 1.0D0 - DLOG(PHI) - DIGAMM(1.0D0/PHI,IFAIL)
      ENDIF

*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      MAXU = 0.0D0
      IDMAXU = 0
      XWZ = 0.0D0
      DO 105 I=0,P
       XWZ(I) = 0.0D0
       SCORE(I) = 0.0D0
       SUMSCR(I) = 0.0D0
       IF (XRC.EQ.1) THEN
        DO 106 J=0,P
         XWX(I,J) = 0.0D0
 106   CONTINUE
       ENDIF
 105  CONTINUE
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF ((MODEL.EQ.10).OR.(MODEL.EQ.11)) SUMSCR(P1) = 0.0D0
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*
*     Decide whether we're going to write residual stuff to scratch 
*     file
*
      IF ((DABS(LLINC).LT.1.0D2*TOL).OR.(RECALC.EQ.2).OR.
     +    ((MAXIT.GE.0).AND.(ITER.GE.MAXIT))) THEN
       RESWRT = 1
      ELSE
       RESWRT = 0
      ENDIF
*
*	Now start on likelihood etc. corresponding to current
*	parameter vector.
*
      TRWRN = 0
      DO 150 I = 1,NOBS
*
*       Compute linear predictor, fitted values etc. 
*       
       CALL FITCAL(BETA,NP,P,MXP,MODEL,UnitNos(90),I,ITER,TRWRN,TRACE,
     +                          MLPHI,LLINVL,X,Y,CaseWt,SITE,YY,MM,DD,
     +                          LOGL,MU,ETA,IFAIL)
     
       IF (IFAIL.EQ.2) GOTO 170
       IF (IFAIL.NE.0) RETURN
*
*	Now finish off contribution to adjusted dependent variate
*	and scores
*
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       IF (MODEL.EQ.1) THEN
*
*	Logistic ... 
*
        W = CaseWt * MU * (1.0D0 - MU)
        U = (Y-MU) / W
*
*       Normal
*
       ELSEIF (MODEL.EQ.10) THEN
        W = CaseWt
        U = (Y - MU)
        SUMSQ = SUMSQ + (CaseWt*U*U)
       ELSEIF (MODEL.EQ.11) THEN
*
*	Gamma. 
*
        W = CaseWt
        U = (Y - MU) / MU
        SUMSQ = SUMSQ + (CaseWt*U*U)
*
*     Need sum of log Y's for estimation of dispersion parameter
*     and calculation of deviance - the Ys may have changed if there 
*     are trace values present.
*
        IF ((ITER.EQ.0).OR.(TRACE.GT.0.0D0)) 
     +                            SUMLNY = SUMLNY + (CaseWt*DLOG(Y))
       ENDIF
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*
*       There's a bit to add on to Z if we're estimating nonlinear parameters
*       (ETA is the linear predictor, but the vector of first derivatives
*       has contributions from the derivatives of the nonlinear params
*       as well).
*
       Z = U + ETA
       IF (NP(7).GT.NP(6)) THEN
        DO 160 J=NP(6)+1,NP(7)
         Z = Z + (BETA(J)*X(J))
 160    CONTINUE
       ENDIF
*
*	Calculate this case's contribution to X'Wz, and to
*	X'WX if this is necessary. (i.e. if the weights have
*	changed on this iteration - which happens for the logistic
*	model and if there are nonlinear transformations of 
*	predictors). Also update score vector. Note that the case 
*       weights have been incorporated into W.
* 
       WZ = W*Z
       WU = W*U
       DO 165 J=0,P
        XWZ(J) = XWZ(J) + (X(J)*WZ)
        SCORE(J) = SCORE(J) + (X(J)*WU)
        IF (XRC.EQ.1) THEN
         DO 166 K=J,P
          XWX(J,K) = XWX(J,K) + ( X(J)*W*X(K) )
 166     CONTINUE
        ENDIF
 165   CONTINUE 
*
*	Write prediction stuff to scratch file for use later on. We only 
*       do this if it looks like we're about to converge (or if we've 
*       already converged but this hasn't been done yet), so as not to 
*       waste a load of time with disk access.
*
       IF (RESWRT.EQ.1) THEN
        WRITE(UnitNos(91),REC=I) SITE,YY,MM,DD,Y,MU,ETA,CaseWt
*
*     Write daily sums of scores to scratch file. Scores for mean components 
*     need to be divided by PHI (NB PHI is stored as a negative number for 
*     models where its value is known). For models with a dispersion 
*     parameter, the relevant scaling is already included in the calculations. 
* 
        IF ((DD.NE.OD).OR.(MM.NE.OM).OR.(YY.NE.OY)) THEN
         DO J=0,P
          SUMSCR(J) = SUMSCR(J) / DABS(PHI)
         END DO         
         WRITE(UFILNO) (SUMSCR(J),J=0,P1)
*
*	Now reset all objects so as to start computations for the new day
*
         OD = DD
         OM = MM
         OY = YY
         NDAYS = NDAYS + 1
         DO 157 J=0,P
          SUMSCR(J) = WU*X(J)
 157     CONTINUE
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*
*       Append score contributions for dispersion parameters, for models
*       that have them. NB these score contributions themselves depend on
*       the dispersion parameter value - use the estimates from the  
*       previous iteration. 
*
         IF (MODEL.EQ.10) THEN
          SUMSCR(P+1) = 
     +            CaseWt*( ((Y-MU)*(Y-MU) / MLPHI) - 1.0D0) / MLPHI
         ELSE IF (MODEL.EQ.11) THEN
          SUMSCR(P+1) = 
     +          -CaseWt*((DLOG(Y/MU)-(Y/MU)) + TMP) / (MLPHI**2)
         ENDIF         
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        ELSE
         DO 158 J=0,P
          SUMSCR(J) = SUMSCR(J) + (WU*X(J))
 158     CONTINUE
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         IF (MODEL.EQ.10) THEN
          SUMSCR(P+1) = SUMSCR(P+1) + 
     +          (CaseWt*( ((Y-MU)*(Y-MU) / MLPHI) - 1.0D0) / MLPHI)
         ELSE IF (MODEL.EQ.11) THEN
          SUMSCR(P+1) = SUMSCR(P+1) -
     +          (CaseWt*((DLOG(Y/MU)-(Y/MU)) + TMP) / (MLPHI**2))
         ENDIF
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         IF (I.EQ.NOBS) THEN
          DO J=0,P
           SUMSCR(J) = SUMSCR(J) / DABS(PHI)
          END DO
          WRITE(UFILNO) (SUMSCR(J),J=0,P1)
          NDAYS = NDAYS + 1
         ENDIF
        ENDIF
       ENDIF

 150  CONTINUE
*
*     Now calculate deviance and finish off likelihood calculations
*     for models that need it
*
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF (MODEL.EQ.1) THEN
*
*     For logistic with cell sizes all equal to 1, log-likelihood
*     for saturated model is zero so deviance is just -2 log L
* 
       DVNCE = -(2.0D0*LOGL)
      ELSEIF (MODEL.EQ.10) THEN
*
*     For Gaussian model, calculate both ML and moment estimates 
*     of dispersion parameter (the first for likelihood calculations
*     in the next iteration, and the second for subsequent use). Messy 
*     but necessary. NB the calculation of LOGL in FITCAL already 
*     includes the effect of dispersion parameter (necessary to allow
*     for normal-heteroscedastic models with non-constant variances),
*     so no need to include it here
*
       DVNCE = SUMSQ
       PHI = SUMSQ / DBLE(NOBS - P - 1)
       MLPHI = SUMSQ / DBLE(NOBS)
       LOGL = LOGL - 
     +        (DBLE(NOBS) * DLOG(8.0D0*DATAN(1.0D0)) / 2.0D0)
      ELSEIF (MODEL.EQ.11) THEN
*
*     Similarly for gamma model
*
       DVNCE = -2.0D0*(LOGL + SUMLNY + DBLE(NOBS))
       PHI = SUMSQ / DBLE(NOBS - P - 1)
       IF(.NOT.(PHI.GT.-1.0D11)) PHI = 0.0D0
       MLPHI = PHI
       CALL MLNUES(MLPHI,SUMLNY,NOBS,LOGL,IFAIL)
      ENDIF
      IF (NP(9).GT.0) BETA(P+1) = MLPHI
      IF (ITER.EQ.0) LLR0 = LOGL
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

*
*       Set RECALC flag to 2 if we've converged and, if we haven't yet
*       written the correct residual information to scratch file, go
*       back and do so (without changing any parameter values).
*
      LLINC = LOGL-OLDLOGL
      IF ( ((DABS(LLINC).LE.TOL).AND.(OVRSHT.EQ.0)).OR.
     +     ((ITER.EQ.MAXIT).AND.(LOGL.GE.LLTOL)).OR.
     +     ((ITER.GT.MAXIT).AND.(MAXIT.GE.0)) ) THEN
       RECALC = 2
       IF (RESWRT.EQ.0) GOTO 100
      ENDIF
*
*	Now finish off X'WX and XWz. The aim is to end up with the Cholesky
*	decomposition of X'WX in XWX. If we didn't have to recalculate
*	XWX on this iteration, it will already contain the decomposition
*	from last time around, so nothing needs doing. Otherwise, just
*       the upper triangle has been calculated - fill in the rest here. 
*
      IF (XRC.EQ.1) THEN
       DO 183 I=0,P
        SCSTD(I) = DSQRT(XWX(I,I))
        DO 184 J=0,I-1
         XWX(I,J) = XWX(J,I)
 184     CONTINUE
 183   CONTINUE
*
*	Fortran 90 syntax in next line - this is the inverse of the naive
*	covariance matrix at the current iteration
*
       CovNaive(0:P,0:P) = XWX(0:P,0:P)
       CALL CHOLDEC(XWX,P+1,MXP+1,IFAIL)
       IF (IFAIL.NE.0) THEN
        WRITE(MESSAGE,10) IFAIL
        CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(3)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(4)),-1,0,0)
        IFAIL = 101
        RETURN
       ENDIF
      ENDIF
*
*     Convergence check: calculate largest standardised element 
*     of score vector
*
      DO 187 I=0,P
       TMP=SCORE(I)/SCSTD(I)
       IF (DABS(TMP).GT.DABS(MAXU)) THEN
        MAXU = TMP
        IDMAXU = I+1
       ENDIF
 187  CONTINUE
*
*	Output summary for this iteration, and dispersion parameter
*	if necessary. We check for NaNs etc. in log-likelihoods, and set
*	them to something huge and negative.
*
 170  IF (.NOT.(LOGL.GT.LLINVL)) THEN
       LOGL = LLINVL
       DVNCE = -LLINVL
      ENDIF
      IFAIL = 0
      If (Verbosity.Gt.1) then
       WRITE(MESSAGE(1),2) ITER,LOGL,MAXU,IDMAXU
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      End if
*
*       Now another iteration if we haven't converged. If the last
*	iteration didn't take the log-likelihood below LLTOL, we proceed
*	in the usual way, by solving (X'WX)b = X'Wz for b. In this 
*	case we calculate X'Wz, then throw it, with the Cholesky factor 
*	of X'W'X, at a standard linear equation solver. Otherwise, we 
*       backtrack till we find an increase in log-likelihood (may slow
*	things down, but it guarantees convergence when looking for
*	nonlinear transformations).
*
      IF (RECALC.NE.2) THEN
       ITER = ITER + 1
       IF ((ITER.EQ.MAXIT).AND.(Verbosity.Gt.1)) THEN
        CALL INTPR('Recalculating best estimates found ...',-1,0,0)
       ENDIF
       IF ( (LOGL.GT.LLTOL).OR.(ITER.LE.1).OR.
     +                         (OVRSHT.EQ.10) ) THEN
        OVRSHT = 0
        LAMBDA(1) = 0.0D0
        LAMBDA(2) = 1.0D0
        LAMBDA(3) = 1.0D0
        LLTAB(1) = LOGL
*
*       Equation solver. Store the adjustments that have been made.
*
        DO 200 I=0,P
         BETADJ(I) = BETA(I)
 200    CONTINUE
        CALL CDSOLV(XWX,BETA,XWZ,P+1,0,MXP+1,IFAIL)
        IF (IFAIL.NE.0) THEN
         IFAIL = 101
         RETURN
        ENDIF
        DO 205 I=0,P
         BETADJ(I) = BETA(I) - BETADJ(I)
 205    CONTINUE
       ELSE
*
*	This is the case when we overshot, so backtrack a bit (except
*       if we're resetting results of a final iteration, in which
*       case we backtrack the whole way). First we reset everything to 
*       original values, then compute reduced step.
*
        OVRSHT = OVRSHT + 1
        IF (OVRSHT.EQ.1) THEN
         LLTAB(2) = LOGL
        ELSE
         LLTAB(3) = LOGL
        ENDIF
        DO 210 I=0,P
         BETA(I) = BETA(I) - BETADJ(I)
         BETADJ(I) = BETADJ(I)/LAMBDA(3)
 210    CONTINUE
        IF ((ITER.LE.MAXIT).OR.(MAXIT.EQ.-1)) THEN
         IF (OVRSHT.EQ.1) THEN
          LAMBDA(3) = 0.5D0
         ELSE
*
*     Here's the backtracking. All being well, we use inverse
*     parabolic interpolation, but constrain ourselves to the
*     line segment within which we know a maximum must lie.
* 
          IF (LOGL.GT.LLINVL) THEN
           CALL INVPAR(LAMBDA,LLTAB,1)
           IF (LAMBDA(3).LT.0.0D0) THEN
            LAMBDA(3) = LAMBDA(2)/1.0D1
           ELSEIF (LAMBDA(3).GT.LAMBDA(2)) THEN
            LAMBDA(3) = LAMBDA(2)/2.0D0
           ENDIF
          ELSE
           LAMBDA(3) = LAMBDA(3)/1.0D2
          ENDIF
         ENDIF
         If (Verbosity.Gt.1) then
          WRITE(MESSAGE(1),
     +        '(''Reducing step size to '',F7.3,''% of original'')')
     +        1.0D2*LAMBDA(3)
          CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
         End if
        ELSE
         LAMBDA(3) = 0.0D0
        ENDIF
        DO 211 I=0,P
         BETADJ(I) = LAMBDA(3)*BETADJ(I)
         BETA(I) = BETA(I) + BETADJ(I)
 211    CONTINUE
       ENDIF

*
*       Update nonlinear parameters if there are any, and if we're
*	iterating (doesn't work if MAXIT=0)
*
       IF ( (NP(7).GT.NP(6)).AND.(MAXIT.NE.0) ) THEN
        CALL NLUPDT(THETA,BETA,NP,COVCODE,MXP,IFAIL)
        IF (IFAIL.NE.0) RETURN
*
*       Also the matrix X. We use the routine COVSET again to recalculate X; 
*       this requires that we set some elements of COVS to their current 
*       values, as RECALC is nonzero and only those elements that have
*       changed will be reset in COVSET - if we don't pass the correct 
*       values across, the interactions will be buggered on exit. Another 
*       tricky feature is that COVSET updates a whole days' responses at 
*       once. This is a pain in the bum for this program but very useful 
*       when the routine is used for simulation. Hence the inelegancies 
*       to do with tracking dates etc! We do it by reading in elements of 
*       the RSPNSE and COVS arrays, modifying as appropriate, then rewriting 
*       them where RSPNSE is non-missing (i.e. skip over sites with missing 
*       data for the current day). As the site attributes may also have 
*       changed, we recalculate them as well (may not be necessary but it 
*       doesn't take any time).
*
        CALL ATTRXFM(SITINF,NSITES,ATTRTXT,NATTR,COVCODE,
     +                   SITXFM,FOUIDX,LEGIDX,THETA,NP,1,IFAIL,MXP)
        IF (IFAIL.NE.0) RETURN
        DO I=1,NSITES
         MISSFL(I) = 1
         RSPNSE(I) = -1.0D101
         Do J=1,NVARS
          IF (PRVRST.EQ.1) THEN
           DO K=0,MXLAG(J)
            DatArray(I,J,K) = -1.0D101
           End Do
          ENDIF
         End Do
        End Do
*
*     If we have to recalculate the DatArray array, we'll have to 
*     rewind and reread the gauge input file
* 
        IF (PRVRST.EQ.1) THEN
         Y0 = 0
         M0 = 0
         D0 = 0
         REWIND(UnitNos(1))
         DONE = 0
        ENDIF

        READ(UnitNos(90),REC=1) SITE,YY,MM,DD
        I = 0
        FORCE = 1
 99     IF (I.LT.NOBS) READ(UnitNos(90),REC=I+1) 
     +                     (CASEID(J),J=1,4),Y,CaseWt,(X(J),J=1,P)
        IF ((CASEID(2).EQ.YY).AND.
     +      (CASEID(3).EQ.MM).AND.
     +      (CASEID(4).EQ.DD).AND.(I.LT.NOBS) ) THEN
         I = I+1
         SITE = CASEID(1)
         MISSFL(SITE) = 0
*
*	NB we are NOT ALLOWED to change year or month effects
*	in COVS between successive calls to COVSET, if COVSET is 
*       going to recalculate them on this pass because of a change
*       in the nonlinear parametrisation. (Possible source of much 
*       headache, already realised). We *DO*, however, have to 
*       explicitly set things here if COVSET's going to skip
*       over them. See COVSET header (file glm_base.f) for details.
*       Here, we set year effects if there's no nonlinear transformation
*       involved; we *always* set all other effects (currently, there
*       are no nonlinear transformations of month effects so COVSET
*       will skip over them all).
*
         DO 310 J=1,NP(1)
          COVS(SITE,J) = X(J)
 310     CONTINUE
         DO 311 J=NP(1)+1,NP(2)
          IF(THETA(J,1).GT.1.0D8) COVS(SITE,J) = X(J)
 311     CONTINUE
         DO 312 J=NP(2)+1,P
          COVS(SITE,J) = X(J)
 312     CONTINUE
         RSPNSE(SITE) = Y
         WtVec(SITE) = CaseWt
        ELSE
*
*       At this point we have populated the array of responses
*       and covariates for the current day, as calculated earlier.
*       Now need to update the covariates associated with 
*       nonlinear transformations. If any of these are associated 
*       with lagged observations, the easiest way to deal with it
*       is by calling the UPDATE routine. 
*
         IF (PRVRST.EQ.1) THEN
 315      CALL UPDATE(UnitNos(1),DatArray,THRESH,THRTYP,MissVal,
     +                RespIdx,NSITES,NVARS,Y0,M0,D0,MXLAG,SCODES,
     +                DONE,IFAIL)
          IF (IFAIL.NE.0) RETURN
          IF ((Y0.NE.YY).OR.(M0.NE.MM).OR.(D0.NE.DD)) GOTO 315
         ENDIF
         CALL COVSET(COVS,UnitNos(3:5),NP,NSITES,NVARS,RespIdx,
     +        COVCODE,TWO,THREE,SITINF,Distance,YY,MM,DD,DatArray,
     +        AllowIncAvge,PWTIDX,MXLAG,BETA,THETA,TRACE,RECALC,FORCE,
     +        MISSFL,1,IFAIL,MXP,0,0)
         IF (IFAIL.NE.0) RETURN
         FORCE = 0
*
*	Write the modified records to the data file (backwards -
*	we know the current record is I so we'll start there and
*	write to the previous line wherever there's a valid RSPNSE
*	value. 
*
         RECNO = I
         DO 320 SITE=NSITES,1,-1
          IF (RSPNSE(SITE).LT.-1.0D99) MISSFL(SITE) = 1
          IF (MISSFL(SITE).EQ.0) THEN
           Y = RSPNSE(SITE)
           CaseWt = WtVec(SITE)
           DO 321 J=1,P
            X(J) = COVS(SITE,J)
 321       CONTINUE
           WRITE(UnitNos(90),REC=RECNO) SITE,YY,MM,DD,
     +                                  Y,CaseWt,(X(J),J=1,P)
           RECNO = RECNO - 1
          ENDIF
          MISSFL(SITE) = 1
 320     CONTINUE
         IF (I.LT.NOBS) THEN
          READ(UnitNos(90),REC=I+1) SITE,YY,MM,DD
         ELSE
          GOTO 330
         ENDIF
        ENDIF
        IF (I.LE.NOBS) GOTO 99       
       ENDIF
*
*     If we're in here we haven't converged, so go back
*
 330   GOTO 100
      ENDIF

***************************************************************************
*	Only reach here when we've converged. Check that we've written
*       all the necessary residual information (if we haven't, it's
*       a programming error).
*
      IF (RESWRT.NE.1) THEN
       CALL INTPR('****ERROR**** Scratch file #2 not written '//
     +'properly in routine IWLS. This is a programming error '//
     +'- contact REC.',-1,0,0)
       IFAIL = 999
       RETURN
      ENDIF
      CNVRGE = 3
*
*	Calculate naive and robust covariance matrices. The naive one
*	is [X'WX]^-1. As it stands, XWX contains the Cholesky factor of 
*	the matrix X'WX (apart from a factor relating to the dispersion
*	parameter). The robust variance estimate is computed in a call
*	to SNDWCH. NB we also stored the variances of the dispersion 
*	parameter estimates here, because these are needed (strictly
*	speaking) when calculating adjusted likelihood ratios.
*
      If (Verbosity.Gt.0) then
       CALL INTPR('Computing covariance matrix of estimates ...',-1,0,0)
      End if
      CALL MATINV(CovNaive,P+1,1,MXP+1,1,IFAIL)
      IF (IFAIL.NE.0) THEN
       IF (IFAIL.GT.0) THEN
        IFAIL = 101
       ELSE
        IFAIL = 997
       ENDIF
       RETURN
      ENDIF
      CovNaive = DABS(PHI)*CovNaive
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF (MODEL.EQ.10) THEN
       CovNaive(P1,P1) = 2.0D0*MLPHI*MLPHI/DBLE(NOBS)
      ELSE IF (MODEL.EQ.11) THEN
       CovNaive(P1,P1) = (MLPHI**4) / 
     +                (DBLE(NOBS)*(TRIGAMM(1.0D0/MLPHI,IFAIL)-MLPHI))
      ENDIF 
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      CovRobust = CovNaive
      XWSWX = CovNaive
c###########################################################################
c###########################################################################
c
c       NEXT TWO LINES ARE A FUDGE; SOMETHING WRONG WITH ROBUST 
c       COVARIANCE MATRIX ESTIMATION WHEN DISPERSION PARAMETERS 
c       ARE INVOLVED, SO JUST DO THIS FOR THE REMAINING PARAMETERS
c       AND THEN SET THE DISPERSION PARAMETER DIAGONAL ENTRY EQUAL 
c       TO THAT FOR THE NAIVE COVARIANCE MATRIX. TO BE FIXED IN 
c       THE FUTURE ...
c
c###########################################################################
c###########################################################################
c      CALL SNDWCH(CovRobust,XWSWX,UFILNO,P1,NDAYS,MXP)
      CALL SNDWCH(CovRobust,XWSWX,UFILNO,P,NDAYS,MXP)
      IF (P1.GT.P) CovRobust(P+1,P+1) = CovNaive(P+1,P+1)
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF ((MODEL.EQ.11).AND.(TRACE.GT.0.0D0)) THEN
       WRITE(MESSAGE,6) 
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      ENDIF 
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (ITER.EQ.MAXIT) THEN
       WRITE(MESSAGE(1),5)
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      ENDIF
      CLOSE(UFILNO)
*
*	Copy THETAs to appropriate positions in BETA
*
      IF (NP(7).GT.NP(6)) THEN
       CALL NLUPDT(THETA,BETA,NP,COVCODE,MXP,IFAIL)
       IF (IFAIL.NE.0) RETURN
      ENDIF
*
*     If TRACE > 0 for gamma models, pass ML estimate of dispersion 
*     parameter back, as small changes can have a big impact on the 
*     initial log-likelihood for future fits, through the replacement
*     of trace values by their conditional expectations
*
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF ((MODEL.EQ.11).AND.(TRACE.GT.0.0D0)) PHI = MLPHI
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
 2    FORMAT(I15, F20.3,F13.4,2X,'(parameter ',I3,')')
 5    FORMAT(5X,'****WARNING**** Reached maximum ',
     +'number of iterations.')
 6    FORMAT('****WARNING**** Log-likelihoods may be inaccurate,',
     +' because trace values',/,'have been replaced by their',
     +' approximate conditional expectations.')
 10   FORMAT('****ERROR**** in routine IWLS, design matrix seems ',
     +'singular.',/5X,'Problem found on row ',I3,'. Have you ',
     +'accidentally duplicated a covariate',/5X,'that was already ',
     +'defined in another way, or specified a covariate that',/5X,
     +'doesn''t vary in your data set?')
 
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE NLUPDT(THETA,BETA,NP,COVCODE,MXP,IFAIL)
*
*       Updates unknown nonlinear parameters if we're estimating
*       them, at the end of each iteration (just copies them across
*       from BETA to appropriate bits of THETA). Variables:
*       THETA   - nonlinear parameters
*       BETA    - coefficients of linear predictors
*       NP      - Indexing for parameters
*       COVCODE - indicates which bits of THETA we're updating
*       MXP,MXN - for dimensioning
*       IFAIL   - Error flags
******************************************************************************
      INTEGER MXP,COVCODE(MXP),NP(10),IFAIL
      DOUBLE PRECISION BETA(0:MXP),THETA(MXP,3)
******************************************************************************
*       Additional INTEGERs
*       ^^^^^^^^^^^^^^^^^^^
*       I       - counter
*       IDX?    - used to locate correct element of THETA from COVCODE
******************************************************************************
      INTEGER I,IDX1,IDX2
******************************************************************************
*       Additional CHARACTERs
*       ^^^^^^^^^^^^^^^^^^^^^
*       MESSAGE - messages to screen
******************************************************************************
      CHARACTER MESSAGE(4)*255
*
*       THETAs are easy to locate because of the way COVCODE was
*       set up in MDLSET (NB integer division). IDX1 is the 
*       covariate, IDX2 the parameter. We check that parameters
*       aren't going off to infinity, and bail out if they are.
*
      IFAIL = 0
      DO 100 I=NP(6)+1,NP(7)
       IDX1 = COVCODE(I)/1000
       IDX2 = COVCODE(I) - (1000*IDX1)
       IF (DABS(BETA(I)).GT.1.0D8) THEN
        WRITE(MESSAGE,1) IDX2,IDX1
        CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(3)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(4)),-1,0,0)
        IFAIL = 998
        RETURN
       ENDIF
       THETA(IDX1,IDX2) = BETA(I)
 100  CONTINUE

 1    FORMAT('****ERROR**** During iteration, the estimate of ',
     +'parameter ',I1,/,'in the nonlinear transformation of ',
     +'covariate no. ',I2,/,'is approaching infinity. This isn''t',
     +' going to work - I''m bailing out.',/,
     +'Try different starting values.')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      BLOCK DATA CorFiles
*
*       To initialise info for common block EMPCOR below
*
      INTEGER OLDRES,OLDDEP,OLDFNO,CORFNO
      COMMON /CorFileInfo/ CORFNO,OLDRES,OLDDEP,OLDFNO
      
      DATA CORFNO,OLDRES,OLDDEP,OLDFNO /4*-1/

      END Block Data CorFiles
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE EMPCOR(FILNO,NSITES,NOBS,MODEL,RESTYP,DEPTYP,
     +                                         CORMAT,CORRN,IFAIL)
*
*   Calculates an empirical spatial correlation matrix or matrix of
*   joint occurrence probabilities, from data held in file connected 
*   to unit FILNO. Arguments:
*	
* FILNO  (INTEGER, input) Channel number of file to read
*	 data from. The file is a direct access file: each 
*	 record contains 4 integers (site, year, month and
*        day), then 4 doubles (observed and expected responses,
*	 Pearson residual and something that doesn't concern
*        us here).
* NSITES Number of sites used.
* NOBS	 Total number of records in file FILNO
* MODEL  Type of model. 1=logistic regression, 11=gamma.
* RESTYP (INTEGER, input). Codes for measure to compute 
*        correlations over, as follows:
*	    0 - raw responses
*	    1 - Pearson residuals (=raw residuals for Gaussian
*               observations)
*	    2 - Anscombe residuals for gamma observations
* DEPTYP (INTEGER, input). Type of dependence (1 = covariance,
* 	 2 = correlation, 3 = joint occurrence probabilities 
*        for binary data)
* CORMAT (DOUBLE, output). Correlation matrix.
* CORRN	 (DOUBLE,output). Ns associated with entries of CORMAT.
******************************************************************************
      INTEGER FILNO,NSITES,NOBS,MODEL,RESTYP,DEPTYP,IFAIL
      DOUBLE PRECISION CORMAT(NSITES,NSITES),CORRN(NSITES,NSITES)
******************************************************************************
*	Additional INTEGERs
*	^^^^^^^^^^^^^^^^^^^
*       I,J,K	- Counters
*       CASEID  - array of identifiers for each case, for purposes of
*                 residual analysis. Col.1 contains site number, 2 
*                 contains year, 3 contains month and 4 day
*	    OLDDAT }- Used to determine when the date changes
*	    NEWDAT }  
*       OLDRES  - Value of RESTYP on previous call
*       OLDDEP  - Value of DEPTYP on previous call
*       OLDFNO  - Value of FILNO on previous call
*       CORFNO  - Number of scratch file used to store the correlations
*                 when we've computed them
***************************************************************************
      INTEGER I,J,K,CASEID(4),OLDDAT,NEWDAT
      INTEGER OLDRES,OLDDEP,OLDFNO,CORFNO
******************************************************************************
*	Additional DOUBLEs
*	^^^^^^^^^^^^^^^^^^
*   Y       Daily response
*   MU      Predicted MEANS
*   ETA     Linear predictors
*   CaseWt  Case weights (interpreted as inverse variances of normalised
*           responses)
*   ONEDAY  Vector of today's values
*   MEAN    Means for calculating correlations. (i,j)th entry is
*           the mean appropriate for site i when calculating
*           correlation with site j.
*   VAR     Variances, similarly
*   TMP     Temporary storage
*   TOL	    Small number such that variances are treated as zero
*		    if they're less than it.
******************************************************************************
      DOUBLE PRECISION Y,MU,ETA,CaseWt
      DOUBLE PRECISION ONEDAY(NSITES)
      DOUBLE PRECISION MEAN(NSITES,NSITES),VAR(NSITES,NSITES)
      DOUBLE PRECISION TMP,TOL
******************************************************************************
*       Additional CHARACTERs
*       ^^^^^^^^^^^^^^^^^^^^^
*       MESSAGE - messages to screen
******************************************************************************
      CHARACTER MESSAGE*255
***************************************************************************
      COMMON /CorFileInfo/ CORFNO,OLDRES,OLDDEP,OLDFNO
      SAVE /CorFileInfo/
            
      TMP = 0.0D0
      IFAIL = 0
******************************************************************************
*     First scenario: we haven't done this before so we have to calculate 
*     everything
******************************************************************************
      IF ((RESTYP.NE.OLDRES).OR.(DEPTYP.NE.OLDDEP).OR.
     +    (FILNO.NE.OLDFNO)) THEN
*
*     Open scratch file for final matrix
*
       IF (CORFNO.NE.-1) THEN
        CLOSE(CORFNO)
        CORFNO = -1
       ENDIF
       OLDRES = RESTYP
       OLDDEP = DEPTYP
       OLDFNO = FILNO

       CALL SCROPN(2,2,0,CORFNO)
*
*       Initialise arrays (Fortran 90 syntax)
*
       TOL = 1.0D-12
       OLDDAT = 0
       ONEDAY = -1.0D12
       CORMAT = 0.0D0
       CORRN = 0.0D0
       MEAN = 0.0D0
       VAR = 0.0D0
*
*       Read data & compute the requested measure. Going to NOBS+1
*       guarantees we include the final day's measurements
*
       DO 100 I=1,NOBS+1
        IF (I.LE.NOBS) THEN
         READ(FILNO,REC=I) (CASEID(J),J=1,4),Y,MU,Eta,CaseWt
         IF (RESTYP.EQ.0) THEN
          TMP = Y
         ELSEIF (RESTYP.EQ.1) THEN
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          IF (MODEL.EQ.1) THEN
           TMP = (Y-MU)/DSQRT(MU*(1.0D0-MU))
          ELSEIF (MODEL.EQ.10) THEN
           TMP = Y-MU
          ELSEIF (MODEL.EQ.11) THEN
           TMP = (Y-MU)/MU
          ELSE
           WRITE(MESSAGE,1) ' MODEL'
           CALL INTPR(TRIM(MESSAGE),-1,0,0)
           IFAIL = 999
           RETURN
          ENDIF
          TMP = TMP * Dsqrt(CaseWt)
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

         ELSEIF (RESTYP.EQ.2) THEN
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          IF (MODEL.EQ.11) THEN
           TMP = ( (Y/MU)**(1.0D0/3.0D0) )
          ELSE
           WRITE(MESSAGE,1) ' MODEL'
           CALL INTPR(TRIM(MESSAGE),-1,0,0)
           IFAIL = 999
           RETURN
          ENDIF
         ELSE
          WRITE(MESSAGE,1) 'RESTYP'
          CALL INTPR(TRIM(MESSAGE),-1,0,0)
          IFAIL = 999
          RETURN
         ENDIF
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        ENDIF
*
*   Now, if the date has changed OR if we're on the last observation, 
*   update yesterday's contributions to the matrix (missing values are 
*   represented by -1.0D12). Always update today's residual vector. This 
*   scheme relies on the fact that the data are in chronological order.
*
        NEWDAT = (10000*CASEID(2))+(100*CASEID(3))+CASEID(4)
        IF ((NEWDAT.NE.OLDDAT).OR.(I.EQ.NOBS+1)) THEN
         OLDDAT=NEWDAT
         DO 200 J=1,NSITES
          IF (ONEDAY(J).GT.-0.9D12) THEN
           DO 205 K=J,NSITES
            IF (ONEDAY(K).GT.-0.9D12) THEN
             CORMAT(J,K) = CORMAT(J,K) + (ONEDAY(J)*ONEDAY(K))
             CORRN(J,K) = CORRN(J,K) + 1.0D0
             MEAN(J,K) = MEAN(J,K) + ONEDAY(J)
             VAR(J,K) = VAR(J,K) + (ONEDAY(J)**2)
             IF (K.NE.J) THEN 
              VAR(K,J) = VAR(K,J) + (ONEDAY(K)**2)
              MEAN(K,J) = MEAN(K,J) + ONEDAY(K)
             ENDIF
            ENDIF
 205       CONTINUE
           ONEDAY(J) = -1.0D12
          ENDIF
 200     CONTINUE
        ENDIF
*
*	And copy the latest entry into today's array
*
        ONEDAY(CASEID(1)) = TMP
 100   CONTINUE
*
*   That's everything collected now, so calculate the covariance 
*   structure and, if necessary, convert to correlations. If variances 
*   are zero, correlations are treated as zero. NB if DEPTYP=3 we're 
*   calculating joint occurrence probabilities rather than covariances,
*   so in this case there's no need to subtract the product of means.
*       
       DO 220 I=1,NSITES
        DO 221 J=I,NSITES
         IF (CORRN(I,J).GT.0.0D0) THEN
          IF (DEPTYP.NE.3) THEN 
           MEAN(I,J) = MEAN(I,J)/CORRN(I,J)
           IF (J.NE.I) MEAN(J,I) = MEAN(J,I)/CORRN(I,J)
           CORMAT(I,J) = CORMAT(I,J)/CORRN(I,J) - (MEAN(I,J)*MEAN(J,I))
          ELSE
           CORMAT(I,J) = CORMAT(I,J)/CORRN(I,J)
          ENDIF
          IF (DEPTYP.EQ.2) THEN
           VAR(I,J) = VAR(I,J)/CORRN(I,J) - (MEAN(I,J)**2)
           IF (J.NE.I) VAR(J,I) = VAR(J,I)/CORRN(I,J) - (MEAN(J,I)**2)
           IF ((VAR(I,J).GT.TOL).AND.(VAR(J,I).GT.TOL)) THEN
            CORMAT(I,J) = CORMAT(I,J)/DSQRT(VAR(I,J)*VAR(J,I))
           ELSE
            CORMAT(I,J) = 0.0D0
           ENDIF
          ENDIF
         ELSE
          CORMAT(I,J) = -9.999D0
         ENDIF
         CORMAT(J,I) = CORMAT(I,J)
         CORRN(J,I) = CORRN(I,J)
         WRITE(CORFNO) CORRN(I,J),CORMAT(I,J)
 221    CONTINUE
 220   CONTINUE

******************************************************************************
*     Second scenario: reread correlations as computed on a previous pass,
*     from scratch file
******************************************************************************
      ELSE
       REWIND(CORFNO)
       DO 320 I=1,NSITES
        DO 321 J=I,NSITES
         READ(CORFNO) CORRN(I,J),CORMAT(I,J)
         CORMAT(J,I) = CORMAT(I,J)
         CORRN(J,I) = CORRN(I,J)
 321    CONTINUE
 320   CONTINUE
      ENDIF
      
 1    FORMAT('****ERROR**** Invalid value of ',A6,' in routine EMPCOR')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE MLNUES(PHIHAT,SUMLNY,NOBS,LOGL,IFAIL)
*
*     To calculate ML estimate of dispersion parameter in a Gamma
*     GLM. Arguments are:
*     PHIHAT    - Estimated dispersion parameter. INPUT/OUTPUT
*     SUMLNY   - The sum of Ln(Y) values for the dataset. INPUT
*     NOBS     - number of observations. INPUT
*     LOGL     - On entry, the part of the log-likelihood corresponding
*                to the means. On exit, the full log-likelihood.
******************************************************************************
      DOUBLE PRECISION PHIHAT,SUMLNY,LOGL
      INTEGER NOBS,IFAIL
******************************************************************************
*     Extra DOUBLEs
*     ^^^^^^^^^^^^^
*     ALPHA    - Shape parameter of distribution (which is 1/PHI)
*     A,B      - limits of current interval (parametrise as log(ALPHA) to 
*                force everything to stay positive)
*     X        - new value
*     F1,F2,FX - objective function values (trying to get them = 0).
*     DIGAMM   - Digamma function
*     GAMMLN   - Log of gamma function
*     RHS      - Right-hand side of equation (bits not involving PHI)
*     TOL      - required precision
******************************************************************************
      DOUBLE PRECISION ALPHA,A,B,X,F1,F2,FX,DIGAMM,GAMMLN,RHS,TOL
******************************************************************************
*     Extra INTEGERs
*     ^^^^^^^^^^^^^^
*     I        - counter
*     IFAIL    - exit code for GAMMLN etc.
*     MAXIT    - maximum no. of iterations
******************************************************************************
      INTEGER I,MAXIT

      MAXIT = 20
      TOL = 1.0D-8
      IFAIL = 0
*
*       Code was originally written to estimate shape parameter 
*       rather than dispersion parameter, and it makes my brain
*       hurt to figure out the changes needed to estimate
*       the dispersion parameter so I've retained the old code
*       and will just convert back to dispersion parameter at the
*       end. Note that there are cases where PHIHAT will be 0 on
*	entry which causes problems (although most systems can 
*	escape from them - not 32-bit Windows installations 
*	though, except with high "verbosity" values in the R
*	call. I dunno)
*
      ALPHA = 1.0D0 / DMAX1(PHIHAT,1.0D-2)
*
*     Do root bracketing. Take as initial guesses half the input value of 
*     PHIHAT, and double it; check that root is bracketed before 
*     proceeding
*
      A = DLOG(ALPHA/2.0D0)
      B = DLOG(ALPHA*2.0D0)
      RHS = 1.0D0 + (LOGL/DBLE(NOBS)) + (SUMLNY/DBLE(NOBS))
      
 5    F1 = DIGAMM(DEXP(A),IFAIL) - A - RHS
      IF (IFAIL.NE.0) RETURN
      IF (F1.GT.0.0D0) THEN
       A = A - 1.0D0
       GOTO 5
      ENDIF
  
 6    F2 = DIGAMM(DEXP(B),IFAIL) - B - RHS
      IF (IFAIL.NE.0) RETURN
      IF (F2.LT.0.0D0) THEN
       B = B + 1.0D0
       GOTO 6
      ENDIF
      DO 10 I = 1,MAXIT
       IF (DABS(B-A).LT.TOL) GOTO 20
       X = ((A*F2)-(B*F1))/(F2-F1)
       FX = DIGAMM(DEXP(X),IFAIL) - X - RHS
       IF (IFAIL.NE.0) RETURN
       B = A
       F2 = F1
       A = X
       F1 = FX
 10   CONTINUE

 20   ALPHA = DEXP((A+B)/2.0D0)
      LOGL = ALPHA*(LOGL + ((DBLE(NOBS)*DLOG(ALPHA)) + SUMLNY) )
      LOGL = LOGL - (DBLE(NOBS)*GAMMLN(ALPHA,IFAIL))
      PHIHAT = 1.0D0 / ALPHA

      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE INVPAR(X,FX,TTYPE)
*
*     Performas inverse parabolic interpolation using 3 points
*     (implementation of Press et al 1992, equation 10.2.1).
*     Given 3 values of X, and corresponding function values FX,
*     fits a quadratic to estimate a turning point of FX. On exit,
*     X(3) contains the estimated turning point; FX(1) and FX(2)
*     contain the largest (if TTYPE=1) or smallest (if TTYPE=-1)
*     two of the input FX values, and X(1) and X(2) contain the 
*     corresponding Xs. 
******************************************************************************
      DOUBLE PRECISION X(3),FX(3)
      INTEGER TTYPE
******************************************************************************
*     Extra DOUBLES
*     ~~~~~~~~~~~~~
*     A,B,C,   }    All as in Press et al 1992, equation 10.2.1
*     FA,FB,FC }
*     NUMER         Numerator in formula
*     DENOM         Denominator
*     FOPT          Largest (smallest) input value of FX 
******************************************************************************
      DOUBLE PRECISION A,B,C,FA,FB,FC,NUMER,DENOM,FOPT
******************************************************************************
*     Extra INTEGERS
*     ~~~~~~~~~~~~~~
*     SMLIDX        Index of smallest (largest) element in FX
*     I             Counter
******************************************************************************
      INTEGER SMLIDX,I
      A = X(1)
      B = X(2)
      C = X(3)
       
      FA = FX(1)
      FB = FX(2)
      FC = FX(3)

*
*     Get rid of smallest element from FX and move elements of X
*     accordingly
*
      SMLIDX = 1
      FOPT = FX(1)
      DO 100 I=2,3
       IF ( (DBLE(TTYPE)*FX(I)).LT.FOPT) THEN
        SMLIDX = I
        FOPT = FX(I)
       ENDIF
 100  CONTINUE

      IF (SMLIDX.NE.3) THEN
       X(SMLIDX) = X(3)
       FX(SMLIDX) = FX(3)
      ENDIF
      FX(3) = 0.0D0

*
*     OK, now compute new value of X
*
      NUMER = ( ((B-A)**2)*(FB-FC) ) - ( ((B-C)**2)*(FB-FA) )
      DENOM = ( (B-A)*(FB-FC) ) - ( (B-C)*(FB-FA) )
      X(3) = B - (0.5D0*NUMER/DENOM)
      

      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE FITCAL(BETA,NP,P,MAXP,MODEL,SCRFNO,RECNO,ITER,
     +               TRWRN,TRACE,DISPAR,LLINVL,X,Y,CaseWt,SITE,
     +               YY,MM,DD,LOGL,MU,ETA,IFAIL)
*
*     To calculate (a) fitted value (b) linear predictor for a 
*     single observation, and increment the log-likelihood  
*     accordingly.
*
*     Arguments:
*
*     BETA      (DOUBLE, input) Coefficient vector for linear predictor
*     NP        (INTEGER, input) Numbers of parameters in each component
*               of BETA
*     P         (INTEGER,input) Number of parameters being estimated
*     MAXP      (INTEGER, input) Storage allocated to X and P
*     MODEL     (INTEGER, input): model code. See header to IWLS
*               routine for values
*     SCRFNO    (INTEGER, input): number of scratch file holding 
*               covariate info etc.
*     RECNO     (INTEGER, input): the record number to read from SCRFNO
*     ITER      (INTEGER, input): iteration number
*     TRWRN     (INTEGER, input): indicates whether we've issued a 
*               warning message about trace values on the current
*               iteration
*     TRACE     (DOUBLE, input): threshold for `trace' values
*     DISPAR    (DOUBLE, input): dispersion parameter for model
*     LLINVL    (DOUBLE, input): 'Invalid' log-likelihood
*     LOGL      (DOUBLE, input/output): on input, log-likelihood before
*               looking at current case; on output, log-likelihood
*               including current case
*     X         (DOUBLE, output) Vector of covariates, read from file
*     Y         (DOUBLE, output) Response, read from file
*     CaseWt    (DOUBLE, output) Weight for this case, read from file
*     SITE      (INTEGER, output) Site code, read from file
*     YY        (INTEGER, output) Year, read from file
*     MM        (INTEGER, output) Month, read from file
*     DD        (INTEGER, output) Day, read from file
*     MU        (DOUBLE, output): fitted value for current case
*     ETA       (DOUBLE, output): linear predictor for current case
*
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*////////////	NB THIS ROUTINE CONTAINS SECTIONS WHICH ARE MODEL-   ///////
*////////////	SPECIFIC. THESE SECTIONS ARE ENCLOSED BY LINES OF    ///////
*////////////   HASHES LIKE THIS, FOR EASY RECOGNITION IF YOU WANT   ///////
*////////////   TO AMEND THEM SUBSEQUENTLY, OR TO ADD EXTRA MODELS.  ///////
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

******************************************************************************
      INTEGER NP(10),P,MAXP,MODEL,SCRFNO,RECNO,TRWRN,IFAIL,ITER
      INTEGER SITE,YY,MM,DD
      DOUBLE PRECISION BETA(0:MAXP),TRACE,DISPAR,LLINVL
      DOUBLE PRECISION X(0:MAXP),Y,CaseWt,LOGL,MU,ETA
******************************************************************************
*     Additional INTEGERs
*     ^^^^^^^^^^^^^^^^^^^
*     J               Counter
******************************************************************************
      INTEGER J
******************************************************************************
*     Additional DOUBLEs
*     ^^^^^^^^^^^^^^^^^^
*     TMP             Temporary storage
*     LLincr          Increment in log-likelihood     
******************************************************************************
      DOUBLE PRECISION TMP, LLincr
******************************************************************************
*       Additional CHARACTERs
*       ^^^^^^^^^^^^^^^^^^^^^
*       MESSAGE - messages to screen
******************************************************************************
      CHARACTER MESSAGE*255

      IFAIL = 0
*
*       Compute linear predictor - NB using parameters up to NP(6) - the
*       rest of BETA is just used in updating nonlinear parameters
*
       ETA = BETA(0)
       READ(SCRFNO,REC=RECNO) SITE,YY,MM,DD,Y,CaseWt,(X(J),J=1,P)
       DO 100 J=1,NP(6)
        ETA = ETA + (BETA(J)*X(J))
 100   CONTINUE
*
*	And fitted values ...
*
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       IF (MODEL.EQ.1) THEN
        TMP = DEXP(ETA)
        MU = TMP/(1.0D0+TMP)
       ELSEIF (MODEL.EQ.10) THEN
        MU = ETA
       ELSEIF (MODEL.EQ.11) THEN
        MU = DEXP(ETA)
       ELSE
        WRITE(MESSAGE,1)
        CALL INTPR(TRIM(MESSAGE),-1,0,0)
        IFAIL = 999
        RETURN
       ENDIF
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*
*	Replace trace values by their approximate conditional
*	expectations as appropriate. We need to modify the value of 
*       Y in the data file as well if this happens.
*
       IF (Y.LT.TRACE) THEN
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*
*	Gamma model. This is rather a crude approximation, but it 
*	avoids having to do numerical integration, which is important 
*	from the point of view of speed! We also need a NaN trap
*	in case we iterate to somewhere peculiar and end up dividing
*	by zero ... this is accomplished via IF(.NOT.(Y.LE.TRACE)).
*       Use ML estimates of dispersion parameter, so that we're
*       effectively doing an EM algorithm kind of thing.
*
        IF (MODEL.EQ.11) THEN
         Y = (MU*DISPAR*TRACE)/( (2.0D0*MU) + (DISPAR*TRACE) )
         IF (.NOT.(Y.LE.TRACE) ) THEN
          IF (TRWRN.LT.10) THEN
*
*       Next line uses R INTPR routine to output warning message
*
           CALL INTPR('****WARNING**** Some trace values have'//
     +     ' approximate conditional expectation > trace threshold',
     +     -1,0,0)
           TRWRN = TRWRN + 10
          ENDIF
          Y = TRACE*0.99D0
         ELSEIF (Y.LE.0.0D0) THEN
          IF (MOD(TRWRN,10).EQ.0) THEN
*
*       Next line uses R INTPR routine to output warning message
*
           CALL INTPR('****WARNING**** Some trace values have'//
     +     ' approximate conditional expectation < 0',-1,0,0)
           TRWRN = TRWRN + 1
          ENDIF
          Y = TRACE/1.0D2
         ENDIF
         WRITE(SCRFNO,REC=RECNO) SITE,YY,MM,DD,Y,CaseWt,(X(J),J=1,P)
        ELSEIF (MODEL.EQ.10) THEN
         CALL INTPR('****ERROR**** Trace values not handled in '//
     +              'Gaussian model',-1,0,0)
         IFAIL = 1
         RETURN
        ENDIF
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       ENDIF
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*
*	For binary outcomes, check that fitted values aren't 0 or 1 
*	to machine precision. Only stop the prog if this happens with
*	starting values - otherwise we'll backtrack on the iteration
*	later on.
*
       IF (MODEL.EQ.1) THEN
        IF ( (MU.EQ.0.0D0).OR.(MU.EQ.1.0D0) ) THEN
         IF (ITER.EQ.0) THEN
          IFAIL = 100
          RETURN
         ELSE
          LOGL = LLINVL
          IFAIL = 2
         ENDIF
        ENDIF
       ENDIF
*
*	Now update log-likelihood. Next line is unnecessary, and 
*       included merely to avoid compiler warnings about LLincr
*       possibly being uninitialised.
*
       LLincr = 0.0d0
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       IF (MODEL.EQ.1) THEN
*
*	Logistic ... (McCullagh & Nelder, eqn 4.13)
*
        LLincr = - DLOG(1.0D0+DEXP(ETA))
        IF (Y.GT.0.0D0) LLincr = LLincr + ETA
       ELSEIF (MODEL.EQ.10) THEN
*
*       Normal. Usually you'd omit terms involving the dispersion parameter
*	here, but need to include them in case we're fitting a heteroscedastic
*	model using CaseWt as the inverse of the variance 
*
        LLincr = (DLOG(CaseWt/DISPAR) - 
     +            (CaseWt*(Y-MU)*(Y-MU)/DISPAR)) / 2.0D0
       ELSEIF (MODEL.EQ.11) THEN
*
*	Gamma. Omit terms that don't depend on MU, because estimation 
*       of NU at each stage screws things up. So this isn't strictly 
*       the log-likelihood, but it's the bit we think is important 
*       (see McCullagh & Nelder, p.290, for example). NB at the end 
*       we need to remember to multiply through by dispersion parameter,
*       but not here! (save on computation). This is dealt with in 
*       routine MLNUES.
*
        LLincr = - (Y/MU) - DLOG(MU)
       ENDIF
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*
*       And finally, add the (weighted) increment to the log-likelihood
*
      LOGL = LOGL + LLincr

 1    FORMAT(/'*****ERROR***** Invalid model code found in FITCAL.')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE RESID(MODEL,PHI,NP,NOBS,FILNO,NSITES,
     +                 YBAR,SY,ERRBAR,MSE,RSQ,PRBAR,SPR,SEPR,
     +                 ARBARO,SARO,ARBARE,SARE,BYPRED,PRMNTH,PRSITE,
     +                 PRYEAR,FIRSTYR,NYEARS,IFAIL)
*
*       To perform a basic residual analysis - Pearson residuals
*       broken down by month, site and year, along with various other
*       model summaries.
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*////////////   NB THIS ROUTINE CONTAINS QUITE A FEW SECTIONS WHICH  ///////
*////////////   ARE MODEL-SPECIFIC. THESE SECTIONS ARE ENCLOSED BY   ///////
*////////////   LINES OF HASHES LIKE THIS, FOR EASY RECOGNITION IF   ///////
*////////////   YOU SUBSEQUENTLY WANT TO AMEND THEM OR ADD EXTRA     ///////
*////////////   MODELS.                                              ///////
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMETERS      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
***************************************************************************
*       Arguments:
*       ^^^^^^^^^^
*       MODEL   Code to indicate which model we're dealing with.
*               Current values are
*               1 : Logistic regression for binary responses
*               10: Normal distribution with identity link
*               11: Gamma model with log link
*	FILNO	File number of scratch file with residual info in
*       NP      No. of parameters estimated in model (=no. of lost
*               degrees of freedom)
*       NOBS    no. of observations
*       NSITES  no. of sites required
*       YBAR    Mean of Ys
*       SY      Std Dev of Ys
*       ERRBAR  Mean prediction error
*       MSE     Mean squared error
*       RSQ     R-squared (1 - Error sum of squares/Data sum of squares)
*       PRBAR   Mean Pearson residual
*       SPR     Std Dev of pearson residuals
*	SEPR	Standard Error of mean Pearson residual
*       PRMNTH  Pearson residuals by month } Column 0 is N, 1 is mean,
*       PRSITE  And by site                } 2 is Std Dev, 3 is S.e. (mean)
*       PRYEAR  And by year                } if model is correct
*       FIRSTYR - First year for which we'll ever have data
*       NYEARS  - No. of years
*       ARBARO  Observed mean Anscombe residual (gamma model)
*       SARO    Observed Std dev of Anscombe residuals (gamma model)
*       ARBARE  Expected mean Anscombe residual (gamma model)
*       SARE    Expected Std dev of Anscombe residuals (gamma model)
*	BYPRED	Table of occurrence frequencies by predicted       } Logistic
*		probability. Row 1 is observed proportion in each  } model
*		category, Row 2 is expected, and Row 3 is N.       } 
*       IFAIL   Error flag
******************************************************************************
      Integer, intent(in) :: FIRSTYR,NYEARS
      Integer, intent(in) :: MODEL,NSITES,NP,NOBS,FILNO
      Integer, intent (out) :: IFAIL
      Double precision, intent(out) :: YBAR,SY,ERRBAR,MSE,RSQ
      Double precision, intent(out) :: PRBAR,SPR,SEPR
      Double precision, intent(out) :: ARBARO,SARO,ARBARE,SARE
      Double precision, intent(out) :: BYPRED(3,0:9),PRMNTH(12,0:3),
     +      PRSITE(NSITES,0:3),PRYEAR(FIRSTYR:FIRSTYR+NYEARS-1,0:3)
******************************************************************************
*               INTEGER variables
*               ^^^^^^^^^^^^^^^^^
*       I,J,K	- Counters
*       CASEID  - array of identifiers for each case, for purposes of
*                 residual analysis. Col.1 contains site number, 2 
*                 contains year, 3 contains month and 4 day
*       SFLAGS  - Flags indicating which sites have data present for the
*                 current day.
*       CURDAT  - Current date
*       OLDDAT  - Date for previous record
*       NADJM   - Numbers of observations from each pair of sites in 
*                 each month (used for calculating adjustments for
*                 spatial dependence)
*       NADJY   - Ditto, for each year 
*       YEAR    - Corresponding to values in SFLAGS
*       MONTH   - Ditto
*	CTABLE	- Classification table in residual analysis (logistic model)
*	BOXNO	- Used to identify which part of the forecast distribution
*		  each observation falls in (logistic model)
***************************************************************************
      INTEGER I,J,K
      INTEGER CASEID(4),SFLAGS(NSITES),CURDAT,OLDDAT      
      INTEGER NADJM(12,NSITES,NSITES),YEAR,MONTH
      INTEGER NADJY(FIRSTYR:FIRSTYR+NYEARS-1,NSITES,NSITES),BOXNO
***************************************************************************
*               DOUBLE PRECISION variables
*               ^^^^^^^^^^^^^^^^^^^^^^^^^^
*       Y       Daily value (0 or 1 for logistic)
*       MU      Predicted MEANS (probs for logistic)
*       ETA     Linear predictors
*       CaseWt  Weight attached to current case (interpreted as inverse
*               of variance)
*       PHI      Fitted dispersion parameter
*       RAWRES  Raw residual (observed - expected)
*	SESCAL	Scaling factor to correct standard errors for spatial
*		dependence
*	CORMAT	- Spatial correlation matrix for Pearson residuals
*	CORRN	- Ns going to make up entries in CORMAT
*       SEADJM  Adjustments to standard errors for monthly means
*       SEADJY  Ditto, annual means 
*       TMP     Temporary storage
******************************************************************************
      DOUBLE PRECISION Y,MU,ETA,CaseWt,PHI,SESCAL,RAWRES,CURPR
      DOUBLE PRECISION CORMAT(NSITES,NSITES),CORRN(NSITES,NSITES)
      DOUBLE PRECISION SEADJM(12),SEADJY(FIRSTYR:FIRSTYR+NYEARS-1)
      DOUBLE PRECISION TMP
***************************************************************************
*	Check dimensioning:
*
*
*       Initialise arrays (easier this way than in a DATA statement)
*
      YBAR = 0.0D0
      SY = 0.0D0
      ERRBAR = 0.0D0
      MSE = 0.0D0
      CURPR = 0.0D0
      PRBAR = 0.0D0
      PRBAR = 0.0D0
      SPR = 0.0D0
      ARBARO = 0.0D0
      SARO = 0.0D0
      DO 100 J=0,3
       DO 101 I=1,NSITES
        PRSITE(I,J) = 0.0D0
        SFLAGS(I) = 0
 101   CONTINUE
       DO 102 I=1,12
        PRMNTH(I,J) = 0.0D0
        SEADJM(I) = 0.0D0
 102   CONTINUE
       DO 103 I=FIRSTYR,FIRSTYR+NYEARS-1
        PRYEAR(I,J) = 0.0D0
        SEADJY(I) = 0.0D0
 103   CONTINUE
 100  CONTINUE

      DO 105 I=1,NSITES
       DO 106 J=I,NSITES
        DO 107 K=1,12
         NADJM(K,I,J) = 0
 107    CONTINUE
        DO 108 K=FIRSTYR,FIRSTYR+NYEARS-1
         NADJY(K,I,J) = 0
 108    CONTINUE
 106   CONTINUE
 105  CONTINUE
      OLDDAT = 0
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*
*     Logistic model: initialise classification table & theoretical
*     std dev of Pearson residuals
*
      IF (MODEL.EQ.1) THEN
       DO 112 I=1,3
        DO 113 J=0,9
         BYPRED(I,J) = 0.0D0
 113    CONTINUE
 112   CONTINUE
       PHI = 1.0D0
*
*     Gamma model: calculate expected Anscombe residual properties
*
      ELSEIF (MODEL.EQ.11) THEN
       CALL ANSPRM(1.0D0/PHI,ARBARE,SARE,IFAIL)
       IF (IFAIL.NE.0) RETURN
       SARE = DSQRT(SARE)
      ENDIF
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*
*       Go through and calculate all required sums, and sums of squares.
*
      DO 200 I=1,NOBS
       READ(FILNO,REC=I) (CASEID(J),J=1,4),Y,MU,ETA,CaseWt
       RAWRES = Y - MU
       MSE = MSE + (RAWRES**2)
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

*
*     Logistic model: produce classification tables, and standard errors. 
*     Rows of CTABLE represent stratification by predicted outcome.
*
       IF (MODEL.EQ.1) THEN
        CURPR = RAWRES/DSQRT(MU*(1.0D0-MU))
*
*	Now more complete stratification - split forecasts into deciles
*
        BOXNO = INT(1.0D1*MU)
        BYPRED(1,BOXNO) = BYPRED(1,BOXNO) + Y
        BYPRED(2,BOXNO) = BYPRED(2,BOXNO) + MU
        BYPRED(3,BOXNO) = BYPRED(3,BOXNO) + 1.0D0
*
*     Gaussian model: summary of observations and R^2
*
       ELSEIF (MODEL.EQ.10) THEN
        CURPR = RAWRES
        YBAR = YBAR + Y
        SY = SY + (Y**2)
        ERRBAR = ERRBAR + RAWRES
*
*     Gamma model: summary of observations, and Anscombe residuals
*
       ELSEIF (MODEL.EQ.11) THEN
        CURPR = RAWRES / MU
        YBAR = YBAR + Y
        SY = SY + (Y**2)
        ERRBAR = ERRBAR + RAWRES
        TMP = (Y/MU)**(1.0D0/3.0D0)
        ARBARO = ARBARO + TMP
        SARO = SARO + ( TMP**2 )
       ENDIF
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*
*	Now Pearson residuals - firstly, overall scores and secondly,
*	stratified by year, month and site. Also calculate contributions
*       to spatial dependence adjustments
*
       CURPR = CURPR * DSQRT(CaseWt)
       PRBAR = PRBAR + CURPR
       SPR = SPR + (CURPR**2)
*
*	0s are N
*
       PRSITE(CASEID(1),0) = PRSITE(CASEID(1),0) + 1.0D0
       PRYEAR(CASEID(2),0) = PRYEAR(CASEID(2),0) + 1.0D0
       PRMNTH(CASEID(3),0) = PRMNTH(CASEID(3),0) + 1.0D0
*
*       1s are means ...
*
       PRSITE(CASEID(1),1) = PRSITE(CASEID(1),1) + CURPR
       PRYEAR(CASEID(2),1) = PRYEAR(CASEID(2),1) + CURPR
       PRMNTH(CASEID(3),1) = PRMNTH(CASEID(3),1) + CURPR
*
*       2s are observed sds (sum of squares for now)
*
       PRSITE(CASEID(1),2) = PRSITE(CASEID(1),2) + (CURPR**2)
       PRYEAR(CASEID(2),2) = PRYEAR(CASEID(2),2) + (CURPR**2)
       PRMNTH(CASEID(3),2) = PRMNTH(CASEID(3),2) + (CURPR**2)
*
*      Here are the contributions to adjustments for spatial dependence -
*      need to know, for each subset of data, how many times each
*      pair of sites occurs
*
       CURDAT = (10000*CASEID(2)) + (100*CASEID(3)) + CASEID(4)
       SFLAGS(CASEID(1)) = 1
       IF ( ((CURDAT.NE.OLDDAT).AND.(OLDDAT.NE.0)) .OR. 
     +      (I.EQ.NOBS) ) THEN
        YEAR = OLDDAT/10000
        MONTH = MOD(OLDDAT,10000) / 100
        DO 120 J=1,NSITES
         DO 121 K=J,NSITES
          IF ( (SFLAGS(J).EQ.1).AND.(SFLAGS(K).EQ.1) ) THEN
           NADJM(MONTH,J,K) = NADJM(MONTH,J,K) + 1
           NADJY(YEAR,J,K) = NADJY(YEAR,J,K) + 1
          ENDIF
 121     CONTINUE
         SFLAGS(J) = 0
 120    CONTINUE
       ENDIF
       OLDDAT = CURDAT
 200  CONTINUE
*
*	That's everything collected now. Next, compute summary statistics.
*	Start by computing inflation factors for standard errors, to 
*	compensate for spatial dependence. This involves calculating the
*	inter-site correlation matrix. NB FOR CONVENIENCE, DIVIDE
*       DIAGONAL ENTRIES BY 2, SINCE WE DON'T USE THIS MATRIX
*       SUBSEQUENTLY
*
       CALL EMPCOR(FILNO,NSITES,NOBS,MODEL,1,2,CORMAT,CORRN,IFAIL)
       IF (IFAIL.NE.0) RETURN
       SESCAL = 0.0D0
       DO 220 I=1,NSITES
        CORMAT(I,I) = 0.5D0 * CORMAT(I,I)
        DO 221 J=I,NSITES
         IF (CORRN(I,J).GT.0.0D0) THEN
          SESCAL = SESCAL + (2.0D0*CORRN(I,J)*CORMAT(I,J))
         ENDIF
         DO 225 K=1,12
          IF (NADJM(K,I,J).GT.0) THEN
           SEADJM(K) = SEADJM(K) + 
     +                     (2.0D0*CORMAT(I,J)*DBLE(NADJM(K,I,J)))
          ENDIF
 225     CONTINUE
         DO 226 K=FIRSTYR,FIRSTYR+NYEARS-1
          IF (NADJY(K,I,J).GT.0) THEN
           SEADJY(K) = SEADJY(K) + 
     +                     (2.0D0*CORMAT(I,J)*DBLE(NADJY(K,I,J)))
          ENDIF
 226     CONTINUE
 221   CONTINUE
 220  CONTINUE
      SESCAL = DSQRT(SESCAL/DBLE(NOBS))
      DO 230 I=1,12
       IF (PRMNTH(I,0).GT.0.0D0) 
     +                   SEADJM(I) = DSQRT(SEADJM(I)/PRMNTH(I,0))
 230  CONTINUE
      DO 231 I=FIRSTYR,FIRSTYR+NYEARS-1
       IF (PRYEAR(I,0).GT.0.0D0) 
     +                   SEADJY(I) = DSQRT(SEADJY(I)/PRYEAR(I,0))
 231  CONTINUE
*
*       Now overall performance (measures of which are model specific). 
*       Note that all Std Err variables currently contain sums of
*       squares. For model std errs we divide by NOBS - NP
*
      PRBAR = PRBAR/DBLE(NOBS)
      SPR = (SPR - (DBLE(NOBS)*(PRBAR**2)))/DBLE(NOBS-NP)
      SPR = DSQRT(SPR)
      SEPR = SESCAL*SPR/DSQRT(DBLE(NOBS))
      MSE = MSE/DBLE(NOBS-NP)
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*
*     Logistic model
*
      IF (MODEL.EQ.1) THEN
       DO 206 I=0,9
        IF (BYPRED(3,I).GT.0.0D0) THEN
         BYPRED(1,I) = BYPRED(1,I)/BYPRED(3,I)
         BYPRED(2,I) = BYPRED(2,I)/BYPRED(3,I)
        ENDIF
 206   CONTINUE
*
*     Gaussian model
*
      ELSEIF (MODEL.EQ.10) THEN
       YBAR = YBAR/DBLE(NOBS)
       SY = (SY - (DBLE(NOBS)*(YBAR**2)))/DBLE(NOBS-1)
       RSQ = 1.0D0 - (MSE/SY)
       SY = DSQRT(SY)
       ERRBAR = ERRBAR/DBLE(NOBS)
*
*     Gamma model
*
      ELSEIF (MODEL.EQ.11) THEN
       YBAR = YBAR/DBLE(NOBS)
       SY = (SY - (DBLE(NOBS)*(YBAR**2)))/DBLE(NOBS-1)
       RSQ = 1.0D0 - (MSE/SY)
       SY = DSQRT(SY)
       ERRBAR = ERRBAR/DBLE(NOBS)
       ARBARO = ARBARO/DBLE(NOBS)
       SARO = (SARO - (DBLE(NOBS)*(ARBARO**2)))/DBLE(NOBS-NP)
       SARO = DSQRT(SARO)
      ENDIF
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*
*       Broken down by month. Use divisor N-1 for observed standard 
*       deviations (they're deviations from observed means) 
*
      DO 210 I=1,12
       IF (PRMNTH(I,0).GT.0.0D0) THEN
        PRMNTH(I,1) = PRMNTH(I,1)/PRMNTH(I,0)
        IF (PRMNTH(I,0).GT.1.0D0) THEN
         PRMNTH(I,2) = (PRMNTH(I,2)-(PRMNTH(I,0)*(PRMNTH(I,1)**2)))
     +                                          /(PRMNTH(I,0)-1.0D0)
         PRMNTH(I,2) = DSQRT(PRMNTH(I,2))
        ELSE
*
*     If there's only 1 observations, no need to correct for spatial
*     dependence so omit SESCAL!
*
         PRMNTH(I,2) = DSQRT(PHI)
        ENDIF
        PRMNTH(I,3) = SEADJM(I)*DSQRT(PHI)/(DSQRT(PRMNTH(I,0)))
       ENDIF
 210  CONTINUE
*
*       Broken down by site
*
      DO 211 I=1,NSITES
       IF (PRSITE(I,0).GT.0.0D0) THEN
        PRSITE(I,1) = PRSITE(I,1)/PRSITE(I,0)
        IF (PRSITE(I,0).GT.1.0D0) THEN
         PRSITE(I,2) = (PRSITE(I,2)-(PRSITE(I,0)*(PRSITE(I,1)**2)))
     +                                          /(PRSITE(I,0)-1.0D0)
         PRSITE(I,2) = DSQRT(PRSITE(I,2))
        ELSE
         PRSITE(I,2) = DSQRT(PHI)
        ENDIF
        PRSITE(I,3) = DSQRT(PHI)/(DSQRT(PRSITE(I,0)))
       ENDIF
 211  CONTINUE
*
*       Broken down by year
*
      DO 212 I=FIRSTYR,FIRSTYR+NYEARS-1
       IF ( PRYEAR(I,0).GT.0.0D0 ) THEN
        PRYEAR(I,1) = PRYEAR(I,1)/PRYEAR(I,0)
        IF (PRYEAR(I,0).GT.1.0D0) THEN
         PRYEAR(I,2) = (PRYEAR(I,2)-(PRYEAR(I,0)*(PRYEAR(I,1)**2)))
     +                                          /(PRYEAR(I,0)-1.0D0)
         PRYEAR(I,2) = DSQRT(PRYEAR(I,2))
        ELSE
         PRYEAR(I,2) = DSQRT(PHI)
        ENDIF
        PRYEAR(I,3) = SEADJY(I)*DSQRT(PHI)/(DSQRT(PRYEAR(I,0)))
       ENDIF
 212  CONTINUE

      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE SNDWCH(XWXINV,XWSWX,FILNO,P,NDAYS,MXP)
*
*     To calculate sandwich variance estimator in the presence of 
*     inter-site dependence. Arguments:
*
*     XWXINV   - On entry, the `independence' variance estimator
*                (i.e. [X'WX]^{-1}). On exit, the corrected 
*                estimator.
*     XWSWX    - On exit, the filling in the sandwich.
*     FILNO    - Number of scratch file containing daily scores
*     P        - Number of parameters in model (excluding constant)
*     NDAYS    - Number of days of data
*     MXP      - Dimensioning for XWXINV and XWSWX
******************************************************************************
      INTEGER FILNO,P,NDAYS,MXP
      DOUBLE PRECISION XWXINV(0:MXP,0:MXP),XWSWX(0:MXP,0:MXP)
******************************************************************************
*     Extra INTEGERs
*     ^^^^^^^^^^^^^^
*     I,J,K    - counters
******************************************************************************
      INTEGER I,J,K
******************************************************************************
*     Extra DOUBLEs
*     ^^^^^^^^^^^^^
*     U        - Score vector read from file
*     TMP      - temporary storage
******************************************************************************
      DOUBLE PRECISION U(0:MXP),TMP(0:MXP,0:MXP)
*
*     Now compute the matrix XWSWX. It's just the sum of squares
*     of rows of FILNO. As it's symmetric, we just compute the upper 
*     triangle.
*
      REWIND(FILNO)
      XWSWX = 0.0D0

      DO 120 K=1,NDAYS
       READ(FILNO) (U(I),I=0,P)
       DO 130 I=0,P
        DO 131 J=I,P
         XWSWX(I,J) = XWSWX(I,J) + (U(I)*U(J))
 131    CONTINUE
 130   CONTINUE
 120  CONTINUE
*
*     Right: at this point XWSWX is correct in its upper triangle.
*     The lower triangle hasn't been set, but as it's symmetric it
*     doesn't matter. The matrix we require is [XWXINV][XWSWX][XWXINV].
*     It's quickest to do this in 2 steps: first compute [XWXINV][XWSWX]
*     and store in TMP; then compute [TMP][XWXINV]. Since the result
*     is symmetric, we only need to compute one half of it, which is 
*     stored in the lower triangle of XWSWX for the moment. The diagonal
*     goes into U.
*
      DO 160 I=0,P
       DO 165 J=0,P
        TMP(I,J) = 0.0D0
        DO 170 K=0,J-1
         TMP(I,J) = TMP(I,J) + (XWXINV(I,K)*XWSWX(K,J))
 170    CONTINUE
        DO 171 K=J,P
         TMP(I,J) = TMP(I,J) + (XWXINV(I,K)*XWSWX(J,K))
 171    CONTINUE
 165   CONTINUE
 160  CONTINUE

      DO 180 I=0,P
       U(I) = 0.0D0
       DO 185 K=0,P
        U(I) = U(I) + (TMP(I,K)*XWXINV(K,I))
        DO 190 J=0,I-1
         XWSWX(I,J) = XWSWX(I,J) + (TMP(I,K)*XWXINV(K,J))
 190    CONTINUE
 185   CONTINUE
 180  CONTINUE
*
*     Finally, copy correct covariance matrix into XWXINV, and tidy up
*     XWSWX
*
      DO 200 I=0,P
       XWXINV(I,I) = U(I)
       DO 201 J=I+1,P
        XWXINV(I,J) = XWSWX(J,I)
        XWXINV(J,I) = XWXINV(I,J)
        XWSWX(J,I) = XWSWX(I,J)
 201   CONTINUE
 200  CONTINUE

      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE CORFIT(SPMOD,CORMAT,CORRN,SITINF,NATTR,
     +                        SCODES,NSITES,RHO,MXP,IFAIL)
***************************************************************************
*       To fit spatial correlation models to a set of empirical correlations.
*       Arguments:
*
*       SPMOD   (integer, input). Code indicating the model to be fitted
*       CORMAT  (double precision, input). A matrix of empirical inter-
*               site correlations
*       CORRN   (double precision, input). Matrix of sample sizes from
*               which the entries in CORMAT have been calculated. 
*       SITINF  (double precision, input). Array containing site information. 
*               See calling program for details
*       NATTR   (integer, input). Number of site attributes that have been 
*               defined
*       SCODES  - array of short site codes
*       NSITES  (integer, input). Number of sites for which empirical 
*               correlations have been calculated.
*       RHO     (double precision, input / output). Array of correlation
*               model parameters. On input this contains an initial 
*               estimate; on output it contains the fitted value
*       MXP     (integer, input). Maximum number of model parameters.
*               Used for dimensioning RHO.
*       IFAIL   Error flag
*
*       The estimation is done by minimising the weighted sum of 
*       squares
*
*       S = \sum n_{ij} (c(x_{ij},y_{ij},\rho)-r_{ij})^2
*
*       where the sum is over all pairs (i,j) of sites; n_{ij} and r_{ij}
*       are the number of observations and empirical correlation for pair 
*       (i,j); and c(x_{ij},y_{ij},rho) is the fitted correlation for
*       that pair of sites which is possibly dependent on the x- and 
*       y-separations x_{ij} and y_{ij}, and the correlation parameter(s)
*       rho. 
*
*       The general form of the correlation functions in this routine
*       is 
*
*       c(x,y) = lambda[alpha + (1-alpha)*cc(x,y)]
*
*       where cc(.,.) is a valid spatial correlation function and
*       alpha and lambda are in the range (0,1). Where there is 
*       no "threshold" parameter alpha, it is taken as zero in the 
*       code below. Where there is no nugget parameter lambda, it
*       is taken as 1.
*
***************************************************************************
      INTEGER SPMOD,NSITES,NATTR,MXP,IFAIL
      DOUBLE PRECISION CORMAT(NSITES,NSITES),CORRN(NSITES,NSITES)
      DOUBLE PRECISION SITINF(NSITES,MXP,0:3)
      DOUBLE PRECISION RHO(MXP)
      CHARACTER SCODES(NSITES)*4
***************************************************************************
*       Extra INTEGERs
*       --------------
*       NPAR    Number of parameters in dependence model
*       I,J,P,Q Counters
*       SPM2    Number of parameters in spatial correlation model cc
*               (see above)
*       ITER    Current iteration number in search for optimum
*               parameter set
*       MAXIT   Maximum number of iterations
***************************************************************************
      INTEGER NPAR,I,J,P,Q,SPM2,ITER,MAXIT
***************************************************************************
*       Extra DOUBLEs
*       -------------
*       GRAD    Gradient of the least-squares objective function when 
*               doing Newton-Raphson. Can't imagine there will be any
*               dependence models with more than 20 parameters. 
*       HESS    Hessian of the least-squares objective function
*       ADJ     Vector of Newton adjustments to the current parameter
*               vector
*       STEP    Step size for Newton adjustments (used to prevent
*               iteration from going into illegal regions of
*               parameter space
*       LL, UL  Lower and upper limits for parameters
*       ADJMAX  Largest adjustment on current iteration
*       XSEP }  Matrices of inter-site separations in the first
*       YSEP }  two attributes defined in SITINF
*       DIST    Matrix of Euclidean inter-site distances
*       RBAR    Mean of observed correlations
*       DBAR    Mean inter-site distance
*       SUMN    Sum of the entries in CORRN (used to normalise the
*               calculations for stability)
*       R       Theoretical correlation matrix corresponding to 
*               current entries of RHO and ALPHA
*       RR      Values of the correlation function cc (see above)
*       ALPHA   Value of alpha (see above)
*       FR,     Values of the least squares objective function at 
*       FROLD   current and previous iteration
*       RRGRAD  Derivatives of the correlation function cc at current 
*               site
*       D2RR    Second derivatives ...
*       TOL     Tolerance used to judge convergence of Newton-Raphson
*       TMP?    Temporary storage
***************************************************************************
      DOUBLE PRECISION GRAD(20),HESS(20,20),ADJ(20),ADJMAX
      DOUBLE PRECISION STEP,LL(20),UL(20)
      DOUBLE PRECISION XSEP(NSITES,NSITES),YSEP(NSITES,NSITES)
      DOUBLE PRECISION R(NSITES,NSITES),RR(NSITES,NSITES)
      Double precision ALPHA, LAMBDA
      DOUBLE PRECISION FR,FROLD,RRGRAD(20),D2RR(20,20)
      DOUBLE PRECISION DIST(NSITES,NSITES),RBAR,DBAR,SUMN,TOL
      Double precision TMP1,TMP2,TMP3
******************************************************************************
*       Additional CHARACTERs
*       ^^^^^^^^^^^^^^^^^^^^^
*       MESSAGE - messages to screen
******************************************************************************
      CHARACTER MESSAGE(4)*255
***************************************************************************
      DATA (LL(I),I=1,20) /20*-1.0D6/
      DATA (UL(I),I=1,20) /20*1.0D6/      

*
*	First calculate mean inter-site correlation and distance
*       (not needed for all models but it's quick)
*
      RBAR =0.0D0
      DBAR = 0.0D0
      SUMN=0.0D0
      MAXIT = 100
      FROLD = -1.0D100
      ITER = 0
      DO 100 I=1,NSITES-1
       DO 101 J=I+1,NSITES
        XSEP(I,J) = 0.0D0
        YSEP(I,J) = 0.0D0
        IF (NATTR.GE.1) XSEP(I,J)=SITINF(I,1,0)-SITINF(J,1,0)
        IF (NATTR.GE.2) YSEP(I,J)=SITINF(I,2,0)-SITINF(J,2,0)
        DIST(I,J) = DSQRT( (XSEP(I,J)*XSEP(I,J)) + 
     +                     (YSEP(I,J)*YSEP(I,J)) )
        IF (CORRN(I,J).GT.0.0D0) THEN
         RBAR = RBAR + (CORRN(I,J)*CORMAT(I,J))
         DBAR = DBAR + (CORRN(I,J)*DIST(I,J))
         SUMN = SUMN + CORRN(I,J)
        ENDIF
 101   CONTINUE
 100  CONTINUE
      RBAR = RBAR / SUMN
      DBAR = DBAR / SUMN
*
*       Constant inter-site correlation
*
      IF (SPMOD.EQ.2) THEN
       RHO(1) = RBAR
*
*       For the remaining correlation functions, use Newton-Raphson
*       to minimise the weighted sum of squares for the nonlinear
*       part of the model. N-R doesn't work if there's a threshold
*       though, because in this case the equations have a root at
*       ALPHA=1. For threshold models therefore, do an inner N-R
*       iteration for the parameters of the fundamental correlation 
*       function; and then calculate the value of ALPHA 
*       corresponding to this.
* 
      ELSE IF ( (SPMOD.GE.3).AND.(SPMOD.LE.8) ) THEN
       TOL = 1.0D-5
       ALPHA = 0.0D0
       LAMBDA = 1.0D0
       SPM2 = SPMOD
       IF ((SPMOD.EQ.4).OR.(SPMOD.EQ.6)) THEN
        SPM2 = SPMOD-1
       ELSEIF ((SPMOD.EQ.7).OR.(SPMOD.EQ.8)) THEN
        SPM2 = 3 + (2*(SPMOD-7))
       ENDIF
*
*       Initial values for the fundamental correlation functions,
*       if they haven't been set by the user. Exponential is initialised 
*       by equating mean observed correlation with fitted value at 
*       average inter-site distance. Powered exponential the same. 
*
       IF (SPM2.EQ.3) THEN
        NPAR = 1
        IF (DABS(RHO(1)).LT.1.0D6) THEN
         IF (RBAR.GT.0.0D0) THEN 
          RHO(1) = -DLOG(RBAR) / DBAR
         ELSE
          RHO(1) = 1.0D1
         ENDIF
         LL(1) = 0.0D0
        ENDIF
       ELSEIF (SPM2.EQ.5) THEN
        NPAR = 2
        IF (DABS(RHO(1)).LT.1.0D6) THEN
         IF (RBAR.GT.0.0D0) THEN 
          RHO(1) = -DLOG(RBAR) / DBAR
         ELSE
          RHO(1) = 1.0D1
         ENDIF
         IF (DABS(RHO(2)).LT.1.0D6) RHO(2) = 1.0D0
         LL(1) = 0.0D0
         LL(2) = 0.0D0
        ENDIF
       ENDIF
*
*       For model with non-zero asymptote or nugget, initialise 
*       the asymptote / nugget parameter using analytical solution. 
*       (NB there is no model with *both* asymptote and nugget at
*       present). For models with neither parameter, calling CORSET 
*       twice is unnecessary but not expensive; to avoid this would 
*       make the code rather cumbersome
*
 130   ITER=ITER+1
       CALL CORSET(0,SPM2,RHO,RR,NSITES,
     +                               SITINF,NATTR,SCODES,MXP,IFAIL)
       IF (IFAIL.NE.0) RETURN
       IF ( (SPMOD.EQ.4).OR.(SPMOD.EQ.6) ) THEN
        TMP1 = 0.0D0
        TMP2 = 0.0D0
        DO 140 I=1,NSITES-1
         DO 141 J=I+1,NSITES
          TMP1 = TMP1 + (CORRN(I,J)*(CORMAT(I,J)-(LAMBDA*RR(I,J)))*
     +                                              (1.0D0-RR(I,J)))
          TMP2 = TMP2 + (CORRN(I,J)*(1.0D0-RR(I,J))*(1.0D0-RR(I,J)))
 141     CONTINUE
 140    CONTINUE
        ALPHA = TMP1 / TMP2
        RHO(NPAR+1) = ALPHA
       ELSE IF ( (SPMOD.EQ.7).OR.(SPMOD.EQ.8) ) THEN
        TMP1 = 0.0D0
        TMP2 = 0.0D0
        DO 145 I=1,NSITES-1
         DO 146 J=I+1,NSITES
          TMP3 = RR(I,J) + (ALPHA*(1.0D0-RR(I,J)))
          TMP1 = TMP1 + (CORRN(I,J)*CORMAT(I,J)*TMP3)
          TMP2 = TMP2 + (CORRN(I,J)*TMP3*TMP3)
 146     CONTINUE
 145    CONTINUE
        LAMBDA = TMP1 / TMP2
        IF (LAMBDA.GT.1.0D0) LAMBDA = 1.0D0
        IF (LAMBDA.LT.0.0D0) LAMBDA = 0.0D0
        RHO(NPAR+1) = LAMBDA
       ENDIF
       CALL CORSET(0,SPMOD,RHO,R,NSITES,
     +                              SITINF,NATTR,SCODES,MXP,IFAIL)
       IF (IFAIL.NE.0) RETURN
*
*       Set up the gradient and hessian vector
*
       FR = 0.0D0
       DO 150 P = 1,NPAR
        GRAD(P) = 0.0D0
        DO 151 Q=1,NPAR
         HESS(P,Q) = 0.0D0
 151    CONTINUE
 150   CONTINUE

       DO 160 I=1,NSITES-1
        DO 161 J=I+1,NSITES
         FR = FR + (CORRN(I,J) * ((CORMAT(I,J)-RR(I,J))**2))
         IF (SPM2.EQ.3) THEN
          RRGRAD(1) = -DIST(I,J)*RR(I,J)
          D2RR(1,1) = DIST(I,J)*DIST(I,J)*RR(I,J)
         ELSE IF (SPM2.EQ.5) THEN
          TMP1 = DIST(I,J)**RHO(2)
          TMP2 = TMP1*RR(I,J)
          RRGRAD(1) = -TMP2
          RRGRAD(2) = -TMP2*RHO(1)*DLOG(DIST(I,J))
          D2RR(1,1) = TMP2*TMP1
          TMP1 = (RHO(1)*TMP1) - 1.0D0
          D2RR(1,2) = TMP2 * TMP1 * DLOG(DIST(I,J))
          D2RR(2,1) = D2RR(1,2)
          D2RR(2,2) = TMP2 * RHO(1) * TMP1 * (DLOG(DIST(I,J))**2)
         ENDIF
         DO 165 P=1,NPAR
          GRAD(P) = GRAD(P) - (CORRN(I,J) * ( CORMAT(I,J)-R(I,J) ) *
     +                                                 RRGRAD(P) )
          DO 166 Q=P,NPAR
           HESS(P,Q) = HESS(P,Q) - ( CORRN(I,J) * (
     +                       ( (CORMAT(I,J)-R(I,J))*D2RR(P,Q) )
     +                       - (LAMBDA * (1.0D0-ALPHA) * 
     +                             RRGRAD(P)*RRGRAD(Q) ) ) )
 166      CONTINUE
 165     CONTINUE
 161    CONTINUE
 160   CONTINUE
*
*       Scaling by LAMBDA and (1-ALPHA) works in all cases because 
*       LAMBDA is 1 and ALPHA is 0 for models where they're not needed
*
       DO 170 P = 1,NPAR
        GRAD(P) = LAMBDA * (1.0D0-ALPHA) * GRAD(P) / SUMN
        DO 171 Q=P,NPAR
         HESS(P,Q) = LAMBDA * (1.0D0-ALPHA) * HESS(P,Q) / SUMN
         IF (Q.NE.P) HESS(Q,P) = HESS(P,Q)
 171    CONTINUE
 170   CONTINUE
*
*       Invert Hessian and compute adjustments. Use STEP to prevent
*       parameters from straying out of legal boundaries
* 
       CALL MATINV(HESS,NPAR,1,20,0,IFAIL)
       IF (IFAIL.NE.0) THEN
        IF (IFAIL.GT.0) THEN
         IFAIL = 67
        ELSE
         IFAIL = 997
        ENDIF
        RETURN
       ENDIF
       STEP = 1.0D0
       ADJMAX = 0.0D0
       DO 180 P=1,NPAR
        ADJ(P) = 0.0D0
        DO 181 Q=1,NPAR
         ADJ(P) = ADJ(P) + ( HESS(P,Q) * GRAD(Q) ) 
         IF (DABS(ADJ(P)).GT.ADJMAX) ADJMAX = DABS(ADJ(P))
 181    CONTINUE
        STEP = DMIN1(STEP,0.5D0*DABS( (RHO(P)-LL(P)) / ADJ(P) ),
     +                     0.5D0*DABS( (UL(P)-RHO(P)) / ADJ(P) ))
 180   CONTINUE
       DO 185 P=1,NPAR
        RHO(P) = RHO(P) - (STEP*ADJ(P))
 185   CONTINUE
*
*       Convergence criterion: stop either when adjustments or relative 
*       changes in objective function are small
*  
       IF ((ADJMAX.GT.TOL).AND.
     +     (ABS((FR-FROLD)/FROLD).GT.TOL)) THEN
        FROLD = FR
        IF (ITER.LT.MAXIT) GOTO 130
        WRITE(MESSAGE,2) ADJMAX
        Do I=1,4
         CALL INTPR(TRIM(MESSAGE(I)),-1,0,0)
        End Do
       ENDIF
      ELSE
       WRITE(MESSAGE(1),1)
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       IFAIL = 999
       RETURN
      ENDIF

 1    FORMAT('****ERROR**** Invalid value of SPMOD in routine CORFIT')
 2    FORMAT('****WARNING**** Iteration limit reached while ',
     +'searching for optimal',/'spatial correlation parameters: ',
     +'largest standardised objective function ',/'gradient is ',
     +F7.5,'. You should check the fit graphically; if inadequate ',
     +/'try different starting values.')

      END
