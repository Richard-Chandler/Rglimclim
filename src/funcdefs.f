**************************************************************************
*       © UCL 2002-2010
**************************************************************************
      BLOCK DATA DSCTXT
**************************************************************************
*       Initializes descriptive text relating to predictors etc. Whenever 
*       a new function or predictor is added, the corresponding entry 
*       should be updated in this block so that program output corresponds 
*       to what the thing is actually doing.
**************************************************************************
      COMMON /DESCRIBE/ SXPTXT,MOTXT,TRNDTXT,TRPTXT,DYTXT,
     +                  SPATXT,SPPTXT,CORTXT,PWTTXT,GLBTXT
      INTEGER MXP,I,J

      PARAMETER (MXP=80)

      CHARACTER*70 MOTXT(MXP),TRNDTXT(MXP)
      CHARACTER*70 DYTXT(0:MXP),PWTTXT(MXP,3)
      CHARACTER*70 TRPTXT(MXP,3),SXPTXT(MXP,3)
      CHARACTER*70 SPATXT(0:MXP),SPPTXT(MXP,4),CORTXT(MXP)
      CHARACTER*70 GLBTXT(MXP)
 
      DATA ((SXPTXT(I,J),J=1,3),I=1,3) /
     +'Parameter of transformation','NOT USED','NOT USED',
     +'Value of a','NOT USED','NOT USED',
     +'Location parameter','Scale parameter','NOT USED'/
      DATA ((SXPTXT(I,J),I=11,30),J=1,3) /
     +20*'Lower limit for Fourier representation',
     +20*'Upper limit for Fourier representation',20*'NOT USED'/
      DATA ((SXPTXT(I,J),I=31,40),J=1,3) /
     +10*'Lower limit for polynomial representation',
     +10*'Upper limit for polynomial representation',10*'NOT USED'/

      DATA ((PWTTXT(I,J),J=1,3),I=2,3) /
     +'Exponential decay rate','NOT USED','NOT USED',
     +'Offset for site attribute 1','Offset for site attribute 2',
     +'Exponential decay rate'/

      DATA (MOTXT(I),I=1,8) /
     +'Monthly annual cycle, cosine component',
     +'Monthly annual cycle, sine component',
     +'First harmonic of monthly annual cycle, cosine component',
     +'First harmonic of monthly annual cycle, sine component',
     +'Second harmonic of monthly annual cycle, cosine component',
     +'Second harmonic of monthly annual cycle, sine component',
     +'Third harmonic of monthly annual cycle, cosine component',
     +'Third harmonic of monthly annual cycle, sine component'/
      DATA (MOTXT(I),I=11,22) /
     +'January indicator','February indicator','March indicator',
     +'April indicator','May indicator','June indicator',
     +'July indicator','August indicator','September indicator',
     +'October indicator','November indicator','December indicator' /

      DATA (TRNDTXT(I),I=1,3) /
     +'Linear (0.1 per year, 0 in 1950)',
     +'Linear (0.1 per year) after year Y',
     +'Cyclical'/

      DATA ((TRPTXT(I,J),J=1,3),I=1,3) /
     +'NOT USED','NOT USED','NOT USED',
     +'Date of changepoint (year)','NOT USED','NOT USED',
     +'Cycle length (years)','Date of cycle minimum','NOT USED'/
*
*       NB ### is used as a generic marker for a variable name - will be
*       substituted by the R code later on.
*
      DATA (DYTXT(I),I=0,10) /
     +'###[t]',
     +'###[t-1]',
     +'###[t-2]',
     +'###[t-3]',
     +'###[t-4]',
     +'###[t-5]',
     +'###[t-6]',
     +'###[t-7]',
     +'###[t-8]',
     +'###[t-9]',
     +'###[t-10]'/
      DATA (DYTXT(I),I=21,28) /
     +'Daily annual cycle, cosine component',
     +'Daily annual cycle, sine component',
     +'First harmonic of daily annual cycle, cosine component',
     +'First harmonic of daily annual cycle, sine component',
     +'Second harmonic of daily annual cycle, cosine component',
     +'Second harmonic of daily annual cycle, sine component',
     +'Third harmonic of daily annual cycle, cosine component',
     +'Third harmonic of daily annual cycle, sine component'/
      DATA (DYTXT(I),I=31,42) /
     +'Smooth January effect','Smooth February effect',
     +'Smooth March effect','Smooth April effect',
     +'Smooth May effect','Smooth June effect',
     +'Smooth July effect','Smooth August effect',
     +'Smooth September effect','Smooth October effect',
     +'Smooth November effect','Smooth December effect' /
*
*	For spatial structure, default equates to independence.
*	Positions 1-20 are used for amounts model (correlation 
*	structure of Anscombe residuals), 21-40 are for wet/dry 
*	model.
*
      DATA (SPATXT(I),I=0,8) /'Independence',
     +'Empirical correlations for ',
     +'Constant inter-site correlations for',
     +'Exponential correlation function: exp[-phi*d] for',
     +'Thresholded exp. corr. fn: a + (1-a)exp[-phi*d] for',
     +'Powered exp. corr. fn: exp[-phi*(d^k)] for',
     +'Thresholded pow. exp. corr. fn: a + (1-a)exp[-phi*(d^k)] for',
     +'Exponential corr. fn. with nugget: L*exp[-phi*(d^k)] for',
     +'Pow. exp. corr. fn with nugget: L*exp[-phi*(d^k)] for'/
      DATA SPPTXT(2,1) /'Inter-site correlation'/
      DATA SPPTXT(3,1) /'Correlation decay rate, phi'/
      DATA (SPPTXT(4,I),I=1,2) /'Correlation decay rate, phi',
     +              'Limiting correlation at large distances, a'/
      DATA (SPPTXT(5,I),I=1,2) /'Correlation decay rate, phi',
     +                          'Power of distance, k'/
      DATA (SPPTXT(6,I),I=1,3) /'Correlation decay rate, phi',
     +              'Power of distance, d',
     +              'Limiting correlation at large distances, a'/
      DATA (SPPTXT(7,I),I=1,2) /'Correlation decay rate, phi',
     +              'Limiting correlation at zero distance, L'/
      DATA (SPPTXT(8,I),I=1,3) /'Correlation decay rate, phi',
     +              'Power of distance, d',
     +              'Limiting correlation at zero distance, L'/
      DATA (SPATXT(I),I=21,22) /
     +'Conditional independence given ''wet/dry'' weather state',
     +'Beta-Binomial distribution for number of wet sites'/
      DATA SPPTXT(21,1) /'Increase in logit on a ''wet'' day'/
      DATA (SPPTXT(22,I),I=1,3) /
     +'Shape parameter of Beta-Binomial distribution',
     +'Emergency shrinkage parameter',
     +'Print errors? (0=no, else yes)'/
*
*	Text describing the quantities to which correlation 
*       structures correspond
*
      DATA CORTXT(1) /'latent Gaussian variables'/
      DATA CORTXT(10) /'residuals'/
      DATA CORTXT(11) /'Anscombe residuals'/
      DATA CORTXT(20) /'Pearson residuals'/
*
*     Text for global quantities
*
      DATA (GLBTXT(I),I=1,1) /'Threshold for small positive values'/

      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE ATTRXFM(SITINF,NSITES,ATTRTXT,NATTR,CODES,XFM,
     +                   FOUIDX,LEGIDX,THETA,NP,ICHECK,IFAIL,MXP)
*
*       Defines nonlinear functions of site attributes. 
*	ARGUMENTS:
*	^^^^^^^^^^
*       SITINF        - array of properties for each site. DOUBLE
*			  PRECISION, input/output.
*       NSITES          - no. of sites. INTEGER, input.
*       ATTRTXT         - text describing attributes at each site.
*			  CHARACTER, input/output.
*       NATTR           - number of available attributes. INTEGER,
*			  input/output.
*	CODES		- vector of covariate codes for the model.
*			  INTEGER, input.
*	XFM		- the choice of transformation. INTEGER, input
*	FOUIDX		- Indices for positions in THETA where lower
*			  limits of windows for Fourier representations
*			  are defined. FOUIDX(J) references the THETA
*			  position corresponding to attribute J. INTEGER,
*			  input.
*	LEGIDX		- Ditto, windows for Legendre representations.
*	THETA		- parameters in nonlinear transformations.
*			  DOUBLE PRECISION, input.
*	NP		- no. of parameters in the model. INTEGER, input.
*       ICHECK          - Indicator for whether parametrisation has been
*                         checked. INTEGER, input.
*       IFAIL           - Exit status: 0 is OK, anything else is an 
*                         error or warning
*	MXP		- Used in dimensioning other arguments. INTEGER,
*			  input.
******************************************************************************
      INTEGER NP(10),ICHECK,MXP
      INTEGER NSITES,NATTR,CODES(MXP),XFM(MXP,2)
      INTEGER FOUIDX(MXP),LEGIDX(MXP)
      INTEGER, intent(inout) :: IFAIL
      DOUBLE PRECISION SITINF(NSITES,MXP,0:3),THETA(MXP,3)
      CHARACTER*70 ATTRTXT(MXP)
******************************************************************************
*	Additional INTEGERs
*	^^^^^^^^^^^^^^^^^^^
*	NXATTR		- number of transformed attributes
*	FREQNO		- Frequency for Fourier transform
*	PIDEF		- Indicator for whether we've calculated PI
*	DEGREE		- Degree of polynomial being calculated
*	I,J,K		- Counters
******************************************************************************
      INTEGER NXATTR,FREQNO,PIDEF,DEGREE
      INTEGER I,J,K
******************************************************************************
*	Additional DOUBLEs
*	^^^^^^^^^^^^^^^^^^
*	PI		- -i log(-1)
*	X		- value of normalized attribute for orthogonal
*			  xfrmations
*	LPK   }		- Legendre polynomial values for degree K, K-1
*	LPKM1 }		  and K-2 (for recursions)
*	LPKM2 }
*	TMP*		- Temporary storage
******************************************************************************
      DOUBLE PRECISION PI,X,LPK,LPKM1,LPKM2,TMP1,TMP2
******************************************************************************
*	Additional CHARACTERs
*	^^^^^^^^^^^^^^^^^^^^^
*       MESSAGE - For error messages
******************************************************************************
      CHARACTER*80 MESSAGE(2)

      SAVE PI,PIDEF
      DATA PIDEF /0/
	  
*
*	Compute PI if we haven't done so already
*
      IF (PIDEF.EQ.0) THEN
       PI = 4.0D0*DATAN(1.0D0)
       PIDEF = 1
      ENDIF
      NXATTR = NATTR
      IFAIL = 0
	  
*
*	We go through the site-specific covariates one at a time, 
*	checking to see if there's a corresponding entry in XFM.
*	If so we calculate the appropriate transformation of that
*	covariate, store as an extra site attribute and re-index
*	the covariate codes so they point to the right place. Also,
*	we compute the derivatives wrt any unknown parameters. Note
*	that XFM(1,.) will always contain the number of the 
*	untransformed site attribute; CODES will point to the
*	transformed attribute. NB we always check that we've
*	defined the right number of parameters for the transformation
*	requested.
*
      DO 500 J=1,NP(1)
       IF (XFM(J,2).EQ.0) THEN
        IF (ICHECK.EQ.0) CALL THTCHK(THETA,J,0,MXP,IFAIL)
        IF (IFAIL.NE.0) RETURN
        GOTO 500
       ENDIF
       NXATTR = NXATTR + 1
*
*       1: Box-Cox power family (1 parameter -> 1 derivative). Have to
*	treat it separately when the parameter is zero.
*
       IF (XFM(J,2).EQ.1) THEN
        IF (ICHECK.EQ.0) CALL THTCHK(THETA,J,1,MXP,IFAIL)
        IF (IFAIL.NE.0) RETURN        
        ATTRTXT(NXATTR) = 'Box-Cox power transform of '//
     +                                     TRIM(ATTRTXT(XFM(J,1)))
        DO 100 I=1,NSITES
         IF (SITINF(I,XFM(J,1),0).LE.0.0D0) GOTO 920
         IF (DABS(THETA(J,1)).GT.0.0D0) THEN
          SITINF(I,NXATTR,0) = 
     +    ( (SITINF(I,XFM(J,1),0)**THETA(J,1)) - 1.0D0 )/THETA(J,1)
          SITINF(I,NXATTR,1) = ( 1.0D0 + (
     +          (SITINF(I,XFM(J,1),0)**THETA(J,1))* (
     +          (THETA(J,1)*DLOG(SITINF(I,XFM(J,1),0))) - 1.0D0
     +                                                )
     +                         )         ) / 
     +          (THETA(J,1)**2)
         ELSE
          SITINF(I,NXATTR,0) = DLOG(SITINF(I,XFM(J,1),0))
          SITINF(I,NXATTR,1) = 0.5D0*
     +                    (DLOG(SITINF(I,XFM(J,1),0))**2)
         ENDIF
 100    CONTINUE
*
*       2: e^{ax} (1 parameter -> 1 derivative)
*
       ELSEIF (XFM(J,2).EQ.2) THEN
        IF (ICHECK.EQ.0) CALL THTCHK(THETA,J,1,MXP,IFAIL)
        IF (IFAIL.NE.0) RETURN
        ATTRTXT(NXATTR) = 'exp[ax] transform of '//
     +                              TRIM(ATTRTXT(XFM(J,1)))
        DO 120 I=1,NSITES
         SITINF(I,NXATTR,0) = 
     +              DEXP(THETA(J,1)*SITINF(I,XFM(J,1),0))
         SITINF(I,NXATTR,1) = SITINF(I,XFM(J,1),0) * 
     +                                 SITINF(I,NXATTR,0)
 120    CONTINUE
*
*       3: arctan transform (2 parameters -> 2 derivatives)
*
       ELSEIF (XFM(J,2).EQ.3) THEN
        IF (ICHECK.EQ.0) CALL THTCHK(THETA,J,2,MXP,IFAIL)
        IF (IFAIL.NE.0) RETURN
        ATTRTXT(NXATTR) = 'Arctan-transformed '//
     +                                      TRIM(ATTRTXT(XFM(J,1)))
        DO 140 I=1,NSITES
         SITINF(I,NXATTR,0) = 
     +           DATAN((SITINF(I,XFM(J,1),0)-THETA(J,1))/THETA(J,2))
         SITINF(I,NXATTR,1) = -THETA(J,2)/
     +    ( (THETA(J,2)**2)+((SITINF(I,XFM(J,1),0)-THETA(J,1))**2) )
         SITINF(I,NXATTR,2) = (THETA(J,1)-SITINF(I,XFM(J,1),0))/
     +    ( (THETA(J,2)**2)+((SITINF(I,XFM(J,1),0)-THETA(J,1))**2) )
 140    CONTINUE
*
*	11-30: Fourier series representation (parameters always
*	assumed known -> 0 derivatives). NB the idea here is for
*	orthogonal series expansions, so we only need 1 parameter
*	set (stored in THETA(FOUIDX(XFM(J,1)),*). Treat odds (sine
*	terms) and evens (cos terms) separately
*
       ELSEIF ( (XFM(J,2).GE.11).AND.(XFM(J,2).LE.30) ) THEN
        IF (FOUIDX(XFM(J,1)).EQ.0) GOTO 910
        IF (ICHECK.EQ.0) THEN
         CALL THTCHK(THETA,FOUIDX(XFM(J,1)),2,MXP,IFAIL)
         IF (IFAIL.NE.0) RETURN
        ENDIF
*
*	Integer arithmetic gives sequence 1,1,2,2,3,3,4,4 etc.
*
        FREQNO = (XFM(J,2)-9)/2
        TMP1 = 2.0D0*PI*DBLE(FREQNO)
        TMP2 = THETA(FOUIDX(XFM(J,1)),2) - THETA(FOUIDX(XFM(J,1)),1)
        IF ((2*FREQNO)+9.EQ.XFM(J,2) ) THEN
         WRITE(ATTRTXT(NXATTR),5) FREQNO,ATTRTXT(XFM(J,1))
         DO 160 I=1,NSITES
          X = SITINF(I,XFM(J,1),0)-THETA(FOUIDX(XFM(J,1)),1)
          X = X/TMP2
          SITINF(I,NXATTR,0) = DSIN( TMP1*X )
 160     CONTINUE
        ELSE
         WRITE(ATTRTXT(NXATTR),6) FREQNO,ATTRTXT(XFM(J,1))
         DO 161 I=1,NSITES
          X = SITINF(I,XFM(J,1),0)-THETA(FOUIDX(XFM(J,1)),1)
          X = X/TMP2
          SITINF(I,NXATTR,0) = DCOS( TMP1*X )
 161     CONTINUE
        ENDIF
*
*	31-40: Legendre polynomial representation (parameters always
*	assumed known -> 0 derivatives). Basically same idea as Fourier,
*	but interval of representation is (-1,1) instead of (0,1). 
*	Polynomials are evaluated using recursions rather than via
*	tables of coefficients.
*
       ELSEIF ( (XFM(J,2).GE.31).AND.(XFM(J,2).LE.40) ) THEN
        IF (LEGIDX(XFM(J,1)).EQ.0) GOTO 910
        IF (ICHECK.EQ.0) THEN
         CALL THTCHK(THETA,LEGIDX(XFM(J,1)),2,MXP,IFAIL)
         IF (IFAIL.NE.0) RETURN
        ENDIF
        TMP1 = THETA(LEGIDX(XFM(J,1)),2) + THETA(LEGIDX(XFM(J,1)),1)
        TMP2 = THETA(LEGIDX(XFM(J,1)),2) - THETA(LEGIDX(XFM(J,1)),1)
        DEGREE = XFM(J,2)-30
        WRITE(ATTRTXT(NXATTR),7) DEGREE,ATTRTXT(XFM(J,1))
        DO 180 I=1,NSITES
         X = ( (2.0D0*SITINF(I,XFM(J,1),0)) - TMP1 ) / TMP2
         IF (DEGREE.EQ.1) THEN         
          SITINF(I,NXATTR,0) = X
         ELSE
          LPKM2 = 1.0D0
          LPKM1 = X
          DO 185 K=2,DEGREE
           LPK = ( DBLE((2*K)-1)*X*LPKM1 ) - ( DBLE(K-1)*LPKM2 )
           LPK = LPK / DBLE(K)
           LPKM2 = LPKM1
           LPKM1 = LPK
 185      CONTINUE
          SITINF(I,NXATTR,0) = LPK
         ENDIF
 180    CONTINUE
       ENDIF

       CODES(J) = NXATTR

 500  CONTINUE

      RETURN
*
*	Error trapping
*
 910  IFAIL = 37
      RETURN
 920  IFAIL = 41
      WRITE(MESSAGE,3) XFM(J,1),I
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)      

 3    FORMAT('****ERROR**** you have asked for a Box-Cox transform ',
     +'of site attribute ',I2,',',/,'but the value of this ',
     +'attribute is not strictly positive at site ',I3,'.') 
 5    FORMAT('Fourier sine component ',I2,' for ',A40)
 6    FORMAT('Fourier cosine component ',I2,' for ',A38)
 7    FORMAT('Legendre polynomial ',I2,' for ',A43)
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE TRENDSB(TRENDFN,YEAR,COVNO,CHOICE,THETA,
     +                                            ICHECK,MXP,IFAIL)
*
*       Returns the value of a long-term trend in year YEAR,
*	parameterised by the appropriate elements of THETA, which
*	are in row COVNO. Also check (as far as possible) that
*	parameters have been correctly specified. The 0th element
*	of TRENDFN contains the value of the function, the
*	1st and second contain values of derivatives wrt any
*	nonlinear parameters we're trying to estimate. TMP
*	is temporary storage for saving time.
*
      INTEGER YEAR,CHOICE,MXP,COVNO,ICHECK,IFAIL
      DOUBLE PRECISION TRENDFN(0:2),PI,THETA(MXP,3),TMP
      CHARACTER*255 MESSAGE
 
      IFAIL = 0
      PI = 4.0D0*DATAN(1.0D0)
      IF (CHOICE.EQ.1) THEN
*
*	Linear trend - no parameters
*
       IF (ICHECK.EQ.0) CALL THTCHK(THETA,COVNO,0,MXP,IFAIL)
       IF (IFAIL.NE.0) RETURN
       TRENDFN(0) = (DBLE(YEAR)-1950.0D0)/1.0D1
       TRENDFN(1) = -1.0D-1
       TRENDFN(2) = 0.0D0
      ELSEIF (CHOICE.EQ.2) THEN
*
*	Piecewise linear trend (1 parameter -> 1 derivative)
*
       IF (ICHECK.EQ.0) CALL THTCHK(THETA,COVNO,1,MXP,IFAIL)
       IF (IFAIL.NE.0) RETURN
       IF (DBLE(YEAR).LE.THETA(COVNO,1)) THEN
        TRENDFN(0) = 0.0D0
        TRENDFN(1) = 0.0D0
       ELSE
        TRENDFN(0) = (DBLE(YEAR) - THETA(COVNO,1))/1.0D1
        TRENDFN(1) = -1.0D-1
       ENDIF
      ELSEIF (CHOICE.EQ.3) THEN
*
*	Cycle (2 parameters -> 2 derivatives)
*
       IF (ICHECK.EQ.0) CALL THTCHK(THETA,COVNO,2,MXP,IFAIL)
       IF (IFAIL.NE.0) RETURN      
       TMP = 2.0D0*PI*(DBLE(YEAR)-THETA(COVNO,2))/THETA(COVNO,1)
       TRENDFN(0) = -DCOS(TMP)
       TRENDFN(1) = -DSIN(TMP)*TMP/THETA(COVNO,1)
       TRENDFN(2) = -2.0D0*PI*DSIN(TMP)/THETA(COVNO,1)
      ELSE
       WRITE(MESSAGE,1)
       CALL INTPR(TRIM(MESSAGE),-1,0,0)
       IFAIL = 300
      ENDIF

 1    FORMAT('****ERROR**** Illegal trend requested in TRENDFN')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE EXTXST(FILNO,TSCALE,PREDNO,LAG,LABEL,IFAIL)
*
*     To set descriptive text for external predictors. Arguments:
*
*     FILNO   - channel number of (open) file to read. Input
*     TSCALE  - Time scale (1=annual,2=minthly,3=daily). Input
*     PREDNO  - Number of predictor in the file. Input
*     LAG     - Lag of predictor value (in units of TSCALE) relative
*               to current observation.
*     LABEL   - Descriptive text. Output.
*     IFAIL   - Error flag
******************************************************************************
      INTEGER FILNO,TSCALE,PREDNO,LAG,IFAIL
      CHARACTER*70 LABEL
******************************************************************************
*     Additional INTEGERs
*     ^^^^^^^^^^^^^^^^^^^
*     NHEAD   - Numbers of header rows to skip in each type of file
*     NPRDEF  - Number of predictors defined in the file
*     CURLIN  - current line number
*     I       - counter
******************************************************************************
      INTEGER NHEAD(3),NPRDEF,CURLIN,I
******************************************************************************
*     Additional CHARACTERs
*     ^^^^^^^^^^^^^^^^^^^^^
*     SCLTXT  - Text describing the time scale - for error message output
*     TMPSTR  - For reading/writing labels
*     MESSAGE - Error messages
******************************************************************************
      CHARACTER*8 SCLTXT(3)
      CHARACTER*70 TMPSTR
      CHARACTER*255 MESSAGE(2)

      DATA (NHEAD(I),I=1,3) /39,40,41/
      DATA (SCLTXT(I),I=1,3) /'annual','monthly','daily'/
      
      IFAIL = 0
*
*     Skip over the file header
*
      DO 50 I=1,NHEAD(TSCALE)
       CURLIN = I
       READ(FILNO,*,ERR=98)
 50   CONTINUE
*
*     Read no. of defined predictors, and check that the one 
*     requested is valid.
*
      CURLIN = CURLIN + 1
      READ(FILNO,*,ERR=98) NPRDEF
      IF (PREDNO.GT.NPRDEF) GOTO 99
      IF (NPRDEF.GT.950) GOTO 97
*
*     Skip to the line containing the text we want
*
      DO 60 I=1,PREDNO-1
       CURLIN = CURLIN + 1
       READ(FILNO,*,ERR=98)
 60   CONTINUE
*
*     Read the text ...
*
      CURLIN = CURLIN + 1
      READ(FILNO,'(A70)',ERR=98) LABEL
*
*     If there's a non-zero lag involved, amend text 
*
      IF (LAG.NE.0) THEN
       WRITE(TMPSTR,'(''Lag '',I3,'': '',A61)') LAG,LABEL
       READ(TMPSTR,'(A70)') LABEL
      ENDIF
*
*     ... and reposition the file (not worth thinking about
*     tracking this to improve efficiency - it only gets done a couple 
*     of times during the entire program)
*
      REWIND(FILNO)

      RETURN
*
*     Error trapping
*
 97   WRITE(MESSAGE,7) SCLTXT(TSCALE)
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
 98   WRITE(MESSAGE,5) SCLTXT(TSCALE),CURLIN
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      IFAIL = 301
      RETURN
 99   WRITE(MESSAGE,6) SCLTXT(TSCALE)
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      IFAIL = 301
      RETURN      

 5    FORMAT('****ERROR**** unexpected error in ',A8,
     +       ' input data file at line ',I4,'.')
 6    FORMAT('****ERROR**** you have requested an external ',
     +A8,' predictor',/,
     +'that is not defined in your input data file.')
 7    FORMAT('****ERROR**** you have defined more than 950 external ',
     +A8,' predictors.',/,'I hope you didn''t mean to do this!')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE EXTSET(FILNO,TSCALE,PREDNO,DD,MM,YY,LAG,
     +                                          XVAL,MISSFL,IFAIL)
*
*     To read the value of an 'external forcing' predictor for
*     the current day/month/year, or at some time lagged relative to
*     the current one. NB the routine FileReset should be called 
*     prior to the first call to this routine in any fitting exercise - 
*     this (re)initialises the elements of the ExtFileState common block
*     Arguments:
*
*     FILNO   - channel number of (open) file to read. Input
*     TSCALE  - Time scale (1=annual,2=monthly,3=daily). Input
*     PREDNO  - Number of predictor in the file. Input
*     DD      - day of current month. Input
*     MM      - current month. Input
*     YY      - current year. Input.
*     LAG     - lag relative to the current day/month/year (measured
*               in units of days,months or years depending on value of
*               TSCALE).
*     XVAL    - Predictor value. Output.
*     MISSFL  - Flag to indicate mssing data. Output.
*     IFAIL   - Error flag. Output.
******************************************************************************
      INTEGER FILNO,TSCALE,PREDNO,DD,MM,YY,LAG,MISSFL,IFAIL
      DOUBLE PRECISION XVAL
******************************************************************************
*     Additional INTEGERs
*     ^^^^^^^^^^^^^^^^^^^
*     MAXNPR  - Maximum number of external predictors
*     NHEAD   - Numbers of header rows to skip in each type of file
*     CHKFIL  - Indicators for whether input files have been checked. 
*     NPRDEF  - Number of predictors defined in the file
*     LIN1ID  - Index numbers of the first data line in each file
*     REQIDX  - Index numbers of required data line
*     BEGLIN  - Line number of first data line in each file
*     ENDLIN  - Last line number in each file
*     OLDIDX  - Index number of last row read
*     CURLIN  - Last line number read, in each file
*     REQLIN  - Required line numbers
*     FY,FM,FD- year, month and day from input file
*     I       - counter
******************************************************************************
      INTEGER MAXNPR
      PARAMETER(MAXNPR=950)
      INTEGER NHEAD(3),CHKFIL(3),NPRDEF(3),LIN1ID(3),REQIDX(3)
      INTEGER CURLIN(3),REQLIN(3),BEGLIN(3),ENDLIN(3),OLDIDX
      INTEGER FY,FM,FD,I
******************************************************************************
*     Additional CHARACTERs
*     ^^^^^^^^^^^^^^^^^^^^^
*     SCLTXT  - Text describing the time scale - for error message output
*     MESSAGE - Error messages
******************************************************************************
      CHARACTER*8 SCLTXT(3)
      CHARACTER*255 MESSAGE(2)
******************************************************************************
*     Additional DOUBLEs
*     ^^^^^^^^^^^^^^^^^^
*     TMPARR  - To contain a row of predictors from the input file 
******************************************************************************
      DOUBLE PRECISION TMPARR(MAXNPR,3)

      DATA (NHEAD(I),I=1,3) /39,40,41/
      DATA (SCLTXT(I),I=1,3) /'annual','monthly','daily'/
      
      COMMON /ExtFileState/ CHKFIL,LIN1ID,CURLIN,BEGLIN,ENDLIN,NPRDEF
      SAVE /ExtFileState/
      SAVE TMPARR

      MISSFL = 0
      IFAIL = 0
*
*     Check that temporary storage is adequate
*
      IF (PREDNO.GT.MAXNPR) GOTO 96
*
*     If this is the first time we've read this file, scan through it once
*     to identify the first and last lines, and check that all the rows 
*     are present and in the right order. Other checks were performed when 
*     setting the text for the predictors.
*
      IF (CHKFIL(TSCALE).EQ.0) THEN
       CHKFIL(TSCALE) = 1
       REWIND(FILNO)
       DO 50 I=1,NHEAD(TSCALE)
        READ(FILNO,*) 
 50    CONTINUE
       CURLIN(TSCALE) = NHEAD(TSCALE) + 1
       READ(FILNO,*,ERR=98) NPRDEF(TSCALE)
       DO 60 I=1,NPRDEF(TSCALE)+1
        CURLIN(TSCALE) = CURLIN(TSCALE) + 1
        READ(FILNO,*,ERR=98)
 60    CONTINUE
*
*     Scan through the file to (a) set line numbers of first and last
*     lines (b) check that all rows are present and in the correct order
*
 70    CURLIN(TSCALE) = CURLIN(TSCALE) + 1
       CALL RDXREC(FILNO,TSCALE,FY,FM,FD,NPRDEF(TSCALE),
     +                                             TMPARR,MAXNPR,IFAIL)
       IF (IFAIL.EQ.1) GOTO 98
       IF (BEGLIN(TSCALE).EQ.0) THEN
        BEGLIN(TSCALE) = CURLIN(TSCALE)
        CALL IDXSET(TSCALE,FY,FM,FD,LIN1ID(TSCALE))
        REQIDX(TSCALE) = LIN1ID(TSCALE)
        OLDIDX = LIN1ID(TSCALE) - 1
       ENDIF
       IF (IFAIL.NE.2) THEN
        CALL IDXSET(TSCALE,FY,FM,FD,REQIDX(TSCALE))
        IF (REQIDX(TSCALE).NE.OLDIDX+1) GOTO 97
        OLDIDX = REQIDX(TSCALE)
        GOTO 70
       ENDIF
*
*     Only get here if IFAIL=2 i.e. we've reached EOF. After EOF the file
*     position is indeterminate, so rewind to be sure that pointers etc.
*     are in the right place. And reset IFAIL to 0 because this is *not*
*     an error condition. 
*
       ENDLIN(TSCALE) = CURLIN(TSCALE) - 1
       IFAIL = 0
       IF (ENDLIN(TSCALE).LT.BEGLIN(TSCALE)) GOTO 99
       REWIND(FILNO)
       CURLIN(TSCALE) = 0
      ENDIF
*
*     Right: now compute the index number of the required date,
*     and hence the line number.
*
      CALL IDXSET(TSCALE,YY,MM,DD,REQIDX(TSCALE))
      REQIDX(TSCALE) = REQIDX(TSCALE) - LAG
      REQLIN(TSCALE) = BEGLIN(TSCALE)+REQIDX(TSCALE)-LIN1ID(TSCALE)
*
*     First option: we already know the required line isn't in the file.
*
      IF ( (REQIDX(TSCALE).LT.LIN1ID(TSCALE)).OR.
     +         (REQLIN(TSCALE).GT.ENDLIN(TSCALE)) ) THEN
       XVAL = -9999.9D0
       MISSFL = 1
       RETURN
*
*     Or: it's the same line as last time. If so, copy value and scarper.
*
      ELSEIF (REQLIN(TSCALE).EQ.CURLIN(TSCALE)) THEN
       XVAL = TMPARR(PREDNO,TSCALE)
       IF (DABS(XVAL+9999.9D0).LT.1.0D-4) MISSFL = 1
       RETURN       
      ENDIF
*
*     OK, now if the required line is beyond the current point in the 
*     file, read up to and including the line before it (if possible):
*
      IF (REQLIN(TSCALE).GT.CURLIN(TSCALE)) THEN
       DO 110 I = CURLIN(TSCALE)+1,REQLIN(TSCALE)-1
        READ(FILNO,*)
 110   CONTINUE
*
*     If the requested line is *earlier* than the last thing we read,
*     EITHER backspace OR rewind, depending on which is quicker (NB for
*     backspace: CURLIN holds number of last row read. E.g. if last row
*     read was 5 and we want to read row 3, we need to backspace row 5,
*     row 4 and then row 3 so that we end up reading row 3 again. 
*     CURLIN then gets set to 2).
*
      ELSE
       IF (2*REQLIN(TSCALE).GT.CURLIN(TSCALE)) THEN
        DO 120 I = CURLIN(TSCALE),REQLIN(TSCALE),-1
         BACKSPACE(FILNO)
 120    CONTINUE
       ELSE
        REWIND(FILNO)
        DO 125 I = 1,REQLIN(TSCALE)-1
         READ(FILNO,*)
 125    CONTINUE
       ENDIF
       CURLIN(TSCALE) = REQLIN(TSCALE)-1
      ENDIF
*
*     No need to check any more - file has already been checked on 
*     the first pass
* 
      CURLIN(TSCALE) = REQLIN(TSCALE)
      CALL RDXREC(FILNO,TSCALE,FY,FM,FD,NPRDEF(TSCALE),
     +                                           TMPARR,MAXNPR,IFAIL)
      XVAL = TMPARR(PREDNO,TSCALE)
      IF (DABS(XVAL+9999.9D0).LT.1.0D-4) then
       MISSFL = 1
      ENDIF
      RETURN

*
*     Error trapping
*
 96   WRITE(MESSAGE,4)
      IFAIL = 997
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      RETURN
 97   WRITE(MESSAGE(1),5) CURLIN(TSCALE),SCLTXT(TSCALE)
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      IFAIL = 210
      RETURN
 98   WRITE(MESSAGE(1),6) SCLTXT(TSCALE),CURLIN(TSCALE)
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      IFAIL = 210
      RETURN
 99   WRITE(MESSAGE(1),7) SCLTXT(TSCALE)
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      IFAIL = 210
      RETURN

 4    FORMAT('****ERROR**** Inadequate storage in routine EXTSET ',
     +'(file funcdefs.f).',/,'Please increase value of MAXNPR ',
     +'and recompile.')
 5    FORMAT('****ERROR**** line ',I5,' of ',A8,' input data',
     +' file carries incorrect date.')
 6    FORMAT('****ERROR**** unexpected error in ',A8,
     +       ' input data file at line ',I5,'.')
 7    FORMAT('****ERROR**** no data in ',A8,' input data file.')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE RDXREC(FILNO,TSCALE,FY,FM,FD,NPRDEF,XVALS,MXNP,IFAIL)
*
*     Reads a record of an `external' data file. Arguments:
*
*     FILNO    - Channel number of open file. Input
*     TSCALE   - 1 for annual, 2 for monthly and 3 for daily data. Input
*     FY       - year read from file. Output
*     FM       - month read from file. Output
*     FD       - day read from file. Output
*     NPRDEF   - Number of predictors defined. Input
*     XVALS    - array of predictors read from file. Output
*     MXNP     - Maximum number of predictors. Input.
*     IFAIL    - error flag. Output
******************************************************************************
      INTEGER, intent(in) :: FILNO, TSCALE, NPRDEF, MXNP
	  Integer, intent(out):: FY, FM, FD, IFAIL
      DOUBLE PRECISION, intent(out) :: XVALS(MXNP,3)
      INTEGER I
	  
      FY = 0
      FM = 0
      FD = 0
      IFAIL = 0

      IF (TSCALE.EQ.1) THEN
       READ(FILNO,1,ERR=99,END=98,IOSTAT=IFAIL) 
     +                         FY,(XVALS(I,TSCALE),I=1,NPRDEF)
      ELSEIF (TSCALE.EQ.2) THEN
       READ(FILNO,2,ERR=99,END=98,IOSTAT=IFAIL) 
     +                         FY,FM,(XVALS(I,TSCALE),I=1,NPRDEF)
      ELSEIF (TSCALE.EQ.3) THEN
       READ(FILNO,3,ERR=99,END=98,IOSTAT=IFAIL) 
     +                         FY,FM,FD,(XVALS(I,TSCALE),I=1,NPRDEF)
      ENDIF
      RETURN
*
*     Error and EOF trapping
*
 98   IFAIL = 2
      RETURN
 99   IFAIL = 1
      RETURN

 1    FORMAT(I4,T10,950(F10.0))
 2    FORMAT(I4,1X,I2,T10,950(F10.0))
 3    FORMAT(I4,2I2,1X,950(F10.0))
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE IDXSET(TSCALE,YY,MM,DD,IDX)
*
*     Allocates unique indices to individual years, months or days.
*     Arguments:
*     TSCALE     1 for annual, 2 for monthly and 3 for daily data. Input
*     YY         Year. Input
*     MM         Month. Input
*     DD         Day. Input
*     IDX        Index number (it's the number of time units since
*                1st January 2000 i.e. 1/1/2000 gets an index of zero
*                if TSCALE = 3, January 2000 gets an index of zero if
*                TSCALE = 2, and 2000 gets an index of zero if
*                TSCALE = 3). Output
******************************************************************************
      INTEGER TSCALE,YY,MM,DD,IDX
******************************************************************************
*     Additional INTEGERs
*     ^^^^^^^^^^^^^^^^^^^
*     JULBL      Julian day number of baseline date (daily data only)
*     JULREQ     Julian day number of required date.
*     IFAIL      Return exit code for Julian conversion routines
******************************************************************************
      INTEGER JULBL,JULREQ,IFAIL

*
*     1st Jan 2000 was Julian day 2451545
*
      DATA JULBL /2451545/
      SAVE JULBL

      IF (TSCALE.EQ.1) THEN
       IDX = YY-2000
      ELSEIF (TSCALE.EQ.2) THEN
       IDX = (12*(YY-2000)) + MM - 1
      ELSEIF (TSCALE.EQ.3) THEN
       CALL DMYTOJ(DD,MM,YY,JULREQ,IFAIL)
       IDX = JULREQ - JULBL
      ENDIF

      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE MONSB(MONFN,MONTH,CHOICE,IFAIL)
*
*       On exit, the variable MONFN contains the value of a 
*	monthly effect for month MONTH
*
      DOUBLE PRECISION MONFN,PI
      INTEGER MONTH,CHOICE,PICALC,IFAIL
      CHARACTER MESSAGE*255

      DATA PICALC /0/
      SAVE PICALC,PI

      IFAIL = 0
      IF (PICALC.EQ.0) THEN
       PI = 4.0D0*DATAN(1.0D0)
       PICALC = 1
      ENDIF
*
*	Sine waves
*
      IF (CHOICE.EQ.1) THEN
       MONFN=DCOS(DBLE(MONTH)*PI/6.0D0)
      ELSEIF (CHOICE.EQ.2) THEN
       MONFN=DSIN(DBLE(MONTH)*PI/6.0D0)
      ELSEIF (CHOICE.EQ.3) THEN
       MONFN=DCOS(DBLE(MONTH)*PI/3.0D0)
      ELSEIF (CHOICE.EQ.4) THEN
       MONFN=DSIN(DBLE(MONTH)*PI/3.0D0)
      ELSEIF (CHOICE.EQ.5) THEN
       MONFN=DCOS(DBLE(MONTH)*PI/2.0D0)
      ELSEIF (CHOICE.EQ.6) THEN
       MONFN=DSIN(DBLE(MONTH)*PI/2.0D0)
      ELSEIF (CHOICE.EQ.7) THEN
       MONFN=DCOS(2.0D0*DBLE(MONTH)*PI/3.0D0)
      ELSEIF (CHOICE.EQ.8) THEN
       MONFN=DSIN(2.0D0*DBLE(MONTH)*PI/3.0D0)
*
*	Indicators for individual months
*
      ELSEIF ( (CHOICE.GE.11).AND.(CHOICE.LE.22) ) THEN
       IF (MONTH.EQ.CHOICE-10) THEN
        MONFN = 1.0D0
       ELSE
        MONFN = 0.0D0
       ENDIF
      ELSE
       WRITE(MESSAGE,1)
       CALL INTPR(TRIM(MESSAGE),-1,0,0)
       IFAIL = 300 
       RETURN
      ENDIF

 1    FORMAT('****ERROR**** Illegal choice of predictor in MONFN')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE PRVTXT(COVCODE,XLABEL,IFAIL)
*
*     Produces a text label for a predictor which is a nonlinear
*     transformation of some previous day's values. Arguments:
*     
*     COVCODE   - coding for the transformation, in the form (1000x)+d,
*                 where x is the transformation to be used, and d is 
*                 the lag
*     XLABEL    - The text label. Input/Output
*     IFAIL     - Error flag
******************************************************************************
      INTEGER COVCODE,IFAIL
      CHARACTER XLABEL*70
******************************************************************************
*     Additional INTEGERs
*     ^^^^^^^^^^^^^^^^^^^
*     LAG       - number of days ago
*     XFRM      - index of the transformation
*     WTSCHM    - Weighting scheme, for use when a weighted average
*                 of previous days' values at different sites is
*                 being used
*     WXORDR    - Order in which averaging and transformation are
*                 performed
*     TXTLEN    - length of text label
*     XFMLEN    - length of transformation text (see XFMTXT below)
*     WTXLEN    - length of weighted average text (see WTXT below)
*     I,J       - counters
******************************************************************************
      INTEGER LAG,XFRM,WTSCHM,WXORDR
      INTEGER TXTLEN,XFMLEN(0:10,2),WTXLEN(0:10,2)
      INTEGER I,J
******************************************************************************
*     Additional CHARACTERs
*     ^^^^^^^^^^^^^^^^^^^^^
*     XFMTXT    - Text for transformations - first element is to go
*                 at the beginning, last at the end
*     WTXT      - Text for weighted averages of previous days' values
*     TMPTX?    - Temporary storage
*     MESSAGE   - For outputting messages
******************************************************************************
      CHARACTER*40 XFMTXT(0:10,2),WTXT(0:10,2)
      CHARACTER*40 TMPTX1,TMPTX2,TMPTX3,TMPTX4
      CHARACTER MESSAGE*255

*
*       NB ### is used as a generic marker for a variable name - will be
*       substituted by the R code later on.
*
      DATA ((XFMTXT(I,J),J=1,2),I=0,5) /'','','Ln(',')',
     +'Ln(1+',')','I(','>0)','Trace indicator: ','',
     +'I(###[t-k]>0: k=1 to ',')'/
      DATA ((XFMLEN(I,J),J=1,2),I=0,5) /2*0,3,1,5,1,2,3,17,0,21,1/
      DATA ((WTXT(I,J),J=1,2),I=0,3) /'','','Mean of ','',
     +'Distance-weighted mean of ','',
     +'Shift/distance-weighted mean of ',''/
      DATA ((WTXLEN(I,J),J=1,2),I=0,3) /2*0,8,0,26,0,32,0/

      IFAIL = 0
*
*     Interpret COVCODE (integer arithmetic)
*
      LAG = MOD(COVCODE,1000)
      XFRM = MOD(COVCODE/1000,10)
      WTSCHM = MOD(COVCODE/10000,10)
      WXORDR = MOD(COVCODE/100000,10)
*
*     Find out how long the text is at present
*
      TXTLEN = INDEX(XLABEL,'  ') - 1
      IF (TXTLEN.EQ.-1) TXTLEN = LEN(XLABEL)
*
*     Now set text (default core label needs changing for
*     persistence indicators). First the case when we're
*     doing averages of transformations ...
*
      IF (WXORDR.EQ.0) THEN
       IF (XFRM.EQ.5) THEN
        WRITE(XLABEL,'(I2)') LAG 
        TXTLEN = INDEX(XLABEL,'  ') - 1
       ENDIF

       IF (XFRM.LE.5) THEN
        TMPTX1 = WTXT(WTSCHM,1)
        TMPTX2 = XFMTXT(XFRM,1)
        TMPTX3 = XFMTXT(XFRM,2)
        TMPTX4 = WTXT(WTSCHM,2)
        XLABEL = TMPTX1(1:WTXLEN(WTSCHM,1)) // 
     +       TMPTX2(1:XFMLEN(XFRM,1)) // 
     +       XLABEL(1:TXTLEN) // 
     +       TMPTX3(1:XFMLEN(XFRM,2)) //
     +       TMPTX4(1:WTXLEN(WTSCHM,2))
       ELSE
        WRITE(MESSAGE,1)
        CALL INTPR(TRIM(MESSAGE),-1,0,0)
        IFAIL = 5
        RETURN
       ENDIF
*
*     And now the case where we're doing transformations of averages ...
*
      ELSEIF (WXORDR.EQ.1) THEN
       IF (XFRM.LE.5) THEN
        TMPTX1 = XFMTXT(XFRM,1)
        TMPTX2 = WTXT(WTSCHM,1)
        TMPTX3 = WTXT(WTSCHM,2)
        TMPTX4 = XFMTXT(XFRM,2)
        IF (XFRM.LE.4) THEN
         XLABEL = TMPTX1(1:XFMLEN(XFRM,1)) // 
     +       TMPTX2(1:WTXLEN(WTSCHM,1)) // 
     +       XLABEL(1:TXTLEN) // 
     +       TMPTX3(1:WTXLEN(WTSCHM,2)) //
     +       TMPTX4(1:XFMLEN(XFRM,2))
*
*     `Persistence' requires special treatment to create a sensible
*     label
*
        ELSE
         WRITE(XLABEL,'(I2)') LAG 
         TXTLEN = INDEX(XLABEL,'  ') - 1
*
*       See comment above about ###
*
         XLABEL = 'I(' // TMPTX2(1:WTXLEN(WTSCHM,1)) //
     +          '###[t-k] > 0: k=1 to ' // XLABEL(1:TXTLEN) // ')'
        ENDIF
       ELSE
        WRITE(MESSAGE,1)
        CALL INTPR(TRIM(MESSAGE),-1,0,0)
        IFAIL = 5
        RETURN
       ENDIF
      ELSE
       WRITE(MESSAGE,1)
       CALL INTPR(TRIM(MESSAGE),-1,0,0)
       IFAIL = 5
       RETURN
      ENDIF

 1    FORMAT('****ERROR**** Illegal transformation found in PRVTXT')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE DYSET(DD,MM,DatArray,SITE,NSITES,NVARS,RespIdx,
     +                 AllowIncAvge,TRACE,COVCODE,THETA,PWTIDX,
     +                 Distance,ICHECK,IFAIL,DYPRED,MXP)
******************************************************************************
*       Returns a predictor that varies on a daily timescale. Arguments
*       (all input):
*
*       DD      - Day of month. 
*       MM      - Month of year. 
*       SITE    - Site we're working on.
*       NSITES  - Total number of sites
*       NVARS   - Number of variables in input data file
*       RespIdx - Index number of response variable
*       AllowIncAvge    - Indicator for whether or not weighted averages
*                 of non-response variables will be considered as 
*                 valid in model fitting if the variable is missing at
*                 the current site 
*       DatArray- Array of current and up to 10 previous values at each site. 
*       TRACE   - Trace threshold.
*       COVCODE - AS for PRVTXT above
*       THETA   - Parameters in nonlinear transofrmations (for computing
*                 distance-dependent weights)
*       PWTIDX  - Index locating row of THETA to use for calculating
*                 distance-dependent weights.
*       Distance- array of inter-site distances. First slice is X-separation,
*                 second slice is Y-separation and third slice is Euclidean
*                 separation
*       ICHECK  - Indicator to tell us whether the model parametrisation
*                 has been checked
*       RECALC  - Indicator to tell us whether to recalculate inter-site
*                 weighting matrices or not. Weirdly, a value of zero means
*                 yes (see COVSET header).
*       IFAIL   - Error flag
*       DYPRED  - OUTPUT: value of predictor, together with any
*                 derivatives wrt nonlinear parameters.
*       MXP     - Maximum number of site attributes
******************************************************************************
      INTEGER, intent(in) :: DD,MM,SITE,NSITES,NVARS,RespIdx
      INTEGER, intent(in) :: AllowIncAvge,COVCODE,MXP,ICHECK
      INTEGER IFAIL,PWTIDX(MXP,0:3,NVARS)
      DOUBLE PRECISION DatArray(NSITES,NVARS,0:10),TRACE,DYPRED(0:3)
      DOUBLE PRECISION THETA(MXP,3),Distance(3,NSITES,NSITES)
******************************************************************************
*     Additional INTEGERs
*     ^^^^^^^^^^^^^^^^^^^
*     CODE       - code used to define predictor
*     XFRM       - transformation, if any
*     WTSCHM     - Weighting scheme for averages of previous days' values
*     WXORDR     - Flag indicating order of averaging and transformation
*     VARNUM     - Index number of variable being referenced
*     NPREQ      - Number of parameters required by various weighting schemes
*     NDERIV     - Number of derivatives we want to calculate
*     IDERIV     - Indicator - do we want to calculate derivatives of
*                  transformations?
*     PICALC     - have we calculated PI yet?
*     CUMDAY     - cumulative no. of days to beginning of 
*                  current month (calculate as though it's a leap
*                  year - in non-leap years this will result in a 
*                  tiny discontinuity in sin and cos terms, but the
*                  effect will be miniscule and it saves 
*                  on computation time)
*     NDAYS      - numbers of days in each month (ditto)
*     CURDAY     - number of current day (1-366)
*     MONADJ     - Month for which an adjustment is to be made 
*     I,J        - counters
******************************************************************************
      INTEGER I,J,CODE,XFRM,PICALC,WTSCHM,WXORDR,VARNUM,NPREQ(10)
      INTEGER NDERIV,IDERIV,CUMDAY(12),NDAYS(12),CURDAY,MONADJ
******************************************************************************
*     Additional DOUBLEs
*     ^^^^^^^^^^^^^^^^^^
*     PI         - the usual
*     TRIGX      - argument to a trigonometric function
*     TMP        - temporary storage
*     PREVEC     - Vector (as opposed to matrix) of previous values
*     SUMWT      - Sum of weights in spatial averages
*     CURWT      - Current weight, together with any necessary 
*                  derivatives wrt unknown parameters
*     SUMWT      - Sum of weights and their derivatives
*     WGTARR     - Array used for derivative calculation (the weights 
*                  themselves are quotients - we need to find the sums of
*                  the raw weights and their derivatives before we can
*                  differentiate the normalised weights)
******************************************************************************
      DOUBLE PRECISION PI,TRIGX,TMP(0:1),PREVEC(10)
      DOUBLE PRECISION CURWT(0:3),SUMWT(0:3),WGTARR(NSITES,0:3)
*******************************************************************************
*     Additional CHARACTERs
*     ^^^^^^^^^^^^^^^^^^^^^
*     MESSAGE   - Error messages
*******************************************************************************
      CHARACTER MESSAGE(2)*255

      DATA PICALC /0/
      DATA (CUMDAY(I),I=1,12) /0,31,60,91,121,152,182,
     +                             213,244,274,305,335/
      DATA (NDAYS(I),I=1,12) /31,29,31,30,31,30,31,31,30,31,30,31/
      DATA (NPREQ(I),I=1,3) /0,1,3/
      SAVE PICALC,PI

      IF (PICALC.EQ.0) THEN
       PI = 4.0D0*DATAN(1.0D0)
       PICALC = 1
      ENDIF
*
*     Interpret COVCODE
*
      CODE = MOD(COVCODE,1000)
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*     Daily effects representing autocorrelation or inter-variable dependence
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF (CODE.LE.10) THEN
       XFRM = MOD(COVCODE/1000,10)
       WTSCHM = MOD(COVCODE/10000,10)
       WXORDR = MOD(COVCODE/100000,10)
       VARNUM = COVCODE/1000000
       IDERIV = 0
*
*       Missing values - set to large negative number so we can catch
*       them. NB this is allowed for in fitting progs anyway, where 
*	there are flags. Any case is defined as valid if it has the
*       required number of previous days' values at the site of interest,
*       OR if the covariate is a weighted average and the user does not
*       insist on the value being present at the site of interest 
*
       IF (DatArray(SITE,VARNUM,CODE).LT.-1.0D99) THEN
        IF ( (WTSCHM.EQ.0).OR.(AllowIncAvge.EQ.0).OR.
     +                            (VARNUM.EQ.RespIdx)) THEN
         DYPRED(0) = -1.0D100
         RETURN
        ENDIF
       ENDIF
       IF (XFRM.EQ.5) THEN
        DO 50 I=1,CODE
         IF (DatArray(SITE,VARNUM,I).LT.-1.0D99) THEN
          DYPRED(0) = -1.0D100
          RETURN
         ENDIF
 50     CONTINUE
       ENDIF
*
*     For non-missing values, calculate transformations. First the
*     case when we're simply looking at each site's history 
*     individually:
*
       IF (WTSCHM.EQ.0) THEN
        IF (XFRM.EQ.5) THEN
         PREVEC(1:CODE) = DatArray(SITE,VARNUM,1:CODE)
        ENDIF
        CALL PRVXFM(DatArray(SITE,VARNUM,CODE),XFRM,TRACE,CODE,PREVEC,
     +                                          IDERIV,TMP,IFAIL)
        IF (IFAIL.NE.0) RETURN
        DYPRED(0) = TMP(0)
*
*     Now the case where we're averaging. First calculate raw weights
*     (and, if necessary, their partial derivatives). Derivatives
*     are only computed if absolutely necessary, since they're
*     expensive. To see if we need them, look at elements of THETA in
*     'nonlinear' section.
*
       ELSE
        NDERIV = 0
        DO 105 I=0,NPREQ(WTSCHM)
         SUMWT(I) = 0.0D0
         DYPRED(I) = 0.0D0
         IF ((I.GT.0).AND.(PWTIDX(WTSCHM,I,VARNUM).NE.0)) 
     +                                NDERIV = NPREQ(WTSCHM)
 105    CONTINUE
        DO 110 J = 1,NSITES
         IF (DatArray(J,VARNUM,CODE).GE.-1.0D99) THEN
          CALL SETWGT(Distance(1:3,SITE,J),WTSCHM,VARNUM,NPREQ,NDERIV,
     +                THETA,PWTIDX,CODE,ICHECK,IFAIL,CURWT,MXP,NVARS)
          IF (IFAIL.NE.0) RETURN
          SUMWT(0:NDERIV) = SUMWT(0:NDERIV) + CURWT(0:NDERIV)
          WGTARR(J,0:NDERIV) = CURWT(0:NDERIV)
         ENDIF
 110    CONTINUE
*
*       It's possible to get here with the sum of the weights as zero;
*       this must mean that there were no data values for the previous
*       day *anywhere*, and hence that the required covariate is 
*       missing
* 
        IF (DABS(SUMWT(0)).LT.1.0D-12) THEN
         DYPRED(0) = -1.0D100
         RETURN
        ENDIF
*
*     Now compute averages. Compute derivatives of normalised weights along 
*     the way. NB need to check that we're not dividing by zero. This
*     will occur if, for example, there's only 1 site. In this case it 
*     gets a weight of 1, and the correct derivative is clearly zero. 
*
        DO 120 J = 1,NSITES
         IF (DatArray(J,VARNUM,CODE).GE.-1.0D99) THEN
          CURWT(0) = WGTARR(J,0)
          DO 112 I=1,NDERIV
           IF (DABS(SUMWT(0)).GT.1.0D-12) THEN
            CURWT(I) = ( (WGTARR(J,I)*SUMWT(0)) - 
     +                   (WGTARR(J,0)*SUMWT(I)) ) / (SUMWT(0)**2)
           ELSE
            CURWT(I) = 0.0D0
           ENDIF
 112      CONTINUE
*     
*     Here are means of transformed values, and their derivatives
*
          IF (WXORDR.EQ.0) THEN
           IF (XFRM.EQ.5) THEN
            DO 114 I=1,CODE
             PREVEC(I) = DatArray(J,VARNUM,I)
 114        CONTINUE
           ENDIF
           CALL PRVXFM(DatArray(J,VARNUM,CODE),XFRM,TRACE,CODE,PREVEC,
     +                                               IDERIV,TMP,IFAIL)
           IF (IFAIL.NE.0) RETURN
           DO 116 I=0,NDERIV
            DYPRED(I) = DYPRED(I) + (CURWT(I) * TMP(0))
 116       CONTINUE
*
*     And here are transformed values of means. 'Persistence' has
*     to be done separately - basically, as soon as we find a zero
*     with non-zero weight, we're out.
*
          ELSEIF (WXORDR.EQ.1) THEN 
           IF (XFRM.NE.5) THEN
            DO 117 I=0,NDERIV
             DYPRED(I) = DYPRED(I) + 
     +                   (CURWT(I) * DatArray(J,VARNUM,CODE))
 117        CONTINUE
           ELSE
            DO 118 I=1,CODE
             IF ((DatArray(J,VARNUM,I).LE.0.0D0).AND.
     +           (CURWT(0).GT.0.0D0)) THEN
              DYPRED(0) = 0.0D0
              RETURN
             ENDIF
 118        CONTINUE
           ENDIF
          ELSE
           WRITE(MESSAGE,1)
           CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
           CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
           IFAIL = 999
           RETURN
          ENDIF
         ENDIF
 120    CONTINUE
*
*     Right, now finish off. Average of transformed value just has
*     to be turned into an average ...
*
        IF (WXORDR.EQ.0) THEN
         IF (DABS(SUMWT(0)).GT.1.0D-12) THEN
          DYPRED(0) = DYPRED(0)/SUMWT(0)
         ELSE
          DYPRED(0) = 0.0D0
         ENDIF
*
*     Transformed value of average has to be calculated, via
*     a call to PRVXFM unless it's a persistence thing (in
*     which case the only way we could have got here is if it's
*      a 1). Also, if we're estimating parameters (which will
*     happen if NERIV > 0), we'll need derivatives.
*
        ELSE
         IF (DABS(SUMWT(0)).GT.1.0D-12) THEN
          DYPRED(0) = DYPRED(0)/SUMWT(0)
         ELSE
          DYPRED(0) = 0.0D0
         ENDIF
         IF (NDERIV.GT.0) IDERIV = 1
         IF (XFRM.NE.5) THEN
          CALL PRVXFM(DYPRED(0),XFRM,TRACE,CODE,PREVEC,
     +                                          IDERIV,TMP,IFAIL)
          IF (IFAIL.NE.0) RETURN
          DYPRED(0) = TMP(0)
          DO 130 I=1,NDERIV
           DYPRED(I) = TMP(1) * DYPRED(I)
 130      CONTINUE
         ELSE
          DYPRED(0) = 1.0D0
         ENDIF
        ENDIF
       ENDIF

*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*     Other effects
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*
*     Daily representation of seasonality (cos/sin of 2 * pi * w * t/366,
*     where w is 1,2,3 or 4)
*
      ELSEIF ( (CODE.GE.21).AND.(CODE.LE.28) ) THEN
       CURDAY = CUMDAY(MM) + DD       
       IF (CODE.LE.22) THEN
        TRIGX = (PI*DBLE(CURDAY)/183.0D0)
       ELSE IF (CODE.LE.24) THEN
        TRIGX = (PI*DBLE(CURDAY)/91.5D0)
       ELSE IF (CODE.LE.26) THEN
        TRIGX = (PI*DBLE(CURDAY)/61.0D0)
       ELSE
        TRIGX = (PI*DBLE(CURDAY)/45.75D0)
       ENDIF
*
*	Next line picks out odd values of CODE (integer arithmetic)
*
       IF (2*((CODE-1)/2).EQ.CODE-1) THEN
         DYPRED(0) = DCOS(TRIGX)
       ELSE
         DYPRED(0) = DSIN(TRIGX)
       ENDIF
*
*     Smooth seasonal adjustment based on a bisquare function -
*     takes value zero on the 0th day of the required month, and on 
*     the (N+1)th day, where N is the number of days in the month
*
      ELSEIF ( (CODE.GE.31).AND.(CODE.LE.42) ) THEN
       MONADJ = CODE-30
       IF (MM.NE.MONADJ) THEN
        DYPRED(0) = 0.0D0
       ELSE
        TMP(0) = DBLE(NDAYS(MM)+1) / 2.0D0
        TMP(0) = (DBLE(DD) - TMP(0)) / TMP(0)
        DYPRED(0) = (1.0D0-(TMP(0)**2))**2
       ENDIF
      ELSE
       WRITE(MESSAGE(1),3)
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       IFAIL = 999
       RETURN
      ENDIF 

 1    FORMAT('****ERROR**** Invalid value of WXORDR in DYSET. This ',
     +'is a ',/,'programming error - contact REC.')
 3    FORMAT('****ERROR**** Undefined effect found in DYSET') 
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE PRVXFM(X,XFRM,TRACE,LAG,PREVEC,IDERIV,FX,IFAIL)
*
*     Calculates transformations of previous days' values. Arguments:
*
*     X (double precision, input): value to be transformed
*     XFRM (integer, input):       transform to be applied
*     TRACE (double, input):       trace threshold
*     LAG (integer,input):         Lag for persistence indicators
*     PREVEC (double, input):      Vector of previous values (only used to 
*				   calculate persistence indicators)
*     IDERIV (integer, input):     Calculate derivatives?
*     FX:                          Transformed value (output), together
*                                  with first derivative if there is one
*     IFAIL                        Error flag
******************************************************************************
      DOUBLE PRECISION X,TRACE,PREVEC(10),FX(0:1)
      INTEGER XFRM,LAG,IDERIV,IFAIL
******************************************************************************
*     Additional INTEGERs
*     ^^^^^^^^^^^^^^^^^^^
*     I          Counter
******************************************************************************
      INTEGER I
******************************************************************************
*     Additional CHARACTERs
*     ^^^^^^^^^^^^^^^^^^^^^
*     MESSAGE   For error messages
******************************************************************************
      CHARACTER*255 MESSAGE(2)

      IFAIL = 0
      
      IF (XFRM.EQ.0) THEN
       FX(0) = X
       IF (IDERIV.EQ.1) FX(1) = 1.0D0
      ELSEIF (XFRM.EQ.1) THEN
       IF (X.GT.0.0D0) THEN
        FX(0) = DLOG(X)
        IF (IDERIV.EQ.1) FX(1) = 1.0D0/X
       ELSE
        WRITE(MESSAGE,1)
        CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
        FX(0) = -1.0D100
        IF (IDERIV.EQ.1) FX(1) = -1.0D100
        IFAIL = 42
       ENDIF
      ELSEIF (XFRM.EQ.2) THEN
       IF (X.GT.-1.0D0) THEN
        FX(0) = DLOG(1.0D0+X)
        IF (IDERIV.EQ.1) FX(1) = 1.0D0/(1.0D0+X)
       ELSE
        WRITE(MESSAGE,1)
        CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
        FX(0) = -1.0D100
        IF (IDERIV.EQ.1) FX(1) = -1.0D100
        IFAIL = 42
       ENDIF
      ELSEIF (XFRM.EQ.3) THEN
       IF (X.GT.0.0D0) THEN
        FX(0) = 1.0D0
       ELSE
        FX(0) = 0.0D0
       ENDIF
       IF (IDERIV.EQ.1) FX(1) = 0.0D0
      ELSEIF (XFRM.EQ.4) THEN
       IF ((X.GT.0.0D0).AND.
     +     (X.LT.TRACE)) THEN
        FX(0) = 1.0D0
       ELSE
        FX(0) = 0.0D0
       ENDIF
       IF (IDERIV.EQ.1) FX(1) = 0.0D0
      ELSEIF (XFRM.EQ.5) THEN
       FX(0) = 1.0D0
       IF (IDERIV.EQ.1) FX(1) = 0.0D0
       DO 50 I=1,LAG
        IF (PREVEC(I).LE.0.0D0) THEN
         FX(0) = 0.0D0
         RETURN
        ENDIF
 50    CONTINUE
      ELSE
       WRITE(MESSAGE(1),2)
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       IFAIL = 300
      ENDIF
 1    FORMAT('****WARNING**** Attempt to log a non-positive ',
     +'number while',/,'calculating predictors from previous days')
 2    FORMAT('****ERROR**** Illegal transformation found in PRVXFM')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE SETWGT(DistVec,WTSCHM,VARNUM,NPREQ,NDERIV,THETA,
     +                  PWTIDX,LAG,ICHECK,IFAIL,WEIGHT,MXP,NVARS)
******************************************************************************
*     Calculates weights associated with neighbours of a site, when
*     averaging previous days' values. Arguments:
*     DistVec    - Vector of distances from site at which we're 
*                  predicting to neighbouring site. Elements are
*                  (x-separation, y-separation, Euclidean distance). 
*                  INPUT
*     WTSCHM     - Weighting scheme to be used. INPUT
*     VARNUM     - Number of the variable we're looking at.
*     NPREQ      - Number of parameters involved in various weighting 
*                  schemes. INPUT.
*     NDERIV     - Number of derivatives we'll actually calculate
*     THETA      - Parameters used in nonlinear transformations. INPUT
*     PWTIDX     - Index pointing to the row of THETA that holds
*                  the parameters used in WTSCHM. INPUT
*     LAG        - Number of days back we're going. Input
*     ICHECK     - Indicator for whether the model parametrisation
*                  has yet been checked. INPUT
*     IFAIL      - Error flag
*     WEIGHT     - Allocated weight, together with necessary derivatives
*                  with respect to nonlinear parameters.
*     MXP        - Maximum number of parameters } Used for dimensioning.
*     NVARS      - Maximum number of variables  } INPUT
******************************************************************************
      INTEGER WTSCHM,NPREQ(10),NDERIV,MXP,NVARS
      INTEGER PWTIDX(MXP,0:3,NVARS),VARNUM,ICHECK,IFAIL,LAG
      DOUBLE PRECISION DistVec(3),THETA(MXP,3),WEIGHT(0:3)
******************************************************************************
*     Additional INTEGERs
*     ^^^^^^^^^^^^^^^^^^^
*     ROW        - Required row of THETA
*     I,J        - counters
******************************************************************************
      INTEGER ROW
******************************************************************************
*     Additional DOUBLEs
*     ^^^^^^^^^^^^^^^^^^
*     DISTNC     - Distance between sites
*     XSEP,YSEP  - Site separation in X & Y directions
******************************************************************************
      DOUBLE PRECISION DISTNC,XSEP,YSEP
******************************************************************************
*     Additional CHARACTERs
*     ^^^^^^^^^^^^^^^^^^^^^
*     MESSAGE   - Error messages
******************************************************************************
      CHARACTER MESSAGE*200
      
      IFAIL = 0
*
*     First up: check that if there ought to be some parameters 
*     associated with this weighting scheme, they've been defined
*
      IF (WTSCHM.NE.1) THEN
       ROW = PWTIDX(WTSCHM,0,VARNUM)
       IF (ROW.EQ.0) GOTO 920
      ENDIF
*
*     Uniform weights:
*
      IF (WTSCHM.EQ.1) THEN
       WEIGHT(0) = 1.0D0
*
*     Exponentially weighted by distance - check that the required
*     parameter has been defined (and that no other parameter
*     has) and if so, compute the weight and its derivative wrt
*     the exponential decay rate.
*  
      ELSEIF (WTSCHM.EQ.2) THEN
       IF (ICHECK.EQ.0) CALL THTCHK(THETA,ROW,NPREQ(WTSCHM),MXP,IFAIL)
       IF (IFAIL.NE.0) RETURN
       WEIGHT(0) = DEXP(-THETA(ROW,1)*DistVec(3))
       IF (NDERIV.GT.0) THEN
        WEIGHT(1) = -(DistVec(3)*WEIGHT(0))
       ENDIF
      ELSEIF (WTSCHM.EQ.3) THEN
*
*     Shift and exponential decay. 3 parameters -> 3 derivatives.
*     Derivatives require special treatment at the origin :-(
*
       IF (ICHECK.EQ.0) CALL THTCHK(THETA,ROW,NPREQ(WTSCHM),MXP,IFAIL)
       IF (IFAIL.NE.0) RETURN
       XSEP = DistVec(1) - (DBLE(LAG)*THETA(ROW,1))
       YSEP = DistVec(2) - (DBLE(LAG)*THETA(ROW,2))
       DISTNC = DSQRT( (XSEP**2) + (YSEP**2) )
       WEIGHT(0) = DEXP(-THETA(ROW,3)*DISTNC)
       IF (NDERIV.GT.0) THEN
        IF (DISTNC.GT.1.0D-12) THEN
         WEIGHT(1) = THETA(ROW,3) * DBLE(LAG) * XSEP *
     +                                 WEIGHT(0) / DISTNC
         WEIGHT(2) = THETA(ROW,3) * DBLE(LAG) * YSEP *  
     +                                 WEIGHT(0) / DISTNC
         WEIGHT(3) = -(DISTNC*WEIGHT(0))
        ELSE
         WEIGHT(1) = -( THETA(ROW,3) * DBLE(LAG) * WEIGHT(0) )
         WEIGHT(2) = WEIGHT(1)
         WEIGHT(3) = 0.0D0
        ENDIF
       ENDIF
      ELSE
       MESSAGE = '***** ERROR **** Illegal weighting scheme '//
     +           'requested in SETWGT.'  
       CALL INTPR(TRIM(MESSAGE),-1,0,0)
       IFAIL = 999
      ENDIF

      RETURN
*
*     Error trapping beyond here
*

 920  IFAIL = 60

      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE THTCHK(THETA,ROW,NPAR,MXP,IFAIL)
******************************************************************************
*     Checks that the correct number of nonlinear parameters have 
*     been defined. Arguments (all input):
*     THETA     - array of nonlinear parameters
*     ROW       - Row of THETA that we're looking at
*     NPAR      - Number of parameters that should have been defined
*     MXP       - Dimensioning for THETA
*     IFAIL     - Error flag
******************************************************************************
      Integer, intent(in) :: ROW,NPAR,MXP
      Integer, intent(out) :: IFAIL
      Double precision, intent(in) :: THETA(MXP,3)
******************************************************************************
*       Additional INTEGERs
*       -------------------
*     I         - counter
******************************************************************************
      INTEGER I
******************************************************************************
*       Additional DOUBLEs
*       ------------------
*     MESSAGE   - Error messages
******************************************************************************
      CHARACTER MESSAGE(2)*255
      
      IFAIL = 0
*
*     First check: all necessary parameters *have* been defined
*
      DO 10 I=1,NPAR
       IF (THETA(ROW,I).GT.1.0D8) GOTO 910
 10   CONTINUE
*
*     Second: check that no extraneous parameters have been defined
*
      DO 20 I=NPAR+1,3
       IF (THETA(ROW,I).LE.1.0D8) GOTO 910
 20   CONTINUE

      RETURN
*
*     Error trapping
*

 910  WRITE(MESSAGE,1) ROW
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      IFAIL = 300

 1    FORMAT('****ERROR**** Wrong number of parameters supplied ',
     +'for transformation ',/,'of covariate ',I2,'. Check model ',
     +'definition file.')
      END

