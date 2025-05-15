******************************************************************************
*       This file contains routines to facilitate the transfer of
*       text strings from Fortran to R. This is necessary because 
*       the R-Fortran interface doesn't allow the transfer of
*       character arrays. To pass a character array back to R
*       therefore, what is done is to write each element to 
*       standard output and, within R, to wrap the .Fortran() call 
*       in capture.output(). Tail wagging the dog or what! 
******************************************************************************
      SUBROUTINE WriteLabels(MODEL,NSITES,MXP,NATTR,SITINF,TmpFlNo,
     +                  NP,COVCODE,NVARS,THETA,SITXFM,FOUIDX,LEGIDX,
     +                  PWTIDX,WTSCHM,SPMOD,GLBCOD,UnitNos)
******************************************************************************
*       To set up labels for model output (without passing character
*       information directly between R and Fortran - instead, the
*       site codes are in a temporary file indexed by TmpFlNo)
*
*      INPUT ARGUMENTS
*      ^^^^^^^^^^^^^^^
*       NSITES  - Number of sites in study region
*       NATTR   - Number of attributes defined for each site
*       SITINF  - Double precision array containing information 
*                 regarding each site (Eastings, Northings etc.). 
*                 Column 0 holds the actual value, 1 and 2 hold 
*                 derivatives wrt any parameters of nonlinear
*                 transformations
*       TmpFlNo - Number of temporary file from which site codes
*                 will be read
*       NP      - Integer array indicating numbers of parameters 
*                 in each model component. See fitting.f for details.
*       IEXT    - Integer array indicating whether `external' predictors 
*                 have been requested at yearly (1), monthly (2) and
*                 daily (3) time scales
*       MXP       Maximum number of parameters in each model (integer)
*       COVCODE - Integer array containing code numbers of covariates. NB
*                 this gets modified by a call to ATTRXFM
*       NVARS   - Number of variables (actually the final dimension of
*                 PWTIDX)
*	THETA	- Double precision array of parameters in nonlinear
*                 functions of basic covariates. See fitting.f for details
*	SITXFM	- Integer array selecting nonlinear transformations of
*                 site attributes. First column chooses the attribute,
*                 second the xfrmation.
*	FOUIDX	- Integer array indexing Fourier representation of 
*                 site effects. See file funcdefs.f for full details
*	LEGIDX	- Integer array indexing Legendre polynomial
*                 representation of site effects.
*       PWTIDX  - Integer array indexing  weighted averages of previous
*                 days' values. See GLMFIT header for details.
*       WTSCHM  - More indexing for weighted averages of previous days'
*                 values (Ith element is the weighting scheme used in Ith
*                 covariate).
*	SPMOD   - Integer, indicating choice of structure for 
*                 inter-site dependence
*       GLBCOD  - Integer array of codes defining `global' quantities 
*                 in GLBVAL
*       UnitNos - array containing numbers of units connected to 
*                 various input files. See GLMFIT header for details. 
***************************************************************************
      INTEGER, INTENT(IN) :: MODEL,NSITES,NATTR,NVARS
      DOUBLE PRECISION SITINF(NSITES,MXP,0:3)
      INTEGER, INTENT(IN) :: TmpFlNo, MXP, NP(10)
      INTEGER, INTENT(INOUT) :: COVCODE(MXP)
      INTEGER, INTENT(IN) :: SPMOD,GLBCOD(MXP)
      INTEGER, INTENT(IN) :: SITXFM(MXP,2),PWTIDX(MXP,0:3,NVARS)
      INTEGER, INTENT(IN) :: FOUIDX(MXP),LEGIDX(MXP)
      INTEGER, INTENT(IN) :: WTSCHM(MXP),UnitNos(100)
      DOUBLE PRECISION, INTENT(IN) :: THETA(MXP,3)
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
******************************************************************************
*               INTEGER variables
*               ^^^^^^^^^^^^^^^^^
*       ICHECK   - Indicator for whether model parametrisation has been 
*                  checked
*       IFAIL    - Error flag
***************************************************************************
      INTEGER I,J,ICHECK,IFAIL
******************************************************************************
***************************************************************************
***************************************************************************
*               CHARACTER variables
*               ^^^^^^^^^^^^^^^^^^^
*       XLABELS - Variable labels for columns of design matrix
*	NLTXT	- Descriptive text for nonlinear parameters
*	TRPTXT	- Describes nonlinear parameters in trend functions
*       SCODES  - array of short site codes (not used here except as
*                 an argument to MakeScodes)
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
*       TempLabel - Temporary label
***************************************************************************
      CHARACTER XLABELS(0:MXP)*70,NLTXT(MXP,3)*70
      CHARACTER SCODES(NSITES)*4
      CHARACTER*70 ATTRTXT(MXP)
      CHARACTER*70 MOTXT(MXP2),TRNDTXT(MXP2),DYTXT(0:MXP2),GLBTXT(MXP2)
      CHARACTER*70 PWTTXT(MXP2,3),TRPTXT(MXP2,3),SXPTXT(MXP2,3)
      CHARACTER*70 SPATXT(0:MXP2),SPPTXT(MXP2,4)
      CHARACTER*70 CORTXT(MXP2)
      Character*200 TempLabel
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>    COMMON BLOCKS   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /DESCRIBE/ SXPTXT,MOTXT,TRNDTXT,TRPTXT,DYTXT,
     +                  SPATXT,SPPTXT,CORTXT,PWTTXT,GLBTXT
******************************************************************************
******************************************************************************
******************************************************************************
*
*       Initializing ...
*
      DATA ICHECK /0/
**************************************************************************
**********		EXECUTABLES START HERE			**********
**************************************************************************
*
*       More initialisation, here using Fortran 90 array
*       capabilities
*
      
      DO I=1,NATTR
       WRITE(ATTRTXT(I),'(''~~~site attribute#'',I4,''#~~~'')') I
      END DO
      XLABELS = 'Undefined'
      NLTXT = 'Undefined'
      IFAIL = 0
*
*       Read site codes from temporary file
*
      CALL GetScodes(TmpFlNo, NSites, Scodes, Ifail)
      IF (IFAIL.NE.0) RETURN
*
*       Generate descriptive text ...
*
      CALL ATTRXFM(SITINF,NSITES,ATTRTXT,NATTR,COVCODE,
     +                          SITXFM,FOUIDX,LEGIDX,THETA,NP,
     +                          ICHECK,IFAIL,MXP)
      IF (IFAIL.NE.0) GOTO 910

      CALL MakeLabels(UnitNos,NP,COVCODE,GLBCOD,SITXFM,PWTIDX,WTSCHM,
     +       ATTRTXT,SXPTXT,TRNDTXT,TRPTXT,MOTXT,DYTXT,PWTTXT,XLABELS,
     +       NLTXT,GLBTXT,MXP,MXP2,NVARS,IFAIL)
      IF (IFAIL.NE.0) GOTO 910
*
*        ... and output (this will be captured by the R call to this
*        routine)
*
      DO I=0,NP(4)
       WRITE(TempLabel,'(''Covariate || '',I4,'' || '',A)') 
     +                                          I+1,TRIM(XLABELS(I))
       CALL INTPR(TRIM(TempLabel),-1,0,0)
      END DO

      DO I=1,NP(4)
       DO J=1,3
        WRITE(TempLabel,
     +        '(''Nonlinear || '',I4,'' || '',I1,'' || '',A)') 
     +                                          I,J,TRIM(NLTXT(I,J))
        CALL INTPR(TRIM(TempLabel),-1,0,0)
       END DO
      END DO
      
      DO I=1,NP(8)
       WRITE(TempLabel,'(''Global    || '',I4,'' || '',A)') 
     +                                          I,TRIM(GLBTXT(I))
       CALL INTPR(TRIM(TempLabel),-1,0,0)
      END DO
      
      TempLabel = 'Corstruct || Structure used: '//SPATXT(SPMOD)
      IF ( .NOT.((SPMOD.GE.20).OR.(SPMOD.EQ.0)) ) THEN
       TempLabel = TRIM(TempLabel)//' '//CORTXT(MODEL)
      ENDIF
      CALL INTPR(TRIM(TempLabel),-1,0,0)
      IF (SPMOD.NE.0) THEN
       DO I=1,NP(10)
        WRITE(TempLabel,'(''Corparams || '',I4,'' || '',A)') 
     +                                          I,TRIM(SPPTXT(SPMOD,I))
        CALL INTPR(TRIM(TempLabel),-1,0,0)
       END DO
      ENDIF

      DO J=1,3
       If (UnitNos(2+J).GT.0) CLOSE(UnitNos(2+J))
      End Do

      RETURN
      
 910  WRITE(TempLabel,'(''Error code '',I3)') IFAIL
      CALL INTPR(TRIM(TempLabel),-1,0,0)
      
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE GetScodes(TmpFlNo, NSites, Scodes, Ifail)
*
*     To read site codes from a temporary file (required as a result of 
*     R no longer supporting the passing of character vectors to Fortran).
*     Arguments:
*
*     TmpFlNo     Number of temporary file from which to read information.
*		  Integer, input.
*     NSites      Number of sites for which to read codes. Integer, input.
*     Scodes      Vector of 4-character site codes (output)
*     Ifail       Error flag. Integer, output. 
*
      INTEGER, INTENT(IN) :: TmpFlNo, NSites
      INTEGER, INTENT(OUT) :: Ifail
      CHARACTER, INTENT(OUT) :: Scodes(NSites)*4
*****************************************************************************
*	Additional INTEGERs
*	^^^^^^^^^^^^^^^^^^^
*       I       - counter
*	UnitNo	- Number of I/O channel to connect to temporary file
*****************************************************************************
      INTEGER I, UnitNo
*****************************************************************************
*	Additional CHARACTERs
*	^^^^^^^^^^^^^^^^^^^^^
*	TmpFlNm	- Name of temporary file
*       Msg     - Error message
*****************************************************************************
      CHARACTER*12 TmpFlNm
      CHARACTER*100 Msg

      IFAIL = 0
*
*     Construct name of temporary file, open it and read first line
*
      WRITE(TmpFlNm, '(''RGLCtmp_'',I4)') TmpFlNo  
      Call GetHandle(UnitNo)
      OPEN(UnitNo,FILE=TmpFlNm)
      READ(UnitNo,*,ERR=99)
*
*     Now read site codes
*
      Do I=1,NSITES
       READ(UnitNo,'(A4)',ERR=99) Scodes(I)
      End Do
      CLOSE(UnitNo)      

      RETURN

 99   Ifail=1
      WRITE(Msg,'(''*ERROR* reading site codes from temporary file'')')
      CALL INTPR(TRIM(Msg),-1,0,0)

      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE MakeLabels(UnitNos,NP,COVCODE,GLBCOD,SITXFM,PWTIDX,
     +     WTSCHM,ATTRTXT,SXPTXT,TRNDTXT,TRPTXT,MOTXT,DYTXT,PWTTXT,
     +     MDLTXT,NLTXT,GLBTXT,MXP,MXP2,NVARS,IFAIL)
*
*       Sets up descriptive text for the model whose covariate
*       structure and codes are held in NP and COVCODE. Variables:
*
*       UnitNos - array containing numbers of units connected to 
*                 various input files. See GLMFIT header for details. 
*       NP      - array indicating how many elements are associated
*                 with each model component. See main program header.
*                 Input.
*       COVCODE - Model specification for each covariate. Input
*       GLBCOD  - Codes for `global' quantities. Input
*	SITXFM	- Chooses nonlinear transformations of site attribs. Input
*       PWTIDX  - Indexing for weighted averages of previous days
*       WTSCHM  - Chooses weighting scheme for averaging previous days'
*                 values. Input.
*       ATTRTXT - Descriptive text for site attributes. Input
*       TRNDTXT - Descriptive text for trend functions. Input
*	TRPTXT	- Describes parameters in trend fns. Input
*       MOTXT   - Descriptive text for monthly functions. Input
*       DYTXT   - Descriptive text for previous days' dependence. Input
*       MDLTXT  - Text corresponding to each element of chosen model.
*                 Output.
*	NLTXT	- Text describing nonlinear parameters in model. Output.
*       MXP     - Max. number of model parameters. Input
*       NVARS   - Number of variables. Input
*       IFAIL   - Error flag
*****************************************************************************
      INTEGER, intent(in) :: UnitNos(100), NP(10), MXP, MXP2, NVARS
      INTEGER, intent(in), dimension(MXP) :: COVCODE,WTSCHM,GLBCOD
      INTEGER, intent(in) :: SITXFM(MXP,2)
      INTEGER, intent(in) :: PWTIDX(MXP,0:3,NVARS)
      CHARACTER (len=*), intent(in) :: ATTRTXT(MXP)
      CHARACTER (len=*), intent(in), dimension(MXP2) :: TRNDTXT,MOTXT
      CHARACTER (len=*), intent(in) :: DYTXT(0:MXP2)
      CHARACTER (len=*), intent(in), dimension(MXP2,3) :: 
     +                                 PWTTXT,TRPTXT,SXPTXT
      CHARACTER (len=*), intent(out) :: GLBTXT(MXP2)
      CHARACTER (len=*), intent(out), dimension(0:MXP) :: MDLTXT
      CHARACTER (len=*), intent(out), dimension(MXP, 3) :: NLTXT
      INTEGER, intent(out) :: IFAIL
*****************************************************************************
*	Additional INTEGERs
*	^^^^^^^^^^^^^^^^^^^
*	NTRND	- Number of current trend function
*       LAG     - Lag currently being considered
*       VARNUM  - Number of lagged variable currently being considered
* 	I	- counter
*****************************************************************************
      INTEGER NTRND,LAG,VARNUM,I
*****************************************************************************
*	Additional CHARACTERs
*	^^^^^^^^^^^^^^^^^^^^^
*	TMPCH?	- Temporary storage
*****************************************************************************
      CHARACTER TMPCH1*11,TMPCH2*70

      IFAIL = 0
      MDLTXT(0) = 'Constant'
      
      DO 10 I=1,NP(1)
       MDLTXT(I) = ATTRTXT(COVCODE(I))
       IF (SITXFM(I,1).NE.0) THEN
        NLTXT(I,1) = SXPTXT(SITXFM(I,2),1)
        NLTXT(I,2) = SXPTXT(SITXFM(I,2),2)
        NLTXT(I,3) = SXPTXT(SITXFM(I,2),3)
       ENDIF
 10   CONTINUE
      
      DO 20 I=NP(1)+1,NP(2)
*
*	For trend functions, we allow the possibility of more than one
*	and identify them with the prefix `Trend 1' etc.
*     
       IF (MOD(IABS(COVCODE(I)),1000).LE.50) THEN
        NTRND = I - NP(1)
        WRITE(TMPCH1,'(''Trend '',I2,'' - '')') NTRND
        TMPCH2 = TRNDTXT(COVCODE(I))
        MDLTXT(I) = TMPCH1//TMPCH2(1:59)
        NLTXT(I,1) = TRPTXT(COVCODE(I),1)
        NLTXT(I,2) = TRPTXT(COVCODE(I),2)
        NLTXT(I,3) = TRPTXT(COVCODE(I),3)
*
*     Set text for `external' annual influences if necessary
*
       ELSE
        CALL EXTXST(UnitNos(3),1,MOD(IABS(COVCODE(I)),1000)-50,
     +                                 COVCODE(I)/1000,MDLTXT(I),IFAIL)
        IF (IFAIL.NE.0) RETURN
       ENDIF
 20   CONTINUE

*
*     Monthly ...
*
      DO 30 I=NP(2)+1,NP(3)
       IF (MOD(IABS(COVCODE(I)),1000).LE.50) THEN
        MDLTXT(I) = MOTXT(COVCODE(I))
       ELSE
        CALL EXTXST(UnitNos(4),2,MOD(IABS(COVCODE(I)),1000)-50,
     +                                 COVCODE(I)/1000,MDLTXT(I),IFAIL)
        IF (IFAIL.NE.0) RETURN
       ENDIF
 30   CONTINUE
*
*     Daily effects. For those involving previous days' values,
*     we may have to change the existing text to process nonlinear
*     transformations.
*
      DO 40 I=NP(3)+1,NP(4)
       IF (MOD(IABS(COVCODE(I)),1000).LE.50) THEN
        LAG = MOD(COVCODE(I),1000)
        VARNUM = COVCODE(I) / 1000000       
        MDLTXT(I) = DYTXT(LAG)
        IF (COVCODE(I).GT.1000) THEN
         CALL PRVTXT(COVCODE(I),MDLTXT(I),IFAIL)
         IF (IFAIL.NE.0) RETURN
         IF (WTSCHM(I).GT.1) THEN
          IF (PWTIDX(WTSCHM(I),0,VARNUM).EQ.I) THEN
           IF (WTSCHM(I).NE.3) THEN
            NLTXT(I,1) = PWTTXT(WTSCHM(I),1)
            NLTXT(I,2) = PWTTXT(WTSCHM(I),2)
           ELSE
            NLTXT(I,1) = 'Offset for '//ATTRTXT(1)
            NLTXT(I,2) = 'Offset for '//ATTRTXT(2)
           ENDIF
           NLTXT(I,3) = PWTTXT(WTSCHM(I),3)
          ENDIF
         ENDIF
        ENDIF
       ELSE
        CALL EXTXST(UnitNos(5),3,MOD(IABS(COVCODE(I)),1000)-50,
     +                                 COVCODE(I)/1000,MDLTXT(I),IFAIL)
        IF (IFAIL.NE.0) RETURN
       ENDIF
 40   CONTINUE

*
*     'Global' quantities. No need to trap `undefined' quantities
*     here, as this has already been done in MDLSET.
*
      DO 80 I=1,NP(8)
       IF (GLBCOD(I)/1000.EQ.1) THEN
        IF (MOD(GLBCOD(I),1000).EQ.1) THEN
         GLBTXT(I) = 'Trace threshold'
        ELSEIF (MOD(GLBCOD(I),1000).EQ.2) THEN
         GLBTXT(I) = "'Soft' threshold for +ve values"
        ELSEIF (MOD(GLBCOD(I),1000).EQ.3) THEN
         GLBTXT(I) = "'Hard' threshold for +ve values"
        ENDIF
       ENDIF
 80   CONTINUE

      END
******************************************************************************
******************************************************************************
******************************************************************************
c      SUBROUTINE MakeScodes(ScodeChars,NSITES,ScodeDigits,SCODES)
c******************************************************************************
c*
c*	ROUTINE DEFUNCT AS OF RGLIMCLIM V 1.4.0
c*
c*       To generate a vector of 4-character site codes from the
c*       characters held in ScodeChars, according to the indices
c*       held in ScodeDigits (this all because the R-Fortran 
c*       interface doesn't allow the passing of character arrays!)
c*       Arguments:
c*
c*       ScodeChars      - A single string containing the characters
c*                         from which the codes will be constructed
c*                         (CHARACTER, input)
c*       NSITES          - The number of sites for which codes
c*                         are being defined (INTEGER, input)
c*       ScodeDigits     - An NSITES*4 integer array, giving the 
c*                         indices of the characters in ScodeChars
c*                         that go to make up the codes for each 
c*                         site (INTEGER, input)
c*       SCODES          - The vector of 4-character site codes
c*                         (CHARACTER, output)
c*
c******************************************************************************
c      CHARACTER (LEN=*), INTENT(IN) :: ScodeChars
c      INTEGER, INTENT(IN) :: NSITES, ScodeDigits(NSITES,4)
c      CHARACTER*4, INTENT(OUT) :: SCODES(NSITES)
c******************************************************************************
c*       Additional INTEGERs
c*       -------------------
c*
c*       I       - counter
c*
c******************************************************************************
c      INTEGER I
c      
c      DO I=1,NSITES
c       SCODES(I) = ScodeChars(ScodeDigits(I,1):ScodeDigits(I,1)) //
c     +             ScodeChars(ScodeDigits(I,2):ScodeDigits(I,2)) //
c     +             ScodeChars(ScodeDigits(I,3):ScodeDigits(I,3)) //
c     +             ScodeChars(ScodeDigits(I,4):ScodeDigits(I,4))
c      END DO
c      
c      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE WriteHeader(HeadType)
******************************************************************************
*
*       To generate header text for definition files. Arguments:
*
*       HeadType        - Integer indicating the type of definition file
*                         being written. Values:
*
*                         1   Model definition file
*                         2   Yearly external predictor file
*                         3   Monthly   "         "       "
*                         4   Daily     "         "       "
*
******************************************************************************
      INTEGER, INTENT(IN) :: HeadType
******************************************************************************
*       Additional INTEGERS
*       --------------------
*       EOLpos          Position of EOL string \n within header text
******************************************************************************
      INTEGER EOLpos
******************************************************************************
*       Additional CHARACTERs
*       ---------------------
*       ModelHead       Header text for model definition file
*       CurLine         Text for current line
******************************************************************************
      CHARACTER ModelHead*3680, CurLine*80

      IF (HeadType.EQ.1) THEN
       ModelHead = 'MODEL SPECIFICATION FILE\n========================'
     +//'\n \nThis file is used to define a generalised linear '//
     +'model for daily climate\nsequences. The row following this '//
     +'header is currently reserved for future\nuse. The next line '//
     +'can be used to define an overall description of the\nmodel '//
     +'(which will appear in summary output files). Following this '//
     +'is a\nrow for every covariate in the model. Each of these '//
     +'rows looks something\nlike this:\n \n'//
     +'COMPONENT       BETA    CODE1   CODE2   CODE3   TEXT\n \n'//
     +'where:\n'//
     +'        COMPONENT is    0 if we''re defining the constant '//
     +'term in the model\n                        1 if we''re '//
     +'dealing with a main site effect\n                        2 '//
     +'if we''re dealing with a main ''year'' effect\n'//
     +'                        3 if we''re dealing with a main '//
     +'''month'' effect\n                        4 if we''re '//
     +'dealing with previous days'' rainfall\n'//
     +'                        5 if we''re defining 2-way '//
     +'interactions\n                        6 if we''re defining '//
     +'3-way interactions\n                        7 if we''re '//
     +'defining nonlinear transformations\n'//
     +'                        8 to define values of ''global'''//
     +'quantities\n                        9 to specify a '//
     +'dispersion parameter\n                       10 for spatial '//
     +'structure\n \n'//
     +'                   and occupies positions 1-5 of the record\n'//
     +'\n        BETA       is the coefficient for this covariate in '
     +//'the model,\n                   occupying positions 6-15 of '//
     +'the record\n \n        CODE1}     are used to define the '//
     +'covariates to the system.\n        CODE2}     The '//
     +'interpretation varies depending on the value\n'//
     +'        CODE3}     of COMPONENT: tables of codes can be '//
     +'found in the\n                   program documentation. '//
     +'CODE1 occupies positions 16-20;\n                   CODE2 '//
     +'occupies positions 21-25 and CODE3 occupies            \n'//
     +'                   positions 26-30. CODE2 and CODE3 are used '//
     +'only for\n                   defining interactions.\n \n'//
     +'        TEXT       contains descriptive text for this '//
     +'covariate, and\n                   appears after position '//
     +'31. It is not used by the program.\n \nIt is vital '//
     +'that this positioning is adhered to, for each record will\n'//
     +'be read using the FORTRAN format I5,F10.6,3I5,A40.\n \n'//
     +'This header is 46 lines long.\n------------------------'//
     +'-----------------------------------------------------\n'
      ELSE IF (HeadType.Eq.2) THEN
       ModelHead='EXTERNAL PREDICTORS - YEARLY DATA FILE\n'//
     +'======================================\n \n'//
     +'        This file is used to define ''external'' predictors, '//
     +'for which\nannual data are available, to a Generalized '//
     +'Linear Model for some\nclimate/weather variable. By '//
     +'''external'', we mean non-deterministic\ntime-varying '//
     +'quantities (i.e. not trends) other than the variable under\n'//
     +'consideration. Excluding this header, there are 2 sections '//
     +'to the file.\nThe first is used to define and label the '//
     +'predictors, and the second to\nprovide the data.\n \n '//
     +'   Definition of predictors:\n    -------------------------\n'//
     +'        The first line of the file, after this header,'//
     +' should\n    contain a single integer NPRED, which is the'//
     +' number of\n    predictors being defined in this file. '//
     +'The subsequent NPRED lines\n    contain text (up to 70 '//
     +'characters) describing each predictor:\n    this will be '//
     +'used in labelling model output.\n \n'//
     +'    Format of data section\n    ----------------------\n'//
     +'        The "data input" section of the file begins on the '//
     +'2nd line\n    after the predictor definition section. Each '//
     +'line contains data for\n    1 month, and looks something '//
     +'like:\n \n     YYYY  PRED(1)  ...  PRED(2)  ...  '//
     +'PRED(NPRED)\n \n    where YYYY is the year\n'//
     +'          PRED(.) are the values of the various predictors '//
     +'in this\n          month (in the order specified in the '//
     +'''predictor definition''\n          section), each of which '//
     +'occupies the first 10 free positions\n          of the line '//
     +'starting with position 10\n \n    These rows are read using '//
     +'the FORTRAN format I4,T10,N(F10.0).\n \n    Missing data '//
     +'should be entered as -9999.9.\n \nThis header is 39 lines '//
     +'long.\n---------------------------------------------------'//
     +'----------------------------\n'
      ELSE IF (HeadType.Eq.3) THEN
       ModelHead='EXTERNAL PREDICTORS - MONTHLY DATA FILE\n'//
     +'=======================================\n \n'//
     +'        This file is used to define ''external'' predictors, '//
     +'for which\nmonthly data are available, to a Generalized '//
     +'Linear Model for some\nclimate/weather variable. By '//
     +'''external'', we mean non-deterministic\ntime-varying '//
     +'quantities (i.e. not trends) other than the variable under\n'//
     +'consideration. Excluding this header, there are 2 sections '//
     +'to the file.\nThe first is used to define and label the '//
     +'predictors, and the second to\nprovide the data.\n \n '//
     +'   Definition of predictors:\n    -------------------------\n'//
     +'        The first line of the file, after this header,'//
     +' should\n    contain a single integer NPRED, which is the'//
     +' number of\n    predictors being defined in this file. '//
     +'The subsequent NPRED lines\n    contain text (up to 70 '//
     +'characters) describing each predictor:\n    this will be '//
     +'used in labelling model output.\n \n'//
     +'    Format of data section\n    ----------------------\n'//
     +'        The "data input" section of the file begins on the '//
     +'2nd line\n    after the predictor definition section. Each '//
     +'line contains data for\n    1 month, and looks something '//
     +'like:\n \n     YYYY MM  PRED(1)  ...  PRED(2)  ...  '//
     +'PRED(NPRED)\n \n    where YYYY is the year\n'//
     +'          MM is the month\n'//
     +'          PRED(.) are the values of the various predictors '//
     +'in this\n          month (in the order specified in the '//
     +'''predictor definition''\n          section), each of which '//
     +'occupies the first 10 free positions\n          of the line '//
     +'starting with position 10\n \n    These rows are read using '//
     +'the FORTRAN format I4,1X,I2,T10,N(F10.0).\n \n    Missing '//
     +'data should be entered as -9999.9.\n \nThis header is 40 '//
     +'lines long.\n----------------------------------------------'//
     +'---------------------------------\n'
      ELSE IF (HeadType.Eq.4) THEN
       ModelHead='EXTERNAL PREDICTORS - DAILY DATA FILE\n'//
     +'=====================================\n \n'//
     +'        This file is used to define ''external'' predictors, '//
     +'for which\ndaily data are available, to a Generalized '//
     +'Linear Model for some\nclimate/weather variable. By '//
     +'''external'', we mean non-deterministic\ntime-varying '//
     +'quantities (i.e. not trends) other than the variable under\n'//
     +'consideration. Excluding this header, there are 2 sections '//
     +'to the file.\nThe first is used to define and label the '//
     +'predictors, and the second to\nprovide the data.\n \n '//
     +'   Definition of predictors:\n    -------------------------\n'//
     +'        The first line of the file, after this header,'//
     +' should\n    contain a single integer NPRED, which is the'//
     +' number of\n    predictors being defined in this file. '//
     +'The subsequent NPRED lines\n    contain text (up to 70 '//
     +'characters) describing each predictor:\n    this will be '//
     +'used in labelling model output.\n \n'//
     +'    Format of data section\n    ----------------------\n'//
     +'        The "data input" section of the file begins on the '//
     +'2nd line\n    after the predictor definition section. Each '//
     +'line contains data for\n    1 month, and looks something '//
     +'like:\n \n     YYYYMMDD  PRED(1)  ...  PRED(2)  ...  '//
     +'PRED(NPRED)\n \n    where YYYY is the year\n'//
     +'          MM is the month\n'//
     +'          DD is the day\n'//
     +'          PRED(.) are the values of the various predictors '//
     +'in this\n          day (in the order specified in the '//
     +'''predictor definition''\n          section), each of which '//
     +'occupies the first 10 free positions\n          of the line '//
     +'starting with position 10\n \n    These rows are read using '//
     +'the FORTRAN format I4,2I2,1X,N(F10.0).\n \n    Missing '//
     +'data should be entered as -9999.9.\n \nThis header is 41 '//
     +'lines long.\n----------------------------------------------'//
     +'---------------------------------\n'
      ENDIF

      EOLpos = INDEX(ModelHead,'\n')
      DO WHILE (EOLpos.gt.0)
       CurLine = ModelHead(1:(EOLpos-1))
       CALL INTPR(TRIM(CurLine),MAX(1,LEN_TRIM(CurLine)),0,0)
       ModelHead = ModelHead(EOLpos+2:)
       EOLpos = INDEX(ModelHead,'\n') 
      END DO 

      END SUBROUTINE WriteHeader
******************************************************************************
******************************************************************************
******************************************************************************
      Subroutine FileConnect(Filename,Varlen,UnitNo,Ifail)
******************************************************************************
*
*       To connect an ASCII file to a unit for input / output. Arguments:
*
*       Filename        Name of file
*       VarLen          Number of characters in Filename
*       UnitNo          Number of unit allocated to this file
*       Ifail           Error flag
*       
******************************************************************************
      Integer, intent(in) :: Varlen
      Character (len=Varlen), intent(in) :: Filename
      Integer, intent(out) :: UnitNo, Ifail
      
      Ifail = 0
      
      Call GetHandle(UnitNo)
      OPEN(UnitNo,FILE=Filename,ERR=990)
      Return
      
 990  IFAIL = 1
      
      End Subroutine FileConnect
******************************************************************************
******************************************************************************
******************************************************************************
      Subroutine RGLCFileConnect(TmpFlNo,Varlen,UnitNo,Ifail)
******************************************************************************
*
*       To connect an ASCII file to a unit for input / output, where the
*       filename is held in temporary file RGLCtmp_# and # is the value of
*       TmpFlNo (this is necessary because from R version 3.6.1 onwards, 
*       the passing of characters between R and Fortran was deprecated).
*       Arguments:
*
*       Filename        Name of file
*       VarLen          Number of characters in Filename
*       UnitNo          Number of unit allocated to this file
*       Ifail           Error flag
*
*	Additional variables:
*
*	TmpFlNm		Name of temporary file from which to read required
*                       filename
*       FlNmFmt         Format string used to read the filename. It will evaluate
*                       to something like (A13) if VarLen is 13
*       
******************************************************************************
      Integer, intent(in) :: TmpFlNo, Varlen
      Character (len=Varlen) :: Filename

      Integer, intent(out) :: UnitNo, Ifail
      Character (len=12) :: TmpFlNm
      Character (len=10) :: FlNmFmt
      
      Ifail = 0
      WRITE(TmpFlNm, '(''RGLCtmp_'',I4)') TmpFlNo  
      WRITE(FlNmFmt, '(''(A'',I7,'')'')') VarLen
*
*       Find an available file handle, attach it to the temporary file and
*       read the second line of this file
*
      Call GetHandle(UnitNo)
      OPEN(UnitNo,FILE=TmpFlNm)
      READ(UnitNo,*)
      READ(UnitNo,FlNmFmt) Filename
      CLOSE(UnitNo)

      OPEN(UnitNo,FILE=Filename,ERR=990)
      Return
      
 990  IFAIL = 1
      
      End Subroutine RGLCFileConnect
******************************************************************************
******************************************************************************
******************************************************************************
      Subroutine CloseFiles(UnitNos,Nunits)
******************************************************************************
*
*       To close all open ASCII file connections. 
*       
******************************************************************************
      Integer, intent(in) :: Nunits, UnitNos(Nunits)
      Integer i
      Logical IsOpen

      Do i=1,Nunits
       If (UnitNos(i).GT.0) then
        Inquire(Unit=UnitNos(i),Opened=IsOpen)
        If (IsOpen) CLOSE (UnitNos(i))
       End If
      End Do      
      
      End Subroutine CloseFiles
******************************************************************************
******************************************************************************
******************************************************************************
      Subroutine FileReset(UnitNos)
******************************************************************************
*       To reset pointers etc. relating to file units that may be 
*       open and connected following an earlier call within the same
*       R session; failure to reset things can lead to serious confusion 
*       on the part of the software! Argument:
*
*       UnitNos - Vector of three integers, indicating the units connected
*                 to files for annual, monthly and daily external information
*                 respectively
*
*       Note: the external units should already be connected on entry. No 
*       check is performed for this. 
*
*       The routine also closes any scratch files that are connected for 
*       calculation of inter-site correlations (these may remain 
*       open following calls to EMPCOR, because it doesn't know when it's
*       supposed to close them). 
******************************************************************************
      Integer, intent(in) :: UnitNos(3)
******************************************************************************
*     Additional Integers
*     -------------------
*     CHKFIL  - Indicators for whether input files have been checked. 
*     NPRDEF  - Number of predictors defined in the file
*     LIN1ID  - Index numbers of the first data line in each file
*     BEGLIN  - Line number of first data line in each file
*     ENDLIN  - Last line number in each file
*     CURLIN  - Last line number read, in each file
*     CORFNO  - Unit number of scratch file for correlations
*     OLDRES,OLDDEP,OLDFNO - see header to EMPCOR
******************************************************************************
      Integer CHKFIL(3),NPRDEF(3),LIN1ID(3),BEGLIN(3),ENDLIN(3)
      Integer CURLIN(3),I
      Integer CORFNO,OLDRES,OLDDEP,OLDFNO

      COMMON /ExtFileState/ CHKFIL,LIN1ID,CURLIN,BEGLIN,ENDLIN,NPRDEF
      COMMON /CorFileInfo/ CORFNO,OLDRES,OLDDEP,OLDFNO
      SAVE /ExtFileState/,/CorFileInfo/
******************************************************************************
*
*       Files containing external covariate information
*
******************************************************************************      
*
*       Ensure all open units are positioned at the start
*
      Do I=1,3
       If (UnitNos(I).GT.0) Rewind(UnitNos(I))
      End Do
*
*       (Re)initialise all variables to zero (Fortran 90 syntax)
*
      CHKFIL = 0
      BEGLIN = 0
      ENDLIN = 0
      CURLIN = 0
      NPRDEF = 0
******************************************************************************
*
*       Scratch files for correlations
*
******************************************************************************      
      If (CORFNO.NE.-1) CLOSE(CORFNO)
      CORFNO = -1
      OLDRES = -1
      OLDDEP = -1
      OLDFNO = -1
      
      End Subroutine FileReset
******************************************************************************
******************************************************************************
******************************************************************************
      Subroutine ScratchUpdate(FromUnit,ToUnit,UpdateType,Nobs,P,Ifail)
******************************************************************************
*       To update scratch files containing responses, case weights and
*       covariate information. Arguments:
*
*       FromUnit        - Unit number of first scratch file, from which
*                         information will be read
*       ToUnit          - Unit number of second scratch file to which
*                         information will be written
*       UpdateType      - 1 if the update involves extracting squared
*                         residuals from FromUnit and writing these
*                         as responses in ToUnit, 2 if it involves
*                         extracting fitted values from FromUnit and
*                         using these as inverse case weights in ToUnit
*       Nobs            - Number of observations (= number of records
*                         in each unit)
*       P               - Number of parameters in model corresponding to 
*                         ToUnit
*       Ifail           - Error flag
******************************************************************************
      Integer, intent(in) :: FromUnit,ToUnit,UpdateType,Nobs,P
      Integer, intent(out) :: Ifail 
******************************************************************************
*       Additional doubles
*       ~~~~~~~~~~~~~~~~~~
*       Y       - Response variables read from the two scratch files
*       CaseWt  - Weight attached to current case
*       X       - Array of covariate values
*       MU      - Fitted values in FromUnit
******************************************************************************
      Double precision Y(2), CaseWt, X(P), MU
******************************************************************************
*       Additional integers
*       ~~~~~~~~~~~~~~~~~~~
*       yy,mm,} - Year,month,day & site - first element corresponds to
*       dd,   }   FromUnit and the second to ToUnit
*       site  }
*       p       - Number of parameters estimated
*       i,j     - Counters
******************************************************************************
      Integer yy(2),mm(2),dd(2),site(2),i,j
******************************************************************************
*       Additional characters
*       ~~~~~~~~~~~~~~~~~~~~~
*       Message - For screen output
******************************************************************************
      Character Message*255
      
      Ifail = 0
      Do i=1,Nobs
       Read(FromUnit,REC=i) SITE(1),YY(1),MM(1),DD(1),Y(1),MU
       Read(ToUnit,REC=i) SITE(2),YY(2),MM(2),DD(2),
     +                                   Y(2),CaseWt,(X(J),J=1,P)
       If ( (Site(1).Ne.Site(2)).Or.(YY(1).Ne.YY(2)).Or.
     +      (MM(1).Ne.MM(2)).Or.(DD(1).Ne.DD(2)) ) then
        Write(Message,1)
        Call Intpr(Trim(Message),-1,0,0)
        Ifail = 999
        Return
       End If
       If (UpdateType.Eq.1) then
        Y(2) = ((Y(1)-MU)*(Y(1)-MU))
       Else If (UpdateType.Eq.2) then 
        CaseWt = 1.0D0 / MU
       Else
        Write(Message,2)
        Call Intpr(Trim(Message),-1,0,0)
        Ifail = 999
        Return
       End If
       Write(ToUnit,REC=i) SITE(2),YY(2),MM(2),DD(2),
     +                                    Y(2),CaseWt,(X(J),J=1,P)      
      End Do

 1    Format('In routine ScratchUpdate, cases in FromUnit and ',
     +       'ToUnit don''t match')
  2   Format('Illegal value of UpdateType in routine ScratchUpdate')
      End Subroutine ScratchUpdate
******************************************************************************
******************************************************************************
******************************************************************************
      Subroutine ScratchSync(UnitNos,Nobs,P,Ifail)
******************************************************************************
*       To synchronise scratch files for joint mean / disperson models
*	when the initial setup leads to a mismatch due to differing
*       availability of covariate information. Arguments:
*
*       UnitNos         - Unit numbers of scratch files for each model
*                         model. UnitNos(i,j) corresponds to the jth
*			  file for model i. A value of zero indicates
*	                  that we shouldn't read from this file. 
*       Nobs            - Numbers of cases available for each model 
*                         (= number of records in each unit). This is 
*                         defined on input, but not used. 
*       P               - Numbers of parameters in the models 
*       Ifail           - Error flag
******************************************************************************
      Integer, intent(in) :: UnitNos(2,2),P(2)
      Integer, intent(out) :: Nobs(2), Ifail 
******************************************************************************
*       Additional doubles
*       ~~~~~~~~~~~~~~~~~~
*       Y       - Response variables read from the scratch files
*       CaseWt  - Weights attached to current case in scratch files
*       X1,X2   - Arrays of covariate values in models
*       MU      - Fitted values in scratch files
*	CaseIDs	- Case identifiers in scratch files, in the form
*                 YYYYMMDDSSSSSS (year, month, day, site). These
*                 are integer-valued but the values exceed 
*                 2^31-1, so may as well store them as double.
*       CurID   - Case identifier for case to be written
******************************************************************************
      Double precision Y(2,2), CaseWt(2), X1(P(1)), X2(P(2)), MU(2)
      Double precision CaseIDs(2,2), CurID
******************************************************************************
*       Additional integers
*       ~~~~~~~~~~~~~~~~~~~
*	ReadRec - Record numbers being read from each unit
*	WriteRec- Record number to write to
*       ReadStatus - Indicates whether a read error has occurred
*       yy,mm,} - Year,month,day & site - first element corresponds to
*       dd,   }   FromUnit and the second to ToUnit
*       site  }
*       i,j     - Counters
******************************************************************************
      Integer ReadRec(2,2),WriteRec,ReadStatus(2,2)
      Integer yy(2,2),mm(2,2),dd(2,2),site(2,2),i,j
******************************************************************************
*       Additional logicals
*       ~~~~~~~~~~~~~~~~~~~
*       GetRec - indicates which files we want to get a record from on this 
*	         pass
******************************************************************************
      Logical GetRec(2,2)

      Data ((ReadRec(i,j),j=1,2),i=1,2) /4*1/
      
      Ifail = 0
      WriteRec = 0
      ReadStatus(1:2,1:2) = 0
      ReadRec(1:2,1:2) = 1
      GetRec(1:2, 1:2) = .True.

      Do while ( all(ReadStatus.Eq.0) )
*
*	Read next record from all files
*
       If ( (UnitNos(1,1).GT.0) .And. GetRec(1,1)) then
        Read(UnitNos(1,1),REC=ReadRec(1,1),IOSTAT=ReadStatus(1,1)) 
     +            SITE(1,1),YY(1,1),MM(1,1),DD(1,1),Y(1,1),
     +            CaseWt(1),(X1(J),J=1,P(1))
        ReadRec(1,1) = ReadRec(1,1) + 1
        CaseIDs(1,1) = Dble(SITE(1,1)) + 
     +      (1.0D6*Dble( (10000*YY(1,1)) + (100*MM(1,1)) + DD(1,1)))
       End If
       If ( (UnitNos(2,1).GT.0) .And. GetRec(2,1)) then
        Read(UnitNos(2,1),REC=ReadRec(2,1),IOSTAT=ReadStatus(2,1)) 
     +            SITE(2,1),YY(2,1),MM(2,1),DD(2,1),Y(2,1),
     +            CaseWt(2),(X2(J),J=1,P(2))
        ReadRec(2,1) = ReadRec(2,1) + 1
        CaseIDs(2,1) = Dble(SITE(2,1)) + 
     +      (1.0D6*Dble( (10000*YY(2,1)) + (100*MM(2,1)) + DD(2,1)))
       End If
       Do j=1,2
        If ( (UnitNos(j,2).GT.0) .And. GetRec(j,2) ) then
         Read(UnitNos(j,2),REC=ReadRec(j,2),IOSTAT=ReadStatus(j,2)) 
     +                  SITE(j,2),YY(j,2),MM(j,2),DD(j,2),Y(j,2),MU(j)
         ReadRec(j,2) = ReadRec(j,2) + 1
         CaseIDs(j,2) = Dble(SITE(j,2)) + 
     +       (1.0D6*Dble( (10000*YY(j,2)) + (100*MM(j,2)) + DD(j,2)))
        End If
       End Do
*
*	Figure out whether records match: "target" corresponds to the 
*	highest case ID (because records are in order). The "Where"
*	statement replaces all "undefined" case IDs with the
*	matching value, to avoid having to deal with them separately.
* 
       CurID = Maxval(CaseIDs,Mask=UnitNos.Gt.0)
       Where (UnitNos.Eq.0) CaseIds=CurId
       GetRec = (CaseIds.ne.CurId)
       If (.Not.(any(GetRec)) .And. (all(ReadStatus.Eq.0))) then
        WriteRec = WriteRec + 1
        If (UnitNos(1,1).GT.0) then
         Write(UnitNos(1,1),REC=WriteRec) 
     +            SITE(1,1),YY(1,1),MM(1,1),DD(1,1),Y(1,1),
     +            CaseWt(1),(X1(J),J=1,P(1))
        End If
        If (UnitNos(2,1).GT.0) then
         Write(UnitNos(2,1),REC=WriteRec) 
     +            SITE(2,1),YY(2,1),MM(2,1),DD(2,1),Y(2,1),
     +            CaseWt(2),(X2(J),J=1,P(2))
        End If
        Do j=1,2
         If (UnitNos(j,2).GT.0) then
          Write(UnitNos(j,2),REC=WriteRec) 
     +                   SITE(j,2),YY(j,2),MM(j,2),DD(j,2),Y(j,2),MU(j)
         End If
        End Do
*
*       If we got to here, we've just written a complete record and
*       definitely want to read from all open files on the next pass
*
        GetRec(1:2,1:2) = .True.
       End If
      End Do
*
*     Get here when there's a read error, probably due to end of file
*
      Nobs(1:2) = WriteRec

      End Subroutine ScratchSync
