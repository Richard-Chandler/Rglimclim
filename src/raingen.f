      SUBROUTINE WDALLOC(WETDRY,NSITES,NKNOWN,NKWET,FLAGS,SPMOD,
     +                   LOGITS,RHO,LGCORR,SIGMA11,SIGMA12,SIGMA22,
     +                   CVAR,CSD,Recalc,IFAIL)
***************************************************************************
*	Allocation of wet/dry sites. Spatial dependence is incorporated
*	via the mechanism requested in SPMOD. Arguments:
*	WETDRY	- INTEGER array (input/output). On input, it
*                 contains rainfall amounts at known sites; on output it 
*                 contains ones at "wet" sites and zeroes at "dry".
*       NSITES	- INTEGER (input): number of sites to generate
*	NKNOWN	- INTEGER (input): number of sites whose values are
*		  known (useful for shortcuts)
*       NKWET   - INTEGER (input): number of sites known to be wet
*	FLAGS	- INTEGER (input): status of each site
*	SPMOD	- INTEGER (input): choice of spatial dependence structure
*	LOGITS	- DOUBLE PRECISION array (input) containing marginal
*		  logits for rain
*	RHO	- DOUBLE PRECISION array (input) of parameters in 
*		  spatial dependence model
*       LGCORR  - DOUBLE PRECISION matrix (input) containing covariance
*                 matrix of latent Gaussian variables (used in some
*                 spatial dependence models)
*       CSD     - Conditional standard deviations needed during 
*                 Gibbs sampling stage
*       CVAR    - Covariance matrix of "unknown" sites given the known ones
*       Recalc  - Indicator for whether to recalculate covariance matrices etc. 
*       IFAIL   - INTEGER (output): error flag
***************************************************************************
      INTEGER NSITES,NKNOWN,NKWET,SPMOD,Recalc,IFAIL
      INTEGER WETDRY(NSITES),FLAGS(NSITES)
      DOUBLE PRECISION LOGITS(NSITES),RHO(4)
      DOUBLE PRECISION LGCORR(NSITES,NSITES)
      DOUBLE PRECISION, intent(inout), dimension(NSITES,NSITES) ::
     +                                  SIGMA11,SIGMA22,SIGMA12,CVAR
      DOUBLE PRECISION, intent(inout) :: CSD(NSITES)
***************************************************************************
*	Extra DOUBLEs
*	^^^^^^^^^^^^^
*	P	- Probability of wet at current site
*	ZBQLU01	- Random number generator
***************************************************************************
      DOUBLE PRECISION P,ZBQLU01
***************************************************************************
*	Extra INTEGERs
*	^^^^^^^^^^^^^^
*	I	- counter
***************************************************************************
      INTEGER I
***************************************************************************
*	Extra CHARACTERs
*	^^^^^^^^^^^^^^^^
*       MESSAGE - Error message
***************************************************************************
      Character MESSAGE(2)*255
***************************************************************************
      IFAIL = 0
*
*	First the easy case - independence model (easier to do it 
*	directly than as a special case of some other model)
*
      IF (SPMOD.EQ.0) THEN
       DO 100 I=1,NSITES
        IF (FLAGS(I).NE.0) THEN
         P = 1.0D0 / (1.0D0+DEXP(-LOGITS(I)))
         WETDRY(I) = 0
         IF (P.GT.ZBQLU01(0.0D0)) THEN
          WETDRY(I) = 1
         ENDIF
        ENDIF
 100   CONTINUE
*
*        Now the models involving latent Gaussian variables. Start by
*        calculating the Cholesky decomposition of the correlation 
*        matrix if it hasn't been done already.
*
      ELSEIF (SPMOD.LT.20) THEN
       CALL WDLG(WETDRY,NSITES,NKNOWN,FLAGS,LOGITS,LGCORR,
     +                  SIGMA11,SIGMA12,SIGMA22,CVAR,CSD,Recalc,IFAIL)
       IF (IFAIL.NE.0) RETURN
*
*        Here's the `Conditional independence given wet/dry weather state'
*        model. 
*
      ELSEIF (SPMOD.EQ.21) THEN
       CALL WD21(WETDRY,NSITES,NKNOWN,FLAGS,LOGITS,RHO(1))
*
*        And `Beta-binomial # of wet sites' model.
*
      ELSEIF (SPMOD.EQ.22) THEN
       CALL WD22(WETDRY,NSITES,NKNOWN,NKWET,FLAGS,LOGITS,RHO)
      ELSE
*
*        Trap any invalid spatial model selections (probably a 
*        programming error if they get this far)
*
       WRITE(MESSAGE,1)
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
       IFAIL = 999
       RETURN
      ENDIF
      
      RETURN

 1    FORMAT('****ERROR**** In routine WDALLOC, an invalid value ',
     +'of variable',
     +/'SPMOD has been found. This is a programming error - please ',
     +'contact REC.')
      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE WDLG(WETDRY,NSITES,NKNOWN,FLAGS,LOGITS,LGCORR,
     +                SIGMA11,SIGMA12,SIGMA22,CVAR,CSD,Recalc,IFAIL)
***************************************************************************
*	Implementation of occurrence model with spatial dependence 
*       structure defined via a correlated vector of latent Gaussian
*       variables. Arguments:
*	WETDRY	- INTEGER array (input/output) containing 1s and 0s
*       NSITES	- INTEGER (input): number of sites to generate
*	NKNOWN	- INTEGER (input): number of sites with known values
*	FLAGS	- INTEGER (input): status of each site
*	LOGITS	- DOUBLE PRECISION array (input) containing marginal
*		  logits for rain
*	LGCORR	- DOUBLE PRECISION (input) containing correlation matrix
*                 of latent Gaussian variables
*       SIGMA11 }
*       SIGMA22 }  Partitions of latent correlation matrix 
*       SIGMA12 }
*       CSD     - Conditional standard deviations needed during 
*                 Gibbs sampling stage
*       CVAR    - Covariance matrix of "unknown" sites given the known ones
*       Recalc  - Indicator for whether to recalculate covariance matrices
*       IFAIL   - Error flag
***************************************************************************
      INTEGER NSITES,NKNOWN,FLAGS(NSITES),WETDRY(NSITES),Recalc,IFAIL
      DOUBLE PRECISION, intent(in) :: LOGITS(NSITES),
     +                                         LGCORR(NSITES,NSITES)
***************************************************************************
*	Extra DOUBLEs
*	^^^^^^^^^^^^^
*	P       - Vector of probabilities of wet at each site
*       TAU     - Vector of thresholds for the latent Gaussian 
*                 variables; values above / below these correspond 
*                 respectively to wet or dry sites.
*       QNORM   - Quantile function for the normal distribution
*       ZBQLNOR - Normal random number generator
*       RTRNOR  - Truncated normal random number generator
*       ZVEC    - Vector of latent Gaussian variables
*       ZNEW    - New variate generated in Gibbs sampler
*       CMEAN   - Conditional mean of latent variables given those
*                 that have already been allocated correctly
*       TMPARR? - Temporary storage
***************************************************************************
      DOUBLE PRECISION P(NSITES),TAU(NSITES)
      DOUBLE PRECISION QNORM,ZBQLNOR,RTRNOR,ZVEC(NSITES)
      DOUBLE PRECISION SIGMA11(NSITES,NSITES),
     +                 SIGMA22(NSITES,NSITES),
     +                 SIGMA12(NSITES,NSITES)
      DOUBLE PRECISION CMEAN(NSITES),CSD(NSITES),CVAR(NSITES,NSITES)
      DOUBLE PRECISION TMPARR1(NSITES)
***************************************************************************
*	Extra INTEGERs
*	^^^^^^^^^^^^^^
*       CONIDX  - Vector containing indices of sites for which 
*                 observations are available
*       GENIDX  - Vector containing indices of sites for which 
*                 data must be generated
*       CURIDX  - Index of current site
*       KSOFAR  - Number of sites found so far with known rainfalls
*       USOFAR  - Number of sites found so far with unknown rainfalls
*       ITER    - Counter for Gibbs sampler iterations
*       NGIBBS  - Number of Gibbs iterations
*	I,J,K   - counters
***************************************************************************
      INTEGER CONIDX(NSITES),GENIDX(NSITES),CURIDX
      INTEGER KSOFAR, USOFAR
      INTEGER ITER, NGIBBS
      INTEGER I,J,K
***************************************************************************
*	Extra CHARACTERs
*	^^^^^^^^^^^^^^^^
*     MESSAGE   - For writing error messages
***************************************************************************
      CHARACTER*72 MESSAGE(5)
      
***************************************************************************
      
      NGIBBS = 10
*
*       If this is the first pass, compute the Cholesky decomposition 
*       of the correlation matrix since this will be used to generate
*       multivariate normal random numbers.
*
c      IF (DABS(RCHOL(1,1)).GT.1.0D11) THEN
c       DO 50 I=1,NSITES
c        DO 51 J=I+1,NSITES
c         IF (DABS(LGCORR(I,J)).GT.1.0D0) THEN
c          WRITE(MESSAGE,1)
c          CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
c          CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
c          CALL INTPR(TRIM(MESSAGE(3)),-1,0,0)
c          CALL INTPR(TRIM(MESSAGE(4)),-1,0,0)
c          IFAIL = 68
c          RETURN
c         ENDIF
c 51     CONTINUE
c 50    CONTINUE
c      ENDIF
*
*       First step: determine the thresholds for the latent Gaussian
*       variables, corresponding to the marginal logits; and order
*       the sites so that those with known values come first. If the 
*       probabilities are really close to 0 or 1, set the threshold to
*       something extremely large
*
      KSOFAR = 0
      USOFAR = 0
      DO 100 I=1,NSITES
       P(I) = 1.0D0/(1.0D0+DEXP(-LOGITS(I)))
       IF (P(I).LT.1.0D-20) THEN
        TAU(I) = -2.0D1
       ELSE IF ((1.0D0-P(I)).LT.1.0D-20) THEN
        TAU(I) = 2.0D1
       ELSE
        TAU(I) = -QNORM(P(I),0.0D0,1.0D0,IFAIL)
       ENDIF
       IF (FLAGS(I).EQ.0) THEN
        KSOFAR = KSOFAR+1
        CONIDX(KSOFAR) = I
       ELSE IF (FLAGS(I).EQ.1) THEN
        USOFAR = USOFAR + 1
        GENIDX(USOFAR) = I
       ENDIF
 100  CONTINUE 
*
*       Next: if there are any known sites then use Gibbs sampling
*       to sample from the distribution of the latent Gaussian 
*       variables (Z) for these sites. The initial configuration is 
*       obtained by sampling independently from the univariate
*       marginal distributions of each Z given the corresponding
*       observation; the conditionals are then just the usual
*       conditionals for a multivariate normal, truncated to the
*       range that is consistent with the observations. Some 
*       experimentation suggests that 10 Gibbs iterations should 
*       be enough (this doesn't have to be massively accurate anyway).
*       Start by generating an initial configuration and inverting
*       the covariance matrix (NB if there's only one site whose 
*       value is known, this initial configuration is exact so we
*       don't need to go and do Gibbs).
*
      IF (NKNOWN.GT.0) THEN
       DO 150 I=1,NKNOWN
        CURIDX = CONIDX(I)
        ZVEC(CURIDX) = RTRNOR(0.0D0,1.0D0,TAU(CURIDX),
     +                                    DBLE(2*WETDRY(CURIDX)-1))
        IF (DABS(ZVEC(CURIDX)).GT.1.0D100) THEN
         IFAIL = 999
         RETURN
        ENDIF

        If (Recalc.eq.1) then 
         SIGMA22(I,I) = LGCORR(CURIDX,CURIDX)
         DO 151 J=I+1,NKNOWN
          SIGMA22(I,J) = LGCORR(CURIDX,CONIDX(J))
          SIGMA22(J,I) = SIGMA22(I,J) 
 151     CONTINUE
        End if
 150    CONTINUE

        If (Recalc.Eq.1) CALL MATINV(SIGMA22,NKNOWN,1,NSITES,1,IFAIL)
       
       IF (NKNOWN.EQ.1) GOTO 299
*
*       For each element of ZVEC, the conditional variance given all
*       other elements is the inverse of the corresponding element of
*       SIGMA22 - this never changes, so store the conditional standard 
*       deviations to save repeated calls to DSQRT later on. 
*
       If (Recalc.Eq.1) then 
        DO 170 I=1,NKNOWN
         CSD(I) = 1.0D0 / DSQRT(SIGMA22(I,I))
 170    CONTINUE
       End if
*
*       Now the Gibbs iterations. For each variate, the conditional
*       mean is a linear combination of all the others; store these
*       conditional means in CMEAN for the moment (it's not being 
*       used for anything else). 
*
       DO 200 ITER=1,NGIBBS
        DO 210 I=1,NKNOWN
         CMEAN(I) = 0.0D0
         CURIDX = CONIDX(I)
         DO 211 J=1,NKNOWN
          IF (J.EQ.I) GOTO 211
          CMEAN(I) = CMEAN(I) - (SIGMA22(I,J)*ZVEC(CONIDX(J)))
 211     CONTINUE
         CMEAN(I) = CMEAN(I) / SIGMA22(I,I)
         ZVEC(CURIDX) = RTRNOR(CMEAN(I),CSD(I),TAU(CURIDX),
     +                                  DBLE(2*WETDRY(CURIDX)-1))
         IF (DABS(ZVEC(CURIDX)).GT.1.0D100) THEN
          IFAIL = 999
          RETURN
         ENDIF
 210    CONTINUE
 200   CONTINUE
      ENDIF
*
*       That's the known sites dealt with; now sample from the conditional
*       distribution of the remainder. Start by calculating the remaining
*       partions of the full covariance matrix i necessary (SIGMA22 already  
*       holds the inverse of the covariance matrix for the "known" sites). 
*
 299  IF (USOFAR.GT.0) THEN
       If (Recalc.Eq.1) then 
        DO 300 I=1,USOFAR
         DO 301 J=1,USOFAR
          SIGMA11(I,J) = LGCORR(GENIDX(I),GENIDX(J))  
 301     CONTINUE
         DO 302 J=1,NKNOWN
          SIGMA12(I,J) = LGCORR(GENIDX(I),CONIDX(J))
 302     CONTINUE
 300    CONTINUE
       End If
*
*       Here are the conditional mean and covariance matrix of the
*       unobserved sites (see comments in subroutine RAINGEN
*       for logic - it's essentially the same as that used there
*       for allocating rainfall amounts).
*
       DO 320 I=1,USOFAR
        CMEAN(I) = 0.0D0
        DO 325 J=1,NKNOWN
         TMPARR1(J) = 0.0D0
         DO 326 K=1,NKNOWN
          TMPARR1(J) = TMPARR1(J) + (SIGMA12(I,K)*SIGMA22(K,J))
 326     CONTINUE
         CMEAN(I) = CMEAN(I) + (ZVEC(CONIDX(J))*TMPARR1(J))
 325    CONTINUE
        If (Recalc.Eq.1) then
         DO 330 J=1,USOFAR
          CVAR(I,J) = SIGMA11(I,J)
          DO 331 K=1,NKNOWN
           CVAR(I,J) = CVAR(I,J) - (SIGMA12(J,K)*TMPARR1(K))
 331      CONTINUE
 330     CONTINUE
        End If
 320   CONTINUE

       If (Recalc.Eq.1) CALL CHOLDEC(CVAR,USOFAR,NSITES,IFAIL)
 
       IF (IFAIL.NE.0) THEN
        WRITE(MESSAGE,1) GENIDX(IFAIL)
        CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(3)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(4)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(5)),-1,0,0)
        IFAIL = 68
        RETURN
       ENDIF

       DO 400 I=1,USOFAR
        TMPARR1(I) = ZBQLNOR(0.0D0,1.0D0)
        ZVEC(GENIDX(I)) = CMEAN(I)
        DO 401 J=1,I
         ZVEC(GENIDX(I)) = ZVEC(GENIDX(I))+(CVAR(J,I)*TMPARR1(J))
 401    CONTINUE
 400   CONTINUE
*
*       ZBQLNOR generates *pairs* of deviates using the Box-Muller
*       algorithm, and saves the unused one for reuse on the next
*       call. This is done for efficiency, but it means that 
*       simulations potentially are not reproducible unless you
*       can guarantee that there's no spare deviate sitting there
*       at the end of the last run. The way to guarantee this is
*       to ensure that you always generate an even number of deviates.
*
       IF (2*(USOFAR/2).NE.USOFAR) 
     +                  TMPARR1(USOFAR) = ZBQLNOR(0.0D0,1.0D0) 
       
      ENDIF
*
*       And finally, convert to 0/1
*     
      DO 500 I=1,NSITES
       IF (FLAGS(I).EQ.1) THEN
        IF (ZVEC(I).GT.TAU(I)) THEN
         WETDRY(I) = 1
        ELSE
         WETDRY(I) = 0
        ENDIF
       ENDIF
 500  CONTINUE
      
 1    FORMAT('****ERROR**** latent correlation matrix isn''t ',
     +'positive definite.',/5X,'Problem found at site ',I3,'. ',
     +'If you are using an empirical',/5X,'correlation matrix, you ',
     +'may need to fit a valid correlation model',/5X,
     +'to get round this problem. Otherwise, check for duplicate ',
     +'sites,',/5X,'or sites that are very close together.')

      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE WD21(WETDRY,NSITES,NKNOWN,FLAGS,LOGITS,LOGA)
***************************************************************************
*	Implementation of occurrence model with spatial dependence 
*       structure 21 (conditional independence given a binary weather
*       state variable). Arguments:
*	WETDRY	- INTEGER array (output) containing 1s and 0s
*       NSITES	- INTEGER (input): number of sites to generate
*	NKNOWN	- INTEGER (input): number of sites with known values
*	FLAGS	- INTEGER (input): status of each site
*	LOGITS	- DOUBLE PRECISION array (input) containing marginal
*		  logits for rain
*	LOGA	- Coefficient of X in updated logit
***************************************************************************
      INTEGER NSITES,NKNOWN,FLAGS(NSITES),WETDRY(NSITES)
      DOUBLE PRECISION LOGITS(NSITES),LOGA
***************************************************************************
*	Extra DOUBLEs
*	^^^^^^^^^^^^^
*	ALPHA	- Overall probability that today is a `wet' day (X=1)
*	XPROB	- Probability that today is `wet' given observations
*	A	- DEXP(LOGA) (obviously!!)
*	B	- Adjustment needed to make the marginal 
*		  probabilities right
*	P	- Overall probability of wet at current site
*	NUMER}	- Numerator and denominator in Bayes' Theorem for
*	DENOM}	  updating probability of today being wet based on
*		  observations. 
*	ZBQLU01	- Random number generator
***************************************************************************
      DOUBLE PRECISION ALPHA,A,B,ZBQLU01
      DOUBLE PRECISION P,TMP,NUMER,DENOM,XPROB
***************************************************************************
*	Extra INTEGERs
*	^^^^^^^^^^^^^^
*	X	- 1 if today is a `wet' day, 0 otherwise
*	I	- counter
***************************************************************************
      INTEGER X,I
***************************************************************************
*     Start by calculating the probability that today is a 
*	`wet' day - ALPHA is mean of predicted probabilities, and this
*	can be adjusted to take into account observed wet/dry patterns
*	at different sites.
*
      A = DEXP(LOGA)
      ALPHA=0.0D0
      DO 200 I=1,NSITES
       TMP = DEXP(LOGITS(I))
       ALPHA=ALPHA+(TMP/(1.0D0+TMP))
 200  CONTINUE
      ALPHA = ALPHA/DBLE(NSITES)
      XPROB = ALPHA
*
*	That's ALPHA sorted - however,
*	if we've observed values at some of the sites then we can use
*	Bayes' Theorem to obtain an updated estimate of the probability
*	that today is wet. The numerator is equal to
*	P(observed pattern|today is wet) * P(today is wet)
*	and the denominator is equal to the numerator plus
*	P(observed pattern|today is dry) * P(today is dry)
*
      IF (NKNOWN.GT.0) THEN
       NUMER = 1.0D0
       DENOM = 1.0D0
       DO 205 I=1,NSITES
        IF (FLAGS(I).EQ.0) THEN
         P = 1.0D0/(1.0D0+DEXP(-LOGITS(I)))
         CALL BCALC(P,ALPHA,A,B)
*
*	Here are the probabilities of seeing the observed value of
*	WETDRY(I) if it's (a) a wet day (b) a dry day
*
         TMP = A*B*P/(1.0D0-(P*(1.0D0-(A*B))))
         IF (WETDRY(I).EQ.0) TMP = 1.0D0-TMP
         NUMER = NUMER*TMP
         TMP = B*P/(1.0D0-(P*(1.0D0-B)))
         IF (WETDRY(I).EQ.0) TMP = 1.0D0-TMP
         DENOM = DENOM*TMP        
        ENDIF
 205   CONTINUE
       DENOM = (DENOM*(1.0D0-ALPHA)) + (NUMER*ALPHA)
       XPROB = NUMER/DENOM
      ENDIF
*
*	Now - is it a wet day or a dry day?
*
      X = 0
      IF (ZBQLU01(0.0D0).LT.XPROB) X=1
*
*	And set things wet or dry conditional on X. 
*
      DO 210 I=1,NSITES
       IF (FLAGS(I).NE.0) THEN
        WETDRY(I) = 0
        P = DEXP(LOGITS(I))/(1.0D0+DEXP(LOGITS(I)))
        CALL BCALC(P,ALPHA,A,B)
*
*        Here's the updated probability
*
        TMP = (A**X)*B*P/(1.0D0-(P*(1.0D0-((A**X)*B))))
         IF ( TMP.GT.ZBQLU01(0.0D0) )THEN
         WETDRY(I) = 1
        ENDIF
       ENDIF
 210  CONTINUE

      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE WD22(WETDRY,NSITES,NKNOWN,NKWET,FLAGS,LOGITS,RHO)
***************************************************************************
*	Implementation of occurrence model with spatial dependence 
*       structure 22 (dependence specified via a beta-binomial
*       distribution for the total number of wet sites). Arguments:
*	WETDRY	- INTEGER array (output) containing 1s and 0s
*       NSITES	- INTEGER (input): number of sites to generate
*	NKNOWN	- INTEGER (input): number of sites with known values
*	NKWET	- INTEGER (input): number of sites known to be wet
*	FLAGS	- INTEGER (input): status of each site
*	LOGITS	- DOUBLE PRECISION array (input) containing marginal
*		  logits for rain
*       RHO     - DOUBLE PRECISION array (input): parameters for 
*                 dependence scheme  
***************************************************************************
      INTEGER NSITES,NKNOWN,NKWET,FLAGS(NSITES),WETDRY(NSITES)
      DOUBLE PRECISION LOGITS(NSITES),RHO(4)
***************************************************************************
*	Extra INTEGERs
*	^^^^^^^^^^^^^^
*     I         - Counter
*     Z         - Number of wet sites
*     IFAIL     - Error flag
*     N1S       - Number of known 1s in vector being generated
*     IWARN     - Flag to warn user if sum distribution is incompatible
*                 with probabilities.
*     PASS1     - Inidcator for first pass through this routine.
*     NSHRNK    - Number of times we've shrunk the probabilities
*                 so far in this routine.
***************************************************************************
      INTEGER I,Z,IFAIL,N1S,IWARN,PASS1,NSHRNK
***************************************************************************
*	Extra DOUBLEs
*	^^^^^^^^^^^^^
*     P         - array of probabilities for the different sites
*     PZ        - probabilities for the Beta-Binomial distribution
*     PHI       - shape parameter of beta-binomial distribution
*     THETA     - Mean of the Beta-Binomial distribution
*     ALPHA }   - Parameters in standard parametrisation of Beta-
*     BETA  }     Binomial
*     TMP       - temporary storage
*     PMIN }    - Minimum and maximum probabilities at individual 
*     PMAX }      sites.
*     SHRFAC    - Amount by which to shrink the probabilities towards
*                 their trimmed mean, if they aren't consistent with
*                 the distribution of the number of wet sites.
***************************************************************************
      DOUBLE PRECISION P(NSITES),PZ(0:NSITES),PHI,THETA,ALPHA,BETA
      DOUBLE PRECISION TMP,PMIN,PMAX,SHRFAC
***************************************************************************
*	Extra CHARACTERs
*	^^^^^^^^^^^^^^^^
*     MESSAGE   - For warning messages
***************************************************************************
      CHARACTER*80 MESSAGE(7)

      DATA SHRFAC,IWARN /-1.0D0,-1/
      SAVE SHRFAC,IWARN
***************************************************************************
*     First up: determine parameters, and user options re warnings etc. 
*
      PASS1 = 0
      NSHRNK = 0
      PHI = RHO(1)
      IF (SHRFAC.LT.0.0D0) THEN
       IF ((RHO(2).GT.1.0D0).OR.(RHO(2).LE.0.0D0)) THEN
        SHRFAC = 0.01D0
        WRITE(MESSAGE,1) SHRFAC
        CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(3)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(4)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(5)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(6)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(7)),-1,0,0)
       ELSE
        SHRFAC = RHO(2)
       ENDIF
      ENDIF

      IF (IWARN.LT.0) THEN
       IF (RHO(3).GT.0.0D0) THEN
        IWARN = 1
       ELSE
        IWARN = 0
       ENDIF
      ENDIF
*
*     Now calculate probabilities of rain at individual sites.
*     These determine the mean of the Beta-Binomial distribution. The
*     IFs in here are needed if subsequently have to adjust the probs
*     for consistency with the sum distribution.
*
      N1S = NKWET
      THETA = 0.0D0
      PMIN = 1.0D0
      PMAX = 0.0D0
      IFAIL = 0

      DO 100 I=1,NSITES
       TMP = DEXP(LOGITS(I))
       P(I)=TMP/(1.0D0+TMP)
       THETA = THETA + P(I)
       IF (P(I).LT.PMIN) PMIN = P(I)
       IF (P(I).GT.PMAX) PMAX = P(I)
 100  CONTINUE
      THETA = THETA/DBLE(NSITES)
*
*     Next, calculate the Beta-Binomial probabilities and generate the
*     number of 1s.
*

 10   ALPHA = THETA*PHI
      BETA = PHI*(1.0D0-THETA)
      CALL BBZGEN(NSITES,ALPHA,BETA,PMIN,PMAX,PZ,IFAIL)
*
*     If all OK, allocate the Z wets to sites.
*
      IF (IFAIL.EQ.0) THEN
       CALL CBVSIM(NSITES,P,PZ,NKNOWN,WETDRY,Z,1,IFAIL)
      ENDIF

      IF (IFAIL.EQ.0) THEN
       IF ((NSHRNK.GT.0).AND.(IWARN.GT.0)) THEN
        WRITE(MESSAGE,3) 1.0D0-((1.0D0-SHRFAC)**NSHRNK)
        CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       ENDIF
       RETURN
      ENDIF
*
*     If either of those failed, shrink the probabilities towards
*     the mean and have another go. 
*
      IF ((IWARN.GT.0).AND.(PASS1.EQ.0)) THEN
       WRITE(MESSAGE,2) THETA,SHRFAC
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
       CALL INTPR(TRIM(MESSAGE(3)),-1,0,0)
       CALL INTPR(TRIM(MESSAGE(4)),-1,0,0)
      ENDIF
      PASS1 = 1
      NSHRNK = NSHRNK + 1
      DO 200 I=1,NSITES
       P(I) = THETA - ((1.0D0-SHRFAC)*(THETA - P(I)))
       IF (FLAGS(I).NE.0) WETDRY(I) = -1
 200  CONTINUE
      PMIN = THETA - ((1.0D0-SHRFAC)*(THETA-PMIN))
      PMAX = THETA - ((1.0D0-SHRFAC)*(THETA-PMAX))
      GOTO 10
 
 1    FORMAT('****WARNING**** You have asked for spatial dependence ',
     +'in rainfall',/,'occurrence to be modelled via a Beta-Binomial',
     +' distribution for the',/,'number of wet sites. This ',
     +'distribution is not guaranteed to be',/,'consistent with the ',
     +'forecast probabilities of rain at each site.',/,
     +'You haven''t specified a shrinkage parameter between 0 and 1 '
     +'for',/,'adjustment of the probabilities in case of ',
     +'inconsistency,',/,'so I''ve set it to ',F5.3,'.')
 2    FORMAT('****ERROR**** Beta-binomial distribution for number ',
     +'of wet sites is inconsistent',/,'with individual ',
     +'probabilities of rain. These will be repeatedly shrunk towards',
     +/,F5.3,' by a factor of ',F6.4,' until a consistent set is ',
     +'found.',/,'To suppress these messages, add the appropriate ',
     +'line to your logistic.def file.')
 3    FORMAT('Final shrinkage factor was ',F6.4,'.')
      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE BBZGEN(S,ALPHA,BETA,PMIN,PMAX,PZ,IFAIL)
***************************************************************************
*     To calculate the probability mass function for a Beta-Binomial 
*     distribution with parameters (S,ALPHA,BETA). The pmf
*     is returned in array PZ. To ensure that these probabilities
*     are compatible with the marginal probabilities at individual
*     sites (which range from PMIN to PMAX - input), we require
*     PZ(0) <= 1 -PMAX and PZ(S) <= PMIN. 
***************************************************************************
      INTEGER S,IFAIL
      DOUBLE PRECISION ALPHA,BETA,PZ(0:S),PMIN,PMAX
***************************************************************************
*     Extra INTEGERs
*     ^^^^^^^^^^^^^^
*     I        - Counter
***************************************************************************
      INTEGER I
***************************************************************************
*     Extra DOUBLEs
*     ^^^^^^^^^^^^^
*     GAMMLN   - Logarithm of gamma function
*     CUMPRB   - Cumulative probability up to and including current
*                value
*     TMP      - temporary storage
*     ZMEAN    - Mean of Z distribution as calculated. For maximum 
*                accuracy, we adjust the probabilities at the end
*                to ensure this is exactly S*ALPHA/(ALPHA+BETA) to 
*                machine precision.
*     MU    - Theoretical mean.
***************************************************************************
      DOUBLE PRECISION GAMMLN,CUMPRB,TMP,ZMEAN,MU
***************************************************************************
*     First step - calculate all the probabilities. 
* 
      IFAIL = 0
      ZMEAN = 0.0D0

      TMP = GAMMLN(ALPHA+BETA,IFAIL) - GAMMLN(BETA,IFAIL)
      PZ(0) = DEXP(TMP + GAMMLN(DBLE(S)+BETA,IFAIL) - 
     +             GAMMLN(DBLE(S)+ALPHA+BETA,IFAIL))
      CUMPRB = PZ(0)
      TMP = TMP - GAMMLN(ALPHA,IFAIL)

      DO 100 I=1,S
*
*     Use a recurrence relationship for the probabilities when 
*     S <= 100, otherwise try and avoid rounding errors by working 
*     directly with gamma functions.
*
       IF (S.LE.100) THEN
        PZ(I) = DBLE(S-I+1)*(ALPHA+DBLE(I-1)) * PZ(I-1) /
     +                               (DBLE(I)*(BETA+DBLE(S-I)))
       ELSE
        PZ(I) = DEXP( TMP + GAMMLN(DBLE(S+1),IFAIL) - 
     +                GAMMLN(DBLE(I+1),IFAIL) -
     +                GAMMLN(DBLE(S-I+1),IFAIL) + 
     +                GAMMLN(ALPHA+DBLE(I),IFAIL) + 
     +                GAMMLN(DBLE(S-I)+BETA,IFAIL) - 
     +                GAMMLN(ALPHA+BETA+DBLE(S),IFAIL))
       ENDIF
       ZMEAN = ZMEAN + (DBLE(I)*PZ(I))
       CUMPRB = CUMPRB + PZ(I)
 100  CONTINUE

*
*     Here are the adjustments to ensure computed probs sum to exactly
*     1 and have the correct mean. Just adjust probs for 0 and S
*     respectively - adjustments are miniscule, but errors may 
*     propagate later.
*
      MU = DBLE(S)*ALPHA/(ALPHA+BETA)
      TMP = (ZMEAN-MU)/DBLE(S)
      CUMPRB = CUMPRB - PZ(S)
      PZ(S) = DMAX1(0.0D0,PZ(S)-TMP)
      CUMPRB = CUMPRB + PZ(S)
      PZ(0) = PZ(0) - (CUMPRB-1.0D0)
*
*     Finally: check a necessary (but not sufficient!) condition for
*     PZ distribution to be compatible with marginal probabilities.
*

      IF ((PZ(0).GT.1.0D0-PMAX).OR.(PZ(S).GT.PMIN)) IFAIL = 1

      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE CBVSIM(S,P,PZ,NK,Y,Z,QUIET,IFAIL)
**********************************************************************
*     To generate a vector of S Bernoulli random variables, Y, such 
*     that P(Y(i) = 1) = P(i), and with sum whose distribution is
*     given in PZ. The sum is retruned in Z. The routine will perform 
*     imputation if necessary - where Y contains a non-negative element 
*     on entry, the corresponding element is assumed known. Only 
*     negative values will be replaced, conditioning on the observed 
*     values of Y. NK (input) is the number of observed values. IFAIL  
*     is an error flag. QUIET takes the value 1 if we wish to suppress 
*     printing of some error messages, 0 otherwise.
**********************************************************************
      INTEGER S,Z,NK,QUIET,IFAIL,Y(S)
      DOUBLE PRECISION P(S),PZ(0:S)
**********************************************************************
*     Additional INTEGERs
*     ^^^^^^^^^^^^^^^^^^^
*     IS     Site counter
*     IZ     Z counter
*     NUK    Number of sites with unknown values
*     IDX    Indices of known sites, with known ones first.
*     N1S    Number of 1s so far.
*     IALL   Value to put into all remaining vacant slots once they're
*            known.
*     Y1     Value at first site.
*     SSTAR  Number of sites left to consider
*     ZSTAR  Number of wet sites among sites left to consider
*     J      COUNTER
**********************************************************************
      INTEGER IS,IZ,NUK,IDX(S),N1S,IALL,Y1,SSTAR,ZSTAR,J
**********************************************************************
*     Additional DOUBLEs
*     ^^^^^^^^^^^^^^^^^^
*     CP      Array whose (s,z)th entry is equal to P(Y(s)=1|Z=z).
*     PZ2     Probability distribution of Z conditioned on observations.
*     P2      Probabilities for unobserved Ys given observations.
*     U       Uniform random number
*     ZBQLU01 Uniform random number generator
*     CUMPRB  Cumulative Z probability
*     TOL     Small number, used for checking rounding errors. NB
*             errors can propagate a bit in here, so this needs to 
*             be small - no scope for approximate inputs!
**********************************************************************
      DOUBLE PRECISION CP(S,0:S),PZ2(0:S),P2(S)
      DOUBLE PRECISION U,ZBQLU01,CUMPRB,TOL
*
*     First step: calculate conditional distributions of the Ys
*     for each value of Z
*
      TOL = 1.0D-12
      IFAIL = 0
      CALL BERCDB(S,P,PZ,S,TOL,CP,QUIET,IFAIL)
      IF (IFAIL.NE.0) RETURN
*
*     Next: which sites have known values, and which unknown? We take 
*     the opportunity to initialise PZ2 and P2 here, ready for use later.
*
      NUK = 0
      N1S = 0
      IALL = -1

      DO 100 IS=1,S
       IF (Y(IS).GE.0.0D0) THEN
        IDX(IS-NUK) = IS
        IF (Y(IS).GT.0.5D0) N1S = N1S + 1
       ELSE
        NUK = NUK + 1
        IDX(NK+NUK) = IS
       ENDIF
       PZ2(IS) = PZ(IS)
 100  CONTINUE
      PZ2(0) = PZ(0)
*
*     Trap anything trivial (shouldn't get here, but someone may
*     use the routine to calculate number of wet sites!)
*
*
      IF (NUK.EQ.0) THEN
       Z = N1S
       RETURN
      ENDIF
*
*     Now we condition on known values. First thing is to re-order the 
*     probabilities in P2 and CP to correspond to the new ordering of sites.
*     We use CP(0) to do the re-ordering (since its value will always
*     be 0, we can reset them on the last pass). 
*
      DO 120 IZ=1,S
       DO 125 IS=1,S
        CP(IS,0) = CP(IDX(IS),IZ)
 125   CONTINUE
       DO 126 IS=1,S
        CP(IS,IZ) = CP(IS,0)
        IF (IZ.EQ.S) CP(IS,0) = 0.0D0
 126   CONTINUE
       P2(IZ) = P(IDX(IZ))
 120  CONTINUE
*
*     Next, go through the known sites one by one, at each stage 
*     updating the conditional probabilities at the 
*     remaining sites and the distribution of Z.
*
      ZSTAR = Z
      SSTAR = S
      DO 150 IS=1,NK
       Y1 = Y(IDX(IS))
       CALL NEWPZ(Y1,P2(1),CP,SSTAR,PZ2,IFAIL)
       IF (IFAIL.NE.0) RETURN
       CALL Y1COND(Y1,CP,SSTAR,0,S+1-IS,TOL,IFAIL)
       IF (IFAIL.NE.0) RETURN
       SSTAR = SSTAR-1
       CALL NEWP(CP,PZ2,SSTAR,P2,TOL)
       ZSTAR = ZSTAR - Y1
 150  CONTINUE
*
*     The conditional distribution of ZSTAR, the number of wet
*     sites remaining, is now held in PZ2. It can take values in 
*     0,...,S-NK.
*
      U = ZBQLU01(0.0D0)
      CUMPRB = 0.0D0
      DO 155 IZ=0,S-NK
       CUMPRB = CUMPRB + PZ2(IZ)
       IF (U.LE.CUMPRB) THEN
        ZSTAR = IZ
        GOTO 159
       ENDIF
 155  CONTINUE
*
*     Finally, sample values at the unknown sites so that the 
*     required number are wet. Again, go through one by one, updating
*     probabilities at the remainder as we go.
*

 159  Z = N1S + ZSTAR
      DO 160 IS=NK+1,S
       Y1 = 0
       IF (ZBQLU01(0.0D0).LT.CP(1,ZSTAR)) Y1 = 1
       Y(IDX(IS)) = Y1
       CALL Y1COND(Y1,CP,SSTAR,ZSTAR,ZSTAR,TOL,IFAIL)
       IF (IFAIL.NE.0) RETURN
*
*     We now know what the *actual* value of ZSTAR is. It may
*     be that now we know what *all* the other values should be
*     (in which case we set them) ... 
*
       ZSTAR = ZSTAR - Y1
       N1S = N1S + Y1
       SSTAR = SSTAR - 1
       IF (SSTAR.EQ.Z-N1S) IALL = 1
       IF (N1S.EQ.Z) IALL = 0
       IF (IALL.GE.0) THEN
        DO 165 J=IS+1,S
         Y(IDX(J)) = IALL
 165    CONTINUE
        RETURN
       ENDIF
*
*     Otherwise, back for another go.
*
 160  CONTINUE
      
      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE BERCDB(S,P,PZ,ZMAX,TOL,CP,QUIET,IFAIL)
**********************************************************************
*     Given a set of S Bernoulli random variables, with success
*     probabilities contained in P, and a distribution for their sum: 
*     P(z) = PZ(z) (z=0,1,...,S), this routine calculates a valid set 
*     of conditional probabilities for each of the Bernoullis given 
*     a value of z. These conditional probabilities are output in the 
*     array CP, whose (s,z)th entry is the conditional probability for
*     the sth Bernoulli variable given the sum is z. P and PZ are of
*     physical dimensions S and 0:S respectively. The argument
*     ZMAX is the largest value of z for which the conditional 
*     probabilities will be calculated - to calculate all of them,
*     set ZMAX >= S. IFAIL is an error flag, with values as 
*     follows:
*                    0    No error
*                    1    Probabilities in PZ do not sum to 1
*                    2    Mean of z distribution does not correspond 
*                         to given probabilities
*                    3    Some input probabilities are outside (0,1)
*                    4    Some values of Z are outside (0,1)
*                    5    Given z distribution is inconsistent with
*                         specified marginal probabilities
*                    6    Ditto (different symptoms)
*     QUIET is a flag taking the value 0 to write all error messages
*     to screen, 1 to suppress some of them. TOL is the tolerance
*     within which sums etc. are expected to match.
**********************************************************************
      INTEGER S,ZMAX,IFAIL,QUIET
      DOUBLE PRECISION P(S),PZ(0:S),CP(S,0:S),TOL
**********************************************************************
*     Additional INTEGERs
*     ^^^^^^^^^^^^^^^^^^^
*     IS,IZ,J  counters
**********************************************************************
*     Additional INTEGERs
*     ^^^^^^^^^^^^^^^^^^^
*     ADJUST Flag indicating whether we're adjusting weights already
*            allocated or trying to find new ones.
*     PASS2  Flag indicating whether we're on a second adjustment pass
*     ADJRSN Reason for making an adjustment. 1: to try and sandwich
*            the desired row total. 2: to try and keep the allocated
*            probabilities away from 0 and 1.
**********************************************************************
      INTEGER IS,IZ
      INTEGER ADJUST,PASS2,ADJRSN
**********************************************************************
*     Additional DOUBLEs
*     ^^^^^^^^^^^^^^^^^^
*     PZSUM    Sum of the PZs
*     ZMEAN    Mean of the distribution PZ
*     PSUM     Sum of the Ps
*     PLEFT    Amount of probability left at each site.
*     LBSUM    Sum of lower bounds for weights in contingency table
*     UBSUM    Sum of upper bounds, ditto
*     RATIO    Ratio in which to split interval (LBSUM,UBSUM) to
*              achieve correct total.
*     SPZ      Array containing `survivor function' for Z - in fact
*              SPZ(z) = P(Z >= z).
*     LB       Array of lower bounds on weights.
*     UB       Array of upper bounds.
*     LBFREE   Array containing, for each value of z, the amount of
*              'slack' available between LBSUM and the desired 
*              row sum
*     ADDADJ   Potential extra weight that can be added to current
*              row to improve configuration
*     SUBADJ   Potential weight that can be subtracted from current
*              row to improve configuration
*     ADDSUM   Sum of ADDADJ
*     SUBSUM   Sum of SUBADJ
*     ADJSUM   Sum of adjustments actually made
*     LIMIT    Proportion of maximum possible adjustment to make
*     SLACK    Amount of probability to try and reserve in each cell
*     ADDPRP   Proportion of ADDSUM that we're actually going to use
*     SUBPRP   Ditto, SUBSUM
*     TMP1     Temporary storage
**********************************************************************
      DOUBLE PRECISION PZSUM,ZMEAN,PSUM,LBSUM,UBSUM,RATIO
      DOUBLE PRECISION SPZ(0:S),PLEFT(S),TMP1
      DOUBLE PRECISION LB(S,0:S),UB(S,0:S),LBFREE(0:S)
      DOUBLE PRECISION ADDADJ(S),SUBADJ(S),LIMIT,SLACK
      DOUBLE PRECISION ADDSUM,SUBSUM,ADJSUM,ADDPRP,SUBPRP
**********************************************************************
*     Additional CHARACTERs
*     ^^^^^^^^^^^^^^^^^^^^^
*     MESSAGE   - Error messages
**********************************************************************
      CHARACTER*255 MESSAGE(5)

      IFAIL = 0
      ADJUST = 0
*
*     First: calculate mean of Z distribution, check it and check 
*     that PZ sums to 1. Take the opportunity to define the joint 
*     probabilities when Z = 0 and Z = 1 as well (they're known!)
*
     
      IF ((PZ(0).LT.0.0D0).OR.(PZ(0).GT.1.0D0)) GOTO 993
      PZSUM = PZ(0)
      PSUM = 0.0D0
      ZMEAN = 0.0D0
      SPZ(0) = 1.0D0
      SPZ(1) = 1.0D0 - PZ(0)
      DO 100 IS = 1,S
       IF ((P(IS).LT.0.0D0).OR.(P(IS).GT.1.0D0)) GOTO 992
       IF ((PZ(IS).LT.0.0D0).OR.(PZ(IS).GT.1.0D0)) GOTO 993
       PSUM = PSUM + P(IS)
       PLEFT(IS) = P(IS)
       IF (IS.LT.S) SPZ(IS+1) = SPZ(IS) - PZ(IS)
       ZMEAN = ZMEAN + (DBLE(IS)*PZ(IS))
       PZSUM = PZSUM + PZ(IS)
       CP(IS,0) = 0.0D0
       CP(IS,S) = 1.0D0
 100  CONTINUE

      IF (DABS(PZSUM-1.0D0).GT.TOL) GOTO 990
      IF (DABS(ZMEAN-PSUM).GT.TOL) GOTO 991
*
*     Next: calculate how much probability needs to be reserved 
*     for high values of Z.
*
      IF (IFAIL.NE.0) RETURN
*
*     Now go through and calculate the required conditional 
*     probabilities. This is done by setting up a contingency table
*     for the joint distribution at each value of the sum IZ, and
*     finding appropriate weights for this table. We go from IZ = 1 to
*     S-1 (already done 0 and S which are trivial. 
 
      IZ = 1
 150  IF (PZ(IZ).LE.TOL**3) THEN
       IZ = IZ+1
       IF (IZ.GT.ZMAX) GOTO 190
       GOTO 150
      ENDIF
*
*     First step: calculate bounds on entries in the table. The 
*     current row total is equal to IZ*PZ(IZ). Store this in PZSUM, 
*     for want of anywhere better. Check for entries where the lower
*     and upper bounds are equal, so we can try and adjust them later.

*
       PZSUM = DBLE(IZ)*PZ(IZ) 
       LBSUM = 0.0D0
       UBSUM = 0.0D0
       ADJRSN = 0
       SLACK = 1.0D-2*PZ(IZ)/DBLE(S)
       DO 160 IS = 1,S
        CALL BOUNDS(PLEFT(IS),PZ(S),PZ(IZ),SPZ(IZ+1),LB(IS,IZ),
     +                                   UB(IS,IZ),LBSUM,UBSUM,0)
        IF (ADJRSN.EQ.0) THEN
         IF (IZ.LT.S-1) THEN
          IF ((UB(IS,IZ)-LB(IS,IZ).LT.SLACK).AND.
     +        ((UB(IS,IZ).GE.PZ(IZ)).OR.(LB(IS,IZ).LE.0.0D0)))
     +                                              ADJRSN = 2
         ENDIF
        ENDIF
*
*     Check bounds are compatible
*
        IF (LB(IS,IZ)-UB(IS,IZ).GT.TOL) GOTO 994
 160   CONTINUE

*
*     We will also try an adjustment if we're going to set all 
*     weights equal to the upper bound
*
       IF ((IZ.LT.S-1).AND.(UBSUM-PZSUM.LT.TOL)) ADJRSN = 2
*
*     Now check that we can make the required row total. If not, go 
*     back to the previous row and re-allocate something. If 
*     there's nothing to re-allocate, back again. If nothing at 
*     all to re-allocate, distribution of sum is incompatible with 
*     margins. NB we use -TOL as a check rather than zero because of 
*     the potential effects of rounding errors in the routine inputs.
*     We will also attempt to make the upper and lower bounds distinct,
*     so as to avoid having conditional probabilities of exactly 0 or
*     1 if possible (just 1 attempt required for this).
*
       ADJUST = IZ
       PASS2 = 0
       LIMIT = 0.05D0
 161   LBFREE(IZ) = PZSUM-LBSUM 
       IF (LBFREE(IZ).LT.-TOL) ADJRSN = 1
       IF (ADJRSN.NE.0) THEN
 169    ADJUST = ADJUST - 1
*
*     If we've made an unsuccessful pass, we either (i) decide
*     it's not critical, if we were just trying to allow some breathing
*     space between upper and lower bounds (ii) throw caution to the 
*     winds and allocate all the possible extra weight (iii) give up, 
*     if that didn't work either.
*
        IF (ADJUST.LE.0) THEN
         IF (ADJRSN.EQ.2) GOTO 164
         PASS2 = PASS2 + 1
         IF ((PASS2.EQ.2).OR.(IZ.EQ.1)) GOTO 996
         IF (PASS2.EQ.1) LIMIT = 0.0D0
         ADJUST = IZ-1
        ENDIF
        ADDSUM = 0.0D0
        SUBSUM = 0.0D0
        DO 170 IS=1,S
*
*     What can be added? Can't go further than the current upper
*     bound, or than the amount of probability left at this site;
*     and there's no point in adding any more probability than will
*     take the current lower bound to zero.
*
         ADDADJ(IS)=DMIN1(UB(IS,ADJUST)-CP(IS,ADJUST),PLEFT(IS),
     +                                              LB(IS,IZ))
         ADDADJ(IS) = DMAX1(ADDADJ(IS),0.0D0)
         ADDSUM = ADDSUM + ADDADJ(IS)
*
*     And what can be subtracted to compensate? Can't go further
*     than the current lower bound
*
         SUBADJ(IS) = DMIN1(CP(IS,ADJUST)-LB(IS,ADJUST),
     +                               SPZ(IZ+1)-PLEFT(IS))
         SUBADJ(IS) = DMAX1(SUBADJ(IS),0.0D0)
         SUBSUM = SUBSUM + SUBADJ(IS)
 170    CONTINUE
*
*     Now we add on a proportion of those weights - not all, because
*     we want to keep away from the boundaries of the feasible region
*     as far as possible.
*
        ADJSUM = (1.0D0-LIMIT)*DMIN1(ADDSUM,SUBSUM)
        IF (ADJSUM.LE.0.0D0) GOTO 169

        ADDPRP = ADJSUM/ADDSUM
        SUBPRP = ADJSUM/SUBSUM
        DO 175 IS = 1,S
         IF (ADDADJ(IS).GT.0.0D0) THEN
          TMP1 = ADDPRP*ADDADJ(IS)
          CP(IS,ADJUST) = CP(IS,ADJUST) + TMP1
          PLEFT(IS) = PLEFT(IS) - TMP1
         ELSEIF (SUBADJ(IS).GT.0.0D0) THEN
          TMP1 = SUBPRP*SUBADJ(IS) 
          CP(IS,ADJUST) = CP(IS,ADJUST) - TMP1
          PLEFT(IS) = PLEFT(IS) + TMP1
         ENDIF
         CALL BOUNDS(PLEFT(IS),PZ(S),PZ(IZ),SPZ(IZ+1),LB(IS,IZ),
     +                                   UB(IS,IZ),LBSUM,UBSUM,1)
 175    CONTINUE
*
*     Reset ADJRSN - if it was 2 and that didn't work, there's nothing
*     to be done. If it was 1 then we'll go back and check again.
* 
        ADJRSN = 0
        GOTO 161
       ENDIF

*
*     If we get here then we've got some feasible limits. Calculate 
*     ratio in which to divide the (LB,UB) interval for each entry
*     and set weights accordingly.
*
 164   IF (DABS(UBSUM-LBSUM).GT.TOL) THEN
        RATIO = LBFREE(IZ)/(UBSUM-LBSUM)
        IF ((PZSUM-UBSUM.GT.TOL).OR.(LBFREE(IZ).LT.-TOL)) GOTO 995
       ELSE
        RATIO = 0.5D0
       ENDIF
       DO 165 IS=1,S
        CP(IS,IZ) = LB(IS,IZ) + (RATIO*(UB(IS,IZ)-LB(IS,IZ)))
        PLEFT(IS) = PLEFT(IS) - CP(IS,IZ)
 165   CONTINUE
       IF (IZ.LT.MIN0(S-1,ZMAX)) THEN
        IZ = IZ+1
        GOTO 150
       ENDIF
*
*     Finally, convert joint probabilities to conditionals (IZ from 1
*     to S is correct - row zero is all 0 anyway). The DMIN1
*     deals with rounding errors that may have crept in and 
*     caused things to go over 1 - typically by a tiny amount, but
*     they have nasty consequences later! In addition, we need to 
*     check that things haven't gone negative from rounding. As if 
*     that wasn't enough, if PZ(IZ) is very small,
*     there's a chance that the row sum will be nowhere near Z due to
*     rounding errors. If so, we rescale. This will make negligible 
*     difference to the column sums (i.e. the Ps) because the rounding
*     will only be serious when the joint probabilities are tiny.
*
 190  DO 200 IZ = 1,S-1
       IF (PZ(IZ).GT.0.0D0) THEN
        PSUM = 0.0D0
        DO 210 IS=1,S
         CP(IS,IZ) = DMIN1(1.0D0,CP(IS,IZ)/PZ(IZ))
         PSUM = PSUM + CP(IS,IZ)
 210    CONTINUE
        TMP1 = DBLE(IZ)/PSUM
        DO 211 IS=1,S
         CP(IS,IZ) = DMIN1(1.0D0,CP(IS,IZ)*TMP1)
         IF (CP(IS,IZ).LT.0.0D0) CP(IS,IZ) = 0.0D0
 211    CONTINUE
       ENDIF
 200  CONTINUE

      RETURN
*
*     Error trapping
*
 990  WRITE(MESSAGE,21) PZSUM
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      IFAIL = 999
      RETURN
 991  WRITE(MESSAGE,22) PSUM,ZMEAN
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      IFAIL = 999
      RETURN
 992  WRITE(MESSAGE,23) 
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      IFAIL = 999
      RETURN
 993  WRITE(MESSAGE,24) 
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      IFAIL = 999
      RETURN
 994  WRITE(MESSAGE,25) LB(IS,IZ),UB(IS,IZ),IZ,S
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(3)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(4)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(5)),-1,0,0)
      IFAIL = 999
      RETURN
 995  WRITE(MESSAGE,26) PZSUM,LBSUM,UBSUM,IZ,S
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(3)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(4)),-1,0,0)
      IFAIL = 999
      RETURN
 996  IFAIL = 7
      IF (QUIET.EQ.0) THEN
       WRITE(MESSAGE,27) 
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      ENDIF
      RETURN

 21   FORMAT('****ERROR**** in routine BERCDB, probabilities in ',
     +     'PZ sum to ',F7.5,/,'instead of 1.')
 22   FORMAT('****ERROR**** in routine BERCDB, sum of individual ',
     +     'probabilities is ',/,F7.4,', but mean of PZ distribution',
     +     ' is ',F7.4)
 23   FORMAT('****ERROR**** in routine BERCDB, some sites have ',
     +     'probabilities',/,'outside the range [0,1].')
 24   FORMAT('****ERROR**** in routine BERCDB, some probabilities ',
     +     'in PZ distribution',/,'are outside the range [0,1].')
 25   FORMAT('****ERROR**** in routine BERCDB, distribution of sum ',
     +     'is incompatible with ',/,'given marginal probabilities -',
     +     ' I''m trying to sandwich something between a',/,
     +     'lower bound of ',F7.5,' and an upper bound of ',F7.5,
     +     '. The problem arose while',/,'trying to calculate ',
     +     'conditional probabilities when sum is ',I3,
     +     ' (maximum ',/,'is ',I3,').')
 26   FORMAT('****ERROR**** in routine BERCDB, distribution of sum ',
     +     'is incompatible with ',/,'given marginal probabilities -',
     +     1X,F7.4,' doesn''t lie in the range ',/,'(',F7.4,',',F7.4,
     +     '). The problem arose while trying to calculate ',/,
     +     'conditional probabilities when sum is ',I3,
     +     ' (maximum is ',I3,').')
 27   FORMAT('****ERROR**** in routine BERCDB, distribution of sum ',
     +     'is incompatible with ',/,'given marginal probabilities.')
      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE BOUNDS(PLEFT,PZS,PZZ,PZGTZ,LB,UB,LBSUM,UBSUM,ADJ)
***************************************************************************
*     Computes lower and upper bounds for joint probabilities in 
*     wet/dry model 22. Arguments:
*             PLEFT     Amount of probility remaing at current site.
*                       DOUBLE, input.
*             PZS       Probability that all sites are wet. DOUBLE, 
*                       input.
*             PZZ       Probability that z sites are wet. DOUBLE, input.
*             PZGTZ     Probability that more than z sites are wet.
*                       DOUBLE, input.
*             LB        Lower bound on joint probability. DOUBLE, input/
*                       output.
*             UB        Upper bound, similarly.
*             LBSUM     Sum of lower bounds for current value of Z.
*                       DOUBLE, in/out.
*             UBSUM     Sum of upper bounds, similarly.
*             ADJ       Flag indicating whether we're adjusting 
*                       values or computing them for the first time.
*                       INTEGER, input.
***************************************************************************
      DOUBLE PRECISION PLEFT,PZS,PZZ,PZGTZ,LB,UB,LBSUM,UBSUM
      INTEGER ADJ
***************************************************************************
*     Extra DOUBLEs
*     ^^^^^^^^^^^^^
*     TMP   Temporary storage
***************************************************************************
      DOUBLE PRECISION TMP
***************************************************************************
*     Lower bound:
*
      TMP = LB
      LB = DMAX1(0.0D0,PLEFT-PZGTZ)
      IF (ADJ.EQ.1) THEN 
       LBSUM = LBSUM + LB - TMP
      ELSE
       LBSUM = LBSUM + LB
      ENDIF
*
*     And upper
*
      IF (ADJ.EQ.1) TMP = UB
      UB = DMIN1(PLEFT-PZS,PZZ)
      IF (ADJ.EQ.1) THEN 
       UBSUM = UBSUM + UB - TMP
      ELSE
       UBSUM = UBSUM + UB
      ENDIF

      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE NEWPZ(Y1,P1,CP,S,PZ,IFAIL)
***************************************************************************
*     Given a set of S Bernoulli random variables such that P(Y1=1) = P1,
*     their sum (Z say) has pmf in PZ and whose conditional probabilities 
*     given Z=z are held in CP(s,z) (s=1,...,S;z=0,...,S), 
*     and given that the first variable is observed to have the value Y1, 
*     this routine modifies PZ to hold the conditional pmf of Z-Y1. IFAIL 
*     is an error flag
***************************************************************************
      INTEGER Y1,S,IFAIL
      DOUBLE PRECISION P1,CP(S,0:S),PZ(0:S)
***************************************************************************
*     Extra INTEGERs
*     ^^^^^^^^^^^^^^
*     IZ     Z counter
***************************************************************************
      INTEGER IZ
***************************************************************************
*     Extra DOUBLEs
*     ^^^^^^^^^^^^^
*     JNTPRB Joint probability derived from conditionals in CP
*     PY1    Probability of Y1 taking its observed value
*     PZSUM  Sum of probabilities in PZ - used to check for rounding 
*            errors
***************************************************************************
      DOUBLE PRECISION JNTPRB,PY1,PZSUM
***************************************************************************
*     Extra CHARACTERs
*     ^^^^^^^^^^^^^^^^
*     MESSAGE   For output of error messages
***************************************************************************
      CHARACTER*80 MESSAGE(2)
***************************************************************************
*     Start by computing probability of observed Y1 value
*
      IFAIL = 0
      IF (Y1.EQ.1) THEN
       PY1 = P1
      ELSE
       PY1 = 1.0D0 - P1
      ENDIF
*
*     At some stage PY1 is bound to be zero!
*
      IF (PY1.LE.0.0D0) THEN
       IFAIL = 1
       WRITE(MESSAGE,1)
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
       RETURN
      ENDIF
*
*     Z-Y1 can only range from 0 to S-1 since there are only S-1
*     sites left. Note that overwriting PZ is OK, since we only ever use 
*     either the current or the next value (neither of which have
*     yet been overwritten).
*
      PZSUM = 0.0D0
      DO 100 IZ=0,S-1
       IF (Y1.EQ.1) THEN
        JNTPRB = CP(1,IZ+Y1)*PZ(IZ+Y1)
       ELSE
        JNTPRB = (1.0D0-CP(1,IZ+Y1))*PZ(IZ+Y1)
       ENDIF
       PZ(IZ) = JNTPRB/PY1
       PZSUM = PZSUM + PZ(IZ)
 100  CONTINUE
*
*     Final correction for rounding errors (I've checked that PZSUM
*     *is* always close to 1!)
*
      DO 110 IZ=0,S-1
       PZ(IZ) = DMIN1(1.0D0,PZ(IZ)/PZSUM)
 110  CONTINUE

 1    FORMAT('****ERROR**** In routine NEWPZ, I''ve found an ',
     +'observation that has zero',/,'probability under the ',
     +'calculated joint distribution for Bernoulli variables.')

      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE NEWP(CP,PZ,S,P,TOL)
***************************************************************************
*     Given a set of S Bernoulli random variables whose sum has 
*     pmf in PZ and whose conditional probabilities given Z=z are in 
*     CP(s,z), this routine calculates the probabilities for the 
*     individual sites, storing them in P. TOL is used for determining 
*     whether the probabilities sum to zero.
***************************************************************************
      INTEGER S
      DOUBLE PRECISION CP(S,0:S),PZ(0:S),P(S),TOL
***************************************************************************
*     Extra INTEGERs
*     ^^^^^^^^^^^^^^
*     IZ     Z counter
*     IS     S counter
***************************************************************************
      INTEGER IS,IZ
***************************************************************************
*     Extra DOUBLEs
*     ^^^^^^^^^^^^^
*     PSUM   Sum of the new Ps    } Used for checking.
*     EZ     Expected value of Z  }
*     TMP    Temporary storage
***************************************************************************
      DOUBLE PRECISION PSUM,EZ,TMP
***************************************************************************
      PSUM = 0.0D0
      EZ = 0.0D0
*
*     Can ignore Z = 0, because that contributes 0 to the probabilities.
*
      DO 100 IS=1,S
       P(IS) = 0.0D0
       DO 110 IZ=1,S
        P(IS) = P(IS) + (CP(IS,IZ)*PZ(IZ))
 110   CONTINUE
       PSUM = PSUM + P(IS)
       EZ = EZ + (DBLE(IS)*PZ(IS))
 100  CONTINUE

*
*     Correct rounding errors (adjustment is small unless many
*     PZs are themselves very small)
*
      IF (PSUM.LT.TOL) RETURN
      TMP = EZ/PSUM
      DO 120 IS=1,S
       P(IS) = DMIN1(1.0D0,P(IS)*TMP)
 120  CONTINUE
      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE Y1COND(Y1,CP,S,FROM,TO,TOL,IFAIL)
***************************************************************************
*     Given a set of S Bernoulli random variables, whose conditional 
*     probabilities given their sum is z are held in CP(s,z) 
*     (s=1,...,S;z=0,...,S), and given that the first variable is 
*     observed to have the value Y1, this routine finds, for each value 
*     of z in the range (FROM,TO), a set of valid conditional probabilities 
*     for the remaining Ys. It replaces the appropriate part of the CP
*     array with the updated version taking into account Y1, in 
*     preparation for a new call. In particular, the updated probabilities
*     for site 2 is stored in column 1 etc., and the updated probabilities
*     for Z=z are stored in row Z-z1.
*
*     TOL is a small number used for checking equality of things. IFAIL 
*     is an error flag taking the value 1 if Y1 has zero probability, 
*     0 otherwise.
*
*     The task of this routine could actually be performed by a call
*     to BERCDB: however, there is no need for any checks here so 
*     it's quicker to program it directly as a special case.
***************************************************************************
      INTEGER S,Y1,FROM,TO,IFAIL
      DOUBLE PRECISION CP(S,0:S),TOL
***************************************************************************
*     Extra INTEGERs
*     ^^^^^^^^^^^^^^
*     IS     Site counter
*     IZ     Z counter
*     ZSTAR  Actual value of the sum of remaining variables
***************************************************************************
      INTEGER IS,IZ,ZSTAR
***************************************************************************
*     Extra DOUBLEs
*     ^^^^^^^^^^^^^
*     PISTAR Probability that Y1 takes its current value given Z
*     LB     Lower bound for entry in current array position
*     UB     Upper bound, similarly
*     LBSUM  Sum of lower bounds
*     UBSUM  Sum of upper bounds
*     ROWSUM Desired row sum
*     RATIO  Ratio in which to split lower and upper bounds to obtain
*            desired row sum
*     PSUM   Sum of probabilities in current row
*     TMP    Temporary storage
***************************************************************************
      DOUBLE PRECISION PISTAR,LB(S),UB(S),LBSUM,UBSUM
      DOUBLE PRECISION ROWSUM,RATIO,PSUM,TMP
***************************************************************************
*     Extra CHARACTERs
*     ^^^^^^^^^^^^^^^^
***************************************************************************
      CHARACTER*72 MESSAGE(3)
***************************************************************************
      DO 100 IZ=MAX0(FROM,1),MIN0(TO,S)
*
*     The desired row sum depends on the value of Y1
*
       ZSTAR=IZ-Y1
       IF (Y1.EQ.1) THEN
        PISTAR = CP(1,IZ)
       ELSE
        PISTAR = 1.0D0 - CP(1,IZ)
       ENDIF
*
*     Check that, if FROM=TO (so presumably we're sampling rather than
*     conditioning), Y1 doesn't have probability zero. If it does, it's
*     probably a programming or algorithm error. If FROM and TO are
*     unequal, the value of Y1 is inconsistent with the value of Z
*     we're currently considering so all the probabilities are zero.
*
       IF (PISTAR.LT.TOL) THEN
        IF (FROM.EQ.TO) GOTO 990
        DO 105 IS=2,S
         CP(IS-1,ZSTAR) = 0.0D0
 105    CONTINUE
*
*     Otherwise compute lower and upper bounds. We assume that 
*     any time LB>UB, it's a rounding error (this has been checked!).
*
       ELSE        
        ROWSUM = DBLE(ZSTAR)*PISTAR
        LBSUM = 0.0D0
        UBSUM = 0.0D0
        DO 110 IS=2,S
         LB(IS) = DMAX1(0.0D0,CP(IS,IZ)+PISTAR-1.0D0)
         UB(IS) = DMIN1(PISTAR,CP(IS,IZ))
         IF(LB(IS).GT.UB(IS)) THEN
          LB(IS) = (LB(IS)+UB(IS))/2.0D0
          UB(IS) = LB(IS)
         ENDIF
         LBSUM = LBSUM + LB(IS)
         UBSUM = UBSUM + UB(IS)
 110    CONTINUE
*
*     Again, if the computed bounds don't sandwich the desired row sum,
*     it's a rounding error. At this stage, UBSUM-LBSUM is guaranteed 
*     greater than 0.
*
        RATIO = UBSUM-LBSUM
        IF (RATIO.GT.TOL) THEN
         RATIO = (ROWSUM-LBSUM)/RATIO
        ELSE
         RATIO = 0.5D0
        ENDIF
        IF (RATIO.LT.0.0D0) THEN
         RATIO = 0.0D0
        ELSEIF (RATIO.GT.1.0D0) THEN
         RATIO = 1.0D0
        ENDIF
*
*     Allocate joint probabilities, convert to conditionals and
*     store in new position in array. This basically amounts to a 
*     shift left and (possibly) up - we've finished with this 
*     portion of the array so it's not a problem. NB correction for
*     rounding here.
*
        PSUM = 0.0D0
        DO 120 IS=2,S
         TMP = (LB(IS) + (RATIO*(UB(IS)-LB(IS))))/PISTAR
         IF (TMP.LT.0.0D0) THEN
          TMP = 0.0D0
         ELSEIF (TMP.GT.1.0D0) THEN
          TMP = 1.0D0
         ENDIF
         CP(IS-1,ZSTAR) = TMP
         PSUM = PSUM + TMP
 120    CONTINUE
*
*     ... and a final rounding correction, to ensure the correct row
*     sum
*
        IF (PSUM.GT.0.0D0) THEN
         TMP = DBLE(ZSTAR)/PSUM
         DO 121 IS=1,S-1
          CP(IS,ZSTAR) = DMIN1(1.0D0,CP(IS,ZSTAR)*TMP)
 121     CONTINUE
        ENDIF
       ENDIF
 100  CONTINUE

      RETURN

*
*     Error trapping beyond here
*
 990  WRITE(MESSAGE,1)
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(3)),-1,0,0)
      IFAIL = 1
      RETURN

  1    FORMAT('****ERROR**** In routine Y1COND, I''ve found an ',
     +'observation',/,'that has zero probability under the ',
     +'calculated joint',/,'distribution for Bernoulli variables.')
      END
***************************************************************************
***************************************************************************
***************************************************************************
      FUNCTION RTRNOR(MU,SIGMA,TAU,Y)
***************************************************************************
*
*     Returns a pseudo-random number from a truncated normal distribution.
*     Inputs (all double precision):
*
*     MU        Mean of the original (non-truncated) distribution
*     SIGMA     Standard deviation of the original distribution
*     TAU       The truncation point
*     Y         If equal to 1, the returned value will be above the
*               truncation point; if -1, it will be below. 
*
*     Algorithm is as follows: if the probability of being above / below 
*     the threshold is not too small, just sample from the unthresholded 
*     distribution until we get a result in the required range. Otherwise, 
*     use rejection sampling with an exponential proposal distribution 
*     whose parameter is chosen to match the gradient of the truncated 
*     normal at the threshold. 
*
***************************************************************************
      DOUBLE PRECISION MU, SIGMA, TAU, Y
      DOUBLE PRECISION RTRNOR
***************************************************************************
*     Additional INTEGERs
*     -------------------
*     NCALLS    Number of calls made to ZBQLNOR
***************************************************************************
      INTEGER NCALLS
***************************************************************************
*     Additional DOUBLEs
*     ------------------
*     LAMBDA    Parameter of exponential proposal distribution for rejection
*               sampling
*     LOGA      Log of scaling factor to match the values of the proposal and
*               target densities at TAU
*     LOGRATIO  Log ratio of target and proposal densities at sampled value
*     LOGF      Log target density at sampled value
*     LOGG      Log proposal density at sampled value
*     TMP       Temporary storage
*     ZBQLNOR   Normal random number generator
*     ZBQLEXP   Exponential random number generator
*     ZBQLU01   Uniform random number generator
*
***************************************************************************
      DOUBLE PRECISION LAMBDA, LOGA, LOGRATIO, LOGF, LOGG, TMP
      DOUBLE PRECISION ZBQLNOR, ZBQLEXP, ZBQLU01
***************************************************************************
*     Additional CHARACTERs
*     ---------------------
*     MESSAGE   For error messages
***************************************************************************
      CHARACTER MESSAGE*255
*
*       The next condition is true if Y=1, TAU < MU+SIGMA or if 
*       Y=-1, TAU>MU-SIGMA - in either of these cases, brute force 
*       and ignorance is fine
*
      NCALLS = 0
      IF (Y*(MU+(SIGMA*Y)-TAU).GT.0.0D0) THEN 
       RTRNOR = ZBQLNOR(MU,SIGMA)
       NCALLS = NCALLS + 1
 10    IF (Y*(RTRNOR-TAU) < 0.0D0) THEN
        RTRNOR = ZBQLNOR(MU,SIGMA)
        NCALLS = NCALLS + 1
        GOTO 10
       ENDIF
      ELSE
*
*       Otherwise use rejection. Proposal distribution is exponential with
*       parameter chosen to maximise the acceptance rate. In some situations, 
*       MU can be extremely large (e.g. when conditioning on known values at 
*       specific sites with very high inter-site correlations), and this can 
*       affect the arithmetic which is therefore better done on a log scale. 
*       The usual rule is: to sample from a density proportional to F, find a 
*       density G such that G(x) > A*F(x) everywhere; then sample X from G 
*       and accept with probability A*F(X)/G(X). As implemented below, we 
*       evaluate log(A), log(F) and log(G)
*
       LAMBDA = Y*(TAU-MU)
       LAMBDA = LAMBDA + DSQRT((LAMBDA*LAMBDA)+(4.0D0*SIGMA*SIGMA))
       LAMBDA = LAMBDA / (2.0D0*SIGMA*SIGMA)
       LOGA = DLOG(LAMBDA) + (0.5D0*((SIGMA*LAMBDA)**2)) - 
     +                   (LAMBDA*(Y*(MU-TAU)+(LAMBDA*SIGMA*SIGMA))) 
 20    RTRNOR = TAU+(Y*ZBQLEXP(1.0D0/LAMBDA))
       LOGF = -0.5D0*((RTRNOR-MU)/SIGMA)**2
       LOGG = DLOG(LAMBDA) - (Y*LAMBDA*(RTRNOR-TAU))
       LOGRATIO = LOGA + LOGF - LOGG
       IF (LOGRATIO.GT.0.0D0) GOTO 99
       IF (DLOG(ZBQLU01(0.0D0)).GT.LOGRATIO) GOTO 20
      ENDIF
*
*       If we've had an odd number of calls to ZBQLNOR, call it again
*       to clean up any spare deviates that are sitting there (see 
*       comments above in routine WDLG). 
*
      IF (2*(NCALLS/2).NE.NCALLS) TMP = ZBQLNOR(0.0D0,1.0D0)
      
      IF (Y*(RTRNOR-TAU).GE.0.0D0) RETURN

 99   WRITE(MESSAGE,1)
      CALL INTPR(TRIM(MESSAGE),-1,0,0)
      RTRNOR = 1.0D101
      
 1    FORMAT('****ERROR**** Programming error in RTRNOR - contact REC')
 
      END
      
