******************************************************************************
*       PROGRAM AUTHOR:         RICHARD E. CHANDLER
*                               (email: r.chandler@.ucl.ac.uk)
*       AFFILIATION:            DEPARTMENT OF STATISTICAL SCIENCE
*                                UNIVERSITY COLLEGE LONDON
*                                 GOWER STREET
*                                  LONDON WC1E 6BT
*                                   UK
*       © UCL 2012-
******************************************************************************
*	This file contains a collection of miscellanous FORTRAN routines 
*       for doing mathematical functions and operations. 
*
*	REVISION HISTORY:
*	=================
*       JULY 2012       Created from routines that were originally in 
*                       rec_math.f, cutting out just the ones that are
*                       needed by GLIMCLIM and removing any write(*,...)
*                       statements for compatibility with Windows R Gui.  
******************************************************************************
******************************************************************************
******************************************************************************
      FUNCTION CDFBVN(H,K,RHO,TOL)
*
*       Calculates the upper right tail area for a bivariate normal
*       distribution:
*       
*       P(Z1 > H, Z2 > K)
*
*       where Z1 and Z2 are standard normal random variables with
*       correlation RHO. The algorithm is from T.G. Donnelly, 
*       Comm. ACM 16(10), p.638 (1973). 
*
*       TOL is the (absolute) precision required
*     
      DOUBLE PRECISION H, K, RHO, TOL, CDFBVN
      DOUBLE PRECISION TWOPI, B, GH, GK, RR, CDFNOR
      DOUBLE PRECISION H2, A2, H4, EX, W2, AP, S2, SP, S1, SN, SQR
      DOUBLE PRECISION CONV, WH, WK, GW, T, G2, CONEX, CN, SGN
      INTEGER IS
      CHARACTER*255 MESSAGE

      IF (DABS(RHO).GT.1.0D0) THEN
       WRITE(MESSAGE,1)
       CDFBVN = -1.0D0
       RETURN
      ENDIF
      IF (.NOT.(TOL.GT.0.0D0)) THEN
       WRITE(MESSAGE,2)
       CDFBVN = -1.0D0
       RETURN
      ENDIF
*
*     May as well code up the "easy" values explicitly.
*
      GH = CDFNOR(H,0.0D0,1.0D0,1)
      GK = CDFNOR(K,0.0D0,1.0D0,1)
      IF (DABS(RHO).LT.TOL) THEN
       CDFBVN = GH*GK
       RETURN
*
*     If the latent correlation is 1, the joint probability is obviously 
*     the smaller of the two marginal probs.
*
      ELSEIF (DABS(1.0D0-RHO).LT.TOL) THEN
       CDFBVN = DMIN1(GH,GK)
       RETURN
*
*     If the correlation is -1 then we have Z2=-Z1; thus 
*     P(Z1>H,Z2>K) = P(Z1>H,Z1<-K) = P(H < Z1 < -K) = 
*     Phi(-K)-Phi(H)=1-Phi(K)-Phi(H) [or zero]. But GH=1-Phi(H) 
*     and GK=1-Phi(K); expression programmed below follows 
*     straightforwardly from this.
*
      ELSEIF (DABS(1.0D0+RHO).LT.TOL) THEN
       CDFBVN =  DMAX1(GH+GK-1.0D0,0.0D0)
       RETURN
      ENDIF
*
*     For all other values we're going to have to do it the hard way.
*
      TWOPI = 8.0D0*DATAN(1.0D0)
      B = 0.0D0
      GH = GH / 2.0D0
      GK = GK / 2.0D0
      RR = 1.0D0 - (RHO*RHO)
      SQR = DSQRT(RR)
      CONV = TWOPI*TOL
      IF (DABS(H).LT.TOL) THEN
       IF (DABS(K).LT.TOL) THEN
        B = 0.25D0 + (DATAN(RHO/SQR) / TWOPI)
        GOTO 350
       ELSE
        B = GK
        IS = -1
        GOTO 340
       ENDIF
      ELSE
       IF (DABS(K).LT.TOL) THEN
        B = GH
       ELSE
        B = GH + GK
        IF (H*K.LT.0.0D0) B = B - 0.5D0
       ENDIF
      ENDIF

      WH = -H
      WK = ( (K/H)-RHO ) / SQR
      GW = 2.0D0*GH
      IS = -1
      
 210  SGN = -1.0D0
      T = 0.0D0
      IF (DABS(WK).GT.TOL) THEN
       IF (DABS(DABS(WK)-1.0D0).LT.TOL) THEN
        T = WK * GW * (1.0D0-GW) / 2.0D0
       ELSE 
        IF ((DABS(WK)-1.0D0).GT.0.0D0) THEN
         SGN = -SGN
         WH = WH*WK
         G2 = CDFNOR(WH,0.0D0,1.0D0,0)
         WK = 1.0D0 / WK
         IF (WK.LT.0.0D0) B = B + 0.5D0
         B = B - ( (GW+G2) / 2.0D0 ) + ( GW*G2 )
        ENDIF
        H2 = WH*WH
        A2 = WK*WK
        H4 = H2 / 2.0D0
        EX = DEXP(-H4)
        W2 = H4*EX
        AP = 1.0D0
        S2 = AP - EX
        SP = AP
        S1 = 0.0D0
        SN = S1
        CONEX = DABS(CONV/WK)
        GOTO 290
 280    SN = SP
        SP = SP + 1.0D0
        S2 = S2 - W2
        W2 = W2 * H4 / SP
        AP = -AP * A2
 290    CN = AP * S2 / ( SN+SP )
        S1 = S1 + CN
        IF (DABS(CN).GT.CONEX) GOTO 280
        T = ( DATAN(WK) - (WK*S1) ) / TWOPI 
       ENDIF
       B = B + (SGN*T)
      ENDIF
 
 340  IF ( (IS.LT.0).AND.(DABS(K).GT.TOL) ) THEN
       WH = -K
       WK = ( (H/K) - RHO ) / SQR
       GW = 2.0D0 * GK
       IS = 1
       GOTO 210
      ENDIF
      
 350  IF (B.LT.0.0D0) B = 0.0D0
      IF (B.GT.1.0D0) B = 1.0D0
      CDFBVN = B
 
 1    FORMAT('****ERROR**** In routine CDFBVN, RHO must be in range ',
     +       '(-1,1)')
 2    FORMAT('****ERROR**** in routine CDFBVN, TOL must be strictly ',
     +       'positive')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      FUNCTION GAMMLN(X,IFAIL)
*
*     Returns the logarithm of the GammFn function evaluated at X (>0).
*     Algorithm taken from NUMERICAL RECIPES in FORTRAN (second
*     edition) by Press, Teukolsky, Vetterling and Flannery (CUP, 1992).
*
      DOUBLE PRECISION X,Z,TMP,TWOPI,C(6),GAMMLN
      INTEGER I,IFAIL
      CHARACTER MESSAGE*80
      DATA C/76.18009172947146D0,-86.50532032941677D0,
     +     24.01409824083091D0,-1.231739572450155D0,
     +     0.1208650973866179D-2,-0.5395239384953D-5/

      IFAIL = 0
*
*     Check input
*
      IF (.NOT.(X.GT.0.0D0)) THEN
       IFAIL = 1
       GAMMLN = -1.0D12
       WRITE(MESSAGE,1)
       CALL INTPR(TRIM(MESSAGE),-1,0,0)
       RETURN
      ENDIF
*
*     Now calculate approximation
*
      TWOPI = 8.0D0*DATAN(1.0D0)
      Z = X-1.0D0
      TMP = 1.000000000190015D0
      DO 10 I=1,6
       TMP = TMP + (C(I)/(Z + DBLE(I)))
 10   CONTINUE
      GAMMLN = DLOG(TMP) + (0.5D0*DLOG(TWOPI))
      TMP = Z + 5.5D0
      GAMMLN = GAMMLN + ((Z+0.5D0)*DLOG(TMP)) - TMP 
      RETURN
 
 1    FORMAT('****ERROR**** In GAMMLN, you''ve asked me to ',
     +'calculate ln[GammFn(x)] for x <= 0.')

      END
******************************************************************************
******************************************************************************
******************************************************************************
      FUNCTION DIGAMM(X,IFAIL)
*
*     DiGammFn function (the derivative of the log of the GammFn function).
*     We do it by numerically differentiating GAMMLN at X, using 
*     Ridders' method for accuracy (See Press et al., Section 5.7)
*
      INTEGER MAXIT
      PARAMETER(MAXIT=20)
      DOUBLE PRECISION X,GAMMLN,DIGAMM

      INTEGER IFAIL,I,J
      DOUBLE PRECISION CON,CON2,BIG,SAFE,DX,TOL
      DOUBLE PRECISION ERR,ERRT,FAC,A(MAXIT,MAXIT),HH
      CHARACTER MESSAGE*255

*
*     Check input
*
      DIGAMM = -1.0D12
      IF (.NOT.(X.GT.0.0D0)) THEN
       IFAIL = 998
       WRITE(MESSAGE,1)
       CALL INTPR(TRIM(MESSAGE),-1,0,0)
       RETURN
      ENDIF
*
*       Otherwise continue
*
      CON = 1.4D0
      CON2 = CON*CON
      DX = DMIN1(X,0.1D0)
      BIG = 1.0D12
      TOL = 1.0D-10
      SAFE = 2.0D0
      IFAIL = 0
*
*     Here's the differentiation, based on polynomial interpolation
*
      HH = DX
      A(1,1) = (GAMMLN(X+HH,IFAIL) - GAMMLN(X-HH,IFAIL))/(2.0D0*HH)
      ERR = BIG
      DO 10 I=2,MAXIT
       HH = HH/CON
       A(1,I) = (GAMMLN(X+HH,IFAIL) - GAMMLN(X-HH,IFAIL))/(2.0D0*HH)
       FAC = CON2
       DO 11 J=2,I
        A(J,I) = (A(J-1,I)*FAC-A(J-1,I-1))/(FAC-1.0D0)
        FAC = CON2*FAC
        ERRT = DMAX1(DABS(A(J,I)-A(J-1,I)),DABS(A(J,I)-A(J-1,I-1)))
        IF (ERRT.LE.ERR) THEN
         ERR = ERRT
         DIGAMM = A(J,I)
        ENDIF
 11    CONTINUE
       IF (DABS(A(I,I)-A(I-1,I-1)).GE.(SAFE*ERR)) GOTO 99
 10   CONTINUE

 99   IF (DABS(ERR).GT.TOL) THEN
       IFAIL = 999
       WRITE(MESSAGE,2) X
       CALL INTPR(TRIM(MESSAGE),-1,0,0)
      ENDIF
      RETURN

 1    FORMAT('****ERROR**** In DIGAMM, you''ve asked me to ',
     +'calculate Psi(x) for x <= 0.')
 2    FORMAT('***WARNING*** Can''t compute digamma(',F8.3,') ',
     +'to required accuracy')
      END
***************************************************************************
***************************************************************************
***************************************************************************
      FUNCTION TRIGAMM(X,IFAIL) 
*
*     Trigamma function (the second derivative of the log of the 
*     gamma function). This is Algorithm AS121 (Applied Stats 27(1), 
*     1978), with amendment from remark ASR88 (Applied Stats 40(3),
*     1991). Also, the constants A and B have been changed to increase
*     the accuracy of the evaluations and hence take advantage of the
*     double precision arithmetic
*
      DOUBLE PRECISION A, B, ONE, HALF, B2, B4, B6, B8, X, Y, Z 
      DOUBLE PRECISION OFLO, C 
      DOUBLE PRECISION TRIGAMM
      INTEGER IFAIL
      
      DATA A, B, ONE, HALF /1.0D-6, 1.0D1, 1.0D0, 0.5D0/ 
      DATA OFLO, C /8.4D37,1.0D10/
C
C B2, B4, B6 AND B8 ARE BERNOULLI NUMBERS 
C 
      DATA B2,B4,B6,B8  /0.166666667D0, -0.03333333333D0,
     +                   0.02380952381D0,-0.03333333333D0/ 
C 
C CHECK FOR POSITIVE VALUE OF X 
C 
      TRIGAMM = 0.0D0
      IFAIL = 1 
      IF (X.LE.0.0D0) RETURN 
      IFAIL = 2 
      IF ((X.LT.ONE/(OFLO ** HALF)).OR.(X.GT.OFLO)) RETURN
      IFAIL = 0 
      Z = X 
C 
C USE SMALL VALUE APPROXIMATION IF X LE. A 
C 
      IF (Z.GT.A) GOTO 10 
      TRIGAMM = 1.0D0 / (Z * Z) 
      RETURN 
C 
C INCREASE ARGUMENT TO (X + I) .GE. B 
C 
 10   IF (Z.GE.B) GOTO 20 
      TRIGAMM = TRIGAMM + 1.0D0 / (Z * Z) 
      Z = Z + 1.0D0 
      GOTO 10 
C 
C APPLY ASYMPTOTIC FORMUIA IF ARGUMENT GE. B 
C 
 20   IF (Z.GE.C) GOTO 30 
      Y = ONE / (Z * Z) 
      TRIGAMM = TRIGAMM + (HALF * Y) +
     +                  (ONE + Y*(B2 + Y*(B4 + Y*(B6 + Y*B8)))) / Z 
      RETURN
 30   TRIGAMM = ONE/Z 
      RETURN 
      END
******************************************************************
******************************************************************
******************************************************************
      SUBROUTINE CHOLDEC(X,P,MXP,IFAIL)
*
*       Replaces the symmetric positive definite matrix X with its
*       Cholesky decomposition (i.e. X = U'.U where U is an upper
*       triangular matrix). On exit X contains U, with zeroes below
*       the diagonal. MXP is the physical dimension of X. If X is
*       not positive definite, on exit IFAIL contains the number 
*       of the first row where a problem was detected; otherwise,
*       IFAIL contains zero on exit.
*
      INTEGER P,I,J,K,MXP,IFAIL
      DOUBLE PRECISION X(MXP,MXP)

      IFAIL = 0
      DO 100 I=1,P
*
*       First the diagonal element ...
*
       DO 101 K=1,I-1
        X(I,I) = X(I,I) - (X(K,I)**2)
 101   CONTINUE
       IF (X(I,I).LT.0.0D0) THEN
        IFAIL = I
        RETURN
       ENDIF
       X(I,I) = DSQRT(X(I,I))
*
*       Then the off-diagonal elements on the Ith row (+ fill in the 
*       column as well)
*
       DO 102 J=I+1,P
        DO 103 K=1,I-1
         X(I,J) = X(I,J) - (X(K,I)*X(K,J))
 103    CONTINUE
        X(I,J) = X(I,J)/X(I,I)
        X(J,I) = 0.0D0
 102   CONTINUE
 100  CONTINUE

      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE CDSOLV(A,X,B,N,DECOMP,MXN,IFAIL)
*
*     Solve the set of equations Ax=b, where A is an NxN positive definite
*     matrix, and B and X are N-vectors. The maximum physical dimensions
*     of A, X and B are MXNxMXN, MXN and MXN respectively. If DECOMP=1
*     on entry, A will be replaced by the upper triangle of its Cholesky 
*     decomposition. If DECOMP=0, then A is taken to contain this upper 
*     triangle already.
*     
*
      INTEGER N,DECOMP,MXN,IFAIL
      DOUBLE PRECISION A(MXN,MXN),X(MXN),B(MXN)
      CHARACTER*255 MESSAGE(2)
*
*     Calculate Cholesky decomposition of A (this overwrites the original
*     matrix)
*
      IF (DECOMP.EQ.1) THEN
       CALL CHOLDEC(A,N,MXN,IFAIL)
       IF (IFAIL.NE.0) THEN
        WRITE(MESSAGE,1) IFAIL
        CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
        RETURN
       ENDIF
      ENDIF
*
*     Equations now read U'Ux=b, where U is upper triangular (& is stored 
*     in A). So first solve for Ux ...
*
      CALL UTSOLV(A,X,B,0,N,MXN,IFAIL)
      IF (IFAIL.NE.0) RETURN
*
*     Now we know what Ux is, so solve for x. According to Numerical 
*     recipes (page 90), it's legal to pass X over twice (once for input
*     and once for output).
*
      CALL UTSOLV(A,X,X,1,N,MXN,IFAIL)
      IF (IFAIL.NE.0) RETURN

 1    FORMAT(5X,'****ERROR**** in CDSOLV, matrix A isn''t positive ',
     +'definite',/5X,'(problem found on row ',I3,').')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE UTSOLV(A,X,B,UPPER,N,MXN,IFAIL)
*
*     Solves either the system of linear equations AX=B (if UPPER=1), or
*     A'X=b (if UPPER=0). A is a square NxN upper triangular matrix.
*     The maximum physical dimensions of A, X and B are MXNxMXN, MXN and 
*     MXN respectively. Only X is modified by the routine.
*
      INTEGER N,MXN,UPPER,IFAIL
      DOUBLE PRECISION A(MXN,MXN),X(MXN),B(MXN),TOL
      INTEGER I,J
      CHARACTER*255 MESSAGE(2)

      TOL = 1.0D-12
      IFAIL = 0
*
*     Here's the case A'x=b
*
      IF (UPPER.EQ.0) THEN
       DO 100 I=1,N
        IF (DABS(A(I,I)).LT.TOL) GOTO 990
        X(I) = B(I)
        DO 110 J=1,I-1
         X(I) = X(I) - (A(J,I)*X(J))
 110    CONTINUE
        X(I) = X(I)/A(I,I)
 100   CONTINUE
*
*     And here's the case Ax=b
*
      ELSEIF(UPPER.EQ.1) THEN
       DO 200 I=N,1,-1
        IF (DABS(A(I,I)).LT.TOL) GOTO 990
        X(I) = B(I)
        DO 210 J=I+1,N
         X(I) = X(I) - (A(I,J)*X(J))
 210    CONTINUE
        X(I) = X(I)/A(I,I)
 200   CONTINUE
      ENDIF
      RETURN

*
*     Error trapping
*
 990  IFAIL = I
      WRITE(MESSAGE,1) I
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      RETURN

 1    FORMAT(5X,'****ERROR**** in routine UTSOLV, the matrix A has a ',
     +'zero',/,5X,'diagonal element in row ',I3,'.') 
      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE MATINV(X,P,DECOMP,MXP,METHOD,IFAIL)
*
*       For matrix inversion. Arguments:
*
*       X      On entry, EITHER the matrix to be inverted (if
*              DECOMP is 1) OR an appropriate decomposition (if
*              DECOMP is 0). On exit, the required inverse.
*              DOUBLE PRECISION, input/output.
*       P      the dimension of X. INTEGER, input.
*       DECOMP The inversion methods are all based on a preliminary
*              `in-place' decomposition of X (the Cholesky, if 
*              METHOD=1 or the LU if otherwise). If, on entry, X already
*              contains this decomposition, put DECOMP=0; otherwise
*              put DECOMP=1. IF DECOMP=0 and METHOD=1 then, on 
*              on entry, the *upper triangle* of X should contain 
*              its Cholesky square root. INTEGER, input.
*       MXP    Physical dimension of X (i.e. the storage allocated).
*              INTEGER, input.
*       METHOD Selects the method for calculating the inverse.
*              If METHOD=1, the Cholesky decomposition is 
*              used. For any other value, inversion is based on the
*              LU decomposition. The Cholesky method is fast, but 
*              requires that the matrix to be inverted is positive 
*              definite (this is true of all covariance matrices,
*              for example). INTEGER, input.
*       IFAIL  Error flag
*
ccc
ccc     Extra dimensioning not needed for shared libary usage
ccc
c      INTEGER MAXDIM
c      PARAMETER (MAXDIM=500)
    
      INTEGER P,I,J,K,MXP,DECOMP,METHOD,IFAIL
      DOUBLE PRECISION X(MXP,MXP)
**********************************************************************
*       Extra INTEGERs
*       ^^^^^^^^^^^^^^
*       INDX    Bookkeeping array for inversion using LU decomposition
*       D       Required by routine LUDCMP
*       IFAIL   Error flag from CHOLDEC
**********************************************************************
      INTEGER INDX(MXP),D
**********************************************************************
*       Extra DOUBLEs 
*       ^^^^^^^^^^^^^
*       U       If METHOD=1, inverse of upper triangular matrix in 
*               Cholesky decomposition of X. Otherwise, temporary
*               storage for use in backsubstitution.
**********************************************************************
      DOUBLE PRECISION U(MXP,MXP)
**********************************************************************
*       Extra CHARACTERs 
*       ^^^^^^^^^^^^^^^^
*       MESSAGE Error messages
**********************************************************************
      CHARACTER*255 MESSAGE(3)

      IFAIL = 0
ccc
ccc     Dimension check not needed for shared library usage
ccc
C      IF (P.GT.MAXDIM) THEN
C       WRITE(MESSAGE,1) MAXDIM
C       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
C       CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
C       IFAIL = -1
C       RETURN
C      ENDIF
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
* METHOD 1: CHOLESKY DECOMPOSITION
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF (METHOD.EQ.1) THEN
*
*       First calculate the decomposition (if necessary)
*
       IF (DECOMP.EQ.1) THEN
        CALL CHOLDEC(X,P,MXP,IFAIL)
        IF (IFAIL.NE.0) THEN
         WRITE(MESSAGE,3) IFAIL
         CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
         CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
         CALL INTPR(TRIM(MESSAGE(3)),-1,0,0)
         RETURN        
        ENDIF
       ENDIF
*
*      Now invert the upper triangle
*         
       CALL UTINV(X,U,P,MXP,MXP)
* 
*       And required inverse is now U.U' - speed things up by ignoring 
*       the zero elements
*
       DO 110 I=1,P
        DO 111 J=I,P
         X(I,J) = 0.0D0
         DO 112 K = J,P
          X(I,J) = X(I,J) + (U(I,K)*U(J,K))
 112     CONTINUE
         X(J,I) = X(I,J)
 111    CONTINUE
 110   CONTINUE
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
* METHOD 2: LU DECOMPOSITION
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ELSE
*
*     Compute LU decomposition of X
*
       IF (DECOMP.EQ.1) CALL LUDCMP(X,P,MXP,INDX,D,IFAIL)
       IF (IFAIL.NE.0) RETURN
*
*     Now do backsubstitution a column at a time (see Press et
*     al., p.40). We store the RHS of the required system of 
*     equations for each column, in the corresponding column of U;
*     then solve to replace this column with the corresponding
*     column of the inverse.
*
       DO 200 J=1,P
        DO 201 I=1,P
         U(I,J) = 0.0D0
 201    CONTINUE
        U(J,J) = 1.0D0
        CALL LUBKSB(X,P,MXP,INDX,U(1,J))
 200   CONTINUE
*
*     And copy U into X
*
       DO 210 I=1,P
        DO 211 J=1,P
         X(I,J) = U(I,J)
 211    CONTINUE
 210   CONTINUE
      ENDIF
      RETURN
c 1    FORMAT('****ERROR**** in MATINV, insufficient local storage ',
c     +'has been allocated.',/,
c     +'Increase the parameter MAXDIM in the source code (currently ',
c     +I4,').')
 3    FORMAT(5X,'****ERROR**** in MATINV, matrix isn''t positive ',
     +'definite so Cholesky',/5X,
     +'decomposition (METHOD=1) won''t work. Problem found on row ',I3,
     +'.',/5X,'If matrix is correct, try METHOD=2.')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE LUDCMP(A,N,NP,INDX,D,IFAIL)
*
*       To find the LU decomposition of the N*N matrix A. This
*       is pretty much copied from 'Numerical Recipes in
*       FORTRAN' (Press et al, 1989, CUP), pp. 35-36, so no 
*       comments.
*
      DOUBLE PRECISION TINY,A,VV,AAMAX,SUM,DUM
      INTEGER N,NP,NMAX,INDX,D,IMAX,IFAIL
      INTEGER I,J,K
      CHARACTER*255 MESSAGE(2)
      PARAMETER (NMAX=500,TINY=1.0D-10)
      DIMENSION A(NP,NP),INDX(NP),VV(NMAX)
      
      IFAIL = 0
*
*      Check there's enough storage in here ...
*
      IF (NMAX.LT.N) GOTO 990
*
*      ... and carry on 
*
      D = 1
      IMAX = 1
      DO 12 I=1,N
       AAMAX = 0.0D0
       DO 11 J=1,N
        IF (DABS(A(I,J)).GT.AAMAX) AAMAX = DABS(A(I,J))
 11    CONTINUE
       IF (DABS(AAMAX).LE.TINY) THEN 
        WRITE(MESSAGE,2) I
        CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
        CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
        IFAIL = I
        RETURN
       ENDIF
       VV(I) = 1.0D0/AAMAX
 12   CONTINUE
      DO 19 J=1,N
       DO 14 I=1,J-1
        SUM = A(I,J)
        DO 13 K=1,I-1
         SUM = SUM - ( A(I,K)*A(K,J) )
 13     CONTINUE
        A(I,J) = SUM
 14    CONTINUE
       AAMAX = 0.0D0
       DO 16 I=J,N
        SUM = A(I,J)
        DO 15 K=1,J-1
         SUM = SUM - ( A(I,K)*A(K,J) )
 15     CONTINUE
        A(I,J) = SUM
        DUM = VV(I)*DABS(SUM)
        IF (DUM.GE.AAMAX) THEN
         IMAX = I
         AAMAX = DUM
        ENDIF
 16    CONTINUE
       IF (J.NE.IMAX) THEN
        DO 17 K=1,N
         DUM = A(IMAX,K)
         A(IMAX,K) = A(J,K)
         A(J,K) = DUM
 17     CONTINUE
        D = -D
        VV(IMAX) = VV(J)
       ENDIF
       INDX(J) = IMAX
*
*       REC-added bit to warn of a singular matrix
*
       IF(DABS(A(J,J)).LE.TINY) THEN
        WRITE(MESSAGE,1) J
        IFAIL = J
        RETURN
       ENDIF
       IF (J.NE.N) THEN
        DUM = 1.0D0/A(J,J)
        DO 18 I=J+1,N
         A(I,J) = A(I,J) * DUM
 18     CONTINUE
       ENDIF
 19   CONTINUE
      RETURN

 990  WRITE(MESSAGE,3)
      CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
      IFAIL = -1
      RETURN

 1    FORMAT('****WARNING**** Singular matrix in LUDCMP (row ',I2,')')
 2    FORMAT(5X,'****ERROR**** All numbers in row ',I2,' of matrix ',
     +'passed ',/5X,
     +'to LUDCMP are too small for me to handle. Try rescaling.')
 3    FORMAT(5X,'****ERROR**** Insufficient local storage in LUDCMP.',
     +/5X,'Increase NMAX and recompile.')
      END
******************************************************************************
******************************************************************************
******************************************************************************
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
*
*       To solve the set of N linear equations A.X = B for X.
*       This is pretty much copied from 'Numerical Recipes in
*       FORTRAN' (Press et al, 1989, CUP), pp. 35-36, so no 
*       comments.
*
      DOUBLE PRECISION A,B,SUM
      INTEGER INDX,LL,II,N,NP,I,J
      DIMENSION A(NP,NP),INDX(NP),B(NP)
      II = 0
      DO 12 I=1,N
       LL = INDX(I)
       SUM = B(LL)
       B(LL) = B(I)
       IF (II.NE.0) THEN
        DO 11 J=II,I-1
         SUM = SUM - ( A(I,J)*B(J) )
 11     CONTINUE
       ELSEIF (SUM.NE.0.0D0) THEN
        II = 1
       ENDIF
       B(I) = SUM
 12   CONTINUE
      DO 14 I=N,1,-1
       SUM = B(I)
       DO 13 J=I+1,N
        SUM = SUM - ( A(I,J)*B(J) )
 13    CONTINUE
       B(I) = SUM/A(I,I)
 14   CONTINUE
      RETURN
      END
***************************************************************************
***************************************************************************
***************************************************************************
      SUBROUTINE UTINV(U,V,P,MXPU,MXPV)
*
*     To invert an upper triangular matrix. Arguments:
*
*     U     The matrix to be inverted. DOUBLE PRECISION, input.
*     V     The inverse. DOUBLE PRECISION, output.
*     P     The dimension of U. INTEGER, input.
*     MXPU  The storage allocated to U. INTEGER, input.
*     MXPV  The storage allocated to V. INTEGER, input. 
*
      INTEGER P,MXPU,MXPV
      DOUBLE PRECISION U(MXPU,MXPU),V(MXPV,MXPV)
***************************************************************************
*     Extra INTEGERs
*     ^^^^^^^^^^^^^^
*     I,J,K   Counters
***************************************************************************
      INTEGER I,J,K
*
*       The trick is to do all the elements in the right order ... 
*       use VU = I and work up columns of V
*
       DO 100 J=1,P
        V(J,J) = 1.0D0/U(J,J)
        DO 101 I=J+1,P
         V(I,J) = 0.0D0
 101    CONTINUE
        DO 102 I = J-1,1,-1
         V(I,J) = 0.0D0
         DO 103 K=I,J-1
          V(I,J) = V(I,J) - (V(I,K)*U(K,J))
 103     CONTINUE
         V(I,J) = V(I,J)/U(J,J)
 102    CONTINUE     
 100   CONTINUE

      END
******************************************************************************
******************************************************************************
******************************************************************************
      FUNCTION CDFNOR(X,MU,SIGMA,UPPER)
******************************************************************
*     Returns the cumulative distribution function (if UPPER = 0),
*     or upper tail area (UPPER=1) for a normal 
*     distribution with mean MU and standard deviation SIGMA, 
*     evaluated at X. The routine uses algorithm AS66 (Applied
*     Statistics 22, 1973, pp.424-7). Comparison with the 15-figure
*     tables in Abramowitz & Stegun (1965, Table 26.1) indicates
*     that the absolute error of this routine is always less than
*     10^-10. For |X| < 3, the relative error is also less than 
*     10^-10; however, for |X| = 5, the relative error in the tail
*     area is around 10^-4. 
******************************************************************
      DOUBLE PRECISION X,MU,SIGMA,CDFNOR
      DOUBLE PRECISION LTONE,UTZERO,HALF,ONE,CON,Z,Y
      INTEGER UPPER,IUP
      CHARACTER MESSAGE*255
*
*     LTONE is the value at which the lower tail area becomes 1 
*     to the accuracy of the machine. UTZERO is the value at 
*     which the upper tail area becomes zero to machine accuracy.
*     The values here are intended to be reasonably portable.
*
      DATA LTONE,UTZERO /12.0D0,25.0D0/
      DATA HALF,ONE,CON /0.5D0,1.0D0,1.28D0/

*
*     Error trapping
*
      IF ((UPPER.NE.0).AND.(UPPER.NE.1)) THEN
       WRITE(MESSAGE,1)
       CDFNOR = -1.0D0
       RETURN
      ENDIF

      IUP = UPPER
      Z = (X - MU)/SIGMA

*
*     Work on positive half-line.
*
      IF (Z.LT.0.0D0) THEN
       IUP = 1 - IUP
       Z = -Z
      ENDIF
*
*     Explicit calculation except in extreme tails
*
      IF ((Z.LT.LTONE).OR.((IUP.EQ.1).AND.(Z.LE.UTZERO))) THEN
       Y = HALF*Z*Z
       IF (Z.LE.CON) THEN
        CDFNOR = HALF - (Z * 
     +                (0.398942280444D0 - (0.399903438504D0 * Y /
     +                (Y + 5.75885480458D0 - (29.8213557808D0 /
     +                (Y + 2.62433121679D0 + (48.6959930692D0 /
     +                (Y + 5.92885724438D0))))))))
       ELSE
        CDFNOR = 0.398942280385D0*DEXP(-Y) /
     +         (Z - 3.8052D-8 + (1.00000615302D0 /
     +         (Z + 3.98064794D-4 + (1.98615381364D0 /
     +         (Z - 0.151679116635D0 + (5.29330324926D0 /
     +         (Z + 4.8385912808D0 - (15.1508972451D0 /
     +         (Z + 0.742380924027D0 + (30.789933034D0 /
     +         (Z + 3.99019417011D0)))))))))))
       ENDIF
*
*     Otherwise give up
*
      ELSE
       CDFNOR = 0.0D0
      ENDIF

      IF(IUP.EQ.0) CDFNOR = ONE - CDFNOR

 1    FORMAT(5X,'****ERROR**** In routine CDFNOR, UPPER ',
     +     'must be 0 or 1')
      END
******************************************************************
******************************************************************
******************************************************************
      FUNCTION QNORM(P,MU,SIGMA,IFAIL)
******************************************************************
*     Returns the Pth quantile of a normal distribution with
*     mean MU and standard deviation SIGMA. IFAIL is an error flag.
*     The routine is based on Algorithm AS241 (Applied Statistics
*     37, 1988, pp. 477-484). The relative error in evaluating 
*     quantiles of the standard normal distribution is less than
*     10^-15. Should be adequate for most purposes ...
******************************************************************
      DOUBLE PRECISION P,MU,SIGMA,QNORM
      DOUBLE PRECISION ONE,SPLIT1,SPLIT2,CONST1,CONST2
      DOUBLE PRECISION A(0:7),B(7),C(0:7),D(7),E(0:7),F(7),Q,R
      INTEGER I,IFAIL
      CHARACTER*255 MESSAGE(2)

      PARAMETER(ONE=1.0D0,SPLIT1=0.425D0,
     +     SPLIT2=5.0D0,CONST1=0.180625D0,CONST2=1.6D0)

*
*     Coefficients for P close to 0.5
*
      DATA (A(I),I=0,7) /3.3871328727963666080D0,
     +     1.3314166789178437745D2,1.9715909503065514427D3,
     +     1.3731693765509461125D4,4.5921953931549871457D4,
     +     6.7265770927008700853D4,3.3430575583588128105D4,
     +     2.5090809287301226727D3/
      DATA (B(I),I=1,7) /4.2313330701600911252D1,
     +     6.8718700749205790830D2,5.3941960214247511077D3,
     +     2.1213794301586595867D4,3.9307895800092710610D4,
     +     2.8729085735721942674D4,5.2264952788528545610D3/
*
*     Coefficients for P close to 0 or 1
*
      DATA (E(I),I=0,7) /6.65790464350110377720D0,
     +     5.46378491116411436990D0,1.78482653991729133580D0,
     +     2.96560571828504891230D-1,2.65321895265761230930D-2,
     +     1.24266094738807843860D-3,2.71155556874348757815D-5,
     +     2.01033439929228813265D-7/
      DATA (F(I),I=1,7) /5.99832206555887937690D-1,
     +     1.36929880922735805310D-1,1.48753612908506148525D-2,
     +     7.86869131145613259100D-4,1.84631831751005468180D-5,
     +     1.42151175831644588870D-7,2.04426310338993978564D-15/
*
*     Coefficients for intermediate values
*
      DATA (C(I),I=0,7) /1.42343711074968357734D0,
     +     4.63033784615654529590D0,5.76949722146069140550D0,
     +     3.64784832476320460504D0,1.27045825245236838258D0,
     +     2.41780725177450611770D-1,2.27238449892691845833D-2,
     +     7.74545014278341407640D-4/
      DATA (D(I),I=1,7) /2.05319162663775882187D0,
     +     1.67638483018380384940D0,6.89767334985100004550D-1,
     +     1.48103976427480074590D-1,1.51986665636164571966D-2,
     +     5.47593808499534494600D-4,1.05075007164441684324D-9/

*
*     Check requested probability
*
      IFAIL = 0
      IF ((P.LE.0.0D0).OR.(P.GE.1.0D0)) THEN
       WRITE(MESSAGE,1) P
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
       CALL INTPR(TRIM(MESSAGE(2)),-1,0,0)
       IF (P.LE.0.0D0) THEN
        QNORM = -1.0D100
       ELSE
        QNORM = 1.0D100
       ENDIF
       RETURN
      ENDIF
*
*     Check SIGMA > 0
*
      IF (SIGMA.LT.0.0D0) THEN
       IFAIL = 1
       WRITE(MESSAGE(1),2)
       CALL INTPR(TRIM(MESSAGE(1)),-1,0,0)
      ENDIF
      Q = P - 0.5D0
*
*     P `close' to 0.5
*
      IF (DABS(Q).LE.SPLIT1) THEN
       R = CONST1 - (Q*Q)
       QNORM = Q * (A(0) + (R * (A(1) + (R * (A(2) + (R * (
     +              A(3) + (R * (A(4) + (R * (A(5) + (R * (
     +              A(6) + (R * A(7))))))))))))))) / 
     +              (ONE + (R * (B(1) + (R * (B(2) + (R * (
     +              B(3) + (R * (B(4) + (R * (B(5) + (R * ( 
     +              B(6) + (R * B(7)))))))))))))))
      ELSE
       IF (Q.LT.0.0D0) THEN
        R = P
       ELSE
        R = ONE - P
       ENDIF
       R = DSQRT(-DLOG(R))
*
*     P `intermediate'
*
       IF (R.LE.SPLIT2) THEN
        R = R - CONST2
        QNORM = (C(0) + (R * (C(1) + (R * (C(2) + (R * (
     +           C(3) + (R * (C(4) + (R * (C(5) + (R * (
     +           C(6) + (R * C(7))))))))))))))) / 
     +           (ONE + (R * (D(1) + (R * (D(2) + (R * (
     +           D(3) + (R * (D(4) + (R * (D(5) + (R * ( 
     +           D(6) + (R * D(7)))))))))))))))
*
*     ... and P close to 0 or 1
*
       ELSE
        R = R - SPLIT2
        QNORM = (E(0) + (R * (E(1) + (R * (E(2) + (R * (
     +           E(3) + (R * (E(4) + (R * (E(5) + (R * (
     +           E(6) + (R * E(7))))))))))))))) / 
     +           (ONE + (R * (F(1) + (R * (F(2) + (R * (
     +           F(3) + (R * (F(4) + (R * (F(5) + (R * ( 
     +           F(6) + (R * F(7)))))))))))))))
       ENDIF
       IF (Q.LT.0.0D0) QNORM = -QNORM
      ENDIF
*
*     Convert to required normal distribution
*
      QNORM = (QNORM*SIGMA) + MU

 1    FORMAT(5X,'****ERROR**** In routine QNORM, probability ',
     +     'must be in open interval (0,1)',/5X,
     +     '- value requested is ',F8.3)
 2    FORMAT(5X,'****WARNING**** negative standard deviation ',
     +     'in routine QNORM')
      END

