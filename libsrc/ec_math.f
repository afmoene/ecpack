
      SUBROUTINE EC_M_ABCForm(a,b,c,Root1, Root2, AllReal)
C     Solves ax^2 + bx + c = 0
      INTEGER Re, Im
      PARAMETER(Re = 1, Im = 2)
      REAL*8 a, b, c, Discrim, Root1(Re:Im), Root2(Re:Im)
      LOGICAL AllReal
      Root1(Re) = 0.D0
      Root1(Im) = 0.D0
      Root2(Re) = 0.D0
      Root2(Im) = 0.D0
      Discrim = b*b - 4.D0*a*c
      AllReal = (Discrim .GE. 0.D0)
      IF (AllReal) THEN
        Root1(Re) = (-b + DSQRT(Discrim)) / (2.D0*a)
        Root2(Re) = (-b - DSQRT(Discrim)) / (2.D0*a)
      ELSE
        Root1(Re) = -b/(2.D0*a)
        Root2(Re) =  Root1(Re)
        Root1(Im) =  DSQRT(-Discrim)/(2.D0*a)
        Root2(Im) = -Root1(Im)
      END IF
      RETURN
      END






C
C From large time-series we extract mean quantities and co-variances.
C For both mean values and for covariances the tolerances are estimated.
C
      SUBROUTINE EC_M_Averag(x,NMax,N,MMax,M,Flag,
     &	Mean,TolMean,Cov,TolCov,MIndep,CIndep,Mok,Cok)
C
C From the given set of calibrated samples, calculate the averages,
C variances, covariances and their tolerances.
C
C input : x : array with physical dimensions (NMax,MMax)
C	      First index counts quantities; second counter
C	      counts samples. Only the first N quantities
C	      and the first M samples are used.
C	  flag : LOGICAL array with physical dimensions NMAx by MMAx,
C	      out of which only N by M are used. If flag(j,i) is true, then
C	      quantity j in sample i is not ok.
C
C output : Mean : array with physical dimension NMax out of
C		   which only N are used. Mean contains the average value
C		   of array x. Only samples with Flag = 0 are used.
C	   TolMean : tolerance of Mean, defined as
C		       2 * sigma / sqrt(NIndep)
C		   where sigma is the standarddeviation of the quantities
C		   and where the number of independent samples is estimated
C		   as twice the number of sign-changes of the fluctuations
C		   of the respective quantities around their means.
C	   Cov : REAL*8 (NMax,NMax) : covariances
C	   TolCov : tolerances of Cov, estimated as tolerances of mean.
C	   MIndep : INTEGER(NMAx) : Number of independent samples
C		   in time series from which means are calculated
C	   CIndep : INTEGER(NMAx,NMAx) : Number of independent samples
C		   in time series from which covariances are calculated
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER NMax,N,MMax,M,i,j,k,NChange(NNMax,NNMax),
     &	MIndep(NMax),CIndep(NMax,NMax),NMin(NNMax,NNMax),
     &  Mok(NMax),Cok(NMax,NMax),NPlus(NNMax,NNMax)
      LOGICAL Flag(NMax,MMax)
      REAL*8 x(NMax,MMax),RawMean(NNMax),Mean(NMax),TolMean(NMax),
     &	Cov(NMax,NMax),TolCov(NMax,NMax),xPrime(NNMax),
     &	PrevPrime(NNMax),dTolCov(NNMax,NNMax),PSwap
C
C Initialise the arrays
C
      DO i=1,NMax
	RawMean(i) = 0.D0
	Mean(i) = 0.D0
	TolMean(i) = 0.D0
	Mok(i) = 0
	DO j=1,NMax
	  Cov(i,j) = 0.D0
	  TolCov(i,j) = 0.D0
	  Cok(i,j) = 0
	ENDDO
      ENDDO
C
C Find a rough estimate for the mean
C
      DO i=1,M
	DO j=1,N
	  IF (.NOT.Flag(j,i)) THEN
	    Mok(j) = Mok(j) + 1
	    RawMean(j) = RawMean(j) + x(j,i)
	  ENDIF
	ENDDO
      ENDDO

      DO j=1,N
	IF (Mok(j).GT.0) RawMean(j) = RawMean(j)/DBLE(Mok(j))
      ENDDO
C
C Find final estimates for mean
C
      DO i=1,M
	DO j=1,N
	  IF (.NOT.Flag(j,i)) THEN
	    Mean(j) = Mean(j) + (x(j,i) - RawMean(j))
	  ENDIF
	ENDDO
      ENDDO

      DO j=1,N
	IF (Mok(j).GT.0) THEN
	  Mean(j) = Mean(j)/DBLE(Mok(j))
	  Mean(j) = Mean(j) + RawMean(j)
	ENDIF
      ENDDO
C
C Find (co-)variances and from them the tolerances of the mean
C
      DO j=1,N
	NChange(j,j) = 0
        NMin(j,j) = 0
        NPlus(j,j) = 0
      ENDDO

      DO i=1,M
	DO j=1,N
	  IF (.NOT.Flag(j,i)) THEN
	    xPrime(j) = x(j,i) - Mean(j)
            IF (xPrime(j).LT.0.D0) NMin(j,j) = NMin(j,j)+1
            IF (xPrime(j).GT.0.D0) NPlus(j,j) = NPlus(j,j)+1
	    IF (i.GT.1) THEN
	      IF (PrevPrime(j)*xPrime(j).LE.0.D0)
     &		NChange(j,j) = NChange(j,j) + 1
	    ENDIF
	    PrevPrime(j) = xPrime(j)
	  ENDIF
	ENDDO


	DO j=1,N
	  DO k=j,N
	    IF (.NOT.(Flag(j,i).OR.Flag(k,i))) THEN
	      Cok(j,k) = Cok(j,k) + 1
	      Cov(j,k) = Cov(j,k) + xPrime(j)*xPrime(k)
	    ENDIF
	  ENDDO
	ENDDO
      ENDDO

      DO j=1,N
	DO k=j,N
	  IF (Cok(j,k).GT.0) THEN
	     Cov(j,k) = Cov(j,k)/DBLE(Cok(j,k))
	  ENDIF
	ENDDO
      ENDDO

      DO j=1,N
	DO k=1,N
	  IF (k.LT.j) Cov(j,k) = Cov(k,j)
	ENDDO
      ENDDO

      DO j=1,N
        PSwap = 2.D0*DBLE(NMin(j,j))*DBLE(NPlus(j,j))/
     &    DBLE(NMin(j,j) + NPlus(j,j))**2.D0
	MIndep(j) = ANINT(DBLE(NChange(j,j))/PSwap) - 1
	MIndep(j) = MAX(MIndep(j),1)
	IF (Cok(j,j) .GT. 0) THEN
	   TolMean(j) = 2.D0*(Cov(j,j)/DBLE(MIndep(j)))**0.5D0
	ENDIF
      ENDDO
C
C Find tolerances for (co-)variances
C
      DO i=1,M
        DO j=1,N
	  IF (.NOT.Flag(j,i)) THEN
	    xPrime(j) = x(j,i) - Mean(j)
	  ENDIF
	ENDDO
	DO j=1,N
	  DO k=j,N
	    IF (.NOT.(Flag(j,i).OR.Flag(k,i))) THEN
	      dTolCov(j,k)=xPrime(j)*xPrime(k)-Cov(j,k)
	      TolCov(j,k)=TolCov(j,k)+(dTolCov(j,k))**2
	    ENDIF
	  ENDDO
	ENDDO
      ENDDO

      DO j=1,N
	DO k=j,N
C
C Calculate the standard deviation of the instantaneous contributions
C to the covariances.
C
	  IF (Cok(j,k).GT.0) THEN
	     TolCov(j,k)=TolCov(j,k)/DBLE(Cok(j,k))
C
C Here we estimate the number of independent contributions to the
C covariance by counting the least number of independent contributions
C to either of the factors.
C
	     CIndep(j,k) = MIN(MIndep(j),MIndep(k))
	     TolCov(j,k) = TolCov(j,k)/DBLE(CIndep(j,k))
C
C Tolerance is defined as 2*sigma, where sigma is standard deviation
C
	     TolCov(j,k) = 2.D0*(TolCov(j,k))**0.5
	  ENDIF
	ENDDO
      ENDDO

      DO j=1,N
	DO k=1,(j-1)
	  TolCov(j,k) = TolCov(k,j)
	  CIndep(j,k) = CIndep(k,j)
	ENDDO
      ENDDO

      RETURN
      END





C
C ########################################################################
C
C A routine to calculate the value of certain simple basefunctions
C
C ########################################################################
C





      REAL*8 FUNCTION EC_M_BaseF(x,FType,Order,C)
C
C Purpose : Calculate simple functions. Currently implemented:
C  - Ordinary polynomials
C  - Polynomials in the natural logarithm of x
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER FType,Order,i
      REAL*8 x,C(0:Order),Dum,LogX

      Dum = 0.D0
      IF (FType.EQ.NormPoly) THEN
	DO i=0,Order
	  Dum = Dum + C(i)*x**i
	ENDDO
      ELSE IF (FType.EQ.LogPoly) THEN
	LogX = LOG(x)
	DO i=0,Order
	  Dum = Dum + C(i)*LogX**i
	ENDDO
      ENDIF

      EC_M_BaseF = Dum

      RETURN
      END

      SUBROUTINE EC_M_Cardano(Poly, Root, AllReal)
C Uses the Cardano solution to solve exactly: ax^3 + bx^2 + cx + d = 0
C See e.g. Abramowitz and Stegun. Here "a" is not allowed to be zero.
      INTEGER Re, Im
      PARAMETER(Re = 1, Im = 2)
      REAL*8 Poly(0:3), a(0:2), q, r, Discrim, t1, t2,
     &  s1(Re:Im), s2(Re:Im), Root(Re:Im,3), Radius, Phi
      LOGICAL AllReal
      INTEGER i
      Root(Im,1) = 0.D0
      Root(Im,2) = 0.D0
      Root(Im,3) = 0.D0
      DO i = 0,2
        a(i) = Poly(i)/Poly(3)
      ENDDO
      Root(Re,1) = -a(2)/3.D0
      Root(Re,2) = Root(Re,1)
      Root(Re,3) = Root(Re,1)
      q = a(1)/3.D0 - (a(2)*a(2))/9.D0
      r = (a(1)*a(2) - 3.D0*a(0))/6.D0 - a(2)*a(2)*a(2)/27.D0
      Discrim = q*q*q + r*r
      AllReal = (Discrim .LE. 0.D0)
      IF (.NOT. AllReal) THEN
        t1 = r + DSQRT(Discrim)
        t2 = r - DSQRT(Discrim)
        s1(Re) = DSIGN(1.D0,t1)*DABS(t1)**(1.D0/3.D0)
        s2(Re) = DSIGN(1.D0,t2)*DABS(t2)**(1.D0/3.D0)
        Root(Re,1) =  Root(Re,1) +  s1(Re) + s2(Re)
        Root(Re,2) =  Root(Re,2) - (s1(Re) + s2(Re))/2.D0
        Root(Im,2) =  (DSQRT(3.D0)/2.D0) * (s1(Re) - s2(Re))
        Root(Re,3) =  Root(Re,2)
        Root(Im,3) = -Root(Im,2)
      ELSE IF ((Discrim .LT. 0.D0) .OR. (DABS(r) .GT. 1.D-15)) THEN
        Radius = DSQRT(r*r - Discrim)
        Phi = DACOS(r/Radius)/3.D0
        Radius = Radius**(1.D0/3.D0)
        s1(Re) = Radius * DCOS(Phi)
        s1(Im) = Radius * DSIN(Phi)
        Root(Re,1) = Root(Re,1) + 2.D0*s1(Re)
        Root(Re,2) = Root(Re,2) - s1(Re) - DSQRT(3.D0)*s1(Im)
        Root(Re,3) = Root(Re,3) - s1(Re) + DSQRT(3.D0)*s1(Im)
      END IF
      RETURN
      END




      REAL*8 FUNCTION EC_M_Det2(x)
C     Give determinant of REAL*8 2*2-matrix
      REAL*8 x(2,2)
      EC_M_Det2 = x(1,1)*x(2,2)-x(2,1)*x(1,2)
      END


      REAL*8 FUNCTION EC_M_Determ(x)
C
C     Give determinant of real 3*3-matrix
C
      REAL*8 x(3,3)
      EC_M_Determ = x(1,1)*(x(2,2)*x(3,3)-x(2,3)*x(3,2))
     &        -x(2,1)*(x(1,2)*x(3,3)-x(1,3)*x(3,2))
     &        +x(3,1)*(x(1,2)*x(2,3)-x(1,3)*x(2,2))
      END





C
C ########################################################################
C
C This routine constructs a trend-corrected time-series from a given
C time-series: The model is an additive model. The straight line from
C least squares regression (linear trend) is subtracted from the dataset.
C
C ########################################################################
C





      SUBROUTINE EC_M_Detren(x,NMax,N,MMAx,M,Mean,Cov,y,RC)
C
C Construct a linearly detrended dataset from a given dataset
C
C input : x : REAL*8(NMax,MMax) : x(i,j) = quantity i in sample j
C		only the first N quantities and the first M samples
C		are used.
C	  Cov : REAL*8 (NMax,NMAx) : covariances of quantities.
C		used to find trend.
C output : y : REAL*8(NMAx,MMax) : detrended timeseries.
C	   RC : REAL*8(NMAx) : Directional coefficients of linear
C		regression trend-lines.
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER NMax,N,i,MMax,M,j
      REAL*8 x(NMax,MMax),Mean(NMax),Cov(NMax,NMax),RC(NMax),
     &	y(NMax,MMax),Trend

      DO j = 1,N
	RC(j) = Cov(TTime,j)/Cov(TTime,TTime)
      ENDDO

      DO i= 1,M
	DO j= 1,N
	  IF (j.NE.TTime) THEN
	    Trend = RC(j)*(x(TTime,i)-Mean(TTime))
	    y(j,i) = x(j,i) - Trend
	  ENDIF
	ENDDO
      ENDDO

      RETURN
      END




      SUBROUTINE EC_M_DSwap(x,y)
C     Interchanges x and y
      REAL*8 x,y,dum
      dum = x
      x = y
      y = dum
      RETURN
      END

      SUBROUTINE EC_M_Ell1Q(phi,alpha,ff,ee)
C     Calculates the elliptic integrals F(Phi\Alpha) and
C     E(Phi\Alpha) using the Arithmetic-Geometric Mean process
C     as described in Abramowitz and Stegun, 17.6 (Numbers in
C     text refer to equations in A&S). Only ok for first quadrant
      REAL*8 phi,alpha,ff,ee,a,b,aprev,bprev,edum,eedum
     &   ,c(0:9),psi(0:9),PI,epsilon
      INTEGER It
      PARAMETER (epsilon = 1.D-15)
      PI = 2.D0*DASIN(1.D0)
      It = 0
C     A&S: 17.6.2:
      aprev = 1.D0
C     A&S: 17.6.2:
      bprev = DCOS(alpha)
C     A&S: 17.6.2:
      c(0) = DSIN(alpha)
C     A&S: 17.6.8:
      psi(0) = phi
      edum = c(0)*c(0)
      eedum = 0.D0
      DO WHILE (DABS(c(It)).GT.epsilon)
         It      = It + 1
C        A&S: 17.6.1:
         a       = (aprev + bprev) / 2.D0
C        A&S: 17.6.1:
         b       = DSQRT(aprev * bprev)
C        A&S: 17.6.1:
         c(It)   = (aprev - bprev) / 2.D0
C        A&S: 17.6.8:
         psi(It) = psi(It-1) + DATAN(TAN(psi(It-1))*bprev/aprev)
         psi(It) = psi(It) + PI*DINT((2.D0*psi(It-1)-psi(It))/PI+.5D0)
C        A&S: 17.6.4:
         edum    = edum  + c(It)*c(It)*2.D0**It
C        A&S: 17.6.10:
         eedum   = eedum + c(It)*DSIN(psi(It))
         aprev   = a
         bprev   = b
      END DO
C     A&S: 17.6.9:
      ff = psi(It)/(a*(2.D0**DBLE(It)))
C     A&S: 17.6.10:
      ee = eedum + ff*(1.D0-edum/2.D0)
      RETURN
      END







      SUBROUTINE EC_M_EllCoords(x,b,y,DyDx)
C Calculate the elliptic coordinates (Lambda, Mu, Nu) plus
C derivatives Dy[i]/Dx[j] corresponding to the Carthesian coordinates
C (x[1], x[2], x[3]) for an ellipsoid with semiaxes (b[1], b[2], b[3])
C with b[1]>b[2]>b[3]. Procedure cannot handle points at coordinateplanes.
C Outside the ellipsoid the elliptic coordinates satisfy:
C -b[1]^2 < Nu < -b[2]^2 < Mu < -b[3]^2 < 0 < Lambda.  See e.g. :
C "Einfuehrung in die Kurven- und Flaechentheorie auf vektorieller
C Grundlage" by C.F. Baeschlin, Orell Fuessli Verlag, Zuerich, 1947,
C but mind that there the elliptic coordinates differ from ours.
       INTEGER Re,Im,Lambda,Mu,Nu
       PARAMETER(Re=1,Im=2,Lambda=1,Mu=2,Nu=3)
      REAL*8 x(3), b(3), y(3), DyDx(3,3),Poly(0:3),
     &  Root(Re:Im,3),SQR
      INTEGER i,j
      LOGICAL AllReal
C First check if the semiaxes are in decreasing order and if the point x
C is not on a coordinate plane
      IF (.NOT. ((b(1) .GE. b(2)) .AND. (b(2) .GE. b(3))))
     &  WRITE(*,*) 'Ellipsoid not properly oriented!!!!!!'
      IF (ABS(x(1)*x(2)*x(3)) .LT. 1.D-10*b(3))
     &  WRITE(*,*) 'Sorry : No elliptic coordinates here!!!!'
C Construct the coefficients of the third order equation for the elliptic
C coordinates
      Poly(3) = -1.D0
      Poly(2) = EC_M_SQR(x(1))+EC_M_SQR(x(2))+EC_M_SQR(x(3))
     &         -EC_M_SQR(b(1))-EC_M_SQR(b(2))-EC_M_SQR(b(3))
      Poly(1) = EC_M_SQR(x(1))*(EC_M_SQR(b(2))+EC_M_SQR(b(3)))
     &         -EC_M_SQR(b(1)*b(2))+ EC_M_SQR(x(2))*(EC_M_SQR(b(1))+
     &          EC_M_SQR(b(3)))-EC_M_SQR(b(1)*b(3))+
     &          EC_M_SQR(x(3))*(EC_M_SQR(b(1))+EC_M_SQR(b(2)))-
     &          EC_M_SQR(b(2)*b(3))
      Poly(0) = EC_M_SQR(x(1)*b (2)*b(3)) + EC_M_SQR(x(2)*b(1)*b(3)) +
     &          EC_M_SQR(x(3)*b(1)*b(2)) - EC_M_SQR(b(1)*b(2)*b(3))
C Solve this cubic equation
      CALL EC_M_Cardano(Poly, Root, AllReal)
      IF (.NOT.(AllReal))
     &  WRITE(*,*) 'Error in finding elliptic coordinates!!!'
      DO i=1,3
        y(i) = Root(Re,i)
      END DO
C Put the elliptic coordinates in correct order and check if they satisfy
C the relation given in the header
      DO i=1,2
        DO j=(i+1),3
          IF (y(j) .GT. y(i)) CALL EC_M_DSwap(y(i),y(j))
        END DO
      END DO
      IF (.NOT.((-EC_M_SQR(b(1)) .LT. y(Nu)).AND.(y(Nu) .LT. -EC_M_SQR(b(2)))
     &     .AND.(-EC_M_SQR(b(2)) .LT. y(Mu)).AND.(y(Mu) .LT. -EC_M_SQR(b(3)))
     &     .AND.(0.D0 .LT. y(Lambda))))
     &  WRITE(*,*) 'ERROR!!! in determination of elliptic coordinates'
C Calculate the derivative Dy[i]/Dx[j] of the elliptic coordinates to the
C Carthesian coordinates
      DO i=1,3
        DO j=1,3
          DyDx(j,i) = 2.D0*x(j)*
     &  (EC_M_SQR(b(MOD(j    ,3) + 1)) + y(i)) *
     &  (EC_M_SQR(b(MOD((j+1),3) + 1)) + y(i)) /
     &  ((y(i) - y(MOD(i,3) + 1)) * (y(i) - y(MOD((i+1),3)+1)))
        END DO
      END DO
      END





      SUBROUTINE EC_M_Ellint(phi,alpha,ff,ee)
C     Calculates the elliptic integrals F(Phi\Alpha) and
C     E(Phi\Alpha) using the Arithmetic-Geometric Mean process
C     as described in Abramowitz and Stegun, 17.6 (Numbers in
C     text refer to equations in A&S). ok for all angles.
      REAL*8 phi,alpha,ff,ee,ffcomp,eecomp,PI
      INTEGER SignArg, Wind
      PI = 2.D0*DASIN(1.D0)
      Wind = IDNINT(phi/PI)
      phi = phi - Wind*PI
      IF (Wind .NE. 0) THEN
        CALL EC_M_Ell1Q(0.5D0*PI,alpha,ffcomp,eecomp)
      ELSE
        ffcomp = 0.D0
        eecomp = 0.D0
      END IF
      SignArg = INT(DSIGN(1.D0,phi))
      CALL EC_M_Ell1Q(DABS(phi),alpha,ff,ee)
C     A&S: 17.4.3:
      ff = SignArg*ff + 2.D0*Wind*ffcomp
C     A&S: 17.4.3:
      ee = SignArg*ee + 2.D0*Wind*eecomp
      RETURN
      END






      SUBROUTINE EC_M_InvM(a, aInv)
C
C     Find the inverse of real 3*3 matrix "a"
C
      REAL*8 a(3,3), aInv(3,3), Dum(3,3), Det, EC_M_Determ
      INTEGER i, j
      Det = EC_M_Determ(a)
      IF (DABS(Det) .LT. 1.D-20)
     &    WRITE(*,*) ' Sorry, cannot invert matrix!'
      Dum(1,1) =  (a(2,2)*a(3,3)-a(2,3)*a(3,2))/Det
      Dum(2,1) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))/Det
      Dum(3,1) =  (a(2,1)*a(3,2)-a(2,2)*a(3,1))/Det
      Dum(1,2) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))/Det
      Dum(2,2) =  (a(1,1)*a(3,3)-a(1,3)*a(3,1))/Det
      Dum(3,2) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))/Det
      Dum(1,3) =  (a(1,2)*a(2,3)-a(1,3)*a(2,2))/Det
      Dum(2,3) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))/Det
      Dum(3,3) =  (a(1,1)*a(2,2)-a(1,2)*a(2,1))/Det
      DO i = 1,3
        DO j = 1,3
          aInv(j,i) = Dum(j,i)
        END DO
      END DO
      RETURN
      END

      SUBROUTINE EC_M_InvM2(a, aInv)
C     Find the inverse of REAL*8 2*2 matrix "a"
      REAL*8 a(2,2), aInv(2,2), Dum(2,2), Det, EC_M_Det2
      INTEGER i, j
      Det = EC_M_Det2(a)
      IF (ABS(Det) .LT. 1.D-10)
     &	WRITE(*,*) ' Sorry, cannot invert matrix!'
      Dum(1,1) =  a(2,2)/Det
      Dum(2,1) = -a(2,1)/Det
      Dum(1,2) = -a(1,2)/Det
      Dum(2,2) =  a(1,1)/Det
      DO i = 1,2
	DO j = 1,2
	  aInv(j,i) = Dum(j,i)
	END DO
      END DO
      RETURN
      END

      SUBROUTINE EC_M_ISwap(x,y)
C     Interchanges x and y
      INTEGER x,y,dum
      dum = x
      x = y
      y = dum
      RETURN
      END

      REAL*8 FUNCTION EC_M_ka(x,b)
C     From equation 8
      REAL*8 x,b(3)
      ka = DSQRT((x+b(1)**2)*(x+b(2)**2)*(x+b(3)**2))
      END

      SUBROUTINE EC_M_Map2Vec(a,x,y)
C
C Calculates the image of "x" under the map "a"; y(i) = a(ij)x(j)
C
C input : a : REAL*8 (2,2) : the mapping matrix
C	  x : REAL*8 (2)   : the vector to be mapped
C output : y : REAL*8 (2)  : the image of the map
C
      IMPLICIT NONE

      REAL*8 a(2,2),x(2),y(2),Dum(2)
      INTEGER I,J

      DO I=1,2
	Dum(I) = 0.D0
      ENDDO

      DO I=1,2
	DO J=1,2
	  Dum(I) = Dum(I) + a(I,J)*x(J)
	ENDDO
      ENDDO

      DO I=1,2
	y(I) = Dum(I)
      ENDDO

      RETURN
      END









      SUBROUTINE EC_M_MapMtx(a,x,y)
C
C Calculates the image of "x" under the map "a";
C y(ji) = a(ki)a(lj)x(lk)
C
C input : a : REAL*8 (3,3) : the mapping matrix
C	  x : REAL*8 (3,3) : the tensor to be mapped
C output : y : REAL*8 (3,3) : the image of the map
C
C Revision: June 21, 2001:
C     - indices i and j in the mapping  have
C       been interchanged. This has also been done in the
C       routines that used EC_M_MapMtx (EC_C_T05)
C
      IMPLICIT NONE

      REAL*8 a(3,3),x(3,3),y(3,3),Dum(3,3)
      INTEGER I,J,K,L

      DO I=1,3
	 DO J=1,3
	    Dum(I,J) = 0.D0
	 ENDDO
      ENDDO

      DO I=1,3
	 DO J=1,3
	    DO K=1,3
	       DO L=1,3
		  Dum(I,J) = Dum(I,J) + a(I,K)*a(J,L)*x(K,L)
	       ENDDO
	    ENDDO
	 ENDDO
      ENDDO

      DO I=1,3
	 DO J=1,3
	    y(I,J) = Dum(I,J)
	 ENDDO
      ENDDO

      RETURN
      END






      SUBROUTINE EC_M_MapVec(a,x,y)
C
C Calculates the image of "x" under the map "a"; y(i) = a(ij)x(j)
C
C input : a : REAL*8 (3,3) : the mapping matrix
C	  x : REAL*8 (3)   : the vector to be mapped
C output : y : REAL*8 (3)  : the image of the map
c
C
      IMPLICIT NONE

      REAL*8 a(3,3),x(3),y(3),Dum(3)
      INTEGER I,J

      DO I=1,3
	Dum(I) = 0.D0
      ENDDO

      DO I=1,3
	DO J=1,3
	  Dum(I) = Dum(I) + a(I,J)*x(J)
	ENDDO
      ENDDO

      DO I=1,3
	y(I) = Dum(I)
      ENDDO

      RETURN
      END



C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C
C Mathematics routines
C
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################







      SUBROUTINE EC_M_MinMax(x,NMax,N,MMax, M,Flag,Mins, Maxs)

      INCLUDE 'parcnst.inc'

      REAL*8 x(NMax,MMax),Mins(NMax),Maxs(NMax)
      LOGICAL Flag(NMax,MMax)
      INTEGER I,J, N, M, NSAMP

      DO I=1,N
         Mins(I) = 1e10
	 Maxs(I) = -1e10
      ENDDO

      DO I=1,N
         NSAMP = 0
         DO J=1,M
	    IF (.NOT. FLAG(I,J)) THEN
	      NSAMP = NSAMP + 1
	      Mins(I) = MIN(Mins(I), x(I,J))
	      Maxs(I) = MAX(Maxs(I), x(I,J))
	    ENDIF
	 ENDDO
	 IF (NSAMP .EQ. 0) THEN
            Mins(I) = DUMMY
            Maxs(I) = DUMMY
	 ENDIF
      ENDDO
      END










      SUBROUTINE EC_M_MMul(a,b,c)
C     Matrix C is product of 3*3-matrices A and B
      INTEGER I, J, K
      REAL*8 A(3,3),B(3,3),C(3,3),Dum(3,3)
      DO I=1,3
	DO J=1,3
	  Dum(i,j)=0.D0
	  DO K=1,3
	     Dum(I,J) = Dum(I,J) + A(I,K)*B(K,J)
	  END DO
	END DO
      END DO
      DO I=1,3
	DO J=1,3
	  c(I,J)=Dum(I,J)
	END DO
      END DO
      RETURN
      END





C
C ########################################################################
C
C This is a bundle of simple mathematics routines
C
C ########################################################################
C





      SUBROUTINE EC_M_MulVec(x,y)
      REAL*8 x(3),y
      INTEGER I
      DO I=1,3
         x(I) = x(I)*y
      END DO
      RETURN
      END

      SUBROUTINE EC_M_SortDecr(x,permutation)
C     Sorts the elements of vector x in decreasing order;
C     permutation needed is returned as well
      REAL*8 x(3)
      INTEGER permutation(3),I,J
      permutation(1) = 1
      permutation(2) = 2
      permutation(3) = 3
      DO I = 1,2
         DO J = (I+1),3
            IF (x(J).GT.x(I)) THEN
               CALL EC_M_DSwap(x(I),x(J))
               CALL EC_M_ISwap(permutation(I),permutation(J))
            END IF
         END DO
      END DO
      RETURN
      END

      SUBROUTINE EC_M_SortUse(x,permutation)
C     Reorders the elements of x according to permutation
      REAL*8 x(3),dum(3)
      INTEGER permutation(3),I
      DO I=1,3
         dum(I) = x(I)
      END DO
      DO I=1,3
         x(I) = dum(permutation(I))
      END DO
      RETURN
      END
C
C All following procedures and functions are implementations of the details
C
      SUBROUTINE EC_M_specint(lambda,b,integral)
C
C     Calculates the integrals at the bottom of page seven
C     \INT_(\lambda)^(\infty) \frac{dq}{{b_{i}^{2}+q}k_{q}}
C     References "G&R" in the code are to Gradshteyn and Ryzhik:
C     Tables of Integrals, Series and Products, 4th ed.,Ac. Press,'65
C
      REAL*8 lambda,b(3),integral(3),phi,alpha,ff,ee,dum1
     &   ,dum2,dum3,dum4,dum5,epsilon,PI
      PARAMETER (epsilon = 1.D-18)
      PARAMETER (PI = 0.31415926535897929170D+01)
      IF (DABS(b(1)-b(2)).LT.epsilon) THEN
         dum1 = ((b(1)**2-b(3)**2)**-(1.5D0))*(PI/2.D0
     &   -DATAN(DSQRT((b(3)**2+lambda)/(b(1)**2-b(3)**2))))
         integral(1) = -DSQRT(b(3)**2+lambda)/((b(1)**2-b(3)**2)*
     &      (b(1)**2+lambda))+dum1
         integral(2) = integral(1)
         integral(3) = 2.D0/((b(1)**2+lambda)*DSQRT(b(3)**2+lambda))
     &      -2.D0*integral(1)
      ELSE IF (DABS(b(2)-b(3)).LT.epsilon) THEN
         dum1 = ((b(1)**2-b(3)**2)**-(1.5D0))*
     &      DLOG((DSQRT(b(1)**2+lambda)
     &      -DSQRT(b(1)**2-b(3)**2))**2/(lambda+b(3)**2))
         integral(1) = -2.D0/((b(1)**2-b(3)**2)*DSQRT(b(1)**2+lambda))
     &      -dum1
         integral(2) = DSQRT(b(1)**2+lambda)/((b(1)**2-b(3)**2)*
     &      (b(3)**2+lambda)) + 0.5D0*dum1
         integral(3) = integral(2)
      ELSE
C        alpha and phi : G&R 3.13
         phi   = DASIN(DSQRT((b(1)**2-b(3)**2)/(b(1)**2+lambda)))
         alpha = DASIN(DSQRT((b(1)**2-b(2)**2)/(b(1)**2-b(3)**2)))
         CALL EC_M_Ellint(phi,alpha,ff,ee)
         dum1 = 2.D0/((b(1)**2-b(2)**2)*DSQRT(b(1)**2-b(3)**2))
         dum2 = 2.D0*DSQRT(b(1)**2-b(3)**2)/
     &          ((b(1)**2-b(2)**2)*(b(2)**2-b(3)**2))
         dum3 = 2.D0/(b(2)**2-b(3)**2)*DSQRT((b(3)**2+lambda)/
     &          ((b(1)**2+lambda)*(b(2)**2+lambda)))
         dum4 = 2.D0/((b(3)**2-b(2)**2)*DSQRT(b(1)**2-b(3)**2))
         dum5 = 2.D0/(b(2)**2-b(3)**2)*DSQRT((b(2)**2+lambda)/
     &          ((b(1)**2+lambda)*(b(3)**2+lambda)))
C      G&R 3.133.1
         integral(1) = dum1 * (-ee +        ff)
C      G&R 3.133.7
         integral(2) = dum2 *   ee - dum1 * ff - dum3
C      G&R 3.133.13
         integral(3) = dum4 *   ee             + dum5
      END IF
      RETURN
      END






      REAL*8 FUNCTION EC_M_SQR(x)
C     Give the square of x
      REAL*8 x
      SQR = x*x
      END

      SUBROUTINE EC_M_UnSort(x,permutation)
C     Unsorts the elements of x originally sorted using permutation
      REAL*8 x(3),dum(3)
      INTEGER permutation(3),I
      DO I=1,3
         dum(I) = x(I)
      END DO
      DO I=1,3
         x(permutation(I)) = dum(I)
      END DO
      RETURN
      END
