C  ec_math.f, part of ecpack
C  
C  Copyright (C) 
C    1998-2000   Arjan van Dijk, Wageningen University, The Netherlands
C    2000-2002   Arjan van Dijk, IMAU, The Netherlands
C    1999-2002   Arnold Moene, Wageningen University, The Netherlands
C 
C  This program is free software; you can redistribute it and/or
C  modify it under the terms of the GNU General Public License
C  as published by the Free Software Foundation; either version 2
C  of the License, or (at your option) any later version.
C 
C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.
C 
C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
C 



      SUBROUTINE EC_M_ABCForm(a,b,c,Root1, Root2, AllReal)
C     ****f* ec_math.f/EC_M_ABCForm
C NAME
C     EC_M_ABCForm
C SYNOPSIS
C     CALL EC_M_ABCForm(a,b,c,Root1, Root2, AllReal)
C FUNCTION
C     Solves ax^2 + bx + c = 0
C INPUTS
C     a      : [REAL*8]
C              first coefficient
C     b      : [REAL*8]
C              second coefficient
C     c      : [REAL*8]
C              third coefficient
C OUTPUTS
C     Root1  : [REAL*8]
C              first root
C     Root2  : [REAL*8]
C              second root
C     AllReal: [LOGICAL]
C              all roots real? 
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
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
C     ****f* ec_math.f/EC_M_Averag
C NAME
C     EC_M_Averag
C SYNOPSIS
C     CALL EC_M_Averag(x,NMax,N,MMax,M,Flag,
C    	               Mean,TolMean,Cov,TolCov,MIndep,CIndep,Mok,Cok)
C FUNCTION
C     From the given set of calibrated samples, calculate the averages,
C     variances, covariances and their tolerances.
C INPUTS
C     x      : [REAL*8(NMax,MMax)]
C              Array with calibrated samples.
C	       First index counts quantities; second counter
C	       counts samples. Only the first N quantities
C	       and the first M samples are used.
C     NMax   : [INTEGER]
C              maximum number of quantities
C     N      : [INTEGER]
C              actual number of quantities
C     MMax   : [INTEGER]
C              maximum number of samples
C     M      : [INTEGER]
C              actual number of samples
C     Flag   : [LOGICAL(NMax, MMax)]
C	       If flag(j,i) is true, then quantity j in sample
C              i is not ok.
C OUTPUTS
C     Mean   : [REAL*8(NMax)]
C	       The average value array x. Only samples with Flag = 0 
C              are used.
C     TolMean: tolerance of Mean, defined as 2 * sigma / sqrt(NIndep)
C              where sigma is the standarddeviation of the quantities
C              and where the number of independent samples is estimated
C              as twice the number of sign-changes of the fluctuations
C              of the respective quantities around their means.
C     Cov    : [REAL*8(NMax,NMax)]  
C              covariances
C     TolCov : [REAL*8(NMax,NMax)]
C              tolerances of Cov, estimated in same way as tolerances 
C              of mean.
C     MIndep : [INTEGER(NMAx)] 
C              Number of independent samples in time series from
C              which means are calculated
C     CIndep : [INTEGER(NMAx,NMAx)]  
C              Number of independent samples in time series from 
C              which covariances are calculated
C     Mok    : [INTEGER(NMax)]
C              number of valid samples for each quantity
C     Cok    : [INTEGER(NMax,NMax)]
C              number of valid samples for each combination of 
C              two quantities
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C USES
C     parcnst.inc
C     ***
      IMPLICIT NONE
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
        PrevPrime(j) = 0
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
	MIndep(j) = INT(DBLE(NChange(j,j))/PSwap) - 1
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
	      TolCov(j,k)=TolCov(j,k)+dTolCov(j,k)**2
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
	     TolCov(j,k) = 2.D0*TolCov(j,k)**0.5D0
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
C     ****f* ec_math.f/EC_M_BaseF
C NAME
C     EC_M_BaseF
C SYNOPSIS
C     Value = EC_M_BaseF(x,FType,Order,C)
C FUNCTION
C     Purpose : Calculate simple functions. Currently implemented:
C     - Ordinary polynomials
C     - Polynomials in the natural logarithm of x
C INPUTS
C     x      : [REAL*8]
C              argument of function 
C     FType  : [INTEGER]
C              function type: NormPoly or LogPoly (defined in
C              parcnst.inc
C     Order  : [INTEGER]
C              order of the polynomial
C     C      : [REAL*8(0:Order)]
C              array with coefficients
C RETURN VALUE
C     return value  : [REAL*8]
C              value of polynomial  
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C USES
C     parcnst.inc
C     ***
      IMPLICIT NONE
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
C     ****f* ec_math.f/EC_M_Cardano
C NAME
C     EC_M_Cardano
C SYNOPSIS
C     CALL EC_M_Cardano(Poly, Root, AllReal)
C FUNCTION
C     Uses the Cardano solution to solve exactly:
C     ax^3 + bx^2 + cx + d = 0
C     "a" is not allowed to be zero.
C INPUTS
C     Poly   : [REAL*8(0:3)]
C              coefficient of third order polynomial
C OUTPUTS
C     Root   : [REAL*8(2,3)]
C              array with roots, first index are 
C              real and imaginary part, respectively,
C              second index for three roots
C     AllReal: [LOGICAL]
C              all roots real? 
C AUTHOR
C     Arjan van Dijk
C SEE ALSO
C     Abramowitz and Stegun
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
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
C     ****f* ec_math.f/EC_M_Det2
C NAME
C     EC_M_Det2
C SYNOPSIS
C     Value =  EC_M_Det2(x)
C FUNCTION
C     Give determinant of REAL*8 2*2-matrix
C INPUTS
C     x      : [REAL*8(2,2)]
C              matrix 
C RETURN VALUE
C     return value  : 
C              [REAL*8]
C              determinant 
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
      REAL*8 x(2,2)
      EC_M_Det2 = x(1,1)*x(2,2)-x(2,1)*x(1,2)
      END


      REAL*8 FUNCTION EC_M_Determ(x)
C     ****f* ec_math.f/EC_M_Determ
C NAME
C     EC_M_Determ
C SYNOPSIS
C     Value =  EC_M_Determ(x)
C FUNCTION
C     Give determinant of real 3*3-matrix
C INPUTS
C     x      : [REAL*8(3,3)]
C              matrix 
C RETURN VALUE
C     return value  : 
C              [REAL*8]
C              determinant 
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
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





      SUBROUTINE EC_M_Detren(x,NMax,N,MMAx,M,Mean,Cov,RC)
C     ****f* ec_math.f/EC_M_Detren
C NAME
C     EC_M_Detren
C SYNOPSIS
C     CALL EC_M_Detren(x,NMax,N,MMAx,M,Mean,Cov,y,RC)
C FUNCTION
C     Construct a linearly detrended dataset from a given dataset
C INPUTS
C     x      : [REAL*8(NMax,MMax)] (in/out)
C              array with samples: x(i,j) = quantity i in sample j
C              only the first N quantities and the first M samples
C              are used. On output: detrended series
C     NMax   : [INTEGER]
C              maximum number of quantities in x
C     N      : [INTEGER]
C              actual number of quantities in x
C     MMax   : [INTEGER]
C              maximum number of samples in x
C     M      : [INTEGER]
C              actual number of samples in x
C     Mean   : [REAL*8(NMax)]
C              mean of all quantities
C     Cov    : [REAL*8 (NMax,NMAx)]
C              covariances of quantities used to find trend.
C OUTPUT 
C     RC     : [REAL*8(NMAx)]  
C              Directional coefficients of linear regression
C              trend-lines.
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C Revision 28-01-2003: remove argument y, since it causes aliassing (and is not
C                      needed).
C USES
C     parcnst.inc
C     ***
      IMPLICIT NONE
      INCLUDE 'parcnst.inc'
      INTEGER NMax,N,i,MMax,M,j
      REAL*8 x(NMax,MMax),Mean(NMax),Cov(NMax,NMax),RC(NMax),
     &	Trend

      DO j = 1,N
	RC(j) = Cov(TTime,j)/Cov(TTime,TTime)
      ENDDO

      DO i= 1,M
	DO j= 1,N
	  IF (j.NE.TTime) THEN
	    Trend = RC(j)*(x(TTime,i)-Mean(TTime))
	    x(j,i) = x(j,i) - Trend
	  ENDIF
	ENDDO
      ENDDO

      RETURN
      END




      SUBROUTINE EC_M_DSwap(x,y)
C     ****f* ec_math.f/EC_M_DSwap
C NAME
C     EC_M_DSwap
C SYNOPSIS
C     CALL EC_M_DSwap(x,y)
C FUNCTION
C     Interchanges x and y
C INPUTS
C     x,y    : [REAL*8] 
C              quantities
C OUTPUTS
C     x,y    : [REAL*8] 
C              quantities
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
      REAL*8 x,y,dum
      dum = x
      x = y
      y = dum
      RETURN
      END

      SUBROUTINE EC_M_Ell1Q(phi,alpha,ff,ee)
C     ****f* ec_math.f/EC_M_Ell1Q
C NAME
C     EC_M_Ell1Q
C SYNOPSIS
C     CALL EC_M_Ell1Q(phi,alpha,ff,ee)
C FUNCTION
C     Calculates the elliptic integrals F(Phi\Alpha) and
C     E(Phi\Alpha) using the Arithmetic-Geometric Mean process
C     as described in Abramowitz and Stegun, 17.6 (Numbers in
C     text refer to equations in A&S). Only ok for first quadrant
C INPUTS
C     phi    : [REAL*8] 
C              argument
C     alpha  : [REAL*8] 
C              argument
C OUTPUTS
C     ff     : [REAL*8] 
C              result
C     ee     : [REAL*8] 
C              result
C AUTHOR
C     Arjan van Dijk
C SEE ALSO
C     Abramowitz and Stegun, 17.6 
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
      REAL*8 phi,alpha,ff,ee,a,b,aprev,bprev,edum,eedum
     &   ,c(0:9),psi(0:9),PI,epsilon
      INTEGER It
      PARAMETER (epsilon = 1.D-15)
      PI = 2.D0*DASIN(1.D0)
      It = 0
      a  = 0
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
C     ****f* ec_math.f/EC_M_EllCoords
C NAME
C     EC_M_EllCoords
C SYNOPSIS
C     CALL EC_M_EllCoords(x,b,y,DyDx)
C FUNCTION
C     Calculate the elliptic coordinates (Lambda, Mu, Nu) plus
C     derivatives Dy[i]/Dx[j] corresponding to the Carthesian
C     coordinates (x[1], x[2], x[3]) for an ellipsoid with 
C     semiaxes (b[1], b[2], b[3]) with b[1]>b[2]>b[3]. 
C     Procedure cannot handle points at coordinateplanes.
C     Outside the ellipsoid the elliptic coordinates satisfy:
C     -b[1]^2 < Nu < -b[2]^2 < Mu < -b[3]^2 < 0 < Lambda.  
C INPUTS
C     x    : [REAL*8(3)] 
C            Carthesian coordinate
C     b    : [REAL*8(3)] 
C            semi-axes of ellipsoid
C     alpha  : [REAL*8] 
C              argument
C OUTPUTS
C     y    : [REAL*8(3)] 
C            new coordinate
C     DyDx : [REAL*8(3,3)] 
C            derivative
C HISTORY
C     $Name$
C     $Id$
C AUTHOR
C     Arjan van Dijk
C USES
C     EC_M_Cardano
C     EC_M_DSwap
C     EC_M_SQR
C     ***
      IMPLICIT NONE
      INTEGER Re,Im,Lambda,Mu,Nu
      PARAMETER(Re=1,Im=2,Lambda=1,Mu=2,Nu=3)
      REAL*8 x(3), b(3), y(3), DyDx(3,3),Poly(0:3),
     &  Root(Re:Im,3), EC_M_SQR
      EXTERNAL EC_M_SQR
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
      IF (.NOT. AllReal)
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
      IF (.NOT.((-EC_M_SQR(b(1)) .LT. y(Nu)) .AND.
     &          (y(Nu) .LT. -EC_M_SQR(b(2))) .AND.
     &          (-EC_M_SQR(b(2)) .LT. y(Mu)) .AND.
     &          (y(Mu) .LT. -EC_M_SQR(b(3))) .AND.
     &          (0.D0 .LT. y(Lambda))))
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
C     ****f* ec_math.f/EC_M_Ellint
C NAME
C     EC_M_Ellint
C SYNOPSIS
C     CALL EC_M_Ellint(phi,alpha,ff,ee)
C FUNCTION
C     Calculates the elliptic integrals F(Phi\Alpha) and
C     E(Phi\Alpha) using the Arithmetic-Geometric Mean process
C     as described in Abramowitz and Stegun, 17.6 (Numbers in
C     text refer to equations in A&S). ok for all angles.
C INPUTS
C     phi    : [REAL*8] 
C              argument
C     alpha  : [REAL*8] 
C              argument
C OUTPUTS
C     ff     : [REAL*8] 
C              result
C     ee     : [REAL*8] 
C              result
C AUTHOR
C     Arjan van Dijk
C SEE ALSO
C     Abramowitz and Stegun, 17.6 
C HISTORY
C     $Name$
C     $Id$
C USES
C     EC_M_Ell1Q
C     ***
      IMPLICIT NONE
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
C     ****f* ec_math.f/EC_M_InvM
C NAME
C     EC_M_InvM
C SYNOPSIS
C     CALL EC_M_InvM(a, aInv)
C FUNCTION
C     Find the inverse of real 3*3 matrix "a"
C INPUTS
C     a    : [REAL*8(3,3)] 
C            matrix
C OUTPUTS
C     inv  : [REAL*8(3,3)] 
C            inverse matrix
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C USES
C     EC_M_Determ
C     ***
      IMPLICIT NONE
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
C     ****f* ec_math.f/EC_M_InvM2
C NAME
C     EC_M_InvM2
C SYNOPSIS
C     CALL EC_M_InvM2(a, aInv)
C FUNCTION
C     Find the inverse of real 2*2 matrix "a"
C INPUTS
C     a    : [REAL*8(2,2)] 
C            matrix
C OUTPUTS
C     inv  : [REAL*8(2,2)] 
C            inverse matrix
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C USES
C     EC_M_Det2
C     ***
      IMPLICIT NONE
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
C     ****f* ec_math.f/EC_M_ISwap
C NAME
C     EC_M_ISwap
C SYNOPSIS
C     CALL EC_M_ISwap(x,y)
C FUNCTION
C     Interchanges x and y
C INPUTS
C     x,y    : [INTEGER] 
C              quantities
C OUTPUTS
C     x,y    : [INTEGER] 
C              quantities
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
      INTEGER x,y,dum
      dum = x
      x = y
      y = dum
      RETURN
      END

      REAL*8 FUNCTION EC_M_ka(x,b)
C     ****f* ec_math.f/EC_M_ka
C NAME
C     EC_M_ka
C SYNOPSIS
C     Value =  EC_M_ka(x,b)
C FUNCTION
C     Returns = DSQRT((x+b(1)**2)*(x+b(2)**2)*(x+b(3)**2))
C     From equation 8
C INPUTS
C     x      : [REAL*8] 
C              argument
C     b      : [REAL*8(3)] 
C              coefficient
C RETURN VALUE
C     return value : 
C              [REAL*8] 
C              result
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
      REAL*8 x,b(3)
      EC_M_ka = DSQRT((x+b(1)**2)*(x+b(2)**2)*(x+b(3)**2))
      END


      SUBROUTINE EC_M_Map2Vec(a,x,y)
C     ****f* ec_math.f/EC_M_Map2Vec
C NAME
C     EC_M_Map2Vec
C SYNOPSIS
C     CALL EC_M_Map2Vec(a,x,y)
C FUNCTION
C     Calculates the image of "x" under the map "a"; y(i) = a(ij)x(j)
C INPUTS
C     a      : [REAL*8(2,2)]
C              the mapping matrix
C     x      : [REAL*8(2)] 
C              the vector to be mapped
C OUTPUT 
C     y      : [REAL*8(2)] 
C              the image of the map
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C     ***
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
C     ****f* ec_math.f/EC_M_MapMtx
C NAME
C     EC_M_MapMtx
C SYNOPSIS
C     CALL EC_M_MapMtx(a,x,y)
C FUNCTION
C     Calculates the image of "x" under the map "a";
C    y(ji) = a(ki)a(lj)x(lk)
C INPUTS
C     a      : [REAL*8(3,3)]
C              the mapping matrix
C     x      : [REAL*8(3,3)] 
C              the tensor to be mapped
C OUTPUT 
C     y      : [REAL*8(3,3)] 
C              the image of the map
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     Revision: June 21, 2001:
C     - indices i and j in the mapping  have
C       been interchanged. This has also been done in the
C       routines that used EC_M_MapMtx (EC_C_T05)
C     $Name$
C     $Id$
C     ***
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
C     ****f* ec_math.f/EC_M_MapVec
C NAME
C     EC_M_MapVec
C SYNOPSIS
C     CALL EC_M_MapVec(a,x,y)
C FUNCTION
C     Calculates the image of "x" under the map "a"; y(i) = a(ij)x(j)
C INPUTS
C     a      : [REAL*8(3,3)]
C              the mapping matrix
C     x      : [REAL*8(3)] 
C              the vector to be mapped
C OUTPUT 
C     y      : [REAL*8(3)] 
C              the image of the map
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C     ***
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
C     ****f* ec_math.f/EC_M_MinMax
C NAME
C     EC_M_MinMax
C SYNOPSIS
C     Call EC_M_MinMax(x,NMax,N,MMax, M,Flag,Mins, Maxs)
C FUNCTION
C     Determines minimum and maximum of quantities
C INPUTS
C     x      : [REAL*8(Nmax,MMax)]
C              the sampled quantities
C     NMax   : [INTEGER] 
C              maximum number of quantities
C     N      : [INTEGER] 
C              actual number of quantities
C     MMax   : [INTEGER] 
C              maximum number of samples
C     M      : [INTEGER] 
C              actual number of samples
C     Flag   : [LOGICAL(NMax, MMax)]
C	       If flag(j,i) is true, then quantity j in sample
C              i is not ok.
C OUTPUT 
C     Mins   : [REAL*8(NMax)] 
C              the minimum of each quantity, only taking into
C              account samples with Flag = 0
C     Maxs   : [REAL*8(NMax)] 
C              the maximum of each quantity, only taking into
C              account samples with Flag = 0
C AUTHOR
C     Arnold Moene 
C USES
C     parcnst.inc
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
      include 'parcnst.inc'
      INTEGER I,J, N, M, NSAMP, NMax, MMax
      REAL*8 x(NMax,MMax),Mins(NMax),Maxs(NMax) 
      LOGICAL Flag(NMax,MMax)

      DO I=1,N
         Mins(I) = 1D10
	 Maxs(I) = -1D10
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
C     ****f* ec_math.f/EC_M_MMul
C NAME
C     EC_M_MMul
C SYNOPSIS
C     CALL EC_M_MMul(a,b,c)
C FUNCTION
C     Matrix C is product of 3*3-matrices A and B
C INPUTS
C     a      : [REAL*8(3,3)]
C              first matrix 
C     b      : [REAL*8(3)] 
C              second matrix 
C OUTPUT 
C     c      : [REAL*8(3)] 
C              matrix product 
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
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
C     ****f* ec_math.f/EC_M_MulVec
C NAME
C     EC_M_MulVec
C SYNOPSIS
C     CALL EC_M_MulVec(x,y)
C FUNCTION
C     Multiply vector with a constant
C INPUTS
C     x      : [REAL*8(3)]
C              vector 
C     y      : [REAL*8] 
C              constant
C OUTPUT 
C     x      : [REAL*8(3)] 
C              vector multiplied with y
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
      REAL*8 x(3),y
      INTEGER I
      DO I=1,3
         x(I) = x(I)*y
      END DO
      RETURN
      END

      SUBROUTINE EC_M_SortDecr(x,permutation)
C     ****f* ec_math.f/EC_M_SortDecr
C NAME
C     EC_M_SortDecr
C SYNOPSIS
C     CALL EC_M_SortDecr(x,permutation)
C FUNCTION
C     Sorts the elements of vector x in decreasing order;
C     permutation needed is returned as well
C INPUTS
C     x      : [REAL*8(3)]
C              vector 
C OUTPUT 
C     x      : [REAL*8(3)]
C              sorted vector 
C     permutation :
C              [INTEGER(3)]
C              permutation to arrive at sorted vector
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C USES
C     EC_M_DSwap
C     EC_M_ISwap
C     ***
      IMPLICIT NONE
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
C     ****f* ec_math.f/EC_M_SortUse
C NAME
C     EC_M_SortUse
C SYNOPSIS
C     CALL EC_M_SortUse(x,permutation)
C FUNCTION
C     Reorders the elements of x according to permutation
C INPUTS
C     x      : [REAL*8(3)]
C              vector 
C     permutation :
C              [INTEGER(3)]
C              permutation to arrive at sorted vector
C OUTPUT 
C     x      : [REAL*8(3)]
C              sorted vector 
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
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
C     ****f* ec_math.f/EC_M_specint
C NAME
C     EC_M_specint
C SYNOPSIS
C     CALL EC_M_specint(lambda,b,integral)
C FUNCTION
C     Calculates the integrals:
C     \INT_(\lambda)^(\infty) \frac{dq}{{b_{i}^{2}+q}k_{q}}
C     References "G&R" in the code are to Gradshteyn and Ryzhik:
C     Tables of Integrals, Series and Products, 4th ed.,Ac. Press,'65
C INPUTS
C     lambda   : [REAL*8]
C                one lambda
C     b        : [REAL*8(3)]
C                vector of b-values
C OUTPUT 
C     integral : [REAL*8(3)]
C                result
C AUTHOR
C     Arjan van Dijk
C SEE ALSO 
C     References "G&R" in the code are to Gradshteyn and Ryzhik:
C     Tables of Integrals, Series and Products, 4th ed.,Ac. Press,'65
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
      REAL*8 lambda,b(3),integral(3),phi,alpha,ff,ee,dum1
     &   ,dum2,dum3,dum4,dum5,epsilon,PI
      PARAMETER (epsilon = 1.D-18)
      PARAMETER (PI = 0.31415926535897929170D+01)
      IF (DABS(b(1)-b(2)).LT.epsilon) THEN
         dum1 = ((b(1)**2-b(3)**2)**(-1.5D0))*(PI/2.D0
     &   -DATAN(DSQRT((b(3)**2+lambda)/(b(1)**2-b(3)**2))))
         integral(1) = -DSQRT(b(3)**2+lambda)/((b(1)**2-b(3)**2)*
     &      (b(1)**2+lambda))+dum1
         integral(2) = integral(1)
         integral(3) = 2.D0/((b(1)**2+lambda)*DSQRT(b(3)**2+lambda))
     &      -2.D0*integral(1)
      ELSE IF (DABS(b(2)-b(3)).LT.epsilon) THEN
         dum1 = ((b(1)**2-b(3)**2)**(-1.5D0))*
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
C     ****f* ec_math.f/EC_M_SQR
C NAME
C     EC_M_SQR
C SYNOPSIS
C     Value =  EC_M_SQR(x)
C FUNCTION
C     Give the square of x
C INPUTS
C     x        : [REAL*8]
C RETURN VALUE 
C     return value : 
C                [REAL*8]
C                result
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
      REAL*8 x
      EC_M_SQR = x*x
      END

      SUBROUTINE EC_M_UnSort(x,permutation)
C     ****f* ec_math.f/EC_M_UnSort
C NAME
C     EC_M_UnSort
C SYNOPSIS
C     CALL EC_M_UnSort(x, permutation)
C FUNCTION
C     Unsorts the elements of x originally sorted using permutation
C INPUTS
C     x      : [REAL*8(3)]
C              vector 
C     permutation :
C              [INTEGER(3)]
C              permutation to arrive at sorted vector
C OUTPUT 
C     x      : [REAL*8(3)]
C              unsorted vector 
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE
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
