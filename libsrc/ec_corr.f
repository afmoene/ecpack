C  ec_corr.f, part of ecpack
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


      

      SUBROUTINE EC_C_D01(x,b,aInv)
C     ****f* ec_corr.f/EC_C_D01
C NAME
C     EC_C_D01
C SYNOPSIS
C     CALL EC_C_D01(x, b, aInv)
C FUNCTION
C     This routines computes and applies the matrix for the
C     correction of turbulent air flow measurements for the
C     presence of small disturbing objects, like a box with
C     electronic apparatus. The approach followed is described in:
C     (see SEE ALSO). The actual matrix is computed 
C     in EC_C_D02.
C INPUTS
C      x : [REAL*8(3)]
C          Position vector of the point where measurements have been taken.
C          The ellipsoid is placed in the origin.
C          A right-handed frame of coordinates is chosen.
C          The flow is supposed to be expressed in this coordinate frame.
C          Therefore, when the flow velocity has positive components,
C          upstream measurement points are selected when by giving
C          vector x negative components!
C      b : [REAL*8(3)]
C          Three ellipsoid semi-axes in meters. The ellipsoid is
C          supposed to be oriented along the coordinate axes.
C          Somehow the algorithm does not seem to like it when two or more
C          semi-axes are equal, or when your point x is in one of the
C          coordinate planes (one component of x equal to zero).
C          To circumvent problems one can take values slightly off
C          the problematic values.
C OUTPUT
C     aInv : [REAL*8(3,3)]
C          The matrix which can be used to correct samples and
C          covariances for flow distortion.
C          Sum_j aInv(i,j)*u(j)     gives the distortion-corrected
C          image of measured velocity u.
C AUTHOR
C     Arjan van Dijk
C
C SEE ALSO
C     Oost, W. (1991).  Flow distortion by an ellipsoid and its 
C     application to the analysis of atmospheric measurements. 
C     J.Atm.Oc.Tech., 8 No 3:331-340. 
C     References in the code are to this article.
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     EC_C_D02
C     EC_M_InvM
C     ***
      IMPLICIT NONE
      REAL*8 x(3),b(3),a(3,3),aInv(3,3)
      CALL EC_C_D02(x,b,a)
      CALL EC_M_InvM(a,aInv)
      RETURN
      END


      SUBROUTINE EC_C_D02(x,b,a)
C     ****f* ec_corr.f/EC_C_D02
C NAME
C     EC_C_D02
C SYNOPSIS
C     CALL EC_C_D02(x, b, a)
C FUNCTION
C     This routines computes the distortion matrix for the
C     correction of turbulent air flow measurements for the
C     presence of small disturbing objects, like a box with
C     electronic apparatus. The approach followed is described in:
C     (see SEE ALSO). 
C     Calculates the distortion matrix according to equation 12.
C INPUTS 
C     b : [REAL*8(3)]
C         a vector containing the three semiaxes of the
C         ellipsoid
C     x : [REAL*8(3)]
C         position where the distortion-matrix will be
C         calculated; coordinates are relative to the
C         center of the ellipsoid
C OUTPUT 
C     a : [REAL*8(3,3)]
C         the distortion matrix, when applied to an
C         undisturbed wind, "a" gives the disturbed wind
C AUTHOR
C     Arjan van Dijk
C SEE ALSO
C     Oost, W. (1991).  Flow distortion by an ellipsoid and its 
C     application to the analysis of atmospheric measurements. 
C     J.Atm.Oc.Tech., 8 No 3:331-340. 
C     References in the code are to this article.
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     EC_M_SortDecr
C     EC_M_SortUse
C     EC_M_EllCoords
C     EC_M_specint
C     EC_M_MulVec
C     EC_M_UnSort
C     ***
      IMPLICIT NONE
      INTEGER Lambda
      PARAMETER(Lambda=1)
      REAL*8 x(3),b(3),a(3,3)
     &   ,alphaint(3),integral(3),EC_M_ka,y(3),DyDx(3,3)
      INTEGER i,j,perm(3)
      CALL EC_M_SortDecr(b,perm)
      CALL EC_M_SortUse(x,perm)
      CALL EC_M_EllCoords(x,b,y,DyDx)
      CALL EC_M_specint(0.D0,b,alphaint)
      CALL EC_M_MulVec(alphaint,b(1)*b(2)*b(3))
      CALL EC_M_specint(y(Lambda),b,integral)
      DO i=1,3
         DO j=1,3
            a(perm(i),perm(j)) = -x(j) * b(1)*b(2)*b(3) *
     &      DyDx(i,Lambda)/
     &      ((2.D0-alphaint(j))*(b(j)**2+y(lambda))*
     &       EC_M_ka(y(lambda),b))
            IF (i.EQ.j) a(perm(i),perm(j)) = a(perm(i),perm(j))
     &      + 1.D0 + b(1)*b(2)*b(3)*integral(i)/(2.D0-alphaint(i))
         END DO
      END DO
      CALL EC_M_UnSort(x,perm)
      CALL EC_M_UnSort(b,perm)
      RETURN
      END





      SUBROUTINE EC_C_F01(Mean,Cov,NMax,NSize,WhichTemp,
     &	NSta,NEnd,NumInt,NS,TauD,TauV,CalSonic,CalTherm,CalHyg,
     &  CalCO2,WXT)
C     ****f* ec_corr.f/EC_C_F01
C NAME
C     EC_C_F01
C SYNOPSIS
C     CALL EC_C_F01(Mean, Cov, NMax, NSize, WhichTemp,
C                   NSta, NEnd, NumInt, NS, TauD, TauV, CalSonic,
C                   CalTherm, CalHyg, CalCO2,WXT)
C FUNCTION
C     Calculate frequency response corrections for sensor response, path
C     length averaging, sensor separation, signal processing and dampening
C     of fluctuations in a tube. Based on publications in SEE ALSO.
C
C INPUTS
C     Mean   : [REAL*8(NMax)]  
C              Array of mean values of the quantities in
C              this experiment (only the first N quantities are used).
C     Cov    : [REAL*8(NMax,NMax)]  
C               covariances of the fluctuations.
C     NMax   : [INTEGER]  
C              Physical dimension of array Mean
C     NSize  : [INTEGER]  
C              Number of quantities actually involved in this
C              experiment.
C     WhichTemp: [INTEGER]
C              Which temperature to use: thermocouple (Tcouple)
C              or sonic (TSOnic): see parcnst.inc
C     NSta   : [REAL*8]  
C              Start frequency numerical integration.
C              Popular value: -5.D0 [unit?].
C     NEND   : [REAL*8]  
C              End frequency numerical integration.
C              Popular value: LOG(5) = 0.69897D0 [unit?].
C     NumINT : [INTEGER]  
C              Number of intervals in integration.
C              Popular value: 19.
C     NS     : [INTEGER]
C     TauD   : [REAL*8]  
C              Interval length for running mean.
C              Popular value: 0.D0 [unit?].
C     TAUV   : [REAL*8]  
C              Low pass filter time constant.
C              Popular value: 0.D0 [unit?]
C     CalSonic : [REAL*8(NQQ)]  
C              Calibration specification array of
C              sonic anemometer.
C     CalTherm : [REAL*8(NQQ)]  
C              Calibration specification array of
C              thermometer.
C     CalHyg : [REAL*8(NQQ)]  
C              Calibration specification array of
C              hygrometer.
C     CalCO2 : [REAL*8(NQQ)]  
C              Calibration specification array of
C              CO2 sensor.
C
C OUTPUT
C     WXT    : [REAL*8(NMax,NMax)]  
C              Correction factors for covariances.
C BUGS
C     In the present configuration, only correction factors are computed
C     for the variances and for the covariances involving the vertical
C     velocity. Other covariances get a factor of 1.
C AUTHOR
C     Arjan van Dijk
C
C SEE ALSO
C	Moore, C.J. (1986): 'Frequency Response Corrections for Eddy
C	Correlation Systems'. Boundary Layer Met. 37: 17-35.
C
C	Philip, J.R. (1963): 'The Damping of Fluctuating Concentration
C	by Continuous Sampling Through a tube' Aust. J. Phys. 16: 454-463.
C
C	Leuning, R. and K.M. King (1991): 'Comparison of Eddy-Covariance
C	Measurements of CO2 Fluxes by open- and closed-path CO2 analysers'
C	(unpublished)
C HISTORY
C       28-05-2001: added info on which temperature should be used
C                   in corrections (Sonic or thermocouple)
C       19-09-2002: added CalCO2 in interface
C       27-01-2003: replaced Nint by NumInt
C     $Name$ 
C     $Id$
C USES
C     Pi
C     EC_Ph_Obukhov
C     ***
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER I1,I2,I3,I4,I,J,NumINT,NMax,NSize, WhichTemp
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ),CALCO2(NQQ),
     &	GAIND,GAINV,GAINT1,GAINT2,GAINUV,GAINW,GAINQ1, GAINCO21,
     &	TRANA,TRANQ1,TRANUV,TRANW,TRANT1,TRANT2,TRANWT1,TRANWT2,TRANWQ1,
     &	TRANWUV,ATT,AWW,AUU,AWT,AUW,BTT,BWW,BUU,BWT,BUW,TAUD,TAUV,TAUT,
     &  TRANCO21, TRANWCO21,
     &	XI,ZETA,CW,CU,WT,UW,TT,UU,WW,F,NS,C,NSTA,NEND,N,LF,X,ZL,
     &	INTTS(NNMax,NNMAX),G(NNMAX,NNMAX),
     &  WXT(NMax,NMax),Mean(NMax),Cov(NMax,NMax),Ustar, Tstar, Qstar

      REAL*8 EC_Ph_Obukhov
      EXTERNAL EC_Ph_Obukhov

C
C-- Declarations of constants ----------------------------------------
C
      TAUT = CalTherm(QQTime)/
     &	(1.D0+4.9D0*(CalTherm(QQTime)**0.5D0*Mean(U))**0.45D0)
      LF   = (NEND-NSTA)/DBLE(NumINT-1)	    ! Interval width num. integr.
      C    = 1.D0/(LOG(10.D0)*LF/3.D0)	    ! Constant of integration, eq.1
      LF   = 10.D0**LF
      N    = 10.D0**NSTA
      DO I=1,NNMAX
         DO J=1,NNMAX
            INTTS(I,J)=0.D0
         ENDDO
      ENDDO

C Make scales for Obukhov length calculation 
      Ustar = (Cov(U,W)**2.D0 + Cov(V,W)**2.D0)**0.25D0
      Tstar = -Cov(W, WhichTemp)/Ustar
C Take care of situation when we don't have a hygrometer
      IF (INT(CalHyg(QQType)) .NE. INT(DUMMY)) THEN
         Qstar = -Cov(W, SpecHum)/Ustar
      ELSE 
         Qstar = 0
      ENDIF
      ZL = CalSonic(QQZ)/EC_Ph_Obukhov(Ustar, Tstar, Qstar,
     &                                 Mean(WhichTemp))

      IF(ZL.GE.0.D0) THEN
C
C-- Stable spectral and co-spectral function factors, eq. 20b, 21b
C
        ATT=0.0961D0+0.644D0*ZL**0.6D0
        AWW=0.838D0+1.172D0*ZL
        AUU=0.2D0*AWW
        AWT=0.284D0*(1.D0+6.4D0*ZL)**0.75D0
        AUW=0.124D0*(1.D0+7.9D0*ZL)**0.75D0
        BTT=3.124D0/ATT**0.667D0
        BWW=3.124D0/AWW**0.667D0
        BUU=3.124D0/AUU**0.667D0
        BWT=2.34D0/AWT**1.1D0
        BUW=2.34D0/AUW**1.1D0
      ELSE
C
C-- Unstable spectral and co-spectral function factors, eq. 23 -----
C
        XI   = (-ZL)**0.667D0
        ZETA = (0.001D0*CalSonic(QQZ))**1.667D0
        CW   = 0.7285D0+1.4115D0*XI
        CU   = 9.546D0+1.235D0*XI/(ZETA**0.4D0)
      END IF
C
C-- Start Simpson's rule numerical integration -----------------------
C     From -5 to log(5) in NumInt steps
C
      I3=0
      DO I1=1,10
        DO I2=2,4,2
          I3=I3+1
          IF (I3.GT.NumInt) GOTO 200
          I4=I2
          IF (I3.EQ.1)	I4=I4-1
          IF (I3.EQ.NumInt) I4=I4-1
C
C-- Skip integration if frequency exceeds Nyquist frequency ------
C
          IF (N.GT.(0.5D0*NS)) GOTO 999

C
C-- ELECTRONIC / DIGITAL FILTERS !!
C-- Auto Regressive Moving Average filter response gain, eq. 19 --
C
          GAIND=1.0D0
          X=2.0D0*PI*N*TAUD
          IF (X.LT.6.D0 .AND. X.GT.0.0D0) GAIND=X*X/(1.0D0+X*X)
C
C-- Butterworth filter frequency response gain, eq. 5 ------------
C
          GAINV = 1.0D0
          GAINV=(2.0D0*PI*N*TAUV)**2.0D0
          GAINV=1.0D0+GAINV*GAINV
C
C-- Data aquisition co-spectral transfer function, eq. 15 (b=3) --
C
          TRANA=1.0D0+(N/(NS-N))**3.0D0
C
C-- LUW THERMOCOUPLE TEMPERATURE !!
C-- Thermocouple frequency response gain, eq. 2 ------------------
C
          GAINT1=1.0D0+(2.0D0*PI*N*TAUT)**2.0D0
C
C-- Thermocouple spatial averaging transfer function, eq. 7 ------
C
          TRANT1=1.0D0	       !spatial averaging negligible
C
C-- W-T lateral separation transfer function, eq. 11 -------------
C
          TRANWT1=1.0D0
          X=N/Mean(U)*(CalTherm(QQX)-CalSonic(QQX))
          IF (X.GT.0.01D0) TRANWT1=EXP(-9.9D0*X**1.5D0)
C
C-- SOLENT SONIC TEMPERATURE !!
C-- Sonic temperature frequency response gain, eq. 2 -------------
C
          GAINT2=1.0D0	       !sonic temperature gain negligible
C
C-- Sonic temperature spatial averaging transfer function, eq.7 --
C
          TRANT2=1.0D0
          X=2.0D0*PI*N/Mean(U)*CalSonic(QQExt4)
          IF (X.GT.0.02D0) 
     &	       TRANT2=(3.0D0+EXP(-X)-4.0D0*(1.0D0-EXP(-X))/X)/X
C
C-- W-T lateral separation transfer function, eq. 11 -------------
C
          TRANWT2=1.0D0	       !there is no lateral separation

C
C-- SOLENT SONIC HORIZONTAL WINDSPEED !!
C-- UV-sensor frequency response gain, eq. 2 ---------------------
C
          GAINUV=1.0D0	       !sonic UV windspeed gain negligible
C
C-- UV-sensor spatial averaging transfer function, eq. 7 ---------
C
          TRANUV=1.0D0
          X=2.0D0*PI*N/Mean(U)*CalSonic(QQPath)
          IF (X.GT.0.02D0) 
     &          TRANUV=(3.0D0+EXP(-X)-4.0D0*(1.0D0-EXP(-X))/X)/X
C
C-- W-UV lateral separation transfer function, eq. 11 ------------
C
          TRANWUV=1.0D0
          X = N/Mean(U)*CalSonic(QQExt3)
          IF (X .GT.0.01D0) TRANWUV = EXP(-9.9D0*X**1.5D0)

C
C-- SOLENT SONIC VERTICAL WINDSPEED !!
C-- W-sensor frequency response gain, eq. 2 ----------------------
C
          GAINW=1.0D0	       !sonic W windspeed gain neglectible
C
C-- W-sensor spatial averaging transfer function, eq. 9 ----------
C
          TRANW=1.0D0
          X=2.0D0*PI*N/Mean(U)*CalSonic(QQPath)
          IF (X.GT.0.04D0)
     &	    TRANW=(1.0D0+(EXP(-X)-3.D0*(1.0D0-EXP(-X))/X)/2.0D0)*4.0D0/X

C
C-- LYMANN-ALPHA OPEN PATH HYGROMETER !!
C-- hygrometer frequency response gain, eq. 2 ------------
C
          GAINQ1=1.0D0	       !hygrometer gain neglectible
C
C-- lymann-alpha spatial averaging transfer function, eq.7 ------------
C
          TRANQ1=1.0D0
          X=2.0D0*PI*N/Mean(U)*CalHyg(QQPath)
          IF (X.GT.0.02D0) 
     &           TRANQ1=(3.0D0+EXP(-X)-4.0D0*(1.0D0-EXP(-X))/X)/X
C
C-- W-Q lateral separation transfer function, eq. 11 -------------
C
          TRANWQ1=1.0D0
          X=N/Mean(U)*CalHyg(QQX)
          IF (X.GT.0.01D0) TRANWQ1=EXP(-9.9D0*X**1.5D0)
C
C-- CO2 sensor !!
C-- CO2 sensor frequency response gain, eq. 2 ------------
C
          GAINCO21=1.0D0	       !CO2 sensor gain neglectible
C
C-- lymann-alpha spatial averaging transfer function, eq.7 ------------
C
          TRANCO21=1.0D0
          X=2.D0*PI*N/Mean(U)*CalCO2(QQPath)
          IF (X.GT.0.02D0) THEN
	     TRANCO21=(3.0D0+EXP(-X)-4.0D0*(1.0D0-EXP(-X))/X)/X
          ENDIF
C
C-- W-Q lateral separation transfer function, eq. 11 -------------
C
          TRANWCO21=1.0D0
          X=N/Mean(U)*CalCO2(QQX)
          IF (X.GT.0.01D0) TRANWCO21=EXP(-9.9D0*X**1.5D0)
C
C-- Composite transfer functions, eq. 28 -------------------------
C
          DO I=1,NNMAX
	     DO J=1,NNMAX
	        G(I,J) = 0.0D0
             ENDDO
	  ENDDO

          G(U,U)= TRANUV / GAINUV    !UU
          G(V,V)= TRANUV / GAINUV    !VV
          G(W,W)= TRANW  / GAINW    !WW
          G(TCouple, Tcouple) = TRANT1 / GAINT1    !TTth
          G(TSonic,TSonic)= TRANT2 / GAINT2    !TTso
          G(SpecHum,SpecHum)= TRANQ1 / GAINQ1    !QQ
          G(Humidity,Humidity)= TRANQ1 / GAINQ1    !rhov rhov
          G(CO2,CO2)= TRANCO21 / GAINCO21    !CO2
          G(SpecCO2,SpecCO2)= TRANCO21 / GAINCO21    !CO2
          G(U,W) = TRANWUV * SQRT (G(U,U) * G(W,W))    !WU
          G(W,U) = G(U,W)
          G(V,W) = TRANWUV * SQRT (G(V,V) * G(W,W))    !WV
          G(W,V) = G(V,W)
          G(W,TCouple) = TRANWT1 * 
     &                   SQRT (G(W,W) * G(TCouple,TCouple))	    !WTth
          G(TCouple,W) = G(W,TCouple)
          G(W,TSonic) = TRANWT2 * SQRT (G(W,W) * G(TSonic,TSonic))    !WTso
          G(TSonic,W) = G(W,TSonic)
          G(W,Humidity) = TRANWQ1 * 
     &                   SQRT (G(W,W) * G(Humidity,Humidity))	    !WQkr
          G(Humidity,W) = G(W,Humidity)
          G(W,SpecHum) = TRANWQ1 * 
     &                   SQRT (G(W,W) * G(Humidity,Humidity))	    !WQkr
          G(SpecHum,W) = G(W,SpecHum)
          G(W,CO2) = TRANWCO21 * 
     &                   SQRT (G(W,W) * G(CO2,CO2))	    !WCO2
          G(CO2,W) = G(W,CO2)
          G(W,SpecCO2) = TRANWCO21 * 
     &                   SQRT (G(W,W) * G(CO2,CO2))	    ! WqCO2
          G(SpecCO2,W) = G(W,SpecCO2)
          DO I=1,NNMax
	     DO J=1,NNMax
	        G(I,J) = G(I,J)*TRANA/GAINV*GAIND
             ENDDO
          ENDDO

          F=N/Mean(U)*CalSonic(QQZ)
          IF (ZL.LT.0.0D0) THEN
C
C-- Unstable normalised spectral and co-spectral forms ---------
C
            UU=210.0D0*F/(1.0D0+33.0D0*F)**1.667D0		       !eq. 23
            UU=UU+F*XI/(ZETA+2.2D0*F**1.667D0)
            UU=UU/CU
            WW=16.0D0*F*XI/(1.0D0+17.0D0*F)**1.667D0		       !eq. 22
            WW=WW+F/(1.0D0+5.3D0*F**1.667D0)
            WW=WW/CW
            TT=14.94D0*F/(1.0D0+24.0D0*F)**1.667D0			     !eq. 24
            IF (F.GE.0.15D0) TT=6.827D0*F/(1.0D0+12.5D0*F)**1.667D0
            WT=12.92D0*F/(1.0D0+26.7D0*F)**1.375D0		       !eq. 25
            IF (F.GE.0.54D0) WT=4.378D0*F/(1.0D0+3.8D0*F)**2.4D0
            UW=20.78D0*F/(1.0D0+31.0D0*F)**1.575D0		       !eq. 26
            IF (F.GE.0.24D0) UW=12.66D0*F/(1.0D0+9.6D0*F)**2.4D0
          ELSE
C
C-- Stable normalised spectral and co-spectral forms, eq. 20a --
C
            UU=F/(AUU+BUU*F**1.667D0)
            WW=F/(AWW+BWW*F**1.667D0)
            TT=F/(ATT+BTT*F**1.667D0)
            WT=F/(AWT+BWT*F**2.1D0)
            UW=F/(AUW+BUW*F**2.1D0)
          ENDIF
C
C-- Integral of co-spectral transfer function * atm. co-spectrum -
C
          INTTS(U,U)  = INTTS(U,U) + I4*G(U,U) *UU  !UU
          INTTS(V,V)  = INTTS(V,V) + I4*G(V,V) *UU  !VV
          INTTS(W,W)  = INTTS(W,W) + I4*G(W,W) *WW  !WW
          INTTS(TCouple,TCouple)  = 
     &           INTTS(TCouple,TCouple) + I4*G(TCouple,TCouple)*TT !TTth
          INTTS(TSonic,TSonic)  = 
     &           INTTS(TSonic,TSonic) + I4*G(TSonic,TSonic) *TT ! TTso
          INTTS(Humidity,Humidity)  =
     &           INTTS(Humidity,Humidity) + I4*G(Humidity,Humidity) *TT  !QQkr
          INTTS(CO2,CO2)  =
     &           INTTS(CO2,CO2) + I4*G(CO2,CO2) *TT  !CO2 CO2
          INTTS(U,W)  = INTTS(U,W) + I4*G(U,W) *UW   !WU
          INTTS(W,U) = INTTS(U,W)
          INTTS(V,W)  = INTTS(V,W) + I4*G(V,W)*UW   !WV
          INTTS(W,V) = INTTS(V,W)
          INTTS(W,TCouple) = INTTS(W,TCouple) + I4*G(W,TCouple)*WT   !WTth
          INTTS(TCouple,W) = INTTS(W,TCouple)
          INTTS(W,TSonic) = INTTS(W,TSonic) + I4*G(W,TSonic)*WT   !WTso
          INTTS(TSonic,W) = INTTS(W,TSonic)
          INTTS(W,Humidity) = INTTS(W,Humidity) + I4*G(W,Humidity)*WT  !WQkr
          INTTS(Humidity,W) = INTTS(W,Humidity)
          INTTS(W,CO2) = INTTS(W,CO2) + I4*G(W,CO2)*WT  !WQCO2
          INTTS(CO2,W) = INTTS(W,CO2)

C
C-- Increase frequency by interval width -------------------------
C
          N=N*LF
        ENDDO
      ENDDO

 200  CONTINUE

C
C-- Final flux loss correction factors !! ----------------------------
C

      DO I=1,NSize
        DO J=1,NSize
           WXT(I,J) = 1.0D0
        ENDDO
      ENDDO

999   CONTINUE
      WXT(U,U) = C/INTTS(U,U)		   !UU
      WXT(V,V) = C/INTTS(V,V)		   !VV
      WXT(U,V) = WXT(U,U)		   !UV
      WXT(V,U) = WXT(U,V)		   !VU
      WXT(W,W) = C/INTTS(W,W)		   !WW
      WXT(W,U) = C/INTTS(U,W)		   !WU
      WXT(U,W) = WXT(W,U)
      WXT(W,V) = C/INTTS(V,W)		   !WV
      WXT(V,W) = WXT(W,V)

      WXT(W,TSonic) = C/INTTS(W,TSonic)	 !WTso
      WXT(TSonic,W) = WXT(W,TSonic)
      WXT(TSonic,TSonic) = C/INTTS(TSonic,TSonic)  !TTso

      WXT(W,TCouple) = C/INTTS(W,TCouple)	  !WTth
      WXT(TCouple,W) = WXT(W,TCouple)
      WXT(TCouple,TCouple) = C/INTTS(TCouple,TCouple)    !TTth

      WXT(W,Humidity) = C/INTTS(W,Humidity)	  !WQly
      WXT(Humidity,W) = WXT(W,Humidity)
      WXT(Humidity,Humidity) = C/INTTS(Humidity,Humidity)  !QQly

      WXT(W,SpecHum)	   = WXT(W,Humidity)
      WXT(SpecHum,W)	   = WXT(Humidity,W)
      WXT(SpecHum,SpecHum) = WXT(Humidity,Humidity)

      WXT(W,CO2) = C/INTTS(W,CO2)	  ! W rhoCO2
      WXT(CO2,W) = WXT(W,CO2)
      WXT(CO2,CO2) = C/INTTS(CO2,CO2)  ! rhoCO2^2

      WXT(W,SpecCO2)	   = WXT(W,CO2)
      WXT(SpecCO2,W)	   = WXT(CO2,W)
      WXT(SpecCO2,SpecCO2) = WXT(CO2,CO2)
      RETURN
      END





      SUBROUTINE EC_C_F02(WXT,NMax,NSize,Lower,Upper,Cov,TolCov)
C     ****f* ec_corr.f/EC_C_F02
C NAME
C     EC_C_F02
C SYNOPSIS
C     CALL EC_C_F02(WXT, NMax, NSize, Lower, Upper, Cov, TolCov)
C FUNCTION
C     Apply frequency response corrections for sensor response, path
C     length averaging, sensor separation, signal processing and dampening
C     of fluctuations in a tube. Based on publications given in SEE ALSO
C INPUTS
C     WXT    : [REAL*8(NMax,NMax)]  
C              Correction factors for covariances.
C     NMax   : [INTEGER]  
C              Physical dimension of array Mean
C     NSize  : [INTEGER]  
C              Number of quantities actually involved in this
C              experiment.
C     Lower  : [REAL*8]
C              Lower acceptance limit for frequency-response
C              factors. Correction factors smaller than Lower are 
C              set to 1.
C     Upper  : [REAL*8]
C              Upper acceptance limit for frequency-response
C              factors. Correction factors larger than Upper are 
C              Which temperature to use: thermocouple (Tcouple)
C              or sonic (TSOnic): see parcnst.inc
C OUTPUT
C     Cov    : [REAL*8(NMax,NMax)]  
C              covariances of the fluctuations.
C     TolCov : [REAL*8(NMax,NMax)]
C              tolerances in covariances 
C
C SEE ALSO
C     Moore, C.J. (1986): 'Frequency Response Corrections for Eddy
C     Correlation Systems'. Boundary Layer Met. 37: 17-35.
C
C     Philip, J.R. (1963): 'The Damping of Fluctuating Concentration
C     by Continuous Sampling Through a tube' Aust. J. Phys. 16: 454-463.
C
C     Leuning, R. and K.M. King (1991): 'Comparison of Eddy-Covariance
C     Measurements of CO2 Fluxes by open- and closed-path CO2 analysers'
C     (unpublished)
C
C HISTORY 
C     28-05-2001: added info on which temperature should be used
C                 in corrections (Sonic or thermocouple)
C     $Name$ 
C     $Id$
C     ***

      IMPLICIT NONE
      INTEGER I,J,NMax,NSize
      REAL*8 WXT(NMax,NMax),Cov(NMax,NMax),Lower,Upper,TolCov(NMax,NMax)

      DO i=1,NSize
        DO j=1,NSize
          IF ((WXT(i,j).GE.Lower).AND.(WXT(i,j).LE.Upper)) THEN
            Cov(i,j) = Cov(i,j)*WXT(i,j)
C
C The tolerance is augmented by half the amount added by frequency-
C response correction.
C
            TolCov(i,j) = SQRT(TolCov(i,j)**2.D0+
     &	      (0.5D0*Cov(i,j)*(WXT(i,j)-1.D0))**2.D0)
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END



C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C
C Correction routines
C
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################






C
C General correction routine that calls all individual corrections
C
      SUBROUTINE EC_C_Main(OutF,
     &	DoPrint,
     &	Mean,NMax,N,TolMean,
     &	Cov,TolCov,
     &  DoCorr, PCorr, ExpVar, 
     &	DirYaw, DirPitch, DirRoll,
     &	SonFactr,
     &	O2Factor,
     &	CalSonic,CalTherm,CalHyg,
     &  CalCo2, FrCor,
     &	WebVel,P, Have_Cal)
C     ****f* ec_corr.f/EC_C_Main
C NAME
C     EC_C_MAIN
C SYNOPSIS
C     CALL EC_C_Main(OutF,
C     	DoPrint, Mean,NMax,N,TolMean, Cov,TolCov,
C     	DoCorr, PCorr, ExpVar,
C     	DirYaw, DirPitch, DirRoll,
C     	SonFactr, O2Factor,
C     	CalSonic,CalTherm,CalHyg,
C       CalCo2,FrCor,
C     	P, Have_Cal)
C FUNCTION
C     Integrated correction routine applying ALL (user-selected)
C     corrections in this library on mean values and covariances.
C     All intermediate results can be output to a file.
C     Moreover they are returned to the calling routine in
C     respective variables.
C INPUTS
C     (outputs and combined inputs/outputs are given as well,
C     but also under OUTPUT)
C     OutF    : [INTEGER]
C               Unit number of file for intermediate results
C     DoPrint : [LOGICAL]
C               Print intermediate results ?
C     Mean    : [REAL*8(NMax)] (in/out)
C               Mean values of all calibrated signals
C     NMax    : [INTEGER]
C               Size of arrays with calibrated signals
C     N       : [INTEGER]
C               Actual number of calibrated signals
C     TolMean : [REAL*8(NMax)] (in/out)
C               Tolerances in mean values of all calibrated signals
C     Cov     : [REAL*8(NMax,NMax)] (in/out)
C               Covariances of all calibrated signals
C     TolCov  : [REAL*8(NMax,NMax)] (in/out)
C               Tolerances in covariances of all calibrated signals
C     DoCorr  : [LOGICAL](NMaxCorr)
C               which corrections to do?
C     PCorr   : [LOGICAL](NMaxCorr)
C               intermediate results of which corrections?
C     ExpVar  : [REAL*8](NMaxExp)
C               array with experimental settings
C     DirYaw  : [REAL*8] (out)
C               Yaw angle (degrees)
C     DirPitch: [REAL*8] (out)
C               Pitch angle (degrees)
C     DirRoll : [REAL*8] (out)
C               Roll angle (degrees)
C     SonFactr: [REAL*8(NMax)] (out)
C               Correction factor due to Schotanus correction for
C               covariance of sonic temperature with each calibrated
C               signal.
C     O2Factor: [REAL*8(NMax)] (out)
C               Correction factor due to oxygen correction for
C               covariance of humidity with each calibrated
C     CalSonic: [REAL*8(NQQ)]
C               Calibration info for sonic anemometer
C     CalTherm: [REAL*8(NQQ)]
C               Calibration info for thermocouple
C     CalHyg  : [REAL*8(NQQ)]
C               Calibration info for hygrometer
C     CalCO2  : [REAL*8(NQQ)]
C               Calibration info for CO2 sensor
C     FrCor   : [REAL*8(NMax,NMax)] (out)
C               Correction factors for covariances for frequency
C               response
C     P       : [REAL*8]
C               Atmospheric pressure (Pa)
C     Have_Cal : [LOGICAL(NMax)]
C               Calibrated signal available for given quantity ?
C OUTPUT
C     Mean    : [REAL*8(NMax)] (in/out)
C               Mean values of all calibrated signals
C     TolMean : [REAL*8(NMax)] (in/out)
C               Tolerances in mean values of all calibrated signals
C     Cov     : [REAL*8(NMax,NMax)] (in/out)
C               Covariances of all calibrated signals
C     TolCov  : [REAL*8(NMax,NMax)] (in/out)
C               Tolerances in covariances of all calibrated signals
C     DirYaw  : [REAL*8] (out)
C               Yaw angle (degrees)
C     DirPitch: [REAL*8] (out)
C               Pitch angle (degrees)
C     DirRoll : [REAL*8] (out)
C               Roll angle (degrees)
C     SonFactr: [REAL*8(NMax)] (out)
C               Correction factor due to Schotanus correction for
C               covariance of sonic temperature with each calibrated
C               signal.
C     O2Factor: [REAL*8(NMax)] (out)
C               Correction factor due to oxygen correction for
C               covariance of humidity with each calibrated
C     FrCor   : [REAL*8(NMax,NMax)] (out)
C               Correction factors for covariances for frequency
C               response
C     WebVel  : [REAL*8] (out)
C               Webb velocity
C HISTORY
C     28-05-2001: added info on whether uncalibrated data are
C                 available for a given variable (mainly important
C                 for sonic and/or Couple temperature since that
C                 is used for various corrections)  
C     07-10-2002: added WebVel to interface to pass
C                 Webb velocity independent from Mean(W)
C     27-01-2003: removed BadTC from interface
C                 replaced list of correction switches by DoCorr and PCorr, ExpVar
C                 replaced Nint by NumInt
C     $Name$ 
C     $Id$
C USES
C     parcnst.inc
C     EC_C_T05
C     EC_C_T06
C     EC_C_T07
C     EC_C_T08
C     EC_C_T09
C     EC_C_T10
C     EC_C_T11
C     EC_C_Schot1
C     EC_C_Schot2
C     EC_C_Oxygen1
C     EC_C_Oxygen2
C     EC_C_F01
C     EC_C_F02
C     EC_C_Webb
C     EC_G_Show
C     EC_G_ShwFrq
C     EC_G_Reset
C     ***
      IMPLICIT NONE
      INCLUDE 'parcnst.inc'

      INTEGER N,NumInt,NMAx,OutF, WhichTemp, I
      LOGICAL DoCorr(NMaxCorr), PCorr(NMaxCorr),
     &	DoPrint,
     &	QYaw,QPitch,QRoll,QFreq,QO2,QWebb,QSonic,QTilt,
     &  Have_Cal(NMax), QSchot
      REAL*8 P,Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),FrCor(NMax,NMax),
     &	DirYaw,O2Factor(NMax),
     &	NSTA,NEND,TAUV,TauD,DirPitch,DirRoll, TSonFact,
     &	SonFactr(NMax),
     &  Yaw(3,3), Roll(3,3), Pitch(3,3), WebVel,
     &  ExpVar(NMaxExp)
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ), CalCO2(NQQ)
C
C
C Only print intermediate results if desired
C
C
      QTilt  = (PCorr(QCTilt)  .AND. DoPrint)
      QYaw   = (PCorr(QCYaw)   .AND. DoPrint)
      QPitch = (PCorr(QCPitch) .AND. DoPrint)
      QRoll  = (PCorr(QCRoll)  .AND. DoPrint)
      QSonic = (PCorr(QCSonic) .AND. DoPrint)
      QO2    = (PCorr(QCO2)    .AND. DoPrint)
      QFreq  = (PCorr(QCFreq)  .AND. DoPrint)
      QWebb  = (PCorr(QCWebb)  .AND. DoPrint)
C A Hack: Just set Qschot to QSonic (name confusion, apparently)
      QSchot = QSonic
C
C
C Perform a tilt-correction of the whole system using KNOWN (!!!!!!) angles.
C Used to compensate for KNOWN miss-alignment!!!!!
C The classic "tilt-corrections", which turn Mean(V), Mean(W) and Cov(W,V)
C to zero, to are carried out later!!!!!!
C If no fixed tilt-angles are known, this subroutine may be skipped.
C
C
      IF (DoCorr(QCTilt)) THEN
        CALL EC_C_T06(ExpVar(QEPreYaw),Yaw)
        CALL EC_C_T05(Mean,NMax,N,Cov,Yaw)
        CALL EC_C_T08(ExpVar(QEPrePitch),Pitch)
        CALL EC_C_T05(Mean,NMax,N,Cov,Pitch)
        CALL EC_C_T10(ExpVar(QEPreRoll),Roll)
        CALL EC_C_T05(Mean,NMax,N,Cov,Roll)
        CALL EC_G_Reset(Have_Cal, Mean, TolMean, Cov, TolCov)

        IF (QTilt) THEN
          CALL EC_G_ShwHead(OutF, 
     &	  'Averages of data after fixed tilt-correction: ')
          CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,NMax,N)
        ENDIF
      ENDIF
C
C
C Perform classic yaw-correction : Mean(V) --> 0
C
C
      IF (DoCorr(QCYaw)) THEN
        CALL EC_C_T07(Mean(U),Mean(V),DirYaw)
        CALL EC_C_T06(DirYaw,Yaw)
        CALL EC_C_T05(Mean,NMax,N,Cov,Yaw)

        CALL EC_G_Reset(Have_Cal, Mean, TolMean, Cov, TolCov)

        IF (QYaw) THEN
          WRITE(OutF,*)
          WRITE(OutF,*) 'Yaw angle = ',DirYaw,' degrees'
          CALL EC_G_ShwHead(OutF,
     &	           'After yaw-correction (Mean(V) -> 0) : ')
          CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,NMax,N)
        ENDIF
      ENDIF
C
C
C Perform classic pitch-correction : Mean(W) --> 0
C
C
      IF (DoCorr(QCPitch)) THEN
        CALL EC_C_T09(Mean(U),Mean(W),DirPitch)
        IF (ABS(DirPitch).LE. ExpVar(QEPitchLim)) THEN
          CALL EC_C_T08(DirPitch,Pitch)
          CALL EC_C_T05(Mean,NMax,N,Cov,Pitch)
        ENDIF
        CALL EC_G_Reset(Have_Cal, Mean, TolMean, Cov, TolCov)

        IF (QPitch) THEN
          WRITE(OutF,*)
          WRITE(OutF,*) 'Pitch angle = ',DirPitch,' degrees'
          CALL EC_G_ShwHead(OutF,
     &         'After pitch-correction (Mean(W) -> 0) : ')
          CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,NMax,N)
        ENDIF
      ENDIF
C
C
C Perform classic roll-correction : Cov(W,V) --> 0
C
C
      IF (DoCorr(QCRoll)) THEN
        CALL EC_C_T11(Cov(V,V),Cov(V,W),Cov(W,W),DirRoll)
        IF (ABS(DirRoll) .LE. ExpVar(QERollLim)) THEN
          CALL EC_C_T10(DirRoll,Roll)
          CALL EC_C_T05(Mean,NMax,N,Cov,Roll)
        ENDIF
        CALL EC_G_Reset(Have_Cal, Mean, TolMean, Cov, TolCov)

        IF (QRoll) THEN
          WRITE(OutF,*)
          WRITE(OutF,*) 'Roll angle = ',DirRoll,' degrees'
          CALL EC_G_ShwHead(OutF,
     &          'After roll-correction (Cov(V,W) -> 0) : ')
          CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,NMax,N)
        ENDIF
      ENDIF
C
C Reset the coveriances for which no data available
C
      CALL EC_G_Reset(Have_Cal, Mean, TolMean, Cov, TolCov)
C
C Correct sonic temperature and all covariances with sonic temperature
C for humidity. This is half the Schotanus-correction. Side-wind
C correction is done at calibration time directly after reading raw data.
C
C Now we need to know if and which temperature we have. Default to
C Sonic temperature
C
      IF (HAVE_CAL(TSonic)) THEN
         WhichTemp = Tsonic
      ELSE IF (Have_CAL(TCouple)) THEN
         WhichTemp = TCouple
      ELSE
         WhichTemp = -1
      ENDIF

      IF (DoCorr(QCSonic)) THEN
        CALL EC_C_Schot1(Mean(SpecHum),Mean(TSonic),NMax,N,Cov,
     &    SonFactr,TSonFact)
        CALL EC_C_Schot2(SonFactr,TSonFact,Mean(TSonic),NMax,N,Cov)
        CALL EC_G_Reset(Have_Cal, Mean, TolMean, Cov, TolCov)

        IF (QSchot) THEN
          WRITE(OutF,*)
            CALL EC_G_ShwHead(OutF,
     &          'Correction factors for H2O-sensitivity of sonic')
          DO i=1,N
	    WRITE(OutF,45) QName(i),SonFactr(i)
          ENDDO
   45	FORMAT('For (',a6,',T(son)) : ',F10.5)
          CALL EC_G_ShwHead(OutF,
     &                      'After H2O-correction for sonic : ')
          CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,NMax,N)
	ENDIF
      ENDIF
C
C
C Perform correction for oxygen-sensitivity of hygrometer
C
C
      IF (DoCorr(QCO2)) THEN
        IF (WhichTemp .GT. 0) THEN
          CALL EC_C_Oxygen1(Mean(WhichTemp),NMax,N,Cov,P,CalHyg(QQType),
     &      WhichTemp,O2Factor)
          CALL EC_C_Oxygen2(O2Factor,NMax,N,Cov)
          CALL EC_G_Reset(Have_Cal, Mean, TolMean, Cov, TolCov)

          IF (QO2) THEN
            CALL EC_G_ShwHead(OutF,
     &        'Correction factors for O2-sensitivity of hygrometer')
            DO i=1,N
	      WRITE(OutF,46) QName(i),O2Factor(i)
            ENDDO
 46         FORMAT('For (',a6,',RhoV or q) : ',F10.5)
            CALL EC_G_ShwHead(OutF,
     &           'After oxygen-correction for hygrometer : ')
            CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,NMax,N)
          ENDIF
        ELSE
            WRITE(*,*) 'ERROR: can not perform O2 correction without ',
     &                 'a temperature'
        ENDIF
      ENDIF
C
C
C Perform correction for poor frequency response and large paths
C
C

C
C Constants for integration routine in correction frequency response
C
      NSTA   = -5.0D0	 ! [?] start frequency numerical integration
      NEND   = 0.69897D0 ! [?] end frequency numerical integration (LOG(5))
      NumINT = 19	 ! [1] number of intervals
      TAUV   = 0.0D0	 ! [?] Low pass filter time constant
      TauD   = 0.0D0	 ! [?] interval length for running mean

      IF (DoCorr(QCFreq)) THEN
        IF (WhichTemp .GT. 0) THEN
          CALL EC_C_F01(Mean,Cov,NMax,N,WhichTemp,
     &         NSta,NEnd,NumInt,ExpVar(QEFreq),TauD,TauV,
     &         CalSonic,CalTherm,CalHyg,
     &         CalCO2, FrCor)
          CALL EC_C_F02(FrCor,NMax,N,ExpVar(QELLimit),
     &                  ExpVar(QEULimit),Cov,TolCov)
     
          CALL EC_G_Reset(Have_Cal, Mean, TolMean, Cov, TolCov)
          IF (QFreq) THEN
            CALL EC_G_ShwFrq(OutF,FrCor,NMax,N)
            WRITE(OutF,*) 'Factors accepted between ',
     &                     ExpVar(QELLimit),
     &	                  ' and ',ExpVar(QEULimit)
            CALL EC_G_ShwHead(OutF, 'After frequency-correction : ')
            CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,NMax,N)
          ENDIF
        ELSE
            WRITE(*,*) 'ERROR: can not perform freq. response ',
     &                 ' correction without a temperature'
        ENDIF
      ENDIF
C
C
C Calculate mean vertical velocity according to Webb
C
C
      WebVel = DUMMY
      IF (DoCorr(QCWebb)) THEN
         IF (WhichTemp .GT. 0) THEN
           CALL EC_C_Webb(Mean,NMax,Cov,P, WhichTemp, WebVel)
           CALL EC_G_Reset(Have_Cal, Mean, TolMean, Cov, TolCov)
           IF (QWebb) THEN
	     WRITE(OutF,*) 'Webb-velocity (vertical) = ',WebVel,' m/s'
             CALL EC_G_ShwHead(OutF, 'After addition of Webb-term: ')
	     CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,
     &              NMax,N)
           ENDIF
         ELSE
           WRITE(*,*) 'ERROR: can not perform Webb correction',
     &                 ' without a temperature'
         ENDIF
      ENDIF

      RETURN
      END



      SUBROUTINE EC_C_Oxygen1(MeanT,NMax,N,Cov,P,HygType,WhichTemp,
     &  Factor)
C     ****f* ec_corr.f/EC_C_Oxygen1
C NAME
C     EC_C_Oxygen1
C SYNOPSIS
C     CALL EC_C_Oxygen1(MeanT, NMax, N, Cov, P, 
C          HygType, WhichTemp, Factor)
C FUNCTION
C     Contribution of other gases especially oxygen absorb
C     some of the radiation at the wavelengths at which the
C     hygrometer works.
C     This routine computes the correction factors. The 
C     factors are applied in EC_C_Oxygen2
C INPUTS
C     MeanT     : [REAL*8]
C                 Mean temperature (Kelvin)
C     NMax      : [INTEGER]
C                 Size of array dimensions
C     N         : [INTEGER]
C                 Actual number of calibrated signals
C     Cov       : [REAL*8(Nmax,Nmax)]
C                 Covariances
C     P         : [REAL*8]
C                 Atmospheric pressure (Pa)
C     HygType   : [INTEGER]
C                 Type of hygrometer (see parcnst.inc for codes)
C     WhichTemp : [INTEGER]
C                 Use thermocouple temperature or sonic temperature ?
C                 (Tcouple or TSonic, codes in parcnst.inc)
C OUTPUT
C     Factor    : [REAL*8(NMax)]
C                 Correction factor for the covariance with each of
C                 calibrated signals
C NOTES
C     Currently this routine knows about three types of hygrometer:
C        ApCampKrypton : the KH20 krypton hygrometer from Campbell Sci.
C        ApMierijLyma  : the Lymann-alpha hygrometer from Mierij Meteo
C        ApLiCor7500   : the LiCor IR  hygrometer (not sensitive to O2)
C
C SEE ALSO
C     EC_C_Oxygen2
C HISTORY
C     28-05-2001: added info on which temperature should be used
C                 in corrections (Sonic or thermocouple)
C     $Name$ 
C     $Id$
C USES
C     parcnst.inc
C     Kok
C     KwK
C     Kola
C     Kwla
C     FracO2
C     MO2
C     RGas
C     ***
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER NMax,N,i, WhichTemp
      REAL*8 Cov(NMax,NMax),Factor(NMax),OXC,P, HygType,
     +       GenKo, GenKw, MeanT

      IF (HygType .EQ. ApCampKrypton) THEN
          GenKo = KoK
          GenKw = KwK
      ELSE IF (HygType .EQ. ApMierijLyma) THEN
          GenKo = KoLa
          GenKw = KwLa
      ELSE IF (HygType .EQ. ApLiCor7500) THEN
          GenKo = 0.0D0
	  GenKw = 1.0D0
      ENDIF
      OXC = FracO2*MO2*P*GenKo/(RGas*MeanT**2.D0*GenKw)

      DO i = 1,N
	Factor(i) = 1.D0 + OXC*Cov(i,WhichTemp)/Cov(i,Humidity)
	IF ((i.EQ.Humidity).OR.(i.EQ.SpecHum))
     &	  Factor(i) = Factor(i)**2.D0
      ENDDO

      RETURN
      END





      SUBROUTINE EC_C_Oxygen2(Factor,NMax,N,Cov)
C     ****f* ec_corr.f/EC_C_Oxygen2
C NAME
C     EC_C_Oxygen2
C SYNOPSIS
C     CALL EC_C_Oxygen2(Factor, Nmax, N, Cov)
C FUNCTION
C     Contribution of other gases especially oxygen absorb
C     some of the radiation at the wavelengths at which the
C     hygrometer works.
C     This routine applies the correction factors that were 
C     computed in EC_C_Oxygen1
C INPUTS
C     Factor    : [REAL*8(NMax)]
C                 Correction factor for the covariance with each of
C                 calibrated signals
C     NMax      : [INTEGER]
C                 Size of array dimensions
C     N         : [INTEGER]
C                 Actual number of calibrated signals
C     Cov       : [REAL*8(Nmax,Nmax)] (in/out)
C                 Covariances
C OUTPUT
C     Cov       : [REAL*8(Nmax,Nmax)] (in/out)
C                 Covariances
C SEE ALSO
C     EC_C_Oxygen1
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     parcnst.inc
C     ***
      IMPLICIT NONE
      INCLUDE 'parcnst.inc'

      INTEGER NMax,N,i
      REAL*8 Cov(NMax,NMax),Factor(NMax)
      DO i = 1,N
	Cov(i,Humidity) = Factor(i)*Cov(i,Humidity)
	Cov(i,SpecHum ) = Factor(i)*Cov(i,SpecHum )
	Cov(Humidity,i) = Cov(i,Humidity)
	Cov(SpecHum ,i) = Cov(i,SpecHum )
      ENDDO

      RETURN
      END




      SUBROUTINE EC_C_Scal(Cal, UDum, VDum, WDum,
     &                  UError, VError, WError)
C     ****f* ec_corr.f/EC_C_Scal
C NAME
C     EC_C_Scal
C SYNOPSIS
C     CALL EC_C_Scal(Cal, UDum, VDum, WDum, 
C                    UError, VError, WError)
C FUNCTION
C     To calibrate a sonic signal according to wind tunnel
C     calibration (for the moment this works for apparatus 6,
C     i.e. a wind tunnel calibrated sonic)
C INPUTS
C     Cal   : [REAL*8(NQQ)]
C             Array of length NQQ with calibration info
C     UDum  : [REAL*8 ] (in/out)
C             One horizontal component (on exit: calibrated)
C     VDum  : [REAL*8 ] (in/out)
C             Another horizontal component (on exit: calibrated)
C     WDum  : [REAL*8 ] (in/out)
C             Vertical component (on exit: calibrated)
C OUTPUT
C     UDum  : [REAL*8 ] (in/out)
C             One horizontal component (on exit: calibrated)
C     VDum  : [REAL*8 ] (in/out)
C             Another horizontal component (on exit: calibrated)
C     WDum  : [REAL*8 ] (in/out)
C             Vertical component (on exit: calibrated)
C     UError: [LOGICAL]  
C             error flag for U (.TRUE. if wrong data)
C     VError: [LOGICAL]  
C             error flag for V (.TRUE. if wrong data)
C     WError: [LOGICAL]  
C             error flag for W (.TRUE. if wrong data)
C AUTHOR
C     Arnold Moene
C CREATION DATE
C     September 26, 2000
C NOTES
C     The method is based on a wind tunnel calibration of the sonic
C     The real velocity components can be derived from the
C     measured components and the real azimuth and elevation angle.
C     But the latter are not known and have to be determined
C     iteratively from the measured components. The relationship
C     between the real components and the measured components is:
C
C     Ureal =  Umeas/(UC1*(1 - 0.5*
C                   ((Azi + (Elev/0.5236)*UC2)*
C                    (1 - UC3*Abs(Elev/0.5236)))**2 ))
C     Vreal =  Vmeas*(1 - VC1*Abs(Elev/0.5236))
C     Wreal =  Wmeas/(WC1*(1 - 0.5*(Azi*WC2)**2))
C     and
C     Azi = arctan(V/U)
C     Elev = arctan(W/sqrt(U**2 + V**2))
C
C     where UC1, UC2, UC3, VC1, WC1, WC2 are fitting coefficients.
C     An azimuth angle of zero is supposed to refer to a wind
C     direction from the most optimal direction (i.e. the 'open'
C     side of a sonic). Samples with an absolute azimuth angle of
C     more than 40 degrees are rejected.
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     parcnst.inc
C     ***

      IMPLICIT NONE
      INCLUDE 'parcnst.inc'

      REAL*8   Cal(NQQ), UDum, VDum, WDum
      LOGICAL  UError, VError, WError

      REAL*8   UCorr, VCorr, WCorr, AziNew, AziOld, ElevNew, ElevOld,
     &         UC1, UC2, UC3, VC1, WC1, WC2
      INTEGER  NITER, ITMAX


      UC1 = Cal(QQExt6)
      UC2 = Cal(QQExt7)
      UC3 = Cal(QQExt8)
      VC1 = Cal(QQExt9)
      WC1 = Cal(QQExt10)
      WC2 = Cal(QQExt11)

      ITMAX = 20
      NITER = 0
      AziOld = 9999D0
      ElevOld = 9999D0
      AziNew = ATAN(VDum/UDum)
      ElevNew = ATAN(WDum/SQRT(UDum**2 + VDum**2))

      UCorr = UDum
      VCorr = VDum
      WCorr = WDum

      IF (ABS(AziNew) .LE. 0.698) THEN
         DO WHILE (((ABS(AziNew-AziOld) .GT. 1./60) .OR.
     &              (ABS(ElevNew-ElevOld) .GT. 1./60)) .AND.
     &           (NITER .LT. ITMAX))
            UCorr =  UDum/(UC1*(1.D0 - 0.5D0*
     &              ((AziNew + (ElevNew/0.5236D0)*UC2)*
     &              (1 - UC3*Abs(ElevNew/0.5236D0))
     &              )**2        )
     &                     )
            VCorr =  VDum*(1 - VC1*Abs(ElevNew/0.5236D0))
            WCorr =  WDum/(WC1*(1.D0 - 0.5D0*(ElevNew*WC2)**2))

            AziOld = AziNew
            ElevOld = ElevNew

            AziNew = ATAN(Vcorr/UCorr)
            ElevNew = ATAN(Wcorr/SQRT(UCorr**2 + VCorr**2))

            NITER = NITER + 1
         ENDDO
      ENDIF

      IF ((NITER .EQ. ITMAX) .OR. (ABS(AziNew) .GT. 0.698)) THEN
         UError = .TRUE.
         VError = .TRUE.
         WError = .TRUE.
      ELSE
         UDum = UCorr
         VDum = VCorr
         WDum = WCorr
      ENDIF

      END





      SUBROUTINE EC_C_Schot1(MeanQ,MeanTSon,NMax,N,Cov,Factor,TSonFact)
C     ****f* ec_corr.f/EC_C_Schot1
C NAME
C     EC_C_Schot1
C SYNOPSIS
C     CALL EC_C_Schot1(MeanQ,MeanTSon,NMax,N,Cov,Factor,TSonFact)
C FUNCTION
C     Compute correction factor for partial Schotanus et al. correction.
C     for humidity of sonic temperature, and of all covariances with 
C     sonic temperature.
C     Sidewind-correction has already been applied in the
C     routine where the sonic signal is calibrated.
C INPUTS
C     MeanQ    : [REAL*8]
C                mean specific humidity  (kg/kg)
C     MeanTSon : [REAL*8]
C                mean sonic temperature (Kelvin)
C     NMax     : [INTEGER]
C                size of arrays
C     N        : [INTEGER]
C                actual number of calibrated signals
C     Cov      : [REAL*8(NMax,NMax)]
C                covariance matrix of calibrated signals
C OUTPUT
C     Factor   : [REAL*8(NMax)]
C                correction factor for the covariances with specific
C                humidity
C     TSonFact : [REAL*8]
C                correction factor for sonic temperature
C SEE ALSO
C     EC_C_Schot1
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     parcnst.inc
C     ***
      IMPLICIT NONE
      INCLUDE 'parcnst.inc'

      INTEGER NMax,N,i
      REAL*8 MeanQ,MeanTSon,Cov(NMax,NMax),Factor(NMax),TSonFact

      DO i = 1,N
	Factor(i) = 1.D0 - 0.51D0*MeanQ
     &	  -0.51D0*MeanTSon*Cov(i,SpecHum)/Cov(i,TSonic)
	IF (i.EQ.TSonic) Factor(i) = Factor(i)**2.D0
      ENDDO

      TSonFact = 1.D0/(1.D0+0.51D0*MeanQ)

      RETURN
      END




      SUBROUTINE EC_C_Schot2(Factor,TSonFact,MeanTSon,NMax,N,Cov)
C     ****f* ec_corr.f/EC_C_Schot2
C NAME
C     EC_C_Schot2
C SYNOPSIS
C     CALL EC_C_Schot2(Factor,TSonFact,MeanTSoc, NMax,N,Cov)
C FUNCTION
C     Apply correction factor for partial Schotanus et al. correction
C     as computed in EC_C_Schot1.
C     for humidity of sonic temperature, and of all covariances with 
C     sonic temperature.
C     Sidewind-correction has already been applied in the
C     routine where the sonic signal is calibrated.
C INPUTS
C     Factor   : [REAL*8(NMax)]
C                correction factor for the covariances with specific
C                humidity
C     TSonFact : [REAL*8]
C                correction factor for sonic temperature
C     MeanTSon : [REAL*8]
C                mean sonic temperature (Kelvin)
C     NMax     : [INTEGER]
C                size of arrays
C     N        : [INTEGER]
C                actual number of calibrated signals
C     Cov      : [REAL*8(NMax,NMax)] (in/out)
C                covariance matrix of calibrated signals
C OUTPUT
C     Cov      : [REAL*8(NMax,NMax)] (in/out)
C                covariance matrix of calibrated signals
C SEE ALSO
C     EC_C_Schot1
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     physcnst.inc
C     ***
      IMPLICIT NONE
      INCLUDE 'parcnst.inc'

      INTEGER NMax,N,i
      REAL*8 Cov(NMax,NMax),Factor(NMax), MeanTSon, TSonFact

      DO i = 1,N
	Cov(i,TSonic) = Factor(i)*Cov(i,TSonic)
	Cov(TSonic,i) = Cov(i,TSonic)
      ENDDO

      MeanTSon = MeanTSon*TSonFact

      RETURN
      END




      REAL*8 FUNCTION EC_C_Schot3(Temp, Rhov, Press)
C     ****f* ec_corr.f/EC_C_Schot3
C NAME
C     EC_C_Schot3
C SYNOPSIS
C     TCORR = EC_C_Schot3(Temp,Rhov, Press)
C FUNCTION
C     To do humidity part of Schotanus et al. correction on a
C     raw sample of sonic temperature, while specific humidity is not
C     yet known: to get specific humidity from the measured 
C     Rhov one needs a temperature: if no thermocouple available
C     the only temperature is the sonic temperature.
C     (TSonic depends on specific humidity, to compute specific 
C     humidity, one needs a temperature, Tsonics depends on 
C     specific humidity ....etc.)
C     Sidewind-correction has already been applied in the
C     routine where the sonic signal is calibrated.
C INPUTS
C     TEMP    : [REAL*8]
C               Sonic temperature (without hyumidity correction) (Kelvin)
C     Rhov    : [REAL*8]
C               absolute humidity (kg/m^3)
C     Press   : [REAL*8]
C               atmospheric pressure (Pa)
C RETURN VALUE
C     return value : [REAL*8]
C               corrected sonic temperature (Kelvin)
C AUTHOR
C     Arnold Moene
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     EC_Ph_Q
C     ***
      IMPLICIT NONE

      REAL*8   Temp, Rhov, Press, SPHUM, SPHUMNEW, NEWTEMP
      REAL*8 EC_Ph_Q ! External function calls

      SPHUM = EC_Ph_Q(Temp, Rhov, Press)
      SPHUMNEW = 1.0D0
      NEWTEMP = Temp
C Does this converge ?
      DO WHILE (ABS((SPHUM - SPHUMNEW)/SPHUM) .GT. 0.0001)
         NEWTEMP = TEMP/(1.D0+0.51D0*SPHUM)
         SPHUMNEW = SPHUM
         SPHUM = EC_Ph_Q(NEWTEMP, Rhov, Press)
      ENDDO
      EC_C_Schot3 = NEWTEMP
      RETURN
      END



C
C ########################################################################
C
C Planar fit method for tiltcorrection
C
C ########################################################################
C



      SUBROUTINE EC_C_T01(DoWBias,SingleRun,uMean,NRuns,Apf,
     &  Alpha,Beta,Gamma,WBias)
C     ****f* ec_corr.f/EC_C_T01
C NAME
C     EC_C_T01
C SYNOPSIS
C     CALL EC_C_T01(DoWBias,SingleRun,uMean,NRuns,Apf, 
C                   Alpha,Beta,Gamma,WBias)
C FUNCTION
C     Subroutine performs some preparations for the actual
C     Planar Fit Method. The actual work is done in routine EC_C_T02.
C
C INPUTS
C     DoWBias:   [LOGICAL]
C                compute bias in mean vertical wind (FALSE 
C                implies that the mean vertical wind over all
C                runs is assumed to be zero
C     SingleRun: [LOGICAL]
C                Determine rotation for a single run 
C     uMean    : [REAL*8(3,Nmax)]
C                matrix of run mean velocity vectors
C     NRuns    : [INTEGER]
C                the number of runs
C OUTPUT
C     Apf      : [REAL*8(3,3)]
C                the planar fit 3*3 untilt-matrix
C     Alpha    : [REAL*8]
C                tiltangle alpha in degrees
C     Beta     : [REAL*8]
C                tiltangle beta in degrees
C     Gamma    : [REAL*8]
C                Fixed yaw-angle in degrees associated with mean over all runs
C     WBias    : [REAL*8]
C                The bias in the vertical velocity
C
C AUTHOR
C     Arjan van Dijk 
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     EC_C_T02
C     ***

      IMPLICIT NONE
      LOGICAL SingleRun,DoWBias
      INTEGER i,j,k,NRuns
      REAL*8 uMean(3,NRuns),Apf(3,3),Alpha,Beta,Gamma,WBias,USum(3),
     &  UUSum(3,3)
C
C Make all sums of mean velocity components and of products of
C mean velocity components in relation W.48
C
      DO i=1,3
        USum(i) = 0.D0
	DO j=1,3
	  UUSum(i,j) = 0.D0
	ENDDO
      ENDDO

      DO k=1,NRuns
	DO i=1,3
          USum(i) = USum(i) + UMean(i,k)
	  DO j=1,3
	    UUSum(i,j) = UUSum(i,j) + UMean(i,k)*UMean(j,k)
	  ENDDO
	ENDDO
      ENDDO
      DO i=1,3
 	USum(i) = USum(i)/NRuns
	DO j=1,3
 	  UUSum(i,j) = UUSum(i,j)/NRuns
	ENDDO
      ENDDO
C
C Call to slave-routine for details
C
      CALL EC_C_T02(DoWBias,SingleRun,uSum,UUSum,Apf,
     &  Alpha,Beta,Gamma,WBias)
      END

      SUBROUTINE EC_C_T02(DoWBias,SingleRun,uSum,UUSum,Apf,
     &  Alpha,Beta,Gamma,WBias)
C     ****f* ec_corr.f/EC_C_T02
C NAME
C     EC_C_T02
C SYNOPSIS
C     CALL EC_C_T02(DoWBias,SingleRun,uMean,NRuns,Apf, 
C                   Alpha,Beta,Gamma,WBias)
C FUNCTION
C     Subroutine computes angles and untilt-matrix needed for
C     tilt correction of Sonic data, using Planar Fit
C     Method, as described in James M. Wilczak et al (2001), 'Sonic
C     Anemometer tilt correction algorithms', Boundary Meteorology 99: 127:150
C     References to formulae are to this article.
C     The planar fit matrix is extended with an additional yaw-correction
C     to turn the first coordinate into the direction of the mean wind
C     over all runs. This extra rotation makes results from different
C     eddy-covariance systems comparable.
C     Furthermore, there is the option to determine a planar fit for
C     a single run (using all individual samples within a run,
C     rather than the mean velocities from a collection of runs as
C     in the classic planar fit method).
C
C INPUTS
C     DoWBias:   [LOGICAL]
C                compute bias in mean vertical wind (FALSE 
C                implies that the mean vertical wind over all
C                runs is assumed to be zero
C     SingleRun: [LOGICAL]
C                Determine rotation for a single run 
C     uMean    : [REAL*8(3,Nmax)]
C                matrix of run mean velocity vectors
C     NRuns    : [INTEGER]
C                the number of runs
C OUTPUT
C     Apf      : [REAL*8(3,3)]
C                the planar fit 3*3 untilt-matrix
C     Alpha    : [REAL*8]
C                tiltangle alpha in degrees
C     Beta     : [REAL*8]
C                tiltangle beta in degrees
C     Gamma    : [REAL*8]
C                Fixed yaw-angle in degrees associated with mean over all runs
C     WBias    : [REAL*8]
C                The bias in the vertical velocity
C
C AUTHOR
C     Arjan van Dijk 
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     EC_M_InvM
C     EC_M_InvM2
C     EC_M_MapVec
C     EC_M_Map2Vec
C     EC_M_MMul
C     ***
C
C Supportive routine for planar fit method for tilt-correction
C
      IMPLICIT NONE
      LOGICAL SingleRun,DoWBias
      REAL*8 b(0:2),S(3,3),SInv(3,3),x(3),Apf(3,3),SS2(2,2),SS2Inv(2,2),
     &  Sqrt1,Sqrt2,Alpha,Beta,SinAlpha,SinBeta,CosAlpha,CosBeta,Pi,
     &  Gamma,WBias,SinGamma,CosGamma,UHor,Yaw(3,3),USum(3),UUSum(3,3),
     &  u,v,w,R11,R12,R13,R22,R23,R33,A,Disc

      Pi = 4.D0*ATAN(1.D0)


      IF (.NOT.SingleRun) THEN
        IF (DoWBias) THEN
C
C Make matrix in relation W.48
C
	  S(1,1) = 1.D0
	  S(1,2) = USum(1)
	  S(1,3) = USum(2)
	  S(2,1) = USum(1)
	  S(2,2) = UUSum(1,1)
	  S(2,3) = UUSum(1,2)
	  S(3,1) = USum(2)
	  S(3,2) = UUSum(1,2)
	  S(3,3) = UUSum(2,2)
C
C Invert this matrix
C
          CALL EC_M_InvM(S,Sinv)
C
C Make RHS of relation W.48
C
	  x(1) = USum(3)
	  x(2) = UUSum(1,3)
	  x(3) = UUSum(2,3)
C
C Calculate coefficients b0, b1 and b2 in relation W.48
C
          CALL EC_M_MapVec(SInv,x,b(0))
C
C Find the bias in the vertical velocity via relation W.39
C
          WBias = b(0)
        ELSE
C
C Make the sub-matrix of relation W.48 for a plane through origin
C
          SS2(1,1) = UUSum(1,2)
          SS2(1,2) = UUSum(2,2)
          SS2(2,1) = UUSum(1,1)
          SS2(2,2) = UUSum(2,1)
C
C Invert this matrix
C
          CALL EC_M_InvM2(SS2,SS2inv)
C
C Make RHS of relation W.48 for this submatrix
C
          x(1) = UUSum(2,3)
          x(2) = UUSum(1,3)
C
C Calculate coefficients b1 and b2 in relation W.48
C
          CALL EC_M_Map2Vec(SS2Inv,x,b(1))
C
C Assume that the calibration of the sonic is ok and has no bias in w
C
          WBias = 0.D0
	ENDIF
      ELSE
	U = USum(1)
	V = USum(2)
	W = USum(3)
	R11 = UUSum(1,1) - U*U
	R12 = UUSum(1,2) - U*V
	R13 = UUSum(1,3) - U*W
	R22 = UUSum(2,2) - V*V
	R23 = UUSum(2,3) - V*W
	R33 = UUSum(3,3) - W*W
C
C Find the plane through both origin and mean velocity such that
C velocity-components in the plane, perpendicular to the mean velocity
C do not correlate with velocity components normal to the plane.
C In practice this means <v'w'>=0 (see roll-correction).
C
        IF ((R23.GT.(1.D-10*(R11+R22+R33))).OR.
     &  ((V*V+W*W).GT.((U*U+V*V+W*W)*1.D-10))) THEN
      A = u**3*R23-w*R22*v*u+v*R11*w*u+v**2*R23*u
     &   -v*R13*u**2+w*R12*v**2-v**3*R13-w*R12*u**2
      Disc = w**4*R11**2+4*w**2*R13**2*u**2+4*w**4*R12**2+v**4*R11**2
     &+4*v
     # **2*R13**2*u**2+8*u**2*R23*v*R11*w-
     # 4*u**2*R23*v*R33*w-4*w*R22*v*u**2*R23-
     # 4*v**3*R11*u*R12+4*v**2*R12**2*w**2+
     # 4*v**3*R23*R11*w-4*v**3
     # *R23*R33*w-2*v**2*R11*w**2*R22-
     # 8*u**3*R23*w*R12+2*u**2*R33*w**2*R11-
     # 2*u**2*R22*w**2*R11+4*v*R12*u**3*R33-
     # 4*w**3*R22*v*R23+4*v*R11*w**3*R23+
     # 8*u*R22*v**2*R13*w+4*u**3*R22*w*R13-
     # 2*v**2*R33*u**2*R22-2*u
     # **2*R33*w**2*R22-4*v**2*R33*u*w*R13-
     # 8*v**3*R12*R13*w+4*v**2*R12**2
     # *u**2-4*v*R12*u*w**2*R11+2*v**2*R11**2*w**2-
     # 8*v**3*R23*u*R13+2*v**
     # 2*R33**2*u**2+
     & 4*v**4*R13**2+4*u**4*R23**2+u**4*R22**2+
     # u**4*R33**2+4*w**3*R22*u*R13-
     # 4*w**3*R11*u*R13-
     # 8*w**3*R12*v*R13+2*v**2*R33*w**2*R22-
     #  2*v**2*R33*w**2*R11-4*u**3*R33*w*R13+
     # 8*u*R33*w**2*R12*v-4*v*R12*u**3*
     # R22-4*v*R12*u*w**2*R22-8*u*R23*w**3*R12-
     # 8*u**3*R23*v*R13-4*v**2*R11*u*w*R13+
     # 4*v**3*R12*u*R33+v**4*R33**2+
     # 2*v**2*R11*u**2*R22-2*v**2*
     # R11*u**2*R33+4*u**2*R23**2*w**2+
     # 4*w**2*R12**2*u**2-2*u**4*R22*
     # R33+w**4*R22**2+4*v**2*R23**2*u**2+
     # 2*u**2*R22**2*w**2-2*v**4*R11*R33-
     # 2*w**4*R22*R11+4*v**2*R23**2*w**2+
     # 4*v**2*R13**2*w**2

      b(1) = -0.5D0*(-v**3*R11+2*u*R12*v**2+2*w**2*R12*u-2*u**2*R23*w
     &  -u**2*R22*v+u**2*v*R33+w**2*R22*v-w**2*R11*v-2*w*R23*v**2
     &  +v**3*R33+v*SQRT(Disc))/a
      b(2) = 0.5D0*(-v**2*R11*u+2*v*R12*u**2+v**2*R33*u
     &  -u**3*R22+u**3*R33-w**2*R22*u+w**2*R11*u
     &  +2*w**2*R12*v-2*v**2*R13*w-2*w*R13*u**2+u*SQRT(Disc))/a
	ELSE
	  b(1) = 0.D0
          b(2) = 0.D0
       ENDIF
      ENDIF
C
C Construct the factors involved in the planar angles
C
      Sqrt1 = SQRT(b(2)*b(2)+1)
      Sqrt2 = SQRT(b(1)*b(1)+b(2)*b(2)+1)
C
C Planar tilt angles alpha and beta in relation W.44
C
      SinAlpha = -b(1)/Sqrt2
      CosAlpha = Sqrt1/Sqrt2
      SinBeta = b(2)/Sqrt1
      CosBeta = 1.D0/Sqrt1

      Alpha = 180.D0/Pi*ATAN2(SinAlpha,CosAlpha)
      Beta = 180.D0/Pi*ATAN2(SinBeta,CosBeta)
C
C Planar (un-)tilt matrix P from relation W.36
C
      Apf(1,1) = CosAlpha
      Apf(1,2) = SinAlpha*SinBeta
      Apf(1,3) = -SinAlpha*CosBeta
      Apf(2,1) = 0.D0
      Apf(2,2) = CosBeta
      Apf(2,3) = SinBeta
      Apf(3,1) = SinAlpha
      Apf(3,2) = -CosAlpha*SinBeta
      Apf(3,3) = CosAlpha*CosBeta
C
C Additional yaw-correction to align the first coordinate axis with
C the mean velocity over all runs according to relation W.45
C
      CALL EC_M_MapVec(Apf,USum,USum)

      UHor = (USum(1)**2+USum(2)**2)**0.5D0
      SinGamma = USum(2)/UHor
      CosGamma = USum(1)/UHor
      Gamma = 180.D0*ACOS(CosGamma)/PI
      IF (SinGamma.LT.0.D0) Gamma = 360.D0-Gamma
      IF (Gamma.GT.180.D0) Gamma = Gamma-360.D0

      Yaw(1,1) = CosGamma
      Yaw(1,2) = SinGamma
      Yaw(1,3) = 0.D0
      Yaw(2,1) = -SinGamma
      Yaw(2,2) = CosGamma
      Yaw(2,3) = 0.D0
      Yaw(3,1) = 0.D0
      Yaw(3,2) = 0.D0
      Yaw(3,3) = 1.D0

      CALL EC_M_MMul(Yaw,Apf,Apf)

      END
C
C ########################################################################
C
C Here comes a set routines that actually implement all corrections
C necessary to process eddy-correlation measurements.
C They do contain interesting models and physics/mathematics!
C All following routines affect mean quantities and co-variances only.
C No raw data-series necessary. Therefore one can use these routines
C on pre-processed data (raw mean values and raw co-variances).
C
C ########################################################################
C



      SUBROUTINE EC_C_T03(Mean,NMax,N,Cov,Speed,Stress,DumVecs,NN)
C     ****f* ec_corr.f/EC_C_T03
C NAME
C     EC_C_T03
C SYNOPSIS
C     CALL EC_C_T03(Mean,NMax,N,Cov,Speed,Stress,DumVecs,NN)
C INPUTS
C     Mean    : [REAL*8(NMax)]
C               means of all variables
C     NMax    : [INTEGER]
C               maximum number of variables (i.e. size of 
C               various matrices)
C     N       : [INTEGER]
C               number of variables in used
C     Cov     : [REAL*8(NMmax, NMax)]
C               covariances of all variables
C     NN      : [INTEGER]
C               maximum size of second axis of DumVecs
C OUTPUT
C     Speed   : [REAL*8(3)]
C               copy of means of all variables 
C     Stress  : [REAL*8(3,3)]
C               copy of covariances of all variables
C     DumVecs : [REAL*8(3,4:NN)]
C               copy of covariances of all variables with velocity components
C FUNCTION
C     Help routine for routines for correction of coordinate system
C     Temporaly stores means and covariances elsewhere 
C AUTHOR
C     Arjan van Dijk 
C HISTORY
C     $Name$ 
C     $Id$
C     ***
C USES 
C     parcnst.inc
      IMPLICIT NONE
      INCLUDE 'parcnst.inc'

      INTEGER i,j,NMax,N,NN
      REAL*8 Stress(3,3),Speed(3),DumVecs(3,4:NN),Mean(NMax),
     &	Cov(NMax,NMax)
C
C Loop : All velocity components
C
      DO i = U,W
	Speed(i) = Mean(i)
	DO j = U,W
	  Stress(i,j) = Cov(i,j)
	ENDDO
	DO j = 4,N
	  DumVecs(i,j) = Cov(i,j)
	ENDDO
      ENDDO

      RETURN
      END


      SUBROUTINE EC_C_T04(Speed,Stress,DumVecs,NN,Mean,NMax,N,Cov)
C     ****f* ec_corr.f/EC_C_T04
C NAME
C     EC_C_T04
C SYNOPSIS
C     CALL EC_C_T04(Speed,Stress,DumVecs,NN, Mean,NMax,N,Cov)
C INPUTS
C     Speed   : [REAL*8(3)]
C               copy of means of all variables 
C     Stress  : [REAL*8(3,3)]
C               copy of covariances of all variables
C     DumVecs : [REAL*8(3,4:NN)]
C               copy of covariances of all variables with velocity components
C     NMax    : [INTEGER]
C               maximum number of variables (i.e. size of 
C               various matrices)
C     N       : [INTEGER]
C               number of variables in used
C     NN      : [INTEGER]
C               maximum size of second axis of DumVecs
C OUTPUT
C     Mean    : [REAL*8(NMax)]
C               means of all variables
C     Cov     : [REAL*8(NMmax, NMax)]
C               covariances of all variables
C FUNCTION
C     Help routine for routines for correction of coordinate system
C     Copies back temporarily stored copies of means and covariances 
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
      INTEGER i,j,NMax,N,NN
      REAL*8 Stress(3,3),Speed(3),DumVecs(3,4:NN),Mean(NMax),
     &	Cov(NMax,NMax)
C
C Loop : All velocity components
C
      DO i = U,W
	Mean(i) = Speed(i)
	DO j = U,W
	  Cov(i,j) = Stress(i,j)
	ENDDO
	DO j = 4,N
	  Cov(i,j) = DumVecs(i,j)
	  Cov(j,i) = Cov(i,j)
	ENDDO
      ENDDO

      RETURN
      END





      SUBROUTINE EC_C_T05(Mean,NMax,N,Cov,Map)
C     ****f* ec_corr.f/EC_C_T05
C NAME
C     EC_C_T05
C SYNOPSIS
C     CALL EC_C_T05(Mean,NMax,N,Cov,Map)
C INPUTS
C     Mean    : [REAL*8(NMax)]
C               means of all variables
C     NMax    : [INTEGER]
C               maximum number of variables (i.e. size of 
C               various matrices)
C     N       : [INTEGER]
C               number of variables in used
C     Cov     : [REAL*8(NMmax, NMax)]
C               covariances of all variables
C     Map     : [REAL*8(3,3)]
C               rotation tensor          
C OUTPUT
C     Mean    : [REAL*8(NMax)]
C               means of all variables
C     Cov     : [REAL*8(NMmax, NMax)]
C               covariances of all variables
C FUNCTION
C     Routine to change coordinate system according to tensor "Map".
C     This routine is called by all tilt-correction procedures.
C     Both the mean velocity and the Reynoldsstresses and the
C     covariances of all velocity components with other quantities
C     are rotated.
C AUTHOR
C     Arjan van Dijk 
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     EC_C_T03
C     EC_C_T04
C     EC_M_MapVec
C     EC_M_MapMtx
C     parcnst.inc
C     ***
      IMPLICIT NONE
      INCLUDE 'parcnst.inc'

      INTEGER NMax,N,j
      REAL*8 Map(3,3),Stress(3,3),Speed(3),
     &  DumVecs(3,4:NNmax),Mean(NMax),Cov(NMax,NMax)

      CALL EC_C_T03(Mean,NMax,N,Cov,Speed,Stress,DumVecs,NNMax)
      CALL EC_M_MapVec(Map,Speed,Speed)
      CALL EC_M_MapMtx(Map,Stress,Stress)
      DO j = 4,N
        CALL EC_M_MapVec(Map,DumVecs(1,j),DumVecs(1,j))
      ENDDO
      CALL EC_C_T04(Speed,Stress,DumVecs,NNMax,Mean,NMax,N,Cov)

      RETURN
      END





      SUBROUTINE EC_C_T06(Direction,Yaw)
C     ****f* ec_corr.f/EC_C_T06
C NAME
C     EC_C_T06
C SYNOPSIS
C     CALL EC_C_T06(Direction,Yaw)
C INPUTS
C     Direction : [REAL*8]
C                 yaw angle 
C OUTPUT
C     Yaw       : [REAL*8(3,3)]
C                 rotation tensor
C FUNCTION
C     Construct rotation matrix for coordinate system 
C     about a KNOWN yaw-angle around the vertical
C AUTHOR
C     Arjan van Dijk 
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     Pi
C     ***
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'

      REAL*8 SinPhi,CosPhi,Yaw(3,3),Direction

      SinPhi = SIN(PI*Direction/180.D0)
      CosPhi = COS(PI*Direction/180.D0)

      Yaw(1,1) = CosPhi
      Yaw(1,2) = SinPhi
      Yaw(1,3) = 0.D0
      Yaw(2,1) = -SinPhi
      Yaw(2,2) = CosPhi
      Yaw(2,3) = 0.D0
      Yaw(3,1) = 0.D0
      Yaw(3,2) = 0.D0
      Yaw(3,3) = 1.D0

      RETURN
      END





      SUBROUTINE EC_C_T07(MeanU,MeanV,Direction)
C     ****f* ec_corr.f/EC_C_T07
C NAME
C     EC_C_T07
C SYNOPSIS
C     CALL EC_C_T07(MeanU,MeanV,Direction)
C INPUTS
C     MeanU     : [REAL*8]
C                 mean u-velocity
C     MeanV     : [REAL*8]
C                 mean v-velocity
C OUTPUT
C     Direction : [REAL*8]
C                 yaw angle (degree)
C FUNCTION
C     Give yaw-angle to transform coordinate system such that
C     v_mean = 0 (no mean lateral horizontal velocity component)
C AUTHOR
C     Arjan van Dijk 
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     Pi
C     ***
      INCLUDE 'physcnst.inc'

      REAL*8 UHor,SinPhi,CosPhi,MeanU,MeanV,Direction

      UHor = (MeanU**2+MeanV**2)**0.5D0
      SinPhi = MeanV/UHor
      CosPhi = MeanU/UHor
      Direction = 180.D0*ACOS(CosPhi)/PI
      IF (SinPhi.LT.0.D0) Direction = 360.D0-Direction

      RETURN
      END





      SUBROUTINE EC_C_T08(Direction,Pitch)
C     ****f* ec_corr.f/EC_C_T08
C NAME
C     EC_C_T08
C SYNOPSIS
C     CALL EC_C_T08(Direction,Yaw)
C INPUTS
C     Direction : [REAL*8]
C                 pitch angle 
C OUTPUT
C     Pitch     : [REAL*8(3,3)]
C                 rotation tensor
C FUNCTION
C     Give matrix for rotation about a KNOWN pitch-angle 
C     around vector (0,1,0).
C AUTHOR
C     Arjan van Dijk 
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     Pi
C     ***
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'

      REAL*8 SinTheta,CosTheta,Pitch(3,3),Direction

      SinTheta = SIN(PI*Direction/180.D0)
      CosTheta = COS(PI*Direction/180.D0)

      Pitch(1,1) = CosTheta
      Pitch(1,2) = 0.D0
      Pitch(1,3) = SinTheta
      Pitch(2,1) = 0.D0
      Pitch(2,2) = 1.D0
      Pitch(2,3) = 0.D0
      Pitch(3,1) = -SinTheta
      Pitch(3,2) = 0.D0
      Pitch(3,3) = CosTheta

      RETURN
      END




      SUBROUTINE EC_C_T09(MeanU,MeanW,Direction)
C     ****f* ec_corr.f/EC_C_T09
C NAME
C     EC_C_T09
C SYNOPSIS
C     CALL EC_C_T09(MeanU,MeanW,Direction)
C INPUTS
C     MeanU     : [REAL*8]
C                 mean u-velocity
C     MeanV     : [REAL*8]
C                 mean W-velocity
C OUTPUT
C     Direction : [REAL*8]
C                 pitch angle (degree)
C FUNCTION
C     Give pitch angle to transform coordinate system such 
c     that w_mean = 0 (no mean vertical velocity component)
C AUTHOR
C     Arjan van Dijk 
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     Pi
C     ***
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'

      REAL*8 SinTheta,MeanU,MeanW,Direction,UTot

      UTot = (MeanU**2+MeanW**2)**0.5D0
      SinTheta = MeanW/UTot
      Direction = 180.D0*ASIN(SinTheta)/PI

      RETURN
      END





      SUBROUTINE EC_C_T10(Direction,Roll)
C     ****f* ec_corr.f/EC_C_T10
C NAME
C     EC_C_T10
C SYNOPSIS
C     CALL EC_C_T10(Direction,Roll)
C INPUTS
C     Direction : [REAL*8]
C                 yaw angle 
C OUTPUT
C     Roll      : [REAL*8(3,3)]
C                 rotation tensor
C FUNCTION
C     Give matrix to rotate coordinate system about a KNOWN 
C     roll-angle around vector (1,0,0).
C AUTHOR
C     Arjan van Dijk 
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     Pi
C     ***
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'

      REAL*8 SinRoll,CosRoll,Roll(3,3),Direction

      SinRoll = SIN(PI*Direction/180.D0)
      CosRoll = COS(PI*Direction/180.D0)

      Roll(1,1) = 1.D0
      Roll(1,2) = 0.D0
      Roll(1,3) = 0.D0
      Roll(2,1) = 0.D0
      Roll(2,2) = CosRoll
      Roll(2,3) = SinRoll
      Roll(3,1) = 0.D0
      Roll(3,2) = -SinRoll
      Roll(3,3) = CosRoll

      RETURN
      END





      SUBROUTINE EC_C_T11(CovVV,CovVW,CovWW,Direction)
C     ****f* ec_corr.f/EC_C_T11
C NAME
C     EC_C_T11
C SYNOPSIS
C     CALL EC_C_T11(MeanU,MeanW,Direction)
C INPUTS
C     CovVV     : [REAL*8]
C                 vv-covariance
C     CovVW     : [REAL*8]
C                 vw-covariance
C     CovWW     : [REAL*8]
C                 ww-covariance
C OUTPUT
C     Direction : [REAL*8]
C                 roll angle (degree)
C FUNCTION
C     Give roll angle to transform coordinate system such 
C     that Cov(V,W) = 0  (vertical velocity fluctuations are
C     independent from horizontal fluctuations)
C AUTHOR
C     Arjan van Dijk 
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     Pi
C     ***
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'

      REAL*8 CovVV,CovVW,CovWW,Direction,RollAngl,Arg1,Arg2

      Arg1 = 2*CovVW
      Arg2 = -CovVV+CovWW
      IF (Arg2.LT.0.D0) THEN
        RollAngl =  0.5D0*ATAN2(Arg1,-Arg2)
      ELSE
        RollAngl = -0.5D0*ATAN2(Arg1, Arg2)
      ENDIF
      Direction = 180.D0*RollAngl/Pi

      RETURN
      END




      SUBROUTINE EC_C_Webb(Mean,NMax,Cov,P,WhichTemp, WebVel)
C     ****f* ec_corr.f/EC_C_Webb
C NAME
C     EC_C_Webb
C SYNOPSIS
C     CALL EC_C_Webb(Mean,NMax,Cov,P,WhichTemp, WebVel)
C INPUTS
C     Mean     : [REAL*8(NMax)]
C                mean of all variables
C     NMax     : [INTEGER]
C                maximum number of variables 
C     Cov      : [REAL*8(NMax,NMax)]
C                covariance of all variables
C     P        : [REAL*8]
C                atmospheric pressure
C     WhichTemp: [INTEGER]
C                Use thermocouple temperature or sonic temperature ?
C                (Tcouple or TSonic, codes in parcnst.inc)
C     Mean     : [REAL*8(NMax)]
C                mean of all variables 
C OUTPUT
C     WebVel    : [REAL*8]
C                Webb velocity (m/s)
C FUNCTION
C     Mean vertical velocity according to Webb, Pearman and Leuning
C AUTHOR
C     Arjan van Dijk 
C HISTORY
C     28-05-2001: added info on which temperature should be used
C                 in corrections (Sonic or thermocouple)
C     07-10-2002: Webb velocity no longer passed through Mean(W). Now
C                 separate variable, WebVel is used for this.
C     $Name$ 
C     $Id$
C USES
C     Mu
C     parcnst.inc
C     EC_Ph_RhoDry
C     ***
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER NMax, WhichTemp
      REAL*8 Mean(NMax),Cov(NMax,NMax),Sigma,EC_Ph_RhoDry,P, WebVel
C
C Ratio of mean densities of water vapour and dry air
C
      Sigma = Mean(Humidity)/EC_Ph_RhoDry(Mean(Humidity),
     &                                    Mean(WhichTemp),P)
      WebVel = (1.D0+Mu*Sigma)*Cov(W,WhichTemp)/Mean(WhichTemp) +
     &	Mu*Sigma*Cov(W,Humidity)/Mean(Humidity)

      RETURN
      END
