


C
C ########################################################################
C
C
C     This set of routines provides the matrix for the correction of
C     turbulent air flow measurements for the presence of
C     small disturbing objects, like a box with electronic
C     apparatus. The approach followed is described in:
C     "Flow distortion by an ellipsoid and its application
C     to atmpspheric measurements" by W.A. Oost, May 1989,
C     KNMI internal report, accepted by J.Atm.Oc.Tech.,'90.
C     References in the code are to this article.
C
C
C ########################################################################
C

C
C The generic call is to the following subroutine:
C
      SUBROUTINE EC_C_D01(x,b,aInv)
C
C Input: x : Position vector of the point where measurements have been taken.
C            The ellipsoid is placed in the origin.
C            A right-handed frame of coordinates is chosen.
C            The flow is supposed to be expressed in this coordinate frame.
C            Therefore, when the flow velocity has positive components,
C            upstream measurement points are selected when by giving
C            vector x negative components!
C        b : Three ellipsoid semi-axes in meters. The ellipsoid is
C            supposed to be oriented along the coordinate axes.
C            Somehow the algorithm does not seem to like it when two or more
C            semi-axes are equal, or when your point x is in one of the
C            coordinate planes (one component of x equal to zero).
C            To circumvent problems one can take values slightly off
C            the problematic values.
C Output: aInv : The matrix which can be used to correct samples and
C            covariances for flow distortion.
C            Sum_j aInv(i,j)*u(j)     gives the distortion-corrected
C            image of measured velocity u.
C
      REAL*8 x(3),b(3),a(3,3),aInv(3,3)
      CALL EC_C_D02(x,b,a)
      CALL EC_M_InvM(a,aInv)
      RETURN
      END






      SUBROUTINE EC_C_D02(x,b,a)
C     Calculates the distortion matrix according to equation 12
C     INPUT : b = a vector containing the three semiaxes of the
C                 ellipsoid
C             x = position where the distortion-matrix will be
C                 calculated; coordinates are relative to the
C                 center of the ellipsoid
C     OUTPUT: a = the distortion matrix, when applied to an
C                 undisturbed wind, "a" gives the disturbed wind
      INTEGER Lambda, Mu, Nu
      PARAMETER(Lambda=1,Mu=2,Nu=3)
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
     &      ((2.D0-alphaint(j))*(b(j)**2+y(lambda))*EC_M_ka(y(lambda),b))
            IF (i.EQ.j) a(perm(i),perm(j)) = a(perm(i),perm(j))
     &      + 1.D0 + b(1)*b(2)*b(3)*integral(i)/(2.D0-alphaint(i))
         END DO
      END DO
      CALL EC_M_UnSort(x,perm)
      CALL EC_M_UnSort(b,perm)
      RETURN
      END





      SUBROUTINE EC_C_F01(Mean,Cov,NMax,NSize,WhichTemp,
     &	NSta,NEnd,NInt,NS,TauD,TauV,CalSonic,CalTherm,CalHyg,WXT)
C
C Calculate frequency response corrections for sensor response, path
C length averaging, sensor separation, signal processing and dampening
C of fluctuations in a tube. Based on publications;
C
C	Moore, C.J. (1986): 'Frequency Response Corrections for Eddy
C	Correlation Systems'. Boundary Layer Met. 37: 17-35.
C
C	Philip, J.R. (1963): 'The Damping of Fluctuating Concentration
C	by Continuous Sampling Through a tube' Aust. J. Phys. 16: 454-463.
C
C	Leuning, R. and K.M. King (1991): 'Comparison of Eddy-Covariance
C	Measurements of CO2 Fluxes by open- and closed-path CO2 analysers'
C	(unpublished)
C
C input : DoPrint [LOGICAL] : Indicator if intermediate results must be
C           printed to file OutF.
C         OutF [INTEGER] : Unit number of output-file for intermediate
C           results. OutF must be an open file.
C         Mean [REAL*8(NMax)] : Array of mean values of the quantities in
C           this experiment (only the first N quantities are used).
C         TolMean [REAL*8(NMax)] : Tolerances of Mean.
C         NMax [INTEGER] : Physical dimension of array Mean
C         N [INTEGER] : Number of quantities actually involved in this
C           experiment.
C         Cov [REAL*8(NMax,NMax)] : covariances of the fluctuations.
C         TolCov [REAL*8(NMax,NMax)] : Tolerances of covariances (2 sigma).
C         QName [CHARACTER*9(NMax)] : Names of the quantities.
C         UName [CHARACTER*9(NMax)] : Names of the units of the quantities.
C         LLimit [REAL*8] : Lower acceptance limit for frequency-response
C           factors. Correction factors smaller than LLimit are set to 1.
C         ULimit [REAL*8] : Upper acceptance limit for frequency-response
C           factors. Correction factors larger than ULimit are set to 1.
C         NSta [REAL*8] : Start frequency numerical integration.
C           Popular value: -5.D0 [unit?].
C         NEND [REAL*8] : End frequency numerical integration.
C           Popular value: LOG(5) = 0.69897D0 [unit?].
C         NINT [INTEGER] : Number of intervals in integration.
C           Popular value: 19.
C         Freq [REAL*8] : Sampling frequency [Hz].
C         TauD [REAL*8] : Interval length for running mean.
C           Popular value: 0.D0 [unit?].
C         TAUV [REAL*8] : Low pass filter time constant.
C           Popular value: 0.D0 [unit?]
C         CalSonic [REAL*8(NQQ)] : Calibration specification array of
C           sonic anemometer.
C         CalTherm [REAL*8(NQQ)] : Calibration specification array of
C           thermometer.
C         CalHyg [REAL*8(NQQ)] : Calibration specification array of
C           hygrometer.
C
C output : FrCor [REAL*8(NMax,NMax)] : Correction factors for covariances.
C
C Revision 28-05-2001: added info on which temperature should be used
C                      in corrections (Sonic or thermocouple)
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER I1,I2,I3,I4,I,J,NINT,NMax,NSize, WhichTemp
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ),
     &	GAIND,GAINV,GAINT1,GAINT2,GAINUV,GAINW,GAINQ1,
     &	TRANA,TRANQ1,TRANUV,TRANW,TRANT1,TRANT2,TRANWT1,TRANWT2,TRANWQ1,
     &	TRANWUV,ATT,AWW,AUU,AWT,AUW,BTT,BWW,BUU,BWT,BUW,TAUD,TAUV,TAUT,
     &	XI,ZETA,CW,CU,WT,UW,TT,UU,WW,F,NS,C,NSTA,NEND,N,LF,X,ZL,
     &	INTTS(15),G(15),WXT(NMax,NMax),Mean(NMax),Cov(NMax,NMax),UStar,
     &	Lower,Upper,TolCov(NMax,NMax)
C
C-- Declarations of constants ----------------------------------------
C
      TAUT = CalTherm(QQTime)/
     &	(1.D0+4.9D0*(CalTherm(QQTime)**0.5*Mean(U))**0.45D0)
      LF   = (NEND-NSTA)/DBLE(NINT-1)	    ! Interval width num. integr.
      C    = 1.D0/(LOG(10.D0)*LF/3.D0)	    ! Constant of integration, eq.1
      LF   = 10.D0**LF
      N    = 10.D0**NSTA
      DO I=1,15
	INTTS(I)=0.D0
      ENDDO

      Ustar = ((Cov(U,W))**2.D0)**0.25D0
      ZL = (CalSonic(QQZ)*Karman*GG*Cov(W,WhichTemp))/
     &	(-Mean(WhichTemp)*Ustar**3.D0)

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
C     From -5 to log(5) in NInt steps
C
      I3=0
      DO I1=1,10
	DO I2=2,4,2
	  I3=I3+1
	  IF (I3.GT.NInt) GOTO 200
	  I4=I2
	  IF (I3.EQ.1)	I4=I4-1
	  IF (I3.EQ.NInt) I4=I4-1
C
C-- Skip integration if frequency exceeds Nyquist frequency ------
C
	  IF (N.GT.(0.5D0*NS)) GOTO 999

C
C-- ELECTRONIC / DIGITAL FILTERS !!
C-- Auto Regressive Moving Average filter response gain, eq. 19 --
C
	  GAIND=1.D0
	  X=2.D0*PI*N*TAUD
	  IF (X.LT.6.D0 .AND. X.GT.0.D0) GAIND=X*X/(1.D0+X*X)
C
C-- Butterworth filter frequency response gain, eq. 5 ------------
C
	  GAINV = 1.D0
	  GAINV=(2.D0*PI*N*TAUV)**2.D0
	  GAINV=1.D0+GAINV*GAINV
C
C-- Data aquisition co-spectral transfer function, eq. 15 (b=3) --
C
	  TRANA=1.D0+(N/(NS-N))**3.D0
C
C-- LUW THERMOCOUPLE TEMPERATURE !!
C-- Thermocouple frequency response gain, eq. 2 ------------------
C
	  GAINT1=1.D0+(2.D0*PI*N*TAUT)**2.D0
C
C-- Thermocouple spatial averaging transfer function, eq. 7 ------
C
	  TRANT1=1.D0	       !spatial averaging negligible
C
C-- W-T lateral separation transfer function, eq. 11 -------------
C
	  TRANWT1=1.D0
	  X=N/Mean(U)*(CalTherm(QQX)-CalSonic(QQX))
	  IF (X.GT.0.01D0) TRANWT1=EXP(-9.9D0*X**1.5D0)
C
C-- SOLENT SONIC TEMPERATURE !!
C-- Sonic temperature frequency response gain, eq. 2 -------------
C
	  GAINT2=1.D0	       !sonic temperature gain negligible
C
C-- Sonic temperature spatial averaging transfer function, eq.7 --
C
	  TRANT2=1.D0
	  X=2.D0*PI*N/Mean(U)*CalSonic(QQExt4)
	  IF (X.GT.0.02D0) TRANT2=(3.D0+EXP(-X)-4.D0*(1.D0-EXP(-X))/X)/X
C
C-- W-T lateral separation transfer function, eq. 11 -------------
C
	  TRANWT2=1.D0	       !there is no lateral separation

C
C-- SOLENT SONIC HORIZONTAL WINDSPEED !!
C-- UV-sensor frequency response gain, eq. 2 ---------------------
C
	  GAINUV=1.D0	       !sonic UV windspeed gain negligible
C
C-- UV-sensor spatial averaging transfer function, eq. 7 ---------
C
	  TRANUV=1.D0
	  X=2.D0*PI*N/Mean(U)*CalSonic(QQPath)
	  IF (X.GT.0.02D0) TRANUV=(3.D0+EXP(-X)-4.D0*(1.D0-EXP(-X))/X)/X
C
C-- W-UV lateral separation transfer function, eq. 11 ------------
C
	  TRANWUV=1.D0
	  X = N/Mean(U)*CalSonic(QQExt3)
	  IF (X .GT.0.01D0) TRANWUV = EXP(-9.9D0*X**1.5D0)

C
C-- SOLENT SONIC VERTICAL WINDSPEED !!
C-- W-sensor frequency response gain, eq. 2 ----------------------
C
	  GAINW=1.D0	       !sonic W windspeed gain neglectible
C
C-- W-sensor spatial averaging transfer function, eq. 9 ----------
C
	  TRANW=1.D0
	  X=2.D0*PI*N/Mean(U)*CalSonic(QQPath)
	  IF (X.GT.0.04D0)
     &	    TRANW=(1.D0+(EXP(-X)-3.D0*(1.D0-EXP(-X))/X)/2.D0)*4.D0/X

C
C-- LYMANN-ALPHA OPEN PATH HYGROMETER !!
C-- hygrometer frequency response gain, eq. 2 ------------
C
	  GAINQ1=1.D0	       !hygrometer gain neglectible
C
C-- lymann-alpha spatial averaging transfer function, eq.7 ------------
C
	  TRANQ1=1.D0
	  X=2.D0*PI*N/Mean(U)*CalHyg(QQPath)
	  IF (X.GT.0.02D0) TRANQ1=(3.D0+EXP(-X)-4.D0*(1.D0-EXP(-X))/X)/X
C
C-- W-Q lateral separation transfer function, eq. 11 -------------
C
	  TRANWQ1=1.D0
	  X=N/Mean(U)*CalHyg(QQX)
	  IF (X.GT.0.01D0) TRANWQ1=EXP(-9.9D0*X**1.5D0)
C
C-- Composite transfer functions, eq. 28 -------------------------
C
	  DO I=1,15
	    G(I) = 0.D0
	  ENDDO

	  G(1)	= TRANUV / GAINUV		    !UU
	  G(2)	= TRANUV / GAINUV		    !VV
	  G(3)	= TRANW  / GAINW		    !WW
	  G(4)	= TRANT1 / GAINT1		    !TTth
	  G(5)	= TRANT2 / GAINT2		    !TTso
	  G(6)	= TRANQ1 / GAINQ1		    !QQkr
	  G(9)	= TRANWUV * SQRT (G(3) * G(1))	    !WU
	  G(10) = TRANWUV * SQRT (G(3) * G(2))	    !WV
	  G(11) = TRANWT1 * SQRT (G(3) * G(4))	    !WTth
	  G(12) = TRANWT2 * SQRT (G(3) * G(5))	    !WTso
	  G(13) = TRANWQ1 * SQRT (G(3) * G(6))	    !WQkr
	  DO I=1,15
	    G(I) = G(I)*TRANA/GAINV*GAIND
	  ENDDO

	  F=N/Mean(U)*CalSonic(QQZ)
	  IF (ZL.LT.0.D0) THEN
C
C-- Unstable normalised spectral and co-spectral forms ---------
C
	    UU=210.0D0*F/(1.D0+33.D0*F)**1.667D0		       !eq. 23
	    UU=UU+F*XI/(ZETA+2.2D0*F**1.667D0)
	    UU=UU/CU
	    WW=16.D0*F*XI/(1.D0+17.D0*F)**1.667D0		       !eq. 22
	    WW=WW+F/(1.D0+5.3D0*F**1.667D0)
	    WW=WW/CW
	    TT=14.94D0*F/(1.+24.D0*F)**1.667D0			     !eq. 24
	    IF (F.GE.0.15D0) TT=6.827D0*F/(1.D0+12.5D0*F)**1.667D0
	    WT=12.92D0*F/(1.D0+26.7D0*F)**1.375D0		       !eq. 25
	    IF (F.GE.0.54D0) WT=4.378D0*F/(1.D0+3.8D0*F)**2.4D0
	    UW=20.78D0*F/(1.D0+31.D0*F)**1.575D0		       !eq. 26
	    IF (F.GE.0.24D0) UW=12.66D0*F/(1.D0+9.6D0*F)**2.4D0
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
	  INTTS(1)  = INTTS(1)	+ I4*G(1) *UU	   !UU
	  INTTS(2)  = INTTS(2)	+ I4*G(2) *UU	   !VV
	  INTTS(3)  = INTTS(3)	+ I4*G(3) *WW	   !WW
	  INTTS(4)  = INTTS(4)	+ I4*G(4) *TT	   !TTth
	  INTTS(5)  = INTTS(5)	+ I4*G(5) *TT	   !TTso
	  INTTS(6)  = INTTS(6)	+ I4*G(6) *TT	   !QQkr
	  INTTS(9)  = INTTS(9)	+ I4*G(9) *UW	   !WU
	  INTTS(10) = INTTS(10) + I4*G(10)*UW	   !WV
	  INTTS(11) = INTTS(11) + I4*G(11)*WT	   !WTth
	  INTTS(12) = INTTS(12) + I4*G(12)*WT	   !WTso
	  INTTS(13) = INTTS(13) + I4*G(13)*WT	   !WQkr

C
C-- Increase frequency by interval width -------------------------
C
	  N=N*LF
	ENDDO
      ENDDO

200   CONTINUE

C
C-- Final flux loss correction factors !! ----------------------------
C

      DO I=1,NSize
	DO J=1,NSize
	  WXT(I,J) = 1.0D0
	ENDDO
      ENDDO

999   CONTINUE
      WXT(U,U) = C/INTTS(1)		   !UU
      WXT(V,V) = C/INTTS(2)		   !VV
      WXT(U,V) = WXT(U,U)		   !UV
      WXT(V,U) = WXT(U,V)		   !VU
      WXT(W,W) = C/INTTS(3)		   !WW
      WXT(W,U) = C/INTTS(9)		   !WU
      WXT(U,W) = WXT(W,U)
      WXT(W,V) = C/INTTS(10)		   !WV
      WXT(V,W) = WXT(W,V)

      WXT(W,TSonic) = C/INTTS(12)	 !WTso
      WXT(TSonic,W) = WXT(W,TSonic)
      WXT(TSonic,TSonic) = C/INTTS(5)	 !TTso

      WXT(W,TCouple) = C/INTTS(11)	  !WTth
      WXT(TCouple,W) = WXT(W,TCouple)
      WXT(TCouple,TCouple) = C/INTTS(4)    !TTth

      WXT(W,Humidity) = C/INTTS(13)	  !WQly
      WXT(Humidity,W) = WXT(W,Humidity)
      WXT(Humidity,Humidity) = C/INTTS(6)  !QQly

      WXT(W,SpecHum)	   = WXT(W,Humidity)
      WXT(SpecHum,W)	   = WXT(Humidity,W)
      WXT(SpecHum,SpecHum) = WXT(Humidity,Humidity)

      RETURN
      END





      SUBROUTINE EC_C_F02(WXT,NMax,NSize,Lower,Upper,Cov,TolCov)
C
C Apply frequency response corrections for sensor response, path
C length averaging, sensor separation, signal processing and dampening
C of fluctuations in a tube. Based on publications;
C
C	Moore, C.J. (1986): 'Frequency Response Corrections for Eddy
C	Correlation Systems'. Boundary Layer Met. 37: 17-35.
C
C	Philip, J.R. (1963): 'The Damping of Fluctuating Concentration
C	by Continuous Sampling Through a tube' Aust. J. Phys. 16: 454-463.
C
C	Leuning, R. and K.M. King (1991): 'Comparison of Eddy-Covariance
C	Measurements of CO2 Fluxes by open- and closed-path CO2 analysers'
C	(unpublished)
C
C Revision 28-05-2001: added info on which temperature should be used
C                      in corrections (Sonic or thermocouple)
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

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
     &	QName,UName,
     &  BadTc,
     &	DoTilt,PTilt,PreYaw,PrePitch,PreRoll,
     &	DoYaw,PYaw,DirYaw,
     &	DoPitch,PPitch,PitchLim,DirPitch,
     &	DoRoll,PRoll,RollLim,DirRoll,
     &	DoSonic,PSonic,SonFactr,
     &	DoO2,PO2,O2Factor,
     &	DoFreq,PFreq,LLimit,ULimit,Freq,CalSonic,CalTherm,CalHyg,FrCor,
     &	DoWebb,PWebb,P, Have_Uncal)
C
C Purpose : Integrated correction routine applying ALL (user-selected)
C	    corrections in this library on mean values and covariances.
C	    All intermediate results can be output to a file.
C	    Moreover they are returned to the calling routine in
C	    respective variables.
C
C Revision 28-05-2001: added info on whether uncalibrated data are
C                      available for a given variable (mainly important
C                      for sonic and/or Couple temperature since that
C                      is used for various corrections)
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER N,NInt,NMAx,OutF, WhichTemp, I, J
      LOGICAL PYaw,PPitch,PRoll,PFreq,PO2,PWebb,
     &	DoYaw,DoPitch,DoRoll,DoFreq,DoO2,DoWebb,PSonic,
     &	DoSonic,DoTilt,PTilt,DoPrint,BadTc,
     &	QYaw,QPitch,QRoll,QFreq,QO2,QWebb,QSonic,QTilt,
     &  Have_Uncal(NMax), QSchot
      REAL*8 P,Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),LLimit,ULimit,FrCor(NMax,NMax),
     &	PitchLim,RollLim,DirYaw,O2Factor(NMax),
     &	NSTA,NEND,TAUV,TauD,Freq,DirPitch,DirRoll,
     &	SonFactr(NMax),PreYaw,PrePitch,PreRoll,TSonFact,
     &  Yaw(3,3), Roll(3,3), Pitch(3,3), Dirs
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ)
      CHARACTER*6 QName(NMax)
      CHARACTER*9 UName(NMax)
C
C
C Only print intermediate results if desired
C
C
      QTilt  = (PTilt  .AND. DoPrint)
      QYaw = (PYaw .AND. DoPrint)
      QPitch   = (PPitch   .AND. DoPrint)
      QRoll  = (PRoll  .AND. DoPrint)
      QSonic = (PSonic .AND. DoPrint)
      QO2    = (PO2    .AND. DoPrint)
      QFreq  = (PFreq  .AND. DoPrint)
      QWebb  = (PWebb  .AND. DoPrint)
C A Hack: Just set Qschot off
      QSchot = .FALSE.
C
C
C Perform a tilt-correction of the whole system using KNOWN (!!!!!!) angles.
C Used to compensate for KNOWN miss-alignment!!!!!
C The classic "tilt-corrections", which turn Mean(V), Mean(W) and Cov(W,V)
C to zero, to are carried out later!!!!!!
C If no fixed tilt-angles are known, this subroutine may be skipped.
C
C
      IF (DoTilt) THEN
	CALL EC_C_T06(PreYaw,Yaw)
	CALL EC_C_T05(Mean,NMax,N,Cov,Yaw)
	CALL EC_C_T08(PrePitch,Pitch)
	CALL EC_C_T05(Mean,NMax,N,Cov,Pitch)
	CALL EC_C_T10(PreRoll,Roll)
	CALL EC_C_T05(Mean,NMax,N,Cov,Roll)
	CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)

	IF (QTilt) THEN
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'Averages of data after fixed ',
     &	  'tilt-correction: '
	  WRITE(OutF,*)
	  CALL EC_G_Show(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
	ENDIF
      ENDIF
C
C
C Perform classic yaw-correction : Mean(V) --> 0
C
C
      IF (DoYaw) THEN
	CALL EC_C_T07(Mean(U),Mean(V),DirYaw)
	CALL EC_C_T06(DirYaw,Yaw)
	CALL EC_C_T05(Mean,NMax,N,Cov,Yaw)

	CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)

	IF (QYaw) THEN
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'Yaw angle = ',DirYaw,' degrees'
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'After yaw-correction (Mean(V) -> 0) : '
	  WRITE(OutF,*)
	  CALL EC_G_Show(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
	ENDIF
      ENDIF
C
C
C Perform classic pitch-correction : Mean(W) --> 0
C
C
      IF (DoPitch) THEN
	CALL EC_C_T09(Mean(U),Mean(W),DirPitch)
	IF (ABS(DirPitch).LE.PitchLim) THEN
          CALL EC_C_T08(DirPitch,Pitch)
          CALL EC_C_T05(Mean,NMax,N,Cov,Pitch)
	ENDIF
	CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)

	IF (QPitch) THEN
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'Pitch angle = ',DirPitch,' degrees'
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'After pitch-correction (Mean(W) -> 0) : '
	  WRITE(OutF,*)
	  CALL EC_G_Show(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
	ENDIF
      ENDIF
C
C
C Perform classic roll-correction : Cov(W,V) --> 0
C
C
      IF (DoRoll) THEN
	CALL EC_C_T11(Cov(V,V),Cov(V,W),Cov(W,W),DirRoll)
	IF (ABS(DirRoll) .LE. RollLim) THEN
          CALL EC_C_T10(DirRoll,Roll)
          CALL EC_C_T05(Mean,NMax,N,Cov,Roll)
	ENDIF
	CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)

	IF (QRoll) THEN
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'Roll angle = ',DirRoll,' degrees'
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'After roll-correction (Cov(V,W) -> 0) : '
	  WRITE(OutF,*)
	  CALL EC_G_Show(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
	ENDIF
      ENDIF
C
C Reset the coveriances for which no data available
C
      DO I=1,NMax
          DO J=1,NMax
	     IF ((.NOT. HAVE_UNCAL(I)) .OR.
     +           (.NOT. HAVE_UNCAL(J))) THEN
                COV(I,J) = DUMMY
	     ENDIF
	  ENDDO
      ENDDO

C
C
C
C Correct sonic temperature and all covariances with sonic temperature
C for humidity. This is half the Schotanus-correction. Side-wind
C correction is done at calibration time directly after reading raw data.
C
C Now we need to know if and which temperature we have. Default to
C Sonic temperature
C
      IF (HAVE_UNCAL(TSonic)) THEN
         WhichTemp = Tsonic
      ELSE IF (Have_UNCAL(TCouple)) THEN
         WhichTemp = TCouple
      ELSE
         WhichTemp = -1
      ENDIF

      IF (DoSonic) THEN
        CALL EC_C_Schot1(Mean(SpecHum),Mean(TSonic),NMax,N,Cov,
     &    SonFactr,TSonFact)
        CALL EC_C_Schot2(SonFactr,TSonFact,Mean(TSonic),NMax,N,Cov)
	CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)

	IF (QSchot) THEN
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'Correction factors for H2O-sensitivity of ',
     &	  'sonic'
	  DO i=1,N
	    WRITE(OutF,45) QName(i),SonFactr(i)
	  ENDDO
   45	FORMAT('For (',a6,',T(son)) : ',F10.5)
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'After H2O-correction for sonic : '
	  WRITE(OutF,*)
	  CALL EC_G_Show(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
	ENDIF
      ENDIF
C
C
C Perform correction for oxygen-sensitivity of hygrometer
C
C
      IF (DoO2) THEN
        IF (WhichTemp .GT. 0) THEN
          CALL EC_C_Oxygen1(Mean(WhichTemp),NMax,N,Cov,P,CalHyg(QQType),
     &      WhichTemp,O2Factor)
          CALL EC_C_Oxygen2(O2Factor,NMax,N,Cov)
	  CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)

	  IF (QO2) THEN
	    WRITE(OutF,*)
	    WRITE(OutF,*) 'Correction factors for O2-sensitivity of ',
     &        'hygrometer'
	    DO i=1,N
	      WRITE(OutF,46) QName(i),O2Factor(i)
	    ENDDO
 46         FORMAT('For (',a6,',RhoV or q) : ',F10.5)
	    WRITE(OutF,*)
	    WRITE(OutF,*) 'After oxygen-correction for hygrometer : '
	    WRITE(OutF,*)
	    CALL EC_G_Show(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
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
      NSTA   = -5.D0	 ! [?] start frequency numerical integration
      NEND   = 0.69897D0 ! [?] end frequency numerical integration (LOG(5))
      NINT   = 19	 ! [1] number of intervals
      TAUV   = 0.D0	 ! [?] Low pass filter time constant
      TauD   = 0.D0	 ! [?] interval length for running mean

      IF (DoFreq) THEN
        IF (WhichTemp .GT. 0) THEN
	  CALL EC_C_F01(Mean,Cov,NMax,N,WhichTemp,
     &         NSta,NEnd,NInt,Freq,TauD,TauV,CalSonic,CalTherm,CalHyg,FrCor)
	  CALL EC_C_F02(FrCor,NMax,N,LLimit,ULimit,Cov,TolCov)
	  CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)

	  IF (QFreq) THEN
	    CALL EC_G_ShwFrq(OutF,QName,FrCor,NMax,N)
	    WRITE(OutF,*)
	    WRITE(OutF,*) 'After frequency-correction : '
	    WRITE(OutF,*) 'Factors accepted between ',LLimit,
     &	                  ' and ',ULimit
	    WRITE(OutF,*)
	    CALL EC_G_Show(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
	  ENDIF

        ELSE
	    WRITE(*,*) 'ERROR: can not perform freq. response correction',
     &                 ' without a temperature'
	ENDIF
      ENDIF
C
C
C Calculate mean vertical velocity according to Webb
C
C
      IF (DoWebb) THEN
         IF (WhichTemp .GT. 0) THEN
	   CALL EC_C_Webb(Mean,NMax,Cov,P, WhichTemp)
	   CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)

	   IF (QWebb) THEN
	     WRITE(OutF,*)
	     WRITE(OutF,*) 'Webb-velocity (vertical) = ',Mean(W),' m/s'
	     WRITE(OutF,*)
	     WRITE(OutF,*) 'After addition of Webb-term : '
	     WRITE(OutF,*)
	     CALL EC_G_Show(OutF,QName,UName,Mean,TolMean,Cov,TolCov,
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
C
C Contribution of other gases especially oxygen absorb
C some of the radiation at the wavelengths at which the
C Krypton works.
C
C Revision 28-05-2001: added info on which temperature should be used
C                      in corrections (Sonic or thermocouple)
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER NMax,N,i, WhichTemp
      REAL*8 Mean(NMax),Cov(NMax,NMax),Factor(NMax),OXC,P, HygType,
     +       GenKo, GenKw, MeanT

      IF (HygType .EQ. ApCampKrypton) THEN
          GenKo = KoK
          GenKw = KwK
      ELSE IF (HygType .EQ. ApMierijLyma) THEN
          GenKo = KoLa
          GenKw = KwLa
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
C
C Contribution of other gases especially oxygen absorb
C some of the radiation at the wavelengths at which the
C Krypton works.
C
C Revision 28-05-2001: added info on which temperature should be used
C                      in corrections (Sonic or thermocouple)
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER NMax,N,i, WhichTemp
      REAL*8 Cov(NMax,NMax),Factor(NMax)
      DO i = 1,N
	Cov(i,Humidity) = Factor(i)*Cov(i,Humidity)
	Cov(i,SpecHum ) = Factor(i)*Cov(i,SpecHum )
	Cov(Humidity,i) = Cov(i,Humidity)
	Cov(SpecHum ,i) = Cov(i,SpecHum )
      ENDDO

      RETURN
      END




C...........................................................................
C Function  : EC_C_SCal
C Purpose   : to calibrate a sonic signal according to wind tunnel
C             calibration (for th moment this works for apparatus 6,
C             i.e. a wind tunnel calibrated sonic)
C Interface : Cal     IN        array of length NQQ with calibation info
C             UDum    IN/OUT    one horizontal component (on exit: calibrated)
C             VDum    IN/OUT    another horizontal component (on exit:
calibrated)
C             WDum    IN/OUT    vertical component (on exit: calibrated)
C             UERROR  IN/OUT    error flag for U (.TRUE. if wrong data)
C             VERROR  IN/OUT    error flag for V (.TRUE. if wrong data)
C             WERROR  IN/OUT    error flag for W (.TRUE. if wrong data)
C Author    : Arnold Moene
C Date      : September 26, 2000
C Remarks   : The method is based on a wind tunnel calibration of the sonic
C             The real velocity components can be derived from the
C             measured components and the real azimuth and elevation angle.
C             But the latter are not known and have to be determined
C             iteratively from the measured components. The relationship
C             between the real components and the measured components is:
C
C               Ureal =  Umeas/(UC1*(1 - 0.5*
C                                    ((Azi + (Elev/0.5236)*UC2)*
C                                     (1 - UC3*Abs(Elev/0.5236)))**2 ))
C               Vreal =  Vmeas*(1 - VC1*Abs(Elev/0.5236))
C               Wreal =  Wmeas/(WC1*(1 - 0.5*(Azi*WC2)**2))
C
C             and
C               Azi = arctan(V/U)
C               Elev = arctan(W/sqrt(U**2 + V**2))
C
C             where UC1, UC2, UC3, VC1, WC1, WC2 are fitting coefficients.
C             An azimuth angle of zero is supposed to refer to a wind
C             direction from the most optimal direction (i.e. the 'open'
C             side of a sonic). Samples with an absolute azimuth angle of
C             more than 40 degrees are rejected.
C...........................................................................
      SUBROUTINE EC_C_Scal(Cal, UDum, VDum, WDum,
     &                  UError, VError, WError)

      INCLUDE 'parcnst.inc'

      REAL*8   Cal(NQQ), UDum, VDum, WDum
      LOGICAL  UError, VError, WError

      REAL*8   UCorr, VCorr, WCorr, AziNew, AziOld, ElevNew, ElevOld,
     &         UC1, UC2, UC3, VC1, WC1, WC2
      INTEGER  I, NITER, ITMAX


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
            UCorr =  UDum/(UC1*(1 - 0.5*
     &              ((AziNew + (ElevNew/0.5236D0)*UC2)*
     &              (1 - UC3*Abs(ElevNew/0.5236D0))
     &              )**2        )
     &                     )
            VCorr =  VDum*(1 - VC1*Abs(ElevNew/0.5236D0))
            WCorr =  WDum/(WC1*(1 - 0.5*(ElevNew*WC2)**2))

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
C
C Partial Schotanus et al. correction: correction for humidity
C of sonic temperature, and of all covariances with sonic temperature.
C Sidewind-correction has already been applied in the
C routine where the sonic signal is calibrated.
C
      INCLUDE 'physcnst.inc'
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
C
C Partial Schotanus et al. correction: correction for humidity
C of sonic temperature, and of all covariances with sonic temperature.
C Sidewind-correction has already been applied in the
C routine where the sonic signal is calibrated.
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER NMax,N,i
      REAL*8 Mean(NMax),Cov(NMax,NMax),Factor(NMax), MeanTSon, TSonFact

      DO i = 1,N
	Cov(i,TSonic) = Factor(i)*Cov(i,TSonic)
	Cov(TSonic,i) = Cov(i,TSonic)
      ENDDO

      MeanTSon = MeanTSon*TSonFact

      RETURN
      END




      REAL*8 FUNCTION EC_C_Schot3(Temp, Rhov, Press)
C
C To do humidity part of Schotanus et al. correction on a
C raw sample of sonic temperature, with Rhov that was not yet
C corrected
C
C Sidewind-correction has already been applied in the
C routine where the sonic signal is calibrated.
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER NMax,N,i
      REAL*8   Temp, Rhov, Press, SPHUM, SPHUMNEW, NEWTEMP
      REAL*8 EC_Ph_Q,EC_M_BaseF ! External function calls

      SPHUM = EC_Ph_Q(Temp, Rhov, Press)
      SPHUMNEW = 1.0
      NEWTEMP = Temp
C Dit convergeert dus helemaal niet !!!
      DO WHILE (ABS((SPHUM - SPHUMNEW)/SPHUM) .GT. 0.0001)

         NEWTEMP = TEMP/(1+0.51*SPHUM)
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



      SUBROUTINE EC_C_T01(TrustWilczak,SingleRun,uMean,NRuns,Apf,
     &  Alpha,Beta,Gamma,WBias)

C
C Subroutine performs tilt correction of Sonic data, using Planar Fit
C Method, as described in James M. Wilczak et al (2001), 'Sonic
C Anemometer tilt correction algorithms', Boundary Meteorology 99: 127:150
C References to formulae are to this article
C The planar fit matrix is extended with an additional yaw-correction
C to turn the first coordinate into the direction of the mean wind
C over all runs. This extra rotation makes results from different
C eddy-covariance systems comparable.
C
C Input: uMean : NRuns x 3 matrix of run mean velocity vectors, where
C        NRuns : the number of runs
C Output: Apf  : The planar fit 3*3 untilt-matrix
C        Alpha : tiltangle alpha in degrees
C        Beta  : tiltangle beta in degrees
C        Gamma : Fixed yaw-angle in degrees associated with mean over all runs
C        WBias : The bias in the vertical velocity
C
      LOGICAL SingleRun,TrustWilczak
      INTEGER i,j,k,NRuns,UmeanMax
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
      CALL EC_C_T02(TrustWilczak,SingleRun,uSum,UUSum,Apf,
     &  Alpha,Beta,Gamma,WBias)


      END






      SUBROUTINE EC_C_T02(TrustWilczak,SingleRun,uSum,UUSum,Apf,
     &  Alpha,Beta,Gamma,WBias)
C
C Supportive routine for planar fit method for tilt-correction
C
      LOGICAL SingleRun,TrustWilczak
      REAL*8 b(0:2),S(3,3),SInv(3,3),x(3),Apf(3,3),SS2(2,2),SS2Inv(2,2),
     &  Sqrt1,Sqrt2,Alpha,Beta,SinAlpha,SinBeta,CosAlpha,CosBeta,Pi,
     &  Gamma,WBias,SinGamma,CosGamma,UHor,Yaw(3,3),USum(3),UUSum(3,3),
     &  u,v,w,R11,R12,R13,R22,R23,R33,A,Disc

      Pi = 4.D0*ATAN(1.D0)


      IF (.NOT.SingleRun) THEN
        IF (TrustWilczak) THEN
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
     #**2*R13**2*u**2+8*u**2*R23*v*R11*w-4*u**2*R23*v*R33*w-4*w*R22*v*u*
     #*2*R23-4*v**3*R11*u*R12+4*v**2*R12**2*w**2+4*v**3*R23*R11*w-4*v**3
     #*R23*R33*w-2*v**2*R11*w**2*R22-8*u**3*R23*w*R12+2*u**2*R33*w**2*R1
     #1-2*u**2*R22*w**2*R11+4*v*R12*u**3*R33-4*w**3*R22*v*R23+4*v*R11*w*
     #*3*R23+8*u*R22*v**2*R13*w+4*u**3*R22*w*R13-2*v**2*R33*u**2*R22-2*u
     #**2*R33*w**2*R22-4*v**2*R33*u*w*R13-8*v**3*R12*R13*w+4*v**2*R12**2
     #*u**2-4*v*R12*u*w**2*R11+2*v**2*R11**2*w**2-8*v**3*R23*u*R13+2*v**
     #2*R33**2*u**2+
     &4*v**4*R13**2+4*u**4*R23**2+u**4*R22**2+u**4*R33**2+4*w**3*R
     #22*u*R13-4*w**3*R11*u*R13-8*w**3*R12*v*R13+2*v**2*R33*w**2*R22-2*v
     #**2*R33*w**2*R11-4*u**3*R33*w*R13+8*u*R33*w**2*R12*v-4*v*R12*u**3*
     #R22-4*v*R12*u*w**2*R22-8*u*R23*w**3*R12-8*u**3*R23*v*R13-4*v**2*R1
     #1*u*w*R13+4*v**3*R12*u*R33+v**4*R33**2+2*v**2*R11*u**2*R22-2*v**2*
     #R11*u**2*R33+4*u**2*R23**2*w**2+4*w**2*R12**2*u**2-2*u**4*R22*
     #R33+w**4*R22**2+4*v**2*R23**2*u**2+2*u**2*R22**2*w**2-2*v**4*R11*R
     #33-2*w**4*R22*R11+4*v**2*R23**2*w**2+4*v**2*R13**2*w**2

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

      UHor = ((USum(1))**2+(USum(2))**2)**0.5
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
C
C Help routine for routines for correction of coordinate system
C
      INCLUDE 'physcnst.inc'
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
C
C Help routine for routines for correction of coordinate system
C
      INCLUDE 'physcnst.inc'
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
C
C Routine to change coordinate system according to tensor "Map".
C This routine is called by all tilt-correction procedures.
C Both the mean velocity and the Reynoldsstresses and the
C covariances of all velocity components with other quantities
C are rotated.
C
      INCLUDE 'physcnst.inc'
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
C
C Construct rotation matrix for coordinate system about a KNOWN yaw-angle
C around the vertical
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

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
C
C Give yaw-angle to transform coordinate system such that v_mean = 0
C (no mean lateral horizontal velocity component)
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      REAL*8 UHor,SinPhi,CosPhi,MeanU,MeanV,Direction

      UHor = (MeanU**2+MeanV**2)**0.5
      SinPhi = MeanV/UHor
      CosPhi = MeanU/UHor
      Direction = 180.D0*ACOS(CosPhi)/PI
      IF (SinPhi.LT.0.D0) Direction = 360.D0-Direction

      RETURN
      END





      SUBROUTINE EC_C_T08(Direction,Pitch)
C
C Give matrix for rotation about a KNOWN pitch-angle around vector (0,1,0).
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

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
C
C Give pitch angle to transform coordinate system such that w_mean = 0
C (no mean vertical velocity component)
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      REAL*8 SinTheta,MeanU,MeanW,Direction,UTot

      UTot = (MeanU**2+MeanW**2)**0.5
      SinTheta = MeanW/UTot
      Direction = 180.D0*ASIN(SinTheta)/PI

      RETURN
      END





      SUBROUTINE EC_C_T10(Direction,Roll)
C
C Give matrix to rotate coordinate system about a KNOWN roll-angle
C around vector (1,0,0).
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

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
C
C Give roll angle to transform coordinate system such that Cov(V,W) = 0
C (vertical velocity fluctuations are independent from horizontal
C fluctuations)
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

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




      SUBROUTINE EC_C_Webb(Mean,NMax,Cov,P,WhichTemp)
C
C Mean vertical velocity according to Webb, Pearman and Leuning
C
C Revision 28-05-2001: added info on which temperature should be used
C                      in corrections (Sonic or thermocouple)
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER NMax, WhichTemp
      REAL*8 Mean(NMax),Cov(NMax,NMax),Sigma,EC_Ph_RhoDry,P
C
C Ratio of mean densities of water vapour and dry air
C
      Sigma = Mean(Humidity)/EC_Ph_RhoDry(Mean(Humidity),Mean(WhichTemp),P)
      Mean(W) = (1.D0+Mu*Sigma)*Cov(W,WhichTemp)/Mean(WhichTemp) +
     &	Mu*Sigma*Cov(W,Humidity)/Mean(Humidity)

      RETURN
      END
