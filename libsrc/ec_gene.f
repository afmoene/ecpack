C  ec_gene.f, part of ecpack
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




C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C
C Main routine that calls calibrations, all corrections and estimates fluxes
C
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################





      SUBROUTINE EC_G_Main(OutF,DoPrint,
     &	RawSampl,MaxChan,Channels,NMax,N,MMax,M,DoCrMean,PCal,PIndep,
     &	Psychro,Freq,CalSonic,CalTherm,CalHyg,P,
     &	Calibr,
     &	Sample,Flag,Mok,Cok,MIndep,CIndep,Rc,BadTc,
     &	QName,UName,
     &	DoDetren,PDetrend,
     &	DoTilt,PTilt,PreYaw,PrePitch,PreRoll,
     &	DoYaw,PYaw,DirYaw,
     &	DoPitch,PPitch,PitchLim,DirPitch,
     &	DoRoll,PRoll,RollLim,DirRoll,
     &	DoSonic,PSonic,SonFactr,
     &	DoO2,PO2,O2Factor,
     &	DoFreq,PFreq,LLimit,ULimit,FrCor,
     &	DoWebb,PWebb,
     &	Mean,TolMean,Cov,TolCov,
     &	HSonic,dHSonic,HTc,dHTc,
     &	LvE,dLvE,LvEWebb,dLvEWebb,
     &	UStar,dUStar,Tau,dTau,
     &  MEANW, TOLMEANW, HAVE_UNCAL)
C
C Purpose : Integrated routine which:
C  - Calibrates raw samples
C  - Estimates mean values and covariances with respective tolerances
C  - Corrects the resulting mean values and covariances for all
C    effects selected by the user
C  - Estimates, from the final mean values and covariances, the surface-
C    fluxes with tolerances.
C
C Revision: 03-04-2001: get mean W and its tolerance before tilt
C                       correction; export via interface (AM)
C Revision: 28-05-2001: added passing of info on whether uncalibrated
C                       data are available (AM)

      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'
      INCLUDE 'calcomm.inc'

      INTEGER NMax,N,i,M,j,MIndep(NMax),CIndep(NMax,NMax),
     &	OutF,MMax,MaxChan,Channels,Mok(NMax),Cok(NMax,NMax),NTcOut,
     & WhichTemp
      LOGICAL PYaw,PPitch,PRoll,PFreq,PO2,PWebb,PCal,PIndep,
     &	DoYaw,DoPitch,DoRoll,DoFreq,DoO2,DoWebb,DoCrMean,PSonic,
     &	DoSonic,DoTilt,PTilt,DoPrint,Flag(NMAx,MMax),DoDetren,
     &	PDetrend,AnyTilt,DumYaw,DumPitch,DumRoll,DumTilt,BadTc
      REAL*8 RawSampl(MaxChan,MMax)
      REAL*8 P,Psychro,Sample(NMax,MMax),
     &	Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),LLimit,ULimit,FrCor(NMax,NMax),
     &	PitchLim,RollLim,DirYaw,RC(NMax),O2Factor(NMax),
     &	Freq,DirPitch,DirRoll,
     &	SonFactr(NMax),CorMean(NNMax),PreYaw,PrePitch,PreRoll,
     &	HSonic,dHSonic,HTc,dHTc,
     &	LvE,dLvE,LvEWebb,dLvEWebb,
     &	UStar,dUStar,Tau,dTau,Speed(3),DumCov(3,3)
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ)
      REAL*8 MEANW, TOLMEANW, MINS(NNMax), MAXS(NNMAX)
      REAL*8 EC_Ph_Q,Yaw(3,3),Pitch(3,3),Roll(3,3)
      CHARACTER*6 QName(NMax)
      CHARACTER*9 UName(NMax)
      LOGICAL HAVE_CAL(NMax)
      LOGICAL HAVE_UNCAL(NMax)
C
C Check whether we have the needed calibration info
C
      DO I=1,NMax
         HAVE_CAL(i) = .FALSE.
      ENDDO
      IF (HAVE_SONCAL) THEN
         HAVE_CAL(U) = .TRUE.
         HAVE_CAL(V) = .TRUE.
         HAVE_CAL(W) = .TRUE.
         HAVE_CAL(TSonic) = .TRUE.
	 WhichTemp = TSonic
      ENDIF
      IF (HAVE_HYGCAL) THEN
         HAVE_CAL(Humidity) = .TRUE.
         HAVE_CAL(SpecHum) = .TRUE.
      ENDIF
      IF (HAVE_TCCAL) THEN
         HAVE_CAL(TCouple) = .TRUE.
	 WhichTemp = TCouple
      ENDIF

C
C Check whether we have data needed for corrections
C
      IF ((DoSonic) .AND. .NOT. (HAVE_CAL(SpecHum))) THEN
         WRITE(*,*) 'WARNING: Do not have a hygrometer: ',
     &                 'can not do ',
     &                 'humidity correction on sonic temperature'
         DoSonic = .FALSE.
      ENDIF
C
C Names of quantities and units
C
      QName(U	    ) = '     U'
      QName(V	    ) = '     V'
      QName(W	    ) = '     W'
      QName(TSonic  ) = 'T(son)'
      QName(TCouple ) = 'T(cpl)'
      QName(Humidity) = '  RhoV'
      QName(SpecHum ) = '     q'
      QName(TTime   ) = '  Time'

      UName(U	    ) = '[m/s]'
      UName(V	    ) = '[m/s]'
      UName(W	    ) = '[m/s]'
      UName(TSonic  ) = '[K]'
      UName(TCouple ) = '[K]'
      UName(Humidity) = '[kg m^-3]'
      UName(SpecHum ) = '[kg/kg]'
      UName(TTime   ) = '[s]'
C
C
C Calibrate the raw samples for the first time ignoring drift in apparatus.
C
C
      DO i=1,N
       	CorMean(i) = 0.D0
      ENDDO

      NTcOut = 0
      BadTc = (.FALSE.)
      Have_Uncal(TTime) = .TRUE.
      IF (Have_Uncal(Humidity)) Have_Uncal(SpecHum) = .TRUE.

      DO i=1,M
         CALL Calibr(RawSampl(1,i),Channels,P,CorMean,
     &	  CalSonic,CalTherm,CalHyg,BadTc,Sample(1,i),N,Flag(1,i),
     &    Have_Uncal)
        IF (Flag(TCouple,i)) NTcOut = NTcOut + 1
      ENDDO
      BadTc = (NTcOut.GE.(M/2))
      IF (BadTc) THEN
        DO i=1,M
          CALL Calibr(RawSampl(1,i),Channels,P,CorMean,
     &      CalSonic,CalTherm,CalHyg,BadTc,Sample(1,i),N,Flag(1,i),
     &      Have_Uncal)
        ENDDO
      ENDIF
C
C
C Calibrate the raw samples for the second time, now correcting mean
C quantities of drift-sensitive apparatus using slow signals.
C
C
      IF (DoCrMean) THEN
C
C Find the shift/drift in humidity of the krypton hygrometer
C
	CALL EC_M_Averag(Sample,NMax,N,MMax,M,Flag,
     &	  Mean,TolMean,Cov,TolCov,MIndep,CIndep,Mok,Cok)
        CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)
        IF (.NOT. Have_Uncal(Humidity)) THEN
	   Mean(Humidity) = Psychro
	   MEAN(SpecHum) = EC_Ph_Q(Mean(Humidity), Mean(WhichTemp), P)
	ENDIF
        CALL EC_M_MinMax(Sample, NMax, N, MMax, M, Flag, MINS, MAXS)
	IF (DoPrint) THEN
	   WRITE(OUTF,*) 'Min/max of samples'
	   DO I=1,N
	      WRITE(OUTF,*) QName(I),': min = ', Mins(I), ', max = ',
     &                   Maxs(I)
	   ENDDO
	ENDIF

	CorMean(Humidity) = Psychro - Mean(Humidity)
	IF (DoPrint.AND.PCal) THEN
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'Added to humidity : ',
     &	  (Psychro - Mean(Humidity)),' [kg m^{-3}]'
	  WRITE(OutF,*)
        ENDIF

	DO i=1,M
	  CALL Calibr(RawSampl(1,i),Channels,P,CorMean,
     &	    CalSonic,CalTherm,CalHyg,BadTc,Sample(1,i),N,Flag(1,i),
     &      Have_Uncal)
	ENDDO
      ENDIF

      IF (DoPrint.AND.PCal) THEN
	DO i=1,M
	  DO j=1,N
	    IF (Flag(j,i) .AND. HAVE_CAL(j) .AND. HAVE_UNCAL(J))
     &       WRITE(OutF,*) 'Error in sample ',i,' = ',
     &	      QName(j),' = ',Sample(j,i)
	  ENDDO
	ENDDO
      ENDIF

C
C
C Estimate mean values, covariances and tolerances of both
C
C
      CALL EC_M_Averag(Sample,NMax,N,MMax,M,Flag,
     &	Mean,TolMean,Cov,TolCov,MIndep,CIndep,Mok,Cok)
      CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)
      IF (.NOT. Have_Uncal(Humidity)) THEN
        Mean(Humidity) = Psychro
        MEAN(SpecHum) = EC_Ph_Q(Mean(Humidity), Mean(WhichTemp), P)
      ENDIF
      CALL EC_M_MinMax(Sample, NMax, N, MMax, M, Flag, MINS, MAXS)
      IF (DoPrint) THEN
	   WRITE(OUTF,*) 'Min/max of samples'
           DO I=1,N
              WRITE(OUTF,*) QName(I),': min = ', Mins(I), ', max = ',
     &                   Maxs(I)
           ENDDO
      ENDIF

C
C Print number of independent contributions
C
      IF (DoPrint.AND.PIndep) THEN
	CALL EC_G_ShwInd(OutF,QName,MIndep,CIndep,NMax,N,M,Freq)
      ENDIF
C
C Print averages of raw, calibrated data
C
      IF (DoPrint.AND.PCal) THEN
	WRITE(OutF,*)
	WRITE(OutF,*) 'For raw calibrated data : '
	WRITE(OutF,*)
	CALL EC_G_Show(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
      ENDIF
C
C
C Subtract a linear trend from all quantities
C
C
      IF (DoDetren) THEN

	CALL EC_M_Detren(Sample,NMax,N,MMAx,M,Mean,Cov,Sample,RC)
C
C
C Estimate mean values, covariances and tolerances of both
C for the detrended dataset
C
C
	CALL EC_M_Averag(Sample,NMax,N,MMax,M,Flag,
     &	  Mean,TolMean,Cov,TolCov,MIndep,CIndep,Mok,Cok)
        CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)
        IF (.NOT. Have_Uncal(Humidity)) THEN
          Mean(Humidity) = Psychro
          MEAN(SpecHum) = EC_Ph_Q(Mean(Humidity), Mean(WhichTemp), P)
        ENDIF
        CALL EC_M_MinMax(Sample, NMax, N, MMax, M, Flag, MINS, MAXS)

	IF (DoPrint.AND.PDetrend) THEN
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'After detrending : '
	  WRITE(OutF,*)

	  WRITE(OUTF,*) 'Min/max of samples'
          DO I=1,N
             WRITE(OUTF,*) QName(I),': min = ', Mins(I), ', max = ',
     &                Maxs(I)
          ENDDO
	  IF (PIndep) THEN
	    CALL EC_G_ShwInd(OutF,QName,MIndep,CIndep,NMax,N,M,Freq)
	  ENDIF

	  CALL EC_G_Show(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
	ENDIF
      ENDIF
C
C Get the mean W before we do tilt correction
C
      MEANW = MEAN(W)
      TOLMEANW = TolMean(W)
C
C Correct mean values and covariances for all thinkable effects
C
      CALL EC_C_Main(OutF,
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
     &	DoWebb,PWebb,P,Have_Uncal)
      CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)
C
C If any transformation of coordinates was required (one of the options
C DoTilt, DoYaw, DoPitch or DoRoll was selected), then the numbers
C of independent samples of the respective quantities has to be
C re-estimated. It is not possible to "tilt the error-bars" on basis of
C the quantities which have been calculated until here.
C Therefore we return to the calibrated time series and
C make the transformations BEFORE averaging. Then averaging and corrections
C are repeated all over again.
C
      AnyTilt = ((DoTilt.OR.DoYaw).OR.(DoPitch.OR.DoRoll))

      IF (AnyTilt) THEN

        DO i=1,3
          DO j=1,3
            DumCov(i,j) = 0.D0
          ENDDO
        ENDDO
C
C Tilt ALL samples
C
        DO i=1,M
          Speed(1) = Sample(U,i)
          Speed(2) = Sample(V,i)
          Speed(3) = Sample(W,i)

          IF (DoTilt) THEN
            CALL EC_C_T06(PreYaw,Yaw)
            CALL EC_C_T05(Speed,3,3,DumCov,Yaw)
            CALL EC_C_T08(PrePitch,Pitch)
            CALL EC_C_T05(Speed,3,3,DumCov,Pitch)
            CALL EC_C_T10(PreRoll,Roll)
            CALL EC_C_T05(Speed,3,3,DumCov,Roll)
          ENDIF

          IF (DoYaw) THEN
            CALL EC_C_T06(DirYaw,Yaw)
            CALL EC_C_T05(Speed,3,3,DumCov,Yaw)
	  ENDIF
          IF (DoPitch) THEN
            CALL EC_C_T08(DirPitch,Pitch)
            CALL EC_C_T05(Speed,3,3,DumCov,Pitch)
	  ENDIF
          IF (DoRoll) THEN
            CALL EC_C_T10(DirRoll,Roll)
            CALL EC_C_T05(Speed,3,3,DumCov,Roll)
	  ENDIF

          Sample(U,i) = Speed(1)
          Sample(V,i) = Speed(2)
          Sample(W,i) = Speed(3)

        ENDDO
C
C Reestablish the averages and covariances in the correct frame of reference.
C
	CALL EC_M_Averag(Sample,NMax,N,MMax,M,Flag,
     &	  Mean,TolMean,Cov,TolCov,MIndep,CIndep,Mok,Cok)
        CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)
        IF (.NOT. Have_Uncal(Humidity)) THEN
          Mean(Humidity) = Psychro
          MEAN(SpecHum) = EC_Ph_Q(Mean(Humidity), Mean(WhichTemp), P)
        ENDIF

	IF (DoPrint) THEN
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'After untilting raw data : '
	  WRITE(OutF,*)

	  IF (PIndep) THEN
	    CALL EC_G_ShwInd(OutF,QName,MIndep,CIndep,NMax,N,M,Freq)
	  ENDIF

	  CALL EC_G_Show(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
	ENDIF

        DumTilt  = (.FALSE.)
        DumYaw = (.FALSE.)
        DumPitch   = (.FALSE.)
        DumRoll  = (.FALSE.)
C
C Perform all necessary corrections on the mean values and (co-)variances.
C
        CALL EC_C_Main(OutF,
     &	  DoPrint,
     &	  Mean,NMax,N,TolMean,
     &	  Cov,TolCov,
     &	  QName,UName,
     &    BadTc,
     &	  DumTilt,PTilt,PreYaw,PrePitch,PreRoll,
     &	  DumYaw,PYaw,DirYaw,
     &	  DumPitch,PPitch,PitchLim,DirPitch,
     &	  DumRoll,PRoll,RollLim,DirRoll,
     &	  DoSonic,PSonic,SonFactr,
     &	  DoO2,PO2,O2Factor,
     &	  DoFreq,PFreq,LLimit,ULimit,Freq,
     &    CalSonic,CalTherm,CalHyg,FrCor,
     &	  DoWebb,PWebb,P, Have_Uncal)
        CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)
      ENDIF
C
C
C Calculate fluxes from processed mean values and covariances.
C
C
      CALL EC_Ph_Flux(Mean,NMax,Cov,TolMean,TolCov,p,BadTc,
     &	HSonic,dHSonic,HTc,dHTc,
     &	LvE,dLvE,LvEWebb,dLvEWebb,
     &	UStar,dUStar,Tau,dTau)

      CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)

      IF ((.NOT. HAVE_Uncal(U)) .OR.
     +    (.NOT. HAVE_Uncal(V)) .OR.
     +    (.NOT. HAVE_Uncal(W))) THEN
         Ustar = DUMMY
         dUstar = DUMMY
      ENDIF
      IF (.NOT. HAVE_Uncal(TSonic)) THEN
	  HSonic = DUMMY
	  dHSonic = DUMMY
      ENDIF
      IF (.NOT. HAVE_Uncal(TCouple)) THEN
	  HTc = DUMMY
	  dHTc = DUMMY
      ENDIF
      IF ((.NOT. HAVE_Uncal(TSonic)) .AND.
     +    (.NOT. HAVE_Uncal(Tcouple))) THEN
        LvEWebb  = DUMMY
        dLvEWebb  = DUMMY
      ENDIF
      IF (.NOT. HAVE_Uncal(Humidity)) THEN
	  Lve = DUMMY
	  dLvE = DUMMY
	  LveWebb = DUMMY
	  dLvEWebb = DUMMY
      ENDIF

      RETURN
      END
      
C...........................................................................
C Routine   : EC_G_CLEARSTR
C Purpose   : To set all characters in a string to CHAR(0)
C Interface : STRING      intent(INOUT)   string to be changed
C Author    : Arnold Moene
C Date      : May 25, 2000
C...........................................................................
      SUBROUTINE EC_G_CLEARSTR(STRING)

      IMPLICIT NONE
      CHARACTER*(*) STRING
      INTEGER       I

      DO 15, I=1,LEN(STRING)
         STRING(I:I) = CHAR(0)
   15 CONTINUE
      END




C.........................................................................
C Routine to reset means and covariances
C.........................................................................
      SUBROUTINE EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)
      INCLUDE 'parcnst.inc'

      LOGICAL Have_Uncal(NNMax)
      REAL*8  Mean(NNmax), TolMean(NNMax), Cov(NNMax, NNMax),
     +        TolCov(NNMax, NNMax)
      INTEGER I,J

      IF ((.NOT. HAVE_Uncal(U)) .OR.
     +    (.NOT. HAVE_Uncal(V)) .OR.
     +    (.NOT. HAVE_Uncal(W))) THEN
         DO I=U,W
	    MEAN(I) = DUMMY
	    TolMean(I) = DUMMY
	    DO J=1,NNMax
	        COV(I,J) = DUMMY
		TolCOV(I,J) = DUMMY
	        COV(J,I) = DUMMY
		TolCOV(J,I) = DUMMY
	    ENDDO
	 ENDDO
      ENDIF
      IF (.NOT. HAVE_Uncal(TSonic)) THEN
          Mean(TSonic) = DUMMY
          TolMean(TSonic) = DUMMY
	  DO I=1,NNMax
	     Cov(Tsonic,I) = DUMMY
	     Cov(I,Tsonic) = DUMMY
	     TolCov(Tsonic,I) = DUMMY
	     TolCov(I,Tsonic) = DUMMY
	  ENDDO
      ENDIF
      IF (.NOT. HAVE_Uncal(TCouple)) THEN
          Mean(TCouple) = DUMMY
          TolMean(TCouple) = DUMMY
	  DO I=1,NNMax
	     Cov(TCouple,I) = DUMMY
	     Cov(I,TCouple) = DUMMY
	     TolCov(TCouple,I) = DUMMY
	     TolCov(I,TCouple) = DUMMY
	  ENDDO
      ENDIF
      IF (.NOT. HAVE_Uncal(Humidity)) THEN
	  DO I=1,NNMax
	    DO J=Humidity, SpecHum
	      Cov(J,I) = DUMMY
	      Cov(I,J) = DUMMY
	      TolCov(J,I) = DUMMY
	      TolCov(I,J) = DUMMY
	    ENDDO
	  ENDDO
      ENDIF
      END


C
C ########################################################################
C
C Some silly routines to print intermediate results to file. Handy
C for checking the origin of unexpected results!
C
C ########################################################################
C





      SUBROUTINE EC_G_Show(OutF,QName,UName,Mean,TolMean,Cov,
     &	TolCov,NMax,N)
C
C Prints mean values and (co-)variances plus respective tolerances
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER i,j,N,NMax,OutF
      REAL*8 Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),Dum(NNMax)
      CHARACTER*6 QName(NMax)
      CHARACTER*9 UName(NMax)

      WRITE(OutF,*)
      WRITE(OutF,*) 'Mean quantities : '
      WRITE(OutF,*)

      DO i=1,N
	WRITE(OutF,10) QName(i),Mean(i),TolMean(i),UName(i)
      ENDDO
 10   FORMAT(a6,' = ',F20.10,' +/- ',F20.10,1X,a9)

      WRITE(OutF,*)
      WRITE(OutF,*) 'Covariances : '
      WRITE(OutF,*) '(2 lines : covariance, tolerance)'
      WRITE(OutF,*)

      WRITE(OutF,30) '      ',(QName(i),i=1,N)
 30   FORMAT(a6,20(a11:,1X))

      DO i=1,N
	WRITE(OutF,20) QName(i),(Cov(i,j),j=1,N)
	WRITE(OutF,20) '      ',(TolCov(i,j),j=1,N)
 20	FORMAT(a6,20(G11.4:,1X))
      ENDDO

      WRITE(OutF,*)
      WRITE(OutF,*) 'Correlation-coefficients : '
      WRITE(OutF,*)

      WRITE(OutF,40) '      ',(QName(i),i=1,N)
 40   FORMAT(a6,20(a10:,1X))

      DO i=1,N
	DO j=1,N
	  IF ((ABS(Cov(i,i)*Cov(j,j)).GT.Epsilon) .AND.
     +        (ABS(Cov(i,i)-DUMMY) .GT. Epsilon) .AND.
     +        (ABS(Cov(j,j)-DUMMY) .GT. Epsilon) .AND.
     +        (ABS(Cov(i,j)-DUMMY) .GT. Epsilon)) THEN
	    Dum(j) = Cov(i,j)/SQRT(Cov(i,i)*Cov(j,j))
	  ELSE
	    Dum(j) = DUMMY
	  ENDIF
	ENDDO
	WRITE(OutF,50) QName(i),(Dum(j),j=1,N)
 50	FORMAT(a6,20(F10.4:,1X))
      ENDDO

      WRITE(OutF,*)
      WRITE(OutF,*) '###############################################'
      WRITE(OutF,*)

      RETURN
      END





      SUBROUTINE EC_G_ShwFrq(OutF,QName,FrCor,NMax,N)
C
C Prints correction factors associated with frequency response
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER i,j,N,NMax,OutF
      REAL*8 FrCor(NMax,NMax)
      CHARACTER*6 QName(NMax)

      WRITE(OutF,*)
      WRITE(OutF,*) 'Frequency respons correction factors'
      WRITE(OutF,*) 'for covariances'
      WRITE(OutF,*)

      WRITE(OutF,30) '  ',(QName(i),i=1,N)
 30   FORMAT(a6,20(1X,a8:))

      DO i=1,N
	IF (i.NE.TTime) WRITE(OutF,20) QName(i),(FrCor(i,j),j=1,N)
 20	FORMAT(a6,20(1X,F8.5:))
      ENDDO

      RETURN
      END





      SUBROUTINE EC_G_ShwInd(OutF,QName,MIndep,
     &	CIndep,NMax,N,M,Freq)
C
C Prints number of independent observations
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER i,j,N,NMax,OutF,MIndep(NMax),CIndep(NMax,NMax),M
      REAL*8 Freq
      CHARACTER*6 QName(NMax)

      WRITE(OutF,*)
      WRITE(OutF,*) ' Number of independent samples: '
      WRITE(OutF,*)

      WRITE(OutF,*) 'For mean quantities : '
      WRITE(OutF,*) '                #    timescale [s]'
      WRITE(OutF,*)

      DO i=1,N
	WRITE(OutF,10) QName(i),MIndep(i),
     &    (DBLE(M)/(Freq*DBLE(MIndep(i))))
      ENDDO
 10   FORMAT(a6,' = ',I8,5X,F12.3)

      WRITE(OutF,*)
      WRITE(OutF,*) 'For covariances : '
      WRITE(OutF,*) '  #'
      WRITE(OutF,*)

      WRITE(OutF,30) '      ',(QName(i),i=1,N)
 30   FORMAT(a6,20(a6:,1X))

      DO i=1,N
	WRITE(OutF,20) QName(i),(CIndep(i,j),j=1,N)
 20	FORMAT(a6,20(I6:,1X))
      ENDDO

      WRITE(OutF,*)
      WRITE(OutF,*) 'timescale [s]'
      WRITE(OutF,*)

      WRITE(OutF,40) '      ',(QName(i),i=1,N)
 40   FORMAT(a6,20(a9:,1X))

      DO i=1,N
	WRITE(OutF,50) QName(i),
     &   (DBLE(M)/(Freq*DBLE(CIndep(i,j))),j=1,N)
 50	FORMAT(a6,20(F9.3:,1X))
      ENDDO

      RETURN
      END
C...........................................................................
C routine   : strcat
C Purpose   : to concatenate a string to another
C Interface : STRING1     intent(INOUT)   string to be added to (supposed
C                                         to be of sufficient length)
C             STRING2     intent(IN)      string to add
C Author    : Arnold Moene
C Date      : May 24, 2000
C...........................................................................
      SUBROUTINE EC_G_STRCAT(STRING1, STRING2)

      IMPLICIT NONE
      CHARACTER*(*) STRING1, STRING2
      INTEGER       I, LEN1
      INTEGER       EC_G_STRLEN
      EXTERNAL      EC_G_STRLEN


      LEN1 = EC_G_STRLEN(STRING1)
      DO 15, I = 1, EC_G_STRLEN(STRING2)
        STRING1(LEN1+I:LEN1+I) = STRING2(I:I)
   15 CONTINUE
      END

C...........................................................................
C Routine   : stripstr
C Purpose   : to replace all trailing spaces in a string by NULL characters
C             and strip all leading spaces. Furthermore, remove
C             opening and closing quotes
C Interface : STRING      intent(INOUT)   string to be checked
C Author    : Arnold Moene
C Date      : May 24, 2000
C...........................................................................
      SUBROUTINE EC_G_STRIPSTR(STRING)

      IMPLICIT NONE
      CHARACTER*(*) STRING
      INTEGER       I, J, K, CHANGED, EC_G_STRLEN
      EXTERNAL EC_G_STRLEN

      CHARACTER     QUOTE(2)

C Remove trailing spaces
      DO 15, I = LEN(STRING), 1, -1
        IF (STRING(I:I) .NE. ' ') GO TO 20
        IF (STRING(I:I) .EQ. ' ') STRING(I:I) = CHAR(0)
   15 CONTINUE
   20 CONTINUE
C Remove leading spaces
      CHANGED = 1      
      I = 1
      DO 30, WHILE (CHANGED .EQ. 1)
        CHANGED = 0
        IF ((STRING(I:I) .EQ. ' ') .OR.
     &      (STRING(I:I) .EQ. CHAR(0))) THEN
           DO 35, J = I+1,LEN(STRING)-1
              STRING(J-1:J-1) = STRING(J:J)
   35      CONTINUE
           CHANGED = 1
        ELSE
           I = I + 1
        ENDIF
   30 CONTINUE
C Remove outer quotes
      QUOTE(1) = "'"
      QUOTE(2) = '"'
      DO 75, K = 1,2
         IF (STRING(1:1) .EQ. QUOTE(K)) THEN
              DO 55, I = LEN(STRING),1,-1
                 IF (STRING(I:I) .EQ. QUOTE(K)) THEN
                    STRING(I:I) = CHAR(0)
                    DO 65, J = I+1,LEN(STRING)-1
                       STRING(J-1:J-1) = STRING(J:J)
   65               CONTINUE
                 ENDIF
   55         CONTINUE
         ENDIF
   75 CONTINUE
      END


C...........................................................................
C Function  : strlen
C Purpose   : to return the length of a string
C Interface : STRING      intent(IN)      string to be checked
C             return value                length of string (integer)
C Author    : after function in book of Clive Page (Professional programmer's
C             guide to Fortran 77)
C Date      : May 24, 2000
C...........................................................................
      INTEGER FUNCTION EC_G_STRLEN(STRING)

      IMPLICIT NONE
      CHARACTER*(*) STRING
      INTEGER       I

      DO 15, I = LEN(STRING), 1, -1
        IF ((STRING(I:I) .NE. ' ') .AND.
     &      (STRING(I:I) .NE. CHAR(0)))
     &     GO TO 20
   15 CONTINUE
   20 EC_G_STRLEN = I
      END

C...........................................................................
C Routine   : EC_G_UPCASE
C Purpose   : To convert a string to upper case
C Interface : STRING      intent(INOUT)   string to be changed
C Author    : Arnold Moene
C Date      : May 24, 2000
C...........................................................................
      SUBROUTINE EC_G_UPCASE(STRING)

      IMPLICIT NONE
      CHARACTER*(*) STRING
      INTEGER       I

      DO 15, I=1,LEN(STRING)
         IF ((ICHAR(STRING(I:I)) .GT. 96) .AND.
     &       (ICHAR(STRING(I:I)) .LT. 124))
     &       STRING(I:I) = CHAR(ICHAR(STRING(I:I))-32)
   15 CONTINUE
      END
