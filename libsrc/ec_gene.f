C  ec_gene.f, part of ecpack
C  
C  Copyright (C) 
C    1998-2000   Arjan van Dijk, Wageningen University, The Netherlands
C    2000-2002   Arjan van Dijk, IMAU, The Netherlands
C    1999-2004   Arnold Moene, Wageningen University, The Netherlands
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
     &	RawSampl,MaxChan,Channels,NMax,N,MMax,M,PCal,PIndep,
     &	Psychro,CalSonic,CalTherm,CalHyg, CalCO2, P,
     &	Calibr,
     &	Sample,Flag,Mok,Cok,MIndep,CIndep,Rc,BadTc,
     &  DoCorr, PCorr, ExpVar, DirYaw, DirPitch, DirRoll,
     &  Apf, SonFactr, O2Factor, FrCor,
     &	Mean,TolMean,Cov,TolCov,
     &  QPhys, dQPhys,
     &  HAVE_UNCAL, HAVE_CAL, DiagFlag, FirstDay)
C     ****f* ec_gene.f/EC_G_Main
C NAME
C     EC_G_Main
C SYNOPSIS
C     CALL EC_G_Main(OutF,DoPrint,
C        RawSampl,MaxChan,Channels,NMax,N,MMax,M, PCal,PIndep,
C        Psychro,CalSonic,CalTherm,CalHyg, CalCO2, P,
C        Calibr,
C        Sample,Flag,Mok,Cok,MIndep,CIndep,Rc,BadTc,
C        DoCorr, PCorr, ExpVar,
C        DirYaw, DirPitch, DirRoll,
C        Apf, SonFactr, O2Factor, FrCor,
C        Mean,TolMean,Cov,TolCov,
C        QPhys, dQPhys,
C        HAVE_UNCAL, Have_cal, DiagFlag, FirstDay)
C FUNCTION
C     Integrated routine which:
C     - Calibrates raw samples
C     - Estimates mean values and covariances with respective tolerances
C     - Corrects the resulting mean values and covariances for all
C       effects selected by the user
C     - Estimates, from the final mean values and covariances, the surface-
C       fluxes with tolerances.
C INPUTS
C     OutF      : [INTEGER]
C                 unit number of file for intermediate results
C     DoPrint   : [LOGICAL]
C                 write intermediate results to file?
C     RawSampl  : [REAL*8(MaxChan, MMax)
C                 raw, uncalibrated samples
C     MaxChan   : [INTEGER]
C                 maximum number of channels in RawSampl
C     Channels  : [INTEGER]
C                 actual number of channels in RawSampl
C     NMax      : [INTEGER]
C                 maximum number of calibrated quantities
C     N         : [INTEGER]
C                 actual number of calibrated quantities  
C     MMax      : [INTEGER]
C                 maximum number of samples
C     M         : [INTEGER]
C                 actual number of samples
C     PCal      : [LOGICAL]
C                 print results of slow sensor correction?
C     PIndep    : [LOGICAL]
C                 print number of independent samples?
C     DoCorr    : [LOGICAL](NMaxCorr)
C                 which corrections to do?
C     PCorr     : [LOGICAL](NMaxCorr)
C                 intermediate results of which corrections?
C     ExpVar    : [REAL*8](NMaxExp)
C                 array with experimental settings
C     Psychro   : [REAL*8]
C                 water vapour density of slow sensor (kg/m^3)
C     CalSonic  : [REAL*8(NNQ)]
C                 array with calibration data of sonic
C     CalTherm  : [REAL*8(NNQ)]
C                 array with calibration data of thermocouple
C     CalHyg    : [REAL*8(NNQ)]
C                 array with calibration data of hygrometer
C     CalCO2    : [REAL*8(NNQ)]
C                 array with calibration data of CO2 sensor
C     P         : [REAL*8]
C                 atmospheric pressure (Pa)
C     Calibr    : [SUBROUTINE]
C                 calibration subroutine
C     Apf       : [REAL*8(3,3)]
C                 planar fit untilt matrix
C     HAVE_UNCAL:   [LOGICAL(NMax)]
C                 switch whether data for uncalibrated data
C                 are available for each channel
C     FirstDay  : [INTEGER]
C                 day number of first sample in array (needed for
C                 detrending data that pass midnight)
C OUTPUTS
C     Sample    : [REAL*8(NMax, MMax)]
C                 array with calibrated samples
C     Flag      : [LOGICAL(NMax, MMax)]
C                 validity flag (true means invalid) for all
C                 calibrated samples
C     Mok       : [INTEGER(NMax)]
C                 number of valid samples for each quantity
C     Cok       : [INTEGER(NMax,NMax)]
C                 number of valid samples for each combination of
C                 two quantities
C     MIndep    : [INTEGER(NMax)]
C                 number of independent samples for each quantity
C     CIndep    : [INTEGER(NMax, NMax)]
C                 number of independent samples for each 
C                 combination of two quantities
C     Rc        : [REAL*8(NMax)]
C                 slope of linear regression in case liner detrending
C                 has been done (for each quanity)
C     BadTc     : [LOGICAL]
C                 if more than half of the thermocouple samples are 
C                 wrong, flag thermocouple as bad
C     DirYaw    : [REAL*8]
C                 yaw angle (degrees)
C     DirPitch  : [REAL*8]
C                 pitch angle (degrees)
C     DirRoll   : [REAL*8]
C                 roll angle (degrees)
C     SonFactr  : [REAL*8(NMax)]
C                 correction factor for the covariances with specific
C                 humidity
C     O2Factor  : [REAL*8(NMax)] 
C                 Correction factor due to oxygen correction for
C                 covariance of humidity with each calibrated
C     FrCor     : [REAL*8(NMax, NMax)]
C                 Correction factors for covariances for frequency
C                 response
C     Mean      : [REAL*8(NMax)] (in/out)
C                 Mean values of all calibrated signals
C     TolMean   : [REAL*8(NMax)] (in/out)
C                 Tolerances in mean values of all calibrated signals
C     Cov       : [REAL*8(NMax,NMax)] (in/out)
C                 Covariances of all calibrated signals
C     TolCov    : [REAL*8(NMax,NMax)] (in/out)
C                 Tolerances in covariances of all calibrated signals
C     QPhys     : [REAL*8](NMaxPhys)] (in/out)
C                 array with physical quantities
C     dQPhys    : [REAL*8](NMaxPhys)] (in/out)
C                 array with tolerances physical quantities
C     HAVE_CAL  : [LOGICAL(NNMax)]
C                 switch whether data for calibrated data
C                 are available for each channel
C     DiagFlag  : [INTEGER](NMaxDiag)
C                 count of flags occuring in diagnostic word of CSAT
C AUTHOR
C     Arjan van Dijk, Arnold Moene
C HISTORY
C     Revision: 03-04-2001: get mean W and its tolerance before tilt
C                           correction; export via interface (AM)
C     Revision: 28-05-2001: added passing of info on whether uncalibrated
C                           data are available (AM)
C     Revision: 18-09-2002: removed calcomm.inc and added FirstDay
C                           to interface (to pass it to calibration routine
C     Revision: 5-12-2002:  added DoPF, PPF, Apf to interface to include planar fit 
C     Revision: 13-01-2003: added vectorwind and dirfrom to interface
C     Revision: 27-01-2003: removed physical quantities from interface
C                           and replace by QPhys
C                           Put all correction info into DoCorr and PCorr, ExpVar
C     REvision: 29-01-2003: added have_cal to interface to
C                           have it available in ec_ncdf
C     $Name$
C     $Id$
C USES
C     EC_C_Main
C     EC_C_T05
C     EC_C_T06
C     EC_C_T08
C     EC_C_T10
C     EC_G_Reset
C     EC_G_ShwInd
C     EC_G_Show
C     EC_M_Averag
C     EC_M_MinMax
C     EC_M_Detren
C     EC_Ph_Q
C     EC_Ph_Flux
C     parcnst.inc
C     ***
      IMPLICIT NONE
      INCLUDE 'parcnst.inc'

      INTEGER NMax,N,i,M,j,MIndep(NMax),CIndep(NMax,NMax),
     &	OutF,MMax,MaxChan,Channels,Mok(NMax),Cok(NMax,NMax),NTcOut,
     &  WhichTemp, DiagFlag(NMaxDiag), FirstDay, DIAG_WORD
      LOGICAL PCal,PIndep,
     &	DoPrint,Flag(NMAx,MMax),
     &	AnyTilt,BadTc,
     &  DoCorr(NMaxCorr), PCorr(NMaxCorr),
     &  DumCorr(NMaxCorr), PDumCorr(NMaxCorr)
      REAL*8 RawSampl(MaxChan,MMax), ExpVar(NMaxExp)
      REAL*8 P,Psychro,Sample(NMax,MMax),
     &	Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),FrCor(NMax,NMax),
     &	DirYaw,RC(NMax),O2Factor(NMax),
     &	DirPitch,DirRoll,
     &	SonFactr(NMax),CorMean(NNMax),
     &	Speed(3),DumCov(3,3), QPhys(NMaxPhys), dQPhys(NMaxPhys)
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ), CalCO2(NQQ)
      REAL*8 MINS(NNMax), MAXS(NNMAX), WebVel
      REAL*8 EC_Ph_Q,Yaw(3,3),Pitch(3,3),Roll(3,3), Apf(3,3), 
     &        DumOut(3) 
      LOGICAL HAVE_CAL(NMax)
      LOGICAL HAVE_UNCAL(NMax)
C
C Check whether we have the needed calibration info
C
      DO I=1,NMax
         HAVE_CAL(i) = .FALSE.
      ENDDO
      IF (HAVE_UNCAL(QUU)) HAVE_CAL(U) = .TRUE.
      IF (HAVE_UNCAL(QUV)) HAVE_CAL(V) = .TRUE.
      IF (HAVE_UNCAL(QUW)) HAVE_CAL(W) = .TRUE.
      IF (CalSonic(QQType) .EQ. ApCSATSonic) THEN
        IF (Have_Uncal(QUDiagnost)) THEN
            DiagFlag(QDDelta) = 0
            DiagFlag(QDLock) = 0
            DiagFlag(QDHigh) = 0
            DiagFlag(QDLow) = 0
            DO I=1,M
               DIAG_WORD = INT(RawSampl(QUDiagnost, I)/4096)
               DiagFlag(QDDelta) = 
     &            DiagFlag(QDDelta) + MOD(INT(DIAG_WORD/8),2)
               DiagFlag(QDLock) = 
     &            DiagFlag(QDLock) + MOD(INT(DIAG_WORD/4),2)
               DiagFlag(QDHigh) = 
     &            DiagFlag(QDHigh) + MOD(INT(DIAG_WORD/2),2)
               DiagFlag(QDLow) = 
     &            DiagFlag(QDLow) + MOD(INT(DIAG_WORD),2)
            ENDDO
        ENDIF
      ENDIF
      WhichTemp = 0
      IF (HAVE_UNCAL(QUTSonic)) THEN
         HAVE_CAL(TSonic) = .TRUE.
	 WhichTemp = TSonic
      ENDIF
      IF (HAVE_UNCAL(QUHumidity)) THEN
         HAVE_CAL(Humidity) = .TRUE.
         HAVE_CAL(SpecHum) = .TRUE.
      ENDIF
      IF (HAVE_UNCAL(QUTCouple)) THEN
         HAVE_CAL(TCouple) = .TRUE.
	 WhichTemp = TCouple
      ENDIF
      IF (HAVE_UNCAL(QUCO2)) THEN
         HAVE_CAL(CO2) = .TRUE.
         HAVE_CAL(SpecCO2) = .TRUE.
      ENDIF
      IF (HAVE_UNCAL(QUDoy) .AND.
     &    HAVE_UNCAL(QUHourMin)) THEN
         Have_Cal(TTime) = .TRUE.
      ENDIF

C
C Check whether we have data needed for corrections
C
      IF (DoCorr(QCSonic) .AND. .NOT. HAVE_CAL(SpecHum)) THEN
         WRITE(*,*) 'WARNING: Do not have a hygrometer: ',
     &                 'can not do ',
     &                 'humidity correction on sonic temperature'
         DoCorr(QCSonic) = .FALSE.
      ENDIF
C
C
C Calibrate the raw samples for the first time ignoring drift in apparatus.
C
C
      DO i=1,N
       	CorMean(i) = 0.0D0
      ENDDO
C
C First try it with thermocouple
C
      IF (HAVE_UNCAL(QUTCouple)) THEN
         NTcOut = 0
         BadTc = (.FALSE.)
         DO i=1,M
            CALL Calibr(RawSampl(1,i),Channels,P,CorMean,
     &               CalSonic,CalTherm,CalHyg,CalCO2,
     &               BadTc,Sample(1,i),N,Flag(1,i),
     &               Have_Uncal, FirstDay, i)
            IF (Flag(TCouple,i)) NTcOut = NTcOut + 1
         ENDDO
         BadTc = (NTcOut.GE.(M/2))
      ELSE
         BadTc = .TRUE.
      ENDIF
C
C If thermocouple too bad, or absent
C
      IF (BadTc) THEN
        DO i=1,M
          CALL Calibr(RawSampl(1,i),Channels,P,CorMean,
     &      CalSonic,CalTherm,CalHyg,CalCO2,
     &      BadTc,Sample(1,i),N,Flag(1,i),
     &      Have_Uncal, FirstDay, i)
        ENDDO
      ENDIF
C
C Calibrate the raw samples for the second time, now correcting mean
C quantities of drift-sensitive apparatus using slow signals.
C
      IF (DoCorr(QCMean)) THEN
C
C Find the shift/drift in humidity of the krypton hygrometer
C
	CALL EC_M_Averag(Sample,NMax,N,MMax,M,Flag,
     &	  Mean,TolMean,Cov,TolCov,MIndep,CIndep,Mok,Cok)
        CALL EC_G_Reset(Have_cal, Mean, TolMean, Cov, TolCov, 
     &                MIndep,  CIndep)
        IF (.NOT. Have_Cal(Humidity)) THEN
	   Mean(Humidity) = Psychro
	   MEAN(SpecHum) = EC_Ph_Q(Mean(Humidity), Mean(WhichTemp), P)
	ENDIF
	IF (DoPrint) THEN
           CALL EC_M_MinMax(Sample, NMax, N, MMax, M, Flag, MINS, MAXS)
           CALL EC_G_ShwMinMax(OutF, N, Mins, Maxs)
	ENDIF

	CorMean(Humidity) = Psychro - Mean(Humidity)
	IF (DoPrint.AND.PCal) THEN
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'Added to humidity : ',
     &	     (Psychro - Mean(Humidity)),' [kg m^{-3}]'
	  WRITE(OutF,*)
        ENDIF

	DO i=1,M
	  CALL Calibr(RawSampl(1,i),Channels,P,CorMean,
     &	    CalSonic,CalTherm,CalHyg,CalCO2,
     &      BadTc,Sample(1,i),N,Flag(1,i),
     &      Have_Uncal, FirstDay, i)
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
     &	               Mean,TolMean,Cov,TolCov,MIndep,
     &                 CIndep,Mok,Cok)
      CALL EC_G_Reset(Have_cal, Mean, TolMean, Cov, TolCov, 
     &                MIndep,  CIndep)
      IF (.NOT. Have_Cal(Humidity)) THEN
        Mean(Humidity) = Psychro
        MEAN(SpecHum) = EC_Ph_Q(Mean(Humidity), Mean(WhichTemp), P)
      ENDIF
      IF (DoPrint) THEN
           CALL EC_G_ShwHead(OutF,'For raw calibrated data : ')
           CALL EC_M_MinMax(Sample, NMax, N, MMax, M, Flag, MINS, MAXS)
           CALL EC_G_ShwMinMax(OutF, N, Mins, Maxs)
           IF (PIndep) THEN
	       CALL EC_G_ShwInd(OutF,MIndep,CIndep,
     &	                         NMax,N,M,ExpVar(QEFreq))
           ENDIF
           IF (PCal) THEN
	       CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,NMax,N)
           ENDIF
      ENDIF

C
C
C Subtract a linear trend from all quantities
C
C
      IF (DoCorr(QCDetrend)) THEN
	CALL EC_M_Detren(Sample,NMax,N,MMAx,M,Mean,Cov,RC)
C
C
C Estimate mean values, covariances and tolerances of both
C for the detrended dataset
C
C
	CALL EC_M_Averag(Sample,NMax,N,MMax,M,Flag,
     &	                 Mean,TolMean,Cov,TolCov,
     &                   MIndep,CIndep,Mok,Cok)
        CALL EC_G_Reset(Have_cal, Mean, TolMean, Cov, TolCov, 
     &                MIndep,  CIndep)
        IF (.NOT. Have_cal(Humidity)) THEN
          Mean(Humidity) = Psychro
          MEAN(SpecHum) = EC_Ph_Q(Mean(Humidity), Mean(WhichTemp), P)
        ENDIF
        CALL EC_M_MinMax(Sample, NMax, N, MMax, M, Flag, MINS, MAXS)

	IF (DoPrint.AND.PCorr(QCDetrend)) THEN
	  CALL EC_G_ShwHead(OutF, 'After detrending : ')
          CALL EC_M_MinMax(Sample, NMax, N, MMax, M, Flag, MINS, MAXS)
          CALL EC_G_ShwMinMax(OutF, N, Mins, Maxs)
	  IF (PIndep) THEN
	    CALL EC_G_ShwInd(OutF,MIndep,CIndep,NMax,N,M,ExpVar(QEFreq))
	  ENDIF
	  CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,NMax,N)
	ENDIF
      ENDIF
C
C Get the mean W before we do tilt correction
C
      QPhys(QPMEANW) = MEAN(W)
      dQPhys(QPMEANW) = TolMean(W)
C
C Correct mean values and covariances for all thinkable effects, first
C time only to get run-based tilt correction (if needed). But 
C before all that, the planar fit untilting needs to be applied.
C
C Do Planar fit tilt correction (only) here 
C
      IF (DoCorr(QCPF)) THEN
C
C Tilt ALL samples
        DO i=1,M
          Speed(1) = Sample(U,i)
          Speed(2) = Sample(V,i)
          Speed(3) = Sample(W,i)
C
C         Speed(3) = Speed3) - WBias

          CALL EC_M_MapVec(Apf,Speed,Dumout)
C
C Feed planarfit rotated sample back into Sample array
C
          Sample(U,i) = Dumout(1)
          Sample(V,i) = Dumout(2)
          Sample(W,i) = Dumout(3)	  
        ENDDO
      ENDIF
      CALL EC_M_Averag(Sample,NMax,N,MMax,M,Flag,
     &	               Mean,TolMean,Cov,TolCov,
     &                 MIndep,CIndep,Mok,Cok)
      CALL EC_G_Reset(Have_cal, Mean, TolMean, Cov, TolCov, 
     &                MIndep,  CIndep)
      IF (.NOT. Have_cal(Humidity)) THEN
          Mean(Humidity) = Psychro
          MEAN(SpecHum) = EC_Ph_Q(Mean(Humidity), Mean(WhichTemp), P)
      ENDIF
      IF (DoPrint.AND. PCorr(QCPF)) THEN
	  CALL EC_G_ShwHead(OUTF, 'After planar fit untilting : ')
	  IF (PIndep) THEN
	    CALL EC_G_ShwInd(OutF,MIndep,CIndep,NMax,N,M,ExpVar(QEFreq))
	  ENDIF
	  CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,NMax,N)
      ENDIF
      
C If 'classic' run-based tilt corrections are not done, 
C the first call to EC_C_Main is redundant and skipped), else
C this run is only done to establish the 'classic' rotation
C angles.
      AnyTilt = ((DoCorr(QCTilt).OR.DoCorr(QCYaw)) .OR.
     &           (DoCorr(QCPitch).OR.DoCorr(QCRoll)))
      IF (AnyTilt) THen
        CALL EC_C_Main(OutF,
     &	  DoPrint,
     &    Mean,NMax,N,TolMean,
     &	  Cov,TolCov,
     &    DoCorr, PCorr, ExpVar,
     &    DirYaw,
     &	  DirPitch,
     &	  DirRoll,
     &	  SonFactr,
     &	  O2Factor,
     &	  CalSonic,CalTherm,CalHyg,
     &    CalCO2, FrCor,
     &	  WebVel, P,Have_cal)
        CALL EC_G_Reset(Have_cal, Mean, TolMean, Cov, TolCov, 
     &                MIndep,  CIndep)
      ENDIF
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
C Estimate mean values, covariances and tolerances of both
C for the plnar fit untilted series
C
C
      IF (AnyTilt) THEN
        DO i=1,3
          DO j=1,3
            DumCov(i,j) = 0.0D0
          ENDDO
        ENDDO
C
C Tilt ALL samples
C     
        DO i=1,M
          Speed(1) = Sample(U,i)
          Speed(2) = Sample(V,i)
          Speed(3) = Sample(W,i)

          IF (DoCorr(QCTilt)) THEN
            CALL EC_C_T06(ExpVar(QEPreYaw),Yaw)
            CALL EC_C_T05(Speed,3,3,DumCov,Yaw)
            CALL EC_C_T08(ExpVar(QEPrePitch),Pitch)
            CALL EC_C_T05(Speed,3,3,DumCov,Pitch)
            CALL EC_C_T10(ExpVar(QEPreRoll),Roll)
            CALL EC_C_T05(Speed,3,3,DumCov,Roll)
          ENDIF

          IF (DoCorr(QCYaw)) THEN
            CALL EC_C_T06(DirYaw,Yaw)
            CALL EC_C_T05(Speed,3,3,DumCov,Yaw)
	  ENDIF
          IF (DoCorr(QCPitch)) THEN
            CALL EC_C_T08(DirPitch,Pitch)
            CALL EC_C_T05(Speed,3,3,DumCov,Pitch)
	  ENDIF
          IF (DoCorr(QCRoll)) THEN
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
     &	                 Mean,TolMean,Cov,TolCov,
     &                   MIndep,CIndep,Mok,Cok)
        CALL EC_G_Reset(Have_cal, Mean, TolMean, Cov, TolCov, 
     &                MIndep,  CIndep)
        IF (.NOT. Have_cal(Humidity)) THEN
          Mean(Humidity) = Psychro
          MEAN(SpecHum) = EC_Ph_Q(Mean(Humidity), Mean(WhichTemp), P)
        ENDIF

	IF (DoPrint) THEN
	  CALL EC_G_ShwHead(OutF, 'After untilting raw data : ')
	  IF (PIndep) THEN
	    CALL EC_G_ShwInd(OutF,MIndep,CIndep,NMax,N,M,ExpVar(QEFreq))
	  ENDIF
	  CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,NMax,N)
	ENDIF

      ENDIF
C
C Make copy of DoCorr and PCorr, and reset tilt info's
C
      DO I=1,NMaxCorr
         DumCorr(I) = DoCorr(I)
         PDumCorr(I) = PCorr(I)
      ENDDO
      DumCorr(QCTilt)  = (.FALSE.)
      DumCorr(QCYaw) = (.FALSE.)
      DumCorr(QCPitch)   = (.FALSE.)
      DumCorr(QCRoll)  = (.FALSE.)

C
C Perform all necessary corrections on the mean values and (co-)variances.
C This is the real thing (irrespective of tilt correction)
C
      CALL EC_C_Main(OutF,
     &	  DoPrint,
     &	  Mean,NMax,N,TolMean,
     &	  Cov,TolCov,
     &    DumCorr, PDumCorr, ExpVar,
     &	  DirYaw,
     &	  DirPitch,
     &	  DirRoll,
     &	  SonFactr,
     &	  O2Factor,
     &	  CalSonic,CalTherm,CalHyg,CalCO2,FrCor,
     &	  WebVel, P, Have_cal)
      CALL EC_G_Reset(Have_cal, Mean, TolMean, Cov, TolCov, 
     &                MIndep,  CIndep)
      
C
C
C Calculate fluxes from processed mean values and covariances.
C
C
C The 180 degrees is to trick ec_Ph_Flux
C
      IF (.NOT. DoCorr(QCYaw)) THEN
         DirYaw = 180D0
      ELSE
         DirYaw = DirYaw + 180D0
      ENDIF
      
      CALL EC_Ph_Flux(Mean,NMax,Cov,TolMean,TolCov,p,BadTc,
     &                QPhys, dQPhys, WebVel, DirYaw)

      CALL EC_G_Reset(Have_cal, Mean, TolMean, Cov, TolCov, 
     &                MIndep,  CIndep)

      RETURN
      END
      




      SUBROUTINE EC_G_Reset(Have_cal, Mean, TolMean, Cov, TolCov,
     &                      Mindep, CIndep)
C     ****f* ec_gene.f/EC_G_Reset
C NAME
C     EC_G_Reset
C SYNOPSIS
C     CALL EC_G_Reset(Have_cal, Mean, TolMean, Cov, TolCov, MIndep,
C     CIndep)
C FUNCTION
C     Routine to reset means and covariances based on availability
C     of the uncalibrated data
C INPUTS
C     Have_cal  : [LOGICAL(NMax)]
C                   switch for each channel whether calibrated
C                   data are available
C     Mean        : [REAL*8(NMax)]
C                   mean of quantities
C     TolMean     : [REAL*8(NMax)]
C                   tolerance in mean of quantities
C     Cov         : [REAL*8(NMax,NMax)]
C                   covariances of quantities
C     TolMean     : [REAL*8(NMax,NMax)]
C                   tolerance in covariances of quantities
C     MIndep      : [REAL*8(NMax)]
C                   number of independent samples
C     CIndep      : [REAL*8(NMax,NMax)]
C                   number of independent samples in covariances 
C OUTPUT
C     Mean        : [REAL*8(NMax)]
C                   mean of quantities
C     TolMean     : [REAL*8(NMax)]
C                   tolerance in mean of quantities
C     Cov         : [REAL*8(NMax,NMax)]
C                   covariances of quantities
C     TolMean     : [REAL*8(NMax,NMax)]
C                   tolerance in covariances of quantities
C     MIndep      : [REAL*8(NMax)]
C                   number of independent samples
C     CIndep      : [REAL*8(NMax,NMax)]
C                   number of independent samples in covariances 
C AUTHOR
C     Arnold Moene
C HISTORY
C     $Name$
C     $Id$
C USES
C     parcnst.inc
C     ***
      IMPLICIT NONE
      INCLUDE 'parcnst.inc'

      LOGICAL Have_cal(NNMax)
      REAL*8  Mean(NNmax), TolMean(NNMax), Cov(NNMax, NNMax),
     +        TolCov(NNMax, NNMax)
      INTEGER I,J, MIndep(NNMax), CIndep(NNMax, NNMax)

      DO I=1,NNMax
         IF (.NOT. Have_CAL(I)) THEN
             MEAN(I) = DUMMY
             TolMean(I) = DUMMY
             MIndep(I) = DUMMY
         ENDIF
         DO J=1,NNMax
            IF ((.NOT. HAVE_CAL(I)).OR.
     +          (.NOT. HAVE_CAL(J))) THEN
                COV(I,J) = DUMMY
                TolCov(I,J) = DUMMY
                CIndep(I,J) = INT(DUMMY)
            ENDIF
         ENDDO
      ENDDO
 
      END


C
C ########################################################################
C
C Some silly routines to print intermediate results to file. Handy
C for checking the origin of unexpected results!
C
C ########################################################################
C





      SUBROUTINE EC_G_Show(OutF,Mean,TolMean,Cov,
     &	TolCov,NMax,N)
C    ****f* ec_gene.f/EC_G_Show
C NAME
C     EC_G_Show
C SYNOPSIS
C     CALL EC_G_Show(OutF,Mean,TolMean,Cov,
C     	             TolCov,NMax,N)
C FUNCTION
C     Prints mean values and (co-)variances plus respective 
C     tolerances
C INPUTS
C     OUTF    : [INTEGER]
C               unit number of file
C     Mean    : [REAL*8(NMax)]
C               means of quantities
C     TolMean : [REAL*8(NMax)]
C               tolerance in means of quantities
C     Cov     : [REAL*8(NMax,NMax)]
C               covariances  of quantities
C     TolCov  : [REAL*8(NMax,NMax)]
C               tolerance in covariances of quantities
C     NMax    : [INTEGER]
C               maximum number of quantities
C     N       : [INTEGER]
C               actual number of quantities
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$
C     $Id$
C USES
C     Epsilon
C     parcnst.inc
C     ***
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER i,j,N,NMax,OutF
      REAL*8 Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),Dum(NNMax)

      WRITE(OutF,*)
      WRITE(OutF,*) 'Mean quantities : '
      WRITE(OutF,*)

      DO i=1,N
	WRITE(OutF,10) QName(i),Mean(i),TolMean(i),UName(i)
      ENDDO
 10   FORMAT(a6,' = ',F20.10,' +/- ',F20.10,1X,a10,1X)

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





      SUBROUTINE EC_G_ShwFrq(OutF,FrCor,NMax,N)
C     ****f* ec_gene.f/EC_G_ShwFrq
C NAME
C     EC_G_ShwFrq
C SYNOPSIS
C     CALL EC_G_ShwFrq(OutF,FrCor,NMax,N)
C FUNCTION
C     Prints correction factors associated with frequency response
C INPUTS
C     OUTF    : [INTEGER]
C               unit number of file
C     FrCor   : [REAL*8(NMax,NMax)]
C               freqency response correction factors for
C               each combination of quantities
C     NMax    : [INTEGER]
C               maximum number of quantities
C     N       : [INTEGER]
C               actual number of quantities
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
      INTEGER i,j,N,NMax,OutF
      REAL*8 FrCor(NMax,NMax)

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





      SUBROUTINE EC_G_ShwInd(OutF,MIndep,
     &	CIndep,NMax,N,M,Freq)
C     ****f* ec_gene.f/EC_G_ShwInd
C NAME
C     EC_G_ShwInd
C SYNOPSIS
C     CALL EC_G_ShwInd(OutF,MIndep,
C                      CIndep,NMax,N,M,Freq)
C FUNCTION
C     Prints number of independent observations
C INPUTS
C     OUTF    : [INTEGER]
C               unit number of file
C     MIndep  : [INTEGER(NMax)]
C               number of independent samples for
C               means
C     CIndep  : [INTEGER(NMax,NMax)]
C               number of independent samples for
C               covariances
C     NMax    : [INTEGER]
C               maximum number of quantities
C     N       : [INTEGER]
C               actual number of quantities
C     M       : [INTEGER]
C               total number of samples
C     Freq    : [REAL*8]
C               sampling frequency (Hz)
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
      INTEGER i,j,N,NMax,OutF,MIndep(NMax),CIndep(NMax,NMax),M
      REAL*8 TmpTime(NMax)
      REAL*8 Freq

      WRITE(OutF,*)
      WRITE(OutF,*) ' Number of independent samples: '
      WRITE(OutF,*)

      WRITE(OutF,*) 'For mean quantities : '
      WRITE(OutF,*) '                #    timescale [s]'
      WRITE(OutF,*)

      DO i=1,N
        IF (MIndep(I).EQ. INT(DUMMY)) THEN
	  WRITE(OutF,10) QName(i),DUMMY,DUMMY
        ELSE
	  WRITE(OutF,10) QName(i),MIndep(i),
     &      (DBLE(M)/(Freq*DBLE(MIndep(i))))
        ENDIF
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
        DO J=1,N
          IF (CIndep(I,J) .EQ. INT(DUMMY)) THEN
             TmpTime(J) = DUMMY
          ELSE
             TmpTime(J) = DBLE(M)/(Freq*DBLE(CIndep(i,j)))
          ENDIF
        ENDDO
	WRITE(OutF,50) QName(i),(TmpTime(J),J=1,N)
 50	FORMAT(a6,20(F9.3:,1X))
      ENDDO

      RETURN
      END

      SUBROUTINE EC_G_ShwMinMax(OutF, N, Mins, Maxs)
C     ****f* ec_gene.f/EC_G_ShwMinMax
C NAME
C     EC_G_ShwMinMax
C SYNOPSIS
C     CALL EC_G_ShwInd(OutF, N, Mins, Maxs)
C FUNCTION
C     Prints min/max of series
C INPUTS
C     OUTF    : [INTEGER]
C               unit number of file
C     N       : [INTEGER]
C               number of series
C     Mins    : [INTEGER(NMax)]
C               min value of series
C     Maxs    : [INTEGER(NMax)]
C               max value of series
C AUTHOR
C     Arnold Moene
C HISTORY
C     $Name$
C     $Id$
C USES
C     parcnst.inc
C     ***
      IMPLICIT NONE

      INCLUDE 'parcnst.inc'
      INTEGER i,N,OutF
      REAL*8   Mins(N),maxs(N)
      
      WRITE(OUTF,*) 'Min/max of samples'
      DO I=1,N
         WRITE(OUTF,*) QName(I),': min = ', Mins(I), ', max = ',
     &                 Maxs(I)
      ENDDO
      RETURN
      END

      SUBROUTINE EC_G_ShwHead(OutF, STRING)
C     ****f* ec_gene.f/EC_G_ShwHead
C NAME
C     EC_G_ShwHead
C SYNOPSIS
C     CALL EC_G_ShwHead(OutF,String)
C FUNCTION
C     Prints header to intermediate results file
C INPUTS
C     OUTF    : [INTEGER]
C               unit number of file
C     N       : [CHARACTER](*)
C               String to write
C AUTHOR
C     Arnold Moene
C HISTORY
C     $Name$
C     $Id$
C     ***
      IMPLICIT NONE

      INTEGER OutF
      CHARACTER*(*) String
      
      WRITE(OUTF,*) 
      WRITE(OUTF,*) String
      WRITE(OUTF,*) 
      RETURN
      END




