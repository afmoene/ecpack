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
     &	Psychro,Freq,CalSonic,CalTherm,CalHyg, CalCO2, P,
     &	Calibr,
     &	Sample,Flag,Mok,Cok,MIndep,CIndep,Rc,BadTc,
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
     & FCO2, dFCO2, FCO2Webb, dFCO2Webb,
     &  MEANW, TOLMEANW, HAVE_UNCAL, FirstDay)
C     ****f* ec_gene.f/EC_G_Main
C NAME
C     EC_G_Main
C SYNOPSIS
C     CALL EC_G_Main(OutF,DoPrint,
C        RawSampl,MaxChan,Channels,NMax,N,MMax,M,DoCrMean,PCal,PIndep,
C        Psychro,Freq,CalSonic,CalTherm,CalHyg, CalCO2, P,
C        Calibr,
C        Sample,Flag,Mok,Cok,MIndep,CIndep,Rc,BadTc,
C        DoDetren,PDetrend,
C        DoTilt,PTilt,PreYaw,PrePitch,PreRoll,
C        DoYaw,PYaw,DirYaw,
C        DoPitch,PPitch,PitchLim,DirPitch,
C        DoRoll,PRoll,RollLim,DirRoll,
C        DoSonic,PSonic,SonFactr,
C        DoO2,PO2,O2Factor,
C        DoFreq,PFreq,LLimit,ULimit,FrCor,
C        DoWebb,PWebb,
C        Mean,TolMean,Cov,TolCov,
C        HSonic,dHSonic,HTc,dHTc,
C        LvE,dLvE,LvEWebb,dLvEWebb,
C        UStar,dUStar,Tau,dTau,
C        FCO2, dFCO2, FCO2Webb, dFCO2Webb,
C        MEANW, TOLMEANW, HAVE_UNCAL, FirstDay)
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
C     DoCrMean  : [LOGICAL]
C                 correct mean quantities of drift-sensitive apparatus 
C                 using slow signal?
C     PCal      : [LOGICAL]
C                 print results of slow sensor correction?
C     PIndep    : [LOGICAL]
C                 print number of independent samples?
C     Psychro   : [REAL*8]
C                 water vapour density of slow sensor (kg/m^3)
C     Freq      : [REAL*8]
C                 sampling frequency
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
C     DoDetren  : [LOGICAL]
C                 do linear deterending?
C     PDetrend  : [LOGICAL]
C                 print results of linear detrending?
C     DoTilt    : [LOGICAL]
C                 do tilt correction with known tilt angles
C     PTilt     : [LOGICAL]
C                 print results of tilting?
C     PreYaw    : [REAL*8]
C                 known yaw angle (degree)
C     PrePitch  : [REAL*8]
C                 known pitch angle (degree)
C     PreRoll   : [REAL*8]
C                 known roll angle (degree)
C     DoYaw     : [LOGICAL]
C                 do yaw correction?
C     PYaw      : [LOGICAL]
C                 print results of yaw correction?
C     DoPitch   : [LOGICAL]
C                 do pitch correction?
C     PPitch    : [LOGICAL]
C                 print results of pitch correction?
C     PitchLim  : [REAL*8]
C                 maximum allowable pith rotation angle (degrees)
C     DoRoll    : [LOGICAL]
C                 do roll correction?
C     PRoll     : [LOGICAL]
C                 print results of roll correction?
C     RollLim   : [REAL*8]
C                 maximum allowable roll rotation angle (degrees)
C     DoSonic   : [LOGICAL]
C                 do correction of sonic temperature (Schotanus)?
C     PSonic    : [LOGICAL]
C                 print results of sonic temperature correction?
C     DoO2      : [LOGICAL]
C                 Correct hygrometer data for oxygen sensitivity ?
C     PO2       : [LOGICAL]
C                 Print intermediate result for oxygen correction ?
C     DoFreq    : [LOGICAL]
c                 Do frequency response correction ?
C     PFreq     : [LOGICAL]
C                 Print intermediate results for frequency response
C                 correction ?
C     LLimit    : [REAL*8]
C                 Lower limit for correction factor for frequency response
C     ULimit    : [REAL*8]
C                 Upper limit for correction factor for frequency response
C     DoWebb    : [LOGICAL]
C                 Do Webb correction for humidity flux ?
C     PWebb     : [LOGICAL]
C                 Print intermediate results for Webb correction ?
C     HAVE_UNCAL: [LOGICAL(NMax)]
C                 switch whether data for uncalibrated data
C                 are available for each channel
C     FirstDay  : [INTEGER]
C                 day number of first sample in array (needed for
C                 detrending data that pass midnight
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
C     HSonic    : [REAL*8]
C                 sensible heat flux with sonic temperature (W/m^2)
C     dHSonic   : [REAL*8]
C                 tolerance in sensible heat flux with sonic temperature (W/m^2)
C     HTc       : [REAL*8]
C                 sensible heat flux with thermocouple temperature (W/m^2)
C     dHTc      : [REAL*8]
C                 tolerance in sensible heat flux with thermocouple temperature (W/m^2)
C     LvE       : [REAL*8]
C                 wq-covariance latent heat flux (W/m^2)
C     dLvE      : [REAL*8]
C                 tolerance in wq-covariance latent heat flux (W/m^2)
C     LvEWebb   : [REAL*8]
C                 Webb term for latent heat  flux (W/m^2)
C     dLvEWebb   : [REAL*8]
C                 tolerance in Webb term for latent heat  flux (W/m^2)
C     UStar     : [REAL*8]
C                 friction velocity (m/s)
C     dUStar    : [REAL*8]
C                 tolerance in friction velocity (m/s)
C     Tau       : [REAL*8]
C                 surface shear stress (N/m^2)
C     dTau      : [REAL*8]
C                 tolerance in surface shear stress (N/m^2)
C     FCO2      : [REAL*8]
C                 CO2 flux from covariance (kg/m^2 s)
C     dFCO2     : [REAL*8]
C                 tolerance in CO2 flux from covariance (kg/m^2 s)
C     FCO2Webb  : [REAL*8]
C                 Webb contribution to CO2 flux (kg/m^2 s)
C     dFCO2Webb : [REAL*8]
C                 tolerance in Webb contribution to CO2 flux (kg/m^2 s)
C     MEANW     : [REAL*8]
C                 mean vertical velocity before tilt correction (m/s)
C     TOLMEANW  : [REAL*8]
C                 tolerance in mean vertical velocity before 
C                 tilt correction (m/s)
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     Revision: 03-04-2001: get mean W and its tolerance before tilt
C                           correction; export via interface (AM)
C     Revision: 28-05-2001: added passing of info on whether uncalibrated
C                           data are available (AM)
C     Revision: 18-09-2002: removed calcomm.inc and added FirstDay
C                           to interface (to pass it to calibration routine
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
     & WhichTemp, FirstDay
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
     &	UStar,dUStar,Tau,dTau,Speed(3),DumCov(3,3),
     &  FCO2, dFCO2, FCO2Webb, dFCO2Webb
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ), CalCO2(NQQ)
      REAL*8 MEANW, TOLMEANW, MINS(NNMax), MAXS(NNMAX)
      REAL*8 EC_Ph_Q,Yaw(3,3),Pitch(3,3),Roll(3,3)
      LOGICAL HAVE_CAL(NMax)
      LOGICAL HAVE_UNCAL(NMax)
C
C Check whether we have the needed calibration info
C
      DO I=1,NMax
         HAVE_CAL(i) = .FALSE.
      ENDDO
      DO I=U,W
         IF (HAVE_UNCAL(I)) HAVE_CAL(I) = .TRUE.
      ENDDO
      IF (HAVE_UNCAL(TSonic)) THEN
         HAVE_CAL(TSonic) = .TRUE.
	 WhichTemp = TSonic
      ENDIF
      IF (HAVE_UNCAL(Humidity)) THEN
         HAVE_CAL(Humidity) = .TRUE.
         HAVE_CAL(SpecHum) = .TRUE.
      ENDIF
      IF (HAVE_UNCAL(TCouple)) THEN
         HAVE_CAL(TCouple) = .TRUE.
	 WhichTemp = TCouple
      ENDIF
      IF (HAVE_UNCAL(CO2)) THEN
         HAVE_CAL(CO2) = .TRUE.
         HAVE_CAL(SpecCO2) = .TRUE.
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
C
C Calibrate the raw samples for the first time ignoring drift in apparatus.
C
C
      DO i=1,N
       	CorMean(i) = 0.D0
      ENDDO
C
C First try it with thermocouple
C
      Have_Uncal(TTime) = .TRUE.
      IF (HAVE_UNCAL(TCouple)) THEN
         NTcOut = 0
         BadTc = (.FALSE.)
         DO i=1,M
            CALL Calibr(RawSampl(1,i),Channels,P,CorMean,
     &               CalSonic,CalTherm,CalHyg,CalCO2,
     &               BadTc,Sample(1,i),N,Flag(1,i),
     &               Have_Uncal, FirstDay)
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
     &      Have_Uncal, FirstDay)
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
     &	    CalSonic,CalTherm,CalHyg,CalCO2,
     &      BadTc,Sample(1,i),N,Flag(1,i),
     &      Have_Uncal, FirstDay)
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
	CALL EC_G_ShwInd(OutF,MIndep,CIndep,NMax,N,M,Freq)
      ENDIF
C
C Print averages of raw, calibrated data
C
      IF (DoPrint.AND.PCal) THEN
	WRITE(OutF,*)
	WRITE(OutF,*) 'For raw calibrated data : '
	WRITE(OutF,*)
	CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,NMax,N)
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
	    CALL EC_G_ShwInd(OutF,MIndep,CIndep,NMax,N,M,Freq)
	  ENDIF

	  CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,NMax,N)
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
     &  BadTc,
     &	DoTilt,PTilt,PreYaw,PrePitch,PreRoll,
     &	DoYaw,PYaw,DirYaw,
     &	DoPitch,PPitch,PitchLim,DirPitch,
     &	DoRoll,PRoll,RollLim,DirRoll,
     &	DoSonic,PSonic,SonFactr,
     &	DoO2,PO2,O2Factor,
     &	DoFreq,PFreq,LLimit,ULimit,Freq,CalSonic,CalTherm,CalHyg,
     &  CalCO2, FrCor,
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
	    CALL EC_G_ShwInd(OutF,MIndep,CIndep,NMax,N,M,Freq)
	  ENDIF

	  CALL EC_G_Show(OutF,Mean,TolMean,Cov,TolCov,NMax,N)
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
     &    BadTc,
     &	  DumTilt,PTilt,PreYaw,PrePitch,PreRoll,
     &	  DumYaw,PYaw,DirYaw,
     &	  DumPitch,PPitch,PitchLim,DirPitch,
     &	  DumRoll,PRoll,RollLim,DirRoll,
     &	  DoSonic,PSonic,SonFactr,
     &	  DoO2,PO2,O2Factor,
     &	  DoFreq,PFreq,LLimit,ULimit,Freq,
     &    CalSonic,CalTherm,CalHyg,CalCO2,FrCor,
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
     &	UStar,dUStar,Tau,dTau, FCO2, dFCO2, FCO2Webb, dFCO2Webb)

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
      IF (.NOT. HAVE_Uncal(CO2)) THEN
	  FCO2 = DUMMY
	  dFCO2 = DUMMY
	  FCO2Webb = DUMMY
	  dFCO2Webb = DUMMY
      ENDIF

      RETURN
      END
      




      SUBROUTINE EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)
C     ****f* ec_gene.f/EC_G_Reset
C NAME
C     EC_G_Reset
C SYNOPSIS
C     CALL EC_G_Reset(Have_Uncal, Mean, TolMean, Cov, TolCov)
C FUNCTION
C     Routine to reset means and covariances based on availability
C     of the uncalibrated data
C INPUTS
C     Have_Uncal  : [LOGICAL(NMax)]
C                   switch for each channel whether uncalibrated
C                   data are available
C     Mean        : [REAL*8(NMax)]
C                   mean of quantities
C     TolMean     : [REAL*8(NMax)]
C                   tolerance in mean of quantities
C     Cov         : [REAL*8(NMax,NMax)]
C                   covariances of quantities
C     TolMean     : [REAL*8(NMax,NMax)]
C                   tolerance in covariances of quantities
C OUTPUT
C     Mean        : [REAL*8(NMax)]
C                   mean of quantities
C     TolMean     : [REAL*8(NMax)]
C                   tolerance in mean of quantities
C     Cov         : [REAL*8(NMax,NMax)]
C                   covariances of quantities
C     TolMean     : [REAL*8(NMax,NMax)]
C                   tolerance in covariances of quantities
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
      REAL*8 Freq

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




