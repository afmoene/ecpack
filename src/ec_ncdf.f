C $Id$
C  ec_ncdf.f, part of ecpack
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



C...........................................................................
C Function  : ec_ncdf
C Purpose   : to do eddy correlation calculations on data from
C             a NetCDF file
C Interface : a configuration file is read from standard input
C Author    : Arjan van Dijk
C             Arnold Moene
C Date      : May 24, 2000
C Updates   : September 27 2000: added option for wind tunnel calibration
C                                of sonic (for EBEX data)
C             October 11 2000:   check for validity of reference temperature
C                                for thermocouple and for validity of
C                                iteration used in thermocouple calibration
C             October 18, 2000:  in calibration of wind tunnel calibrated
C                                sonic: check for valid wind direction
C                                before interation (change in ec_pack)
C             October 30, 2000:  check for validity of reference temperature
C                                for thermocouple: Tref is in Celcius (not
C                                Kelvin !)
C             November 1, 2000:  fixed error in sonic calibration (azimuth
C                                angle replaced elevation angle)
C             February 28, 2001: fixed error in computation of density of
C                                dry air (it was correct in the original
C                                version of ECpack, but incorrectly changed
C                                by me (AM)
C             April 3, 2001:     get mean W and its tolerance (before tilt
C                                correction) from ECMain and write it flux
C                                file in two added columns
C             Sept. 18, 2002:    removed calcomm.inc and added HAVE_UNCAL
C                                to store which input data we have
C Current version: 1.11
C...........................................................................
      PROGRAM EC_NCDF

      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'
      INCLUDE 'version.inc'

      CHARACTER*255 FNAME, DatDir, OutDir, ParmDir, FluxName,
     &              ParmName,InterName, PlfName
      CHARACTER*40 OutName,DumName1,DumName2
      LOGICAL PRaw,PYaw,PPitch,PRoll,PFreq,PO2,PWebb,PCal,PIndep,
     &  DoYaw,DoPitch,DoRoll,DoFreq,DoO2,DoWebb,DoCrMean,PSonic,
     &  DoSonic,DoTilt,PTilt,DoPrint,Flag(NNMAx,MMMax),DoDetren,
     &  PDetrend,DoStruct,BadTc
      INTEGER N,i,M,MIndep(NNMax),CIndep(NNMax,NNMax),FOO,
     &  Channels,Delay(NNNMax),Mok(NNMax),Cok(NNMax,NNMax), FirstDay
      REAL*8 RawSampl(NNNMax,MMMax),Sample(NNMax,MMMax),P,Psychro,
     &  Mean(NNMax),TolMean(NNMax),Cov(NNMax,NNMax),
     &  TolCov(NNMax,NNMax),LLimit,ULimit,FrCor(NNMax,NNMax),
     &  PitchLim,RollLim,DirYaw,RC(NNMax),O2Factor(NNMax),
     &  Freq,DirPitch,DirRoll,
     &  SonFactr(NNMax),PreYaw,PrePitch,PreRoll,
     &  HSonic,dHSonic,HTc,dHTc,
     &  LvE,dLvE,LvEWebb,dLvEWebb,
     &  UStar,dUStar,Tau,dTau,CorMean(NNMax),DirFrom,
     &  Gain(NNNMax),Offset(NNNMax),StartTime(3),StopTime(3),
     &  StructSep,
     &  FCO2,dFCO2, FCO2Webb, dFCO2Webb
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ), CalCO2(NQQ),
     &  R,dR,CTSon2,CTCop2,Cq2,CTSonq,CTCopq,
     &  dCTSon2,dCTCop2,dCq2,dCTSonq,dCTCopq
      REAL*8 MeanW, TolMeanW,
     &       StdTs, dStdTs,
     &       StdTc, dStdTc,
     &       Stdq, dStdq,
     &       StdU, dStdU,
     &       StdV, dStdV,
     &       StdW, dStdW,
     &       StdspecCO2, dStdspecCO2,
     &       StdCO2, dStdCO2,
     &       SumLvE, dSumLvE,
     &       SumFCO2, dSumFCO2
      CHARACTER*255 SonName,CoupName,HygName, CO2Name, DUMSTRING
      CHARACTER*255 NCvarname(NNNMax)
      LOGICAL       HAVE_UNCAL(NNNMax)

      INTEGER EC_T_STRLEN
C
C The name of the calibration routine specified at the end of this file
C
      EXTERNAL Calibrat, EC_T_STRLEN

C...........................................................................
C Start of executable statements
C...........................................................................
C
C Open and read configuration file
C

      IF (ECConfFile .NE. 5) THEN
         OPEN(ECConfFile, FILE='ecconf')
      ENDIF

C Report version number (VERSION is defined in version.inc)
      WRITE(*,*) 'EC_NCDF Version: ', VERSION
C Give some RCS info (do not edit this!!, RCS does it for us)
      WRITE(*,*) '$Name$'
      WRITE(*,*) '$Date$'
      WRITE(*,*) '$Revision$'

      CALL EC_F_GetConf(ECConfFile,
     &             DatDir, OutDir, ParmDir,
     &             FluxName, ParmName, InterName,
     &             PlfName,
     &             SonName, CoupName, HygName, CO2Name,
     &             NCVarname, NNNMax)
C
C Assume first we have no uncalibrated samples at all
C
      DO I=1,NNNMax
         HAVE_UNCAL(I) = .FALSE.
      ENDDO
C
C Check whether we have a reference temperature for a thermocouple
      IF (EC_T_STRLEN(NCVarName(Tref)) .GT. 0) THEN
         HAVE_UNCAL(Tref) = .TRUE.
      ENDIF
C
C Check which calibration files we have
C
      IF (EC_T_STRLEN(SonName) .GT. 0) THEN
C For now assume it is always a 3-D sonic
         HAVE_UNCAL(U) = .TRUE.
         HAVE_UNCAL(V) = .TRUE.
         HAVE_UNCAL(W) = .TRUE.
         HAVE_UNCAL(TSonic) = .TRUE.
      ENDIF
      IF (EC_T_STRLEN(HygName) .GT. 0) THEN
         HAVE_UNCAL(Humidity) = .TRUE.
      ENDIF
      IF (EC_T_STRLEN(CoupName) .GT. 0) THEN
         HAVE_UNCAL(TCouple) = .TRUE.
      ENDIF
      IF (EC_T_STRLEN(CO2Name) .GT. 0) THEN
         HAVE_UNCAL(CO2) = .TRUE.
      ENDIF
C
C Open an output file for the fluxes
C
      OPEN(FluxFile,FILE=FluxName,FORM='FORMATTED')
      WRITE(FluxFile,56)
 56   FORMAT(
     &  'DOY Hr Mn ',
     &  'DOY Hr Mn ',
     &  '#samples ',
     &  '  #U     #V     #W  #TSon  #TCop  #RhoV     #q  #time   ',
     &  'Dir d(dir)   ',
     &  'Mean(U)        dMean(U)        ',
     &  'Mean(TSon)     dMean(TSon)     ',
     &  'Mean(TCop)     dMean(TCop)     ',
     &  'Mean(q)        dMean(q)        ',
     &  'sigma(TSon)    dsigma(TSon)    ',
     &  'sigma(TCop)    dsigma(TCop)    ',
     &  'sigma(q)       dsigma(q)       ',
     &  'sigma(u)       dsigma(u)       ',
     &  'sigma(v)       dsigma(v)       ',
     &  'sigma(w)       dsigma(w)       ',
     &  'cov(TSon,q)    dcov(TSon,q)    ',
     &  'cov(TCop,q)    dcov(TCop,q)    ',
     &  'cov(TSon,U)    dcov(TSon,U)    ',
     &  'cov(TCop,U)    dcov(TCop,U)    ',
     &  'cov(q,U)       dcov(q,U)       ',
     &  'H(Sonic)       Tol(HSonic)     ',
     &  'H(TCouple)     Tol(HTCouple    ',
     &  'LvECov         dLvECov         ',
     &  'LvEWebb        dLvEWebb        ',
     &  'LvE            dLvE            ',
     &  'UStar          dUStar          ',
     &  'Tau            dTau            ',
     &  'R(delay)       dR              ',
     &  'CTSon2         dCTSon2         ',
     &  'CTCop2         dCTCop2         ',
     &  'Cq2            dCq2            ',
     &  'CTSonq         dCTSonq         ',
     &  'CTCopq         dCTCopq         ',
     &  'MeanW          dMeanW          ',
     &  '#CO2          ',
     &  'MeanCO2        dMeanCO2        ',
     &  'MeanspecCO2    dMeanspecCO2    ',
     &  'stdCO2         dstdCO2         ',
     &  'stdspecCO2     dstdspecCO2     ',
     &  'FCO2Cov        dFCO2Cov        ',
     &  'FCO2Webb       dFCO2Webb       ',
     &  'FCO2           dFCO2           '
     &  )
C
C Number of quantities involved in this experiment
C
      N = 9  ! Number of calibrated quantities you will follow
C
C Read parameters from file with name ParmName
C
      CALL EC_F_Params(ParmName,
     &  Freq,PitchLim,RollLim,PreYaw,PrePitch,PreRoll,
     &  LLimit,ULimit,DoCrMean,DoDetren,DoSonic,
     &  DoTilt,DoYaw,DoPitch,DoRoll,DoFreq,DoO2,DoWebb,DoStruct,DoPrint,
     &  PRaw,PCal,PDetrend,PIndep,PTilt,PYaw,PPitch,
     &  PRoll,PSonic,PO2,PFreq,PWebb, StructSep)
C
C Read calibration data from files :
C
      DO i=1,NNNMax
        Delay( i) = 0
        Gain(  i) = 1.D0
        Offset(i) = 0.D0
      ENDDO
      IF (EC_T_STRLEN(SonName) .GT. 0) THEN
           CALL EC_F_ReadAp(SonName,CalSonic) ! The sonic
           IF ((CalSonic(QQType) .NE. ApGillSonic) .AND.
     &         (CalSonic(QQType) .NE. ApSon3DCal)  .AND.
     &         (CalSonic(QQType) .NE. ApKaijoTR90)  .AND.
     &         (CalSonic(QQType) .NE. ApKaijoTR61)) THEN
               WRITE(*,*) 'ERROR: Calibration file ',SonName,
     &                    'does not contain sonic info'
               STOP 'Rewrite your calibration file'
           ENDIF
      ENDIF
      IF (EC_T_STRLEN(CoupName) .GT. 0) THEN
           CALL EC_F_ReadAp(CoupName,CalTherm) ! A thermo-couple
           IF (CalTherm(QQType) .NE. ApTCouple) THEN
               WRITE(*,*) 'ERROR: Calibration file ', CoupName,
     &                    'does not contain thermocouple info'
               STOP 'Rewrite your calibration file'
           ENDIF
      ENDIF
      IF (EC_T_STRLEN(HygName) .GT. 0) THEN
           CALL EC_F_ReadAp(HygName,CalHyg) ! A hygrometer
           IF ((CalHyg(QQType) .NE. ApCampKrypton) .AND.
     &         (CalHyg(QQType) .NE. ApMierijLyma) .AND.
     &         (CalHyg(QQType) .NE. ApLiCor7500)) THEN
               WRITE(*,*) 'ERROR: Calibration file ', HygName,
     &                    'does not contain hygrometer info'
               STOP 'Rewrite your calibration file'
           ENDIF
      ENDIF
      IF (EC_T_STRLEN(CO2Name) .GT. 0) THEN
           CALL EC_F_ReadAp(CO2Name,CalCO2) ! A CO2 sensor 
           IF (CalCO2(QQType) .NE. ApLiCor7500) THEN
               WRITE(*,*) 'ERROR: Calibration file ', CO2Name,
     &                    'does not contain CO2 sensor info'
               STOP 'Rewrite your calibration file'
           ENDIF
      ENDIF
C
C Set delays, gains and offset of all channels in raw datafile
C At reading time all channels are corrected by mapping:
C                x --> (x/Gain) + Offset
C

      Delay(U)        = NINT(CalSonic(QQDelay)*0.001D0*Freq)
      Delay(V)        = NINT(CalSonic(QQDelay)*0.001D0*Freq)
      Delay(W)        = NINT(CalSonic(QQDelay)*0.001D0*Freq)
      Delay(Humidity) = NINT(CalHyg(QQDelay)*0.001D0*Freq)
      Delay(Tcouple)  = NINT(CalTherm(QQDelay)*0.001D0*Freq)
      Delay(Tsonic)   = NINT(CalSonic(QQDelay)*0.001D0*Freq)
      Delay(CO2)      = NINT(CalCO2(QQDelay)*0.001D0*Freq)

      Gain(U)        = CalSonic(QQGain)
      Gain(V)        = CalSonic(QQGain)
      Gain(W)        = CalSonic(QQGain)
      Gain(Humidity) = CalHyg(QQGain)
      Gain(Tcouple)  = CalTherm(QQGain)
      Gain(Tsonic)   = CalSonic(QQGain)
      Gain(CO2)      = CalCO2(QQGain)

      Offset(U)        = CalSonic(QQOffset)
      Offset(V)        = CalSonic(QQOffset)
      Offset(W)        = CalSonic(QQOffset)
      Offset(Humidity) = CalHyg(QQOffset)
      Offset(Tcouple)  = CalTherm(QQOffset)
      Offset(Tsonic)   = CalSonic(QQOffset)
      Offset(CO2)      = CalCO2(QQOffset)

C
C Get names of files
C
      OPEN(IntervalFile,FILE = InterName,IOSTAT=FOO,STATUS='OLD')
      IF (FOO.NE.0) THEN
        WRITE(*,*) 'Could not open interval-file ',InterName
        STOP'- Fatal: no time info could be read -'
      ENDIF
C Read one line of comments (blindly, it is supposed to be a comment line)
      READ(IntervalFile,*)
C
C######################################################################
C######################################################################
C
C Here starts the loop over the different files
C
C######################################################################
C######################################################################
C
C
C Read time-interval from file, plus value of psychrometer RhoV and atm.  pressure
C
 36   READ(IntervalFile,*,END=37) (StartTime(i),i=1,3),
     &  (StopTime(i),i=1,3),Psychro,p,DumName1,DumName2
      FName   = DumName1
      OutName = DumName2
C
      DUMSTRING=DatDir
      CALL EC_T_STRCAT(DUMSTRING, '/')
      CALL EC_T_STRCAT(DUMSTRING, Fname)
      Fname=DUMSTRING

      DUMSTRING=OutDir
      CALL EC_T_STRCAT(DUMSTRING, '/')
      CALL EC_T_STRCAT(DUMSTRING, OutName)
      OutName=DUMSTRING

C
C Show which file is currently being analysed
C
      WRITE(*,951) FName(:EC_T_STRLEN(Fname)),(NINT(StartTime(i)),i=1,3)
 951  FORMAT(a,1X,3(I5,1X))
C
C If wanted, make a file for printing intermediate results
C
      IF (DoPrint) THEN
        OPEN(OutFile,FILE=OutName,FORM='FORMATTED')
        WRITE(OutFile,54) FName(:EC_T_STRLEN(Fname))
 54     FORMAT('File = ',A)
        WRITE(OutFile,*) 'From (day,hour,minute)',
     &    (StartTime(i),' ',i=1,3),' to ',(StopTime(i),i=1,3)
        WRITE(OutFile,*) 'At RhoV(slow) = ',Psychro,' [kg m^{-3}]',
     &    ' and atm. pressure = ',p,' [Pascal]'
      ENDIF
C
C For subroutine calibrat below find the first day; transfer via interface 
C via EC_G_Main
C 
C
      FirstDay = StartTime(1)
C      
C Set Number of uncalibrated signals
C
      Channels = NNNMax
C      
C Set Number of calibrated signals
C
      N = NNMax
C
C Reset Have_Uncal
C
      DO I=1,Channels
          HAVE_UNCAL(I) = .FALSE.
      ENDDO
C
C Read raw data from file
C
      CALL EC_F_ReadNCDF(FName,StartTime,StopTime,
     &  Delay,Gain,Offset,RawSampl,NNNMax,MMMax,M,
     &  NCVarName, HAVE_UNCAL)
      IF (DoPrint.AND.PRaw) THEN
        WRITE(OutFile,*)
        WRITE(OutFile,*) 'Received number of samples = ',M
        WRITE(OutFile,*)
      ENDIF
C
C The mean humidity as given by psychrometer. Will be used to correct
C krypton for drift and dirt. Replace this statement by a read-statement
C or so.
C
      CorMean(Humidity) = Psychro ! [kg m^{-3}]
      IF (M.GT.0) THEN
C
C Call to main routine
C
        CALL EC_G_Main(OutFile,DoPrint,
     &    RawSampl,NNNMax,Channels,NNMax,N,MMMax,M,
     &    DoCrMean,PCal,PIndep,
     &    Psychro,Freq,CalSonic,CalTherm,CalHyg,
     &    CalCO2, P,
     &    Calibrat,
     &    Sample,Flag,Mok,Cok,MIndep,CIndep,Rc,BadTc,
     &    DoDetren,PDetrend,
     &    DoTilt,PTilt,PreYaw,PrePitch,PreRoll,
     &    DoYaw,PYaw,DirYaw,
     &    DoPitch,PPitch,PitchLim,DirPitch,
     &    DoRoll,PRoll,RollLim,DirRoll,
     &    DoSonic,PSonic,SonFactr,
     &    DoO2,PO2,O2Factor,
     &    DoFreq,PFreq,LLimit,ULimit,FrCor,
     &    DoWebb,PWebb,
     &    Mean,TolMean,Cov,TolCov,
     &    HSonic,dHSonic,HTc,dHTc,
     &    LvE,dLvE,LvEWebb,dLvEWebb,
     &    UStar,dUStar,Tau,dTau,
     &    FCO2, dFCO2, FCO2Webb, dFCO2Webb,
     &    MeanW, TolMeanW, HAVE_UNCAL, FirstDay)
C
C Calculate structure parameters
C
        IF (DoStruct) THEN
          R = StructSep
          CALL EC_Ph_Struct(Sample,NNMax,N,MMMax,M,Flag,TSonic ,TSonic ,
     &      R,dR,Freq,CIndep,CTSon2,dCTSon2)
          R = StructSep
          CALL EC_Ph_Struct(Sample,NNMax,N,MMMax,M,Flag,TCouple,TCouple,
     &      R,dR,Freq,CIndep,CTCop2,dCTCop2)
          R = StructSep
          CALL EC_Ph_Struct(Sample,NNMax,N,MMMax,M,Flag,SpecHum,SpecHum,
     &      R,dR,Freq,CIndep,Cq2   ,dCq2   )
          R = StructSep
          CALL EC_Ph_Struct(Sample,NNMax,N,MMMax,M,Flag,TSonic ,SpecHum,
     &      R,dR,Freq,CIndep,CTSonq,dCTSonq)
          R = StructSep
          CALL EC_Ph_Struct(Sample,NNMax,N,MMMax,M,Flag,TCouple,SpecHum,
     &      R,dR,Freq,CIndep,CTCopq,dCTCopq)
        ELSE
           CTSon2 = -9999.D0
          dCTSon2 = -9999.D0
           CTCop2 = -9999.D0
          dCTCop2 = -9999.D0
           Cq2    = -9999.D0
          dCq2    = -9999.D0
           CTSonq = -9999.D0
          dCTSonq = -9999.D0
           CTCopq = -9999.D0
          dCTCopq = -9999.D0
        ENDIF
C
C
C Print fluxes
C
C
        IF (DirYaw.GE.180.D0) THEN
          DirFrom = DirYaw-180.D0
        ELSE
          DirFrom = DirYaw+180.D0
        ENDIF
	IF (Cov(Tsonic,TSonic) .EQ. DUMMY) THEN
	   StdTs = DUMMY
	   dStdTs = DUMMY
	ELSE
	   StdTs = SQRT(Cov(TSonic,TSonic))
	   dStdTs = 0.5D0*TolCov(TSonic,TSonic)/
     &              SQRT(Cov(TSonic,TSonic))
        ENDIF
	IF (Cov(TCouple,TCouple) .EQ. DUMMY) THEN
	   StdTc = DUMMY
	   dStdTc = DUMMY
	ELSE
	   StdTc = SQRT(Cov(TCouple,TCouple))
	   dStdTc = 0.5D0*TolCov(TCouple,TCouple)/
     &              SQRT(Cov(TCouple,TCouple))
        ENDIF
	IF (Cov(SpecHum,SpecHum) .EQ. DUMMY) THEN
	   Stdq = DUMMY
	   dStdq = DUMMY
	ELSE
	   Stdq = SQRT(Cov(SpecHum,SpecHum))
	   dStdq = 0.5D0*TolCov(SpecHum,SpecHum)/
     &              SQRT(Cov(SpecHum,SpecHum))
        ENDIF
	IF (Cov(U,U) .EQ. DUMMY) THEN
	   StdU = DUMMY
	   dStdU = DUMMY
	ELSE
	   StdU = SQRT(Cov(U,U))
	   dStdU = 0.5D0*TolCov(U,U)/
     &              SQRT(Cov(U,U))
        ENDIF
	IF (Cov(V,V) .EQ. DUMMY) THEN
	   StdV = DUMMY
	   dStdV = DUMMY
	ELSE
	   StdV = SQRT(Cov(V,V))
	   dStdV = 0.5D0*TolCov(V,V)/
     &              SQRT(Cov(V,V))
        ENDIF
	IF (Cov(W,W) .EQ. DUMMY) THEN
	   StdW = DUMMY
	   dStdW = DUMMY
	ELSE
	   StdW = SQRT(Cov(W,W))
	   dStdW = 0.5D0*TolCov(W,W)/
     &              SQRT(Cov(W,W))
        ENDIF
        IF (Cov(SpecCO2,SpecCO2) .EQ. DUMMY) THEN
	   StdspecCO2 = DUMMY
	   dStdspecCO2 = DUMMY
	ELSE
	   StdspecCO2 = SQRT(Cov(SpecCO2,SpecCO2))
	   dStdspecCO2 = 0.5D0*TolCov(SpecCO2,SpecCO2)/
     &              SQRT(Cov(SpecCO2,SpecCO2))
        ENDIF
        IF (Cov(CO2,CO2) .EQ. DUMMY) THEN
	   StdCO2 = DUMMY
	   dStdCO2 = DUMMY
	ELSE
	   StdCO2 = SQRT(Cov(CO2,CO2))
	   dStdCO2 = 0.5D0*TolCov(CO2,CO2)/
     &              SQRT(Cov(CO2,CO2))
        ENDIF
	IF ((LvE .EQ. DUMMY) .OR. (LvEWebb .EQ. DUMMY)) THEN
	   SumLvE = DUMMY
	   dSumLvE = DUMMY
	ELSE
	   SumLvE = LvE + LvEWebb
           dSumLvE =  SQRT(dLvE**2.D0 + dLvEWebb**2.D0)
	ENDIF
	IF ((FCO2 .EQ. DUMMY) .OR. (FCO2Webb .EQ. DUMMY)) THEN
	   SumFCO2 = DUMMY
	   dSumFCO2 = DUMMY
	ELSE
	   SumFCO2 = FCO2 + FCO2Webb
           dSumFCO2 =  SQRT(dFCO2**2.D0 + dFCO2Webb**2.D0)
	ENDIF
        WRITE(FluxFile,55)
     &    (NINT(StartTime(i)),i=1,3),
     &    (NINT(StopTime(i)),i=1,3),
     &    M,(Mok(i),i=1,7), Mok(TTime),
     &    NINT(DirFrom),
     &    NINT(2.D0*180.D0*ATAN(SQRT(Cov(V,V))/Mean(U))/Pi),
     &    Mean(U),TolMean(U),
     &    Mean(TSonic),TolMean(TSonic),
     &    Mean(TCouple),TolMean(TCouple),
     &    Mean(SpecHum),TolMean(SpecHum),
     &    StdTs,
     &    dStdTs,
     &    StdTc,
     &    dStdTc,
     &    Stdq,
     &    dStdq,
     &    StdU,
     &    dStdU,
     &    StdV,
     &    dStdV,
     &    StdW,
     &    dStdW,
     &    Cov(TSonic,SpecHum),TolCov(TSonic,SpecHum),
     &    Cov(TCouple,SpecHum),TolCov(TCouple,SpecHum),
     &    Cov(TSonic,U),TolCov(TSonic,U),
     &    Cov(TCouple,U),TolCov(TCouple,U),
     &    Cov(SpecHum,U),TolCov(SpecHum,U),
     &    HSonic,dHSonic,HTc,dHTc,
     &    LvE,dLvE,LvEWebb,dLvEWebb,
     &    SumLvE, dSumLvE,
     &    UStar,dUStar,Tau,dTau,
     &    R,dR,
     &    CTSon2,dCTSon2,CTCop2,dCTCop2,Cq2,dCq2,CTSonq,dCTSonq,
     &    CTCopq,dCTCopq,
     &    MeanW, TolMeanW,
     &    Mok(CO2),
     &    Mean(CO2), TolMean(CO2),
     &    Mean(specCO2), TolMean(specCO2),
     &    StdCO2, dStdCO2,
     &    StdspecCO2, dStdspecCO2,
     &    FCO2,dFCO2,FCO2Webb,dFCO2Webb,
     &    SumFCO2, dSumFCO2


 55     FORMAT(2(I3,1X,2(I2,1X)),9(I6,1X),2(I5,1X),29(2(G14.5:,1X)),
     &        I14,  1X,
     &        7(2(G14.5:,1X)), 3X)
      ELSE
        WRITE(FluxFile,55)
     &    (NINT(StartTime(i)),i=1,3),
     &    (NINT(StopTime(i)),i=1,3),
     &    M,(0,i=1,8),-9999,-9999,(-9999.0,i=1,73)
      ENDIF

      IF (DoPrint) CLOSE(OutFile)

      GOTO 36

 37   CONTINUE
C
C######################################################################
C######################################################################
C
C Here ends the loop over the files
C
C######################################################################
C######################################################################
C
      CLOSE(IntervalFile)
      CLOSE(FluxFile)

      END





C
C ########################################################################
C
C Please check ALL LINES in this procedure and modify where applicable!
C
C Make sure the header is not modified!!!!!!!!
C EC-Pack expects it this way!!!!!!!!!
C
C Input : RawSampl : REAL*8(Channels) : Raw, uncalibrated sample from file
C         P        : Atmospheric pressure [Pa]
C         CorMean  : REAL*8(N) : Correction terms for calibrated quantities
C                    especially implemented to allow for correction of
C                    signal of dirty krypton hygrometer with psychrometer
C         CalSonic : Calibration array REAL*8(NQQ) of the anemometer
C         CalTherm : Calibration array REAL*8(NQQ) of the thermometer
C         CalHyg   : Calibration array REAL*8(NQQ) of the hygrometer
C         HAVE_TREF: do we have a reference temperature for thermocouple
C                    (if so apply calibration curve and add reference
C                    temperature)
C OutPut: Sample   : REAL*8(N) : Calibrated sample
C         Error    : LOGICAL(N) : Indicator if quantities are valid
C                    after calibration
C
C Revision 28-05-2001: added info on which variables are available
C                      (mainly important for either sonic or
C                      thermocouple temperature)
C Revision 18-09-2002: added FirstDay to interface, removed calccomm.inc
C                      added CO2 calibration
C ########################################################################
C
      SUBROUTINE Calibrat(RawSampl,Channels,P,CorMean,
     &  CalSonic,CalTherm,CalHyg,CalCO2, BadTc,Sample,N,Error,
     &  Have_Uncal,FirstDay)

      IMPLICIT NONE

      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      INTEGER N,i, FirstDay, Channels
      LOGICAL Error(N),BadTc, Have_Uncal(N)
      REAL*8 RawSampl(Channels),Sample(N),P,UDum,VDum,Hook,Dum,
     &  CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ), CalCO2(NQQ),
     &  CorMean(N),
     &  Hours,Minutes,Days,Secnds, TsCorr, WDum
      REAL*8 EC_PH_Q,EC_M_BaseF, EC_Ph_QCO2 ! External function calls

      REAL*8 C(0:NMaxOrder)
      INTEGER IFLG, DUMORDER, DUMTYP
      REAL*8 ABSERROR, RELERROR, ESTIM, DNLIMIT, UPLIMIT

      REAL*8 DUMFUNC, EC_C_Schot3
      EXTERNAL DUMFUNC, EC_C_Schot3
C
C To communicate with function of which zero is searched
C
      COMMON /DUMCOM/ C, DUM, DUMORDER, DUMTYP
C
C In principle no apparatus has shown to give a valid signal.
C
      DO i=1,N
        Error(i) = (.TRUE.)
      ENDDO
C
C Calculate time from the columns in the raw data
C
      Days = RawSampl(Doy) - DBLE(FirstDay)
      Hours = INT(0.01D0*RawSampl(HourMin))
      Minutes = RawSampl(HourMin) - 100.D0*Hours
      Secnds = RawSampl(Sec)
      Sample(TTime) =
     &  Secnds + 60.D0*(Minutes +60.D0*(Hours + 24.D0*Days))

      Error(TTime) = (.FALSE.)
C
C Take the velocity from the raw data (again assume it's a 3D sonic)
C
      IF (HAVE_UNCAL(U) .AND. 
     &    HAVE_UNCAL(V) .AND. 
     &    HAVE_UNCAL(W)) THEN
         UDum  =  RawSampl(U) ! U [m/s]
         VDum  =  RawSampl(V) ! V [m/s]
         Sample(W) = RawSampl(W) ! W [m/s]

         Error(U) = (ABS(UDum).GT.40.D0)
         Error(V) = (ABS(VDum).GT.40.D0)
         Error(W) = (ABS(Sample(W)).GT.20.D0)
C
C If this is a wind tunnel calibrated sonic,
C we can apply the KNMI calibrations
C
         IF (CalSonic(QQType) .EQ. ApSon3Dcal) THEN
            WDum = Sample(W)
            CALL EC_C_Scal(CalSonic,
     &                  UDUM, VDUM, WDum,
     &                  ERROR(U), ERROR(V), ERROR(W))
            Sample(W) = WDum
         ENDIF

C This is the construction when we have a diagnostic variable
         IF (Have_Uncal(Diagnost)) THEN
            IF ((Error(U).OR.Error(V).OR.Error(W)).OR.
     &          ((RawSampl(Diagnost).LT.0.D0).OR.
     &           (RawSampl(Diagnost).GT.100.D0))) THEN
               Error(U) = (.TRUE.)
               Error(V) = (.TRUE.)
               Error(W) = (.TRUE.)
            ENDIF
C This is the construction when we do not have a diagnostic variable
         ELSE 
	    IF ((Error(U).OR.Error(V)).OR.Error(W)) THEN
               Error(U) = (.TRUE.)
               Error(V) = (.TRUE.)
               Error(W) = (.TRUE.)
             ENDIF
	 ENDIF

         IF (.NOT.((Error(U).OR.Error(V)).OR.Error(W))) THEN

C
C Rotate velocity according to angle of setup's north relative real north
C Rotation angle as given in calibration file should be taken negative
C
C Sonic takes flow that blows into the sonic as positive. To get
C the wind speed with wind to north as positive substract 180 degrees:
C
C Apparently V is defined other way around (explains -VDum)
C
           Hook = PI*(-CalSonic(QQYaw)+180)/180.D0
           Sample(U) =  COS(Hook)*UDum + SIN(Hook)*(-VDum)
           Sample(V) = -SIN(Hook)*UDum + COS(Hook)*(-VDum)
C
C Kaijos have different coordinate system
C
           IF ((CalSonic(QQType) .EQ. ApKaijoTR90) .OR.
     &         (CalSonic(QQType) .EQ. ApKaijoTR61)) THEN
               Hook = -90.D0 * PI/180.D0
               UDum  =  Sample(U) ! U [m/s]
               VDum  =  Sample(V) ! V [m/s]
               Sample(U) =  COS(Hook)*UDum + SIN(Hook)*VDum
               Sample(V) = -SIN(Hook)*UDum + COS(Hook)*VDum
           ENDIF

           Sample(TSonic) = RawSampl(TSonic) + Kelvin
C
C Here the Schotanus et al. correction for sensitivity of the sonic for
C lateral velocity is applied directly to the raw data. This is only
C half of the Schotanus-correction. Humidity-correction comes later.
C
           Sample(TSonic) = Sample(TSonic)
     &         + ((Sample(U))**2 + (Sample(V))**2)/GammaR
C
C Following correction is implemented to correct for the true length of the
C sonic, which is estimated from the ratio of mean sonic temperatures and
C mean thermocouple temperatures.
C
C
C 5 august: this part is made inactive by setting factor to 1
C
           Error(Tsonic) = ((Sample(TSonic).GT.(MaxT))
     &       .OR. (Sample(TSonic).LT.(MinT)))
         ENDIF
	 Error(U) = (Error(U) .OR. (.NOT. Have_Uncal(U)))
	 Error(V) = (Error(V) .OR. (.NOT. Have_Uncal(V)))
	 Error(W) = (Error(W) .OR. (.NOT. Have_Uncal(W)))
	 Error(TSonic) = (Error(TSonic) .OR. (.NOT. Have_Uncal(TSonic)))
      ELSE
         Error(U) = .TRUE.
         Error(V) = .TRUE.
         Error(W) = .TRUE.
         Error(Tsonic) = .TRUE.
      ENDIF
C
C Calibrate thermocouple
C
      IF (HAVE_UNCAL(TCouple)) THEN
         IF (HAVE_UNCAL(Tref)) THEN

C
C This is calibration according to Campbells P14 instruction
C      T =  Calibration(voltage +
C                       inverse calibration(reference temperature))
C Reference temperature is supposed to be in Celcius !!
C
            DO i=0,CalTherm(QQOrder)
              c(i) = CalTherm(QQC0+i)
            ENDDO
            Dum = RawSampl(Tref) 
            Error(TCouple) =
     &           (((DUM) .GT.(MaxT))
     &             .OR. ((DUM) .LT.(MinT)))
            IF (.NOT. Error(TCouple)) THEN
               UPLIMIT = (DUM + 3.0 - CalTherm(QQC0))/CalTherm(QQC1)
               DNLIMIT = (DUM - 3.0 - CalTherm(QQC0))/CalTherm(QQC1)
               ESTIM = (Dum - CalTherm(QQC0))/CalTherm(QQC1)
               RELERROR = 1e-4
               ABSERROR = (0.001)/CalTherm(QQC1)
               DUMORDER = CalTherm(QQOrder)
               DUMTYP = CalTherm(QQFunc)
               CALL DFZERO(DUMFUNC, DNLIMIT, UPLIMIT, ESTIM, RELERROR,
     &                    ABSERROR, IFLG)
               IF (IFLG .GT. 3) THEN
	                ERROR(TCouple) = .TRUE.
	            ENDIF
               Sample(Tcouple) = Kelvin +
     &             EC_M_BaseF(DNLIMIT + RawSampl(Tcouple),
     &              NINT(CalTherm(QQFunc)),
     &              NINT(CalTherm(QQOrder)),c)
            ENDIF
         ELSE
C
C Suppose that sample is already temperature (now in Kelvin, Sept. 18, 2002!)
C
             Sample(TCouple) = RawSampl(TCouple) + Kelvin
         ENDIF

         Error(TCouple) = ((Sample(TCouple).GT.(Kelvin+MaxT))
     &     .OR. (Sample(TCouple).LT.(Kelvin+MinT)))
         Error(TCouple) = (Error(TCouple) .OR.
     &                        (.NOT. Have_Uncal(TCouple)))
      ELSE
         Error(Tcouple) = .TRUE.
      ENDIF
C
C Calibrate hygrometer
C
      IF (HAVE_UNCAL(Humidity)) THEN
         Dum = RawSampl(Humidity)
         IF (Dum .GE. Epsilon) THEN
C
C Add an optional correction to the krypton's humidity to compensate for
C drift in the hygrometer (e.g. dirt). Compensation may be provided by
C the difference between the mean humidity by a (slow) psychrometer's
C output and the mean of the humidity estimated by the hygrometer.
C
            DO i=0,CalHyg(QQOrder)
              c(i) = CalHyg(QQC0+i)
            ENDDO
C Now in kg/m^3 (Sept. 18, 2002)
            Sample(Humidity) = CorMean(Humidity) +
     &         EC_M_BaseF(Dum,NINT(CalHyg(QQFunc)),
     &           NINT(CalHyg(QQOrder)),c)

            Error(Humidity) = ((Sample(Humidity).GT.MaxRhoV)
     &          .OR. (Sample(Humidity).LT.MinRhoV))
            Error(SpecHum) = Error(Humidity)
C
C Calculate specific humidity associated with this absolute humidity
C
            IF (.NOT.Error(Humidity)) THEN
              IF (.NOT.BadTc) THEN
                IF (.NOT.Error(TCouple)) THEN
                  Sample(SpecHum)=EC_PH_Q(Sample(Humidity),
     &                                Sample(TCouple),P)
                ELSE IF (.NOT.Error(TSonic)) THEN
                  TsCorr = EC_C_Schot3(Sample(TSonic),
     &                                Sample(Humidity), P)
                  Sample(SpecHum)=EC_PH_Q(Sample(Humidity),
     &                                TsCorr,P)
                ELSE
                  Error(SpecHum) = (.TRUE.)
                ENDIF
              ELSE
                IF (.NOT.Error(TSonic)) THEN
                  TsCorr = EC_C_Schot3(Sample(TSonic),
     &                                Sample(Humidity), P)
                  Sample(SpecHum)=EC_PH_Q(Sample(Humidity),
     &                                TsCorr,P)
                ELSE
                  Error(SpecHum) = (.TRUE.)
                ENDIF
              ENDIF
            ENDIF
         ENDIF
         Error(Humidity) = (Error(Humidity) .OR.
     &                        (.NOT. Have_Uncal(Humidity)))
         Error(SpecHum) = (Error(SpecHum) .OR. Error(Humidity))
      ELSE
         Error(Humidity) = .TRUE.
         Error(SpecHum) = .TRUE.
      ENDIF
C
C Get CO2 sample
C
      IF (HAVE_UNCAL(CO2)) THEN
         DO i=0,CalCO2(QQOrder)
              c(i) = CalCO2(QQC0+i)
         ENDDO
	 Sample(CO2) = CorMean(CO2) + 
     &         EC_M_BaseF(RawSampl(CO2),NINT(CalCO2(QQFunc)),
     &           NINT(CalCO2(QQOrder)),c)
         Error(CO2) = ((Sample(CO2).GT.MaxRhoCO2)
     &          .OR. (Sample(CO2).LT.MinRhoCO2))
C
C Calculate specific CO2 associated with this absolute CO2 concentration
C
         IF (.NOT.Error(CO2)) THEN
              IF (.NOT.BadTc) THEN
                IF (.NOT.Error(TCouple) .AND.
     &              .NOT.Error(Humidity)) THEN
                  Sample(SpecCO2)=EC_PH_QCO2(
     &                                Sample(CO2),
     &                                Sample(Humidity),
     &                                Sample(TCouple),P)
                ELSE IF (.NOT.Error(TSonic) .AND.
     &                   .NOT.Error(Humidity)) THEN
                  TsCorr = EC_C_Schot3(Sample(TSonic),
     &                                Sample(Humidity), P)
                  Sample(SpecCO2)=EC_PH_QCO2(
     &                                Sample(CO2),
     &                                Sample(Humidity),
     &                                TsCorr,P)
                  Error(SpecCO2) = .FALSE.
                ELSE
                  Error(SpecCO2) = (.TRUE.)
                ENDIF
              ELSE
                IF (.NOT.Error(TSonic) .AND.
     &              .NOT.Error(Humidity)) THEN
                  TsCorr = EC_C_Schot3(Sample(TSonic),
     &                                Sample(Humidity), P)
                  Sample(SpecCO2)=EC_PH_QCO2(
     &                                Sample(CO2),
     &                                Sample(Humidity),
     &                                TsCorr,P)
                  Error(SpecCO2) = .FALSE.
                ELSE
                  Error(SpecCO2) = (.TRUE.)
                ENDIF
              ENDIF
         ENDIF
         Error(CO2) = (Error(CO2) .OR.
     &                        (.NOT. Have_Uncal(CO2)))
         Error(SpecCO2) = (Error(SpecCO2) .OR. Error(CO2))
      ELSE
         Error(CO2) = .TRUE.
         Error(SpecCO2) = .TRUE.
      ENDIF
C
C Time is used for detrending. Therefore we take out all samples with
C a defect in one of the velocities.
C I wouldnt know why so, test has been removed AM
C
      Error(TTime) = (.FALSE.)
C      DO i=U,W
C        Error(TTime) = (Error(TTime).OR.Error(i))
C      ENDDO

      RETURN
      END


C
C Function: DUMFUNC
C Purpose : to return value of calibration function minus the y-value
C           in order to be used in DFZERO (root finding routine)
C
      REAL*8 FUNCTION DUMFUNC(X)

      IMPLICIT NONE
      INCLUDE 'parcnst.inc'

      REAL*8 X
      REAL*8 C(0:NMaxOrder), DUM
      INTEGER DUMORDER, DUMTYP

      REAL*8 EC_M_BaseF
      EXTERNAL EC_M_BaseF

C
C To communicate with function of which zero is searched
C
      COMMON /DUMCOM/ C, DUM, DUMORDER, DUMTYP

      DUMFUNC = EC_M_BaseF(X, DUMTYP, DUMORDER, C) - DUM
      END




