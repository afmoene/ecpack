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
     &              ParmName,InterName, PlfName, PlfIntName
      CHARACTER*255 OutName,DumName1,DumName2
      LOGICAL PRaw,PYaw,PPitch,PRoll,PFreq,PO2,PWebb,PCal,PIndep,
     &  DoYaw,DoPitch,DoRoll,DoFreq,DoO2,DoWebb,DoCrMean,PSonic,
     &  DoSonic,DoTilt,PTilt,DoPrint,Flag(NNMAx,MMMax),DoDetren,
     &  PDetrend,DoStruct,BadTc, DoPF, PPf
      INTEGER N,i,M,MIndep(NNMax),CIndep(NNMax,NNMax),FOO,
     &  Channels,Delay(NNNMax),Mok(NNMax),Cok(NNMax,NNMax), FirstDay,
     &  NDelta, NLock, NHigh, NLow, StartTime(3),StopTime(3)
      REAL*8 RawSampl(NNNMax,MMMax),Sample(NNMax,MMMax),P,Psychro,
     &  Mean(NNMax),TolMean(NNMax),Cov(NNMax,NNMax),
     &  TolCov(NNMax,NNMax),LLimit,ULimit,FrCor(NNMax,NNMax),
     &  PitchLim,RollLim,DirYaw,RC(NNMax),O2Factor(NNMax),
     &  Freq,DirPitch,DirRoll,
     &  SonFactr(NNMax),PreYaw,PrePitch,PreRoll,
     &  HSonic,dHSonic,HTc,dHTc,
     &  LvE,dLvE,LvEWebb,dLvEWebb,
     &  UStar,dUStar,Tau,dTau,CorMean(NNMax),DirFrom,
     &  Gain(NNNMax),Offset(NNNMax),
     &  StructSep,
     &  FCO2,dFCO2, FCO2Webb, dFCO2Webb,
     &  Apf(3,3), Alpha, Beta, Gamma, WBias,
     &  CosGamma, SinGamma, Yaw(3,3), InvYaw(3,3), PFValid
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
      LOGICAL       HAVE_UNCAL(NNNMax), AnglesFound

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
     &             FluxName, ParmName, InterName, PlfIntName,
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
      WRITE(FluxFile,57)
 56   FORMAT(
     &  ' DOY Hr Mn ',
     &  'DOY Hr Mn ',
     &  '#samples',
     &  '  #U     #V     #W  #TSon  #TCop  #RhoV     #q  #time   ',
     &  'Dir d(dir)   ',
     &  'Mean(U)        dMean(U)         ',
     &  'Mean(TSon)     dMean(TSon)      ',
     &  'Mean(TCop)     dMean(TCop)      ',
     &  'Mean(q)        dMean(q)         ',
     &  'sigma(TSon)    dsigma(TSon)     ',
     &  'sigma(TCop)    dsigma(TCop)     ',
     &  'sigma(q)       dsigma(q)        ',
     &  'sigma(u)       dsigma(u)        ',
     &  'sigma(v)       dsigma(v)        ',
     &  'sigma(w)       dsigma(w)        ',
     &  'cov(TSon,q)    dcov(TSon,q)     ',
     &  'cov(TCop,q)    dcov(TCop,q)     ',
     &  'cov(TSon,U)    dcov(TSon,U)     ',
     &  'cov(TCop,U)    dcov(TCop,U)     ',
     &  'cov(q,U)       dcov(q,U)        ',
     &  'H(Sonic)       Tol(HSonic)      ',
     &  'H(TCouple)     Tol(HTCouple     ',
     &  'LvECov         dLvECov          ',
     &  'LvEWebb        dLvEWebb         ',
     &  'LvE            dLvE             ',
     &  'UStar          dUStar           ',
     &  'Tau            dTau             ',
     &  'R(delay)       dR               ',
     &  'CTSon2         dCTSon2          ',
     &  'CTCop2         dCTCop2          ',
     &  'Cq2            dCq2             ',
     &  'CTSonq         dCTSonq          ',
     &  'CTCopq         dCTCopq          ',
     &  'MeanW          dMeanW           ',
     &  '#CO2           ',
     &  'MeanCO2        dMeanCO2         ',
     &  'MeanspecCO2    dMeanspecCO2     ',
     &  'stdCO2         dstdCO2          ',
     &  'stdspecCO2     dstdspecCO2      ',
     &  'FCO2Cov        dFCO2Cov         ',
     &  'FCO2Webb       dFCO2Webb        ',
     &  'FCO2           dFCO2            ',
     &  '#DeltaTemp     #PoorSignalLock  ',
     &  '#AmplHigh      #AmplLow        '
     &  )
  57   FORMAT(     
     &  '  -  -  - ',
     &  ' -  -  - ',    
     &  '[-]     ',
     &  '  [-]    [-]    [-] [-]    [-]    [-]       [-] [-]     ',
     &  '[degr] [degr]',
     &  ' [m/s]         [m/s]            ',
     &  '[K]            [K]              ',
     &  '[K]            [K]              ',
     &  '[kg/kg]        [kg/kg]          ',
     &  '[K]            [K]              ',
     &  '[K]            [K]              ',
     &  '[kg/kg]        [kg/kg]          ',
     &  '[m/s]          [m/s]            ',
     &  '[m/s]          [m/s]            ',
     &  '[m/s]          [m/s]            ',
     &  '[K*kg/kg]      [K*kg/kg]        ',
     &  '[K*kg/kg]      [K*kg/kg]        ',
     &  '[K*m/s]        [K*m/s]          ',
     &  '[K*m/s]        [K*m/s]          ',
     &  '[kg/kg*m/s]    [kg/kg*m/s]      ',
     &  '[W/m^2]        [W/m^2]          ',
     &  '[W/m^2]        [W/m^2]          ',
     &  '[W/m^2]        [W/m^2]          ',
     &  '[W/m^2]        [W/m^2]          ',
     &  '[W/m^2]        [W/m^2]          ',
     &  '[m/s]          [m/s]            ',
     &  '[N/m^2]        [N/m^2]          ',
     &  '[m]            [m]              ',
     &  '[K^2/m^2/3]    [K^2/m^2/3]      ',
     &  '[K^2/m^2/3]    [K^2/m^2/3]      ',
     &  '[1/m^2/3]      [1/m^2/3]        ',
     &  '[K*kg/kg/m^2/3] [K*kg/kg/m^2/3] ',
     &  '[K*kg/kg/m^2/3] [K*kg/kg/m^2/3] ',
     &  '[m/s]          [m/s]            ',
     &  '[-]            ',
     &  '[kg/m^3]       [kg/m^3]         ',
     &  '[kg/kg]        [kg/kg]          ',
     &  '[kg/m^3]       [kg/m^3]         ',
     &  '[kg/kg]        [kg/kg]          ',
     &  '[kg/m^2/s^1]   [kg/m^2/s^1]     ',
     &  '[kg/m^2/s^1]   [kg/m^2/s^1]     ',
     &  '[kg/m^2/s^1]   [kg/m^2/s^1]     ',
     &  '[-]            [-]              ',
     &  '[-]            [-]              ')     
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
     &  DoTilt,DoYaw,DoPitch,DoRoll,DoFreq,DoO2,DoWebb,DoStruct, DoPF,
     &  DoPrint,
     &  PRaw,PCal,PDetrend,PIndep,PTilt,PYaw,PPitch,
     &  PRoll,PSonic,PO2,PFreq,PWebb, PPF, PFValid, StructSep)
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
           IF ((CalSonic(QQType) .NE. ApCSATSonic) .AND.
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
      WRITE(*,951) FName(:EC_T_STRLEN(Fname)),(StartTime(i),i=1,3)
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
C
C If planar fit tilt correction , try to get rotation matrix from file
C
      IF (DoPF) THEN
        CALL EC_F_GetPF(PlfName, StartTime, StopTime, 
     &                  Alpha, Beta, Gamma, 
     &                  WBias, Apf, AnglesFound)
        IF (AnglesFound) THEN
C Undo Yaw angle of planar fit untilt matrix
           CosGamma = COS(GAMMA*pi/180)
           SinGamma = SIN(GAMMA*pi/180)
           Yaw(1,1) = CosGamma
           Yaw(1,2) = SinGamma
           Yaw(1,3) = 0.D0
           Yaw(2,1) = -SinGamma
           Yaw(2,2) = CosGamma
           Yaw(2,3) = 0.D0
           Yaw(3,1) = 0.D0
           Yaw(3,2) = 0.D0
           Yaw(3,3) = 1.D0
	   CALL EC_M_InvM(Yaw,InvYaw)
	   CALL EC_M_MMul(InvYaw, Apf, Apf)
        ELSE
           M = 0
        ENDIF
      ENDIF
      
      
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
     &    DoPF, PPF, Apf, 
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
           IF (.NOT. DoWebb) THEN
	     SumLvE = LvE
	     dSumLvE = dLvE
	   ELSE
	     SumLvE = DUMMY
	     dSumLvE = DUMMY
           ENDIF
	ELSE
	   SumLvE = LvE + LvEWebb
           dSumLvE =  SQRT(dLvE**2.D0 + dLvEWebb**2.D0)
	ENDIF
        IF ((FCO2 .EQ. DUMMY) .OR. (FCO2Webb .EQ. DUMMY)) THEN
          IF (.NOT. DoWebb) THEN
	     SumFCO2 = FCO2
	     dSumFCO2 = dFCO2
          ELSE
             SumFCO2 = DUMMY
             dSumFCO2 = DUMMY
          ENDIF
	ELSE
	   SumFCO2 = FCO2 + FCO2Webb
           dSumFCO2 =  SQRT(dFCO2**2.D0 + dFCO2Webb**2.D0)
	ENDIF
C
C Count error codes of CSAT sonic        
C
        IF (CalSonic(QQType) .EQ. ApCSATSonic) THEN
           NDelta = 0
           NLock = 0
           NHigh = 0
           NLow = 0
           DO I=1,M
              NDelta = NDelta + INT(Sample(DiagDelta, I))
              NLock = NLock + INT(Sample(DiagLock, I))
              NHigh = NHigh + INT(Sample(DiagHigh, I))
              NLow = NLow + INT(Sample(DiagLow, I))
           ENDDO
        ENDIF
        
        WRITE(FluxFile,55)
     &    (StartTime(i),i=1,3),
     &    (StopTime(i),i=1,3),
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
     &    SumFCO2, dSumFCO2,
     &    NDelta, NLock, NHigh, NLow
        
  
 55     FORMAT(2(I3,1X,2(I2,1X)),9(I6,1X),2(I5,1X),29(2(G15.5:,1X)),
     &        I13,  1X,
     &        7(2(G15.5:,1X)), 4(I15,1X))
      ELSE
        WRITE(FluxFile,55)
     &    (StartTime(i),i=1,3),
     &    (StopTime(i),i=1,3),
     &    M,(0,i=1,8), INT(DUMMY),INT(DUMMY),(DUMMY,i=1,58),
     &    0, (DUMMY,i=1,14), (INT(DUMMY), i=1,4)
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








