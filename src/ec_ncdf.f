C $Id$
C  ec_ncdf.f, part of ecpack
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
C             For more recent changes, see Changelog
C Current version: see version.inc
C...........................................................................
      PROGRAM EC_NCDF

      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'
      INCLUDE 'version.inc'

      CHARACTER*255 FNAME, DatDir, OutDir, ParmDir, FluxName,
     &              ParmName,InterName, PlfName, PlfIntName,
     &              OutName,DumName1,DumName2,
     &              SonName,CoupName,HygName, CO2Name, DUMSTRING,
     &              NCvarname(NNNMax)
      LOGICAL PRaw, PCal,PIndep,
     &  DoPrint,Flag(NNMAx,MMMax),
     &  DoStruct,BadTc,
     &  HAVE_UNCAL(NNNMax), HAVE_CAL(NNMax),
     &  AnglesFound, 
     &  DoCorr(NMaxCorr), PCorr(NMaxCorr),
     &  OutMean(NNMax), OutCov(NNMax, NNMax), 
     &  OutStd(NNMax), OutNum(NNMax), OutStr(NNMax, NNMax), 
     &  OutPh(NMaxPhys), Outputs(NMaxOS)     
      INTEGER N,i,J,M,MIndep(NNMax),CIndep(NNMax,NNMax),FOO,
     &  Channels,Delay(NNNMax),Mok(NNMax),Cok(NNMax,NNMax), 
     &  NDelta, NLock, NHigh, NLow, StartTime(3),StopTime(3),
     &  DiagFlag(NMaxDiag)
      REAL*8 RawSampl(NNNMax,MMMax),Sample(NNMax,MMMax),P,Psychro,
     &  Mean(NNMax),TolMean(NNMax),Cov(NNMax,NNMax),
     &  TolCov(NNMax,NNMax),FrCor(NNMax,NNMax),
     &  DirYaw,RC(NNMax),O2Factor(NNMax),
     &  DirPitch,DirRoll,
     &  SonFactr(NNMax),
     &  Gain(NNNMax),Offset(NNNMax),
     &  Apf(3,3), Alpha, Beta, Gamma, WBias,
     &  CosGamma, SinGamma, Yaw(3,3), InvYaw(3,3),
     &  ExpVar(NMaxExp), CorrPar(NMaxCorrPar),
     &  CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ), CalCO2(NQQ),
     &  R,dR,CTSon2,CTCop2,Cq2,CTSonq,CTCopq,
     &  dCTSon2,dCTCop2,dCq2,dCTSonq,dCTCopq,
     &  Qphys(NMaxPhys), dQPhys(NMaxPHys),
     &  Std(NNMax), dStd(NNMax),
     &  Struct(NNMax,NNMax),
     &  dStruct(NNMax,NNMax)

      INTEGER EC_T_STRLEN
C
C This may not be portable: works on gnu-fortran
C
C      CHARACTER*(*) FDATE
C
C The name of the calibration routine specified at the end of this file
C
      EXTERNAL Calibrat, EC_T_STRLEN
C
C This may not be portable: works on gnu-fortran
C
      INTRINSIC FDATE

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

C Get defaults for configuration 
      Outputs(OSOld) = .TRUE.
      Outputs(OSTol) = .TRUE.
      Outputs(OScsv) = .FALSE.
      Outputs(OSDiag) = .FALSE.
      Outputs(OSVers) = .TRUE.
      DO I=1,NNMax
         OutMean(I) = .FALSE.
         OutStd(I) = .FALSE.
         OutNum(I) = .FALSE.
         DO J=1,NNMax
            OutCov(I,J) = .FALSE.
            OutStr(I,J) = .FALSE.
         ENDDO
      ENDDO
      DO I=1,NMaxPhys
         OutPh(I) = .FALSE.
      ENDDO
C Get configuration (file names etc.)
      CALL EC_F_GetConf(ECConfFile,
     &             DatDir, OutDir, ParmDir,
     &             FluxName, ParmName, InterName, PlfIntName,
     &             PlfName,
     &             SonName, CoupName, HygName, CO2Name,
     &             NCVarname, NNNMax,
     &             OutMean, OutCov, OutPh, OutStd, OutNum, OutStr,
     &             Outputs, DoCorr, ExpVar, CorrPar,
     &             DoStruct, DoPrint,
     &             PCorr, PRaw, PCal, PIndep)

C
C Assume first we have no uncalibrated samples at all
C
      DO I=1,NNNMax
         HAVE_UNCAL(I) = .FALSE.
      ENDDO
C
C Reset calibration information
C
      DO I=1,NQQ
        CalSonic(I) = DUMMY
        CalHyg(I) = DUMMY
        CalTherm(I) = DUMMY
        CalCO2(I) = DUMMY
      ENDDO
C
C Check whether we have a reference temperature for a thermocouple
      IF (EC_T_STRLEN(NCVarName(QUTref)) .GT. 0) THEN
         HAVE_UNCAL(QUTref) = .TRUE.
      ENDIF
C
C Check which calibration files we have
C
      IF (EC_T_STRLEN(SonName) .GT. 0) THEN
C For now assume it is always a 3-D sonic
         HAVE_UNCAL(QUU) = .TRUE.
         HAVE_UNCAL(QUV) = .TRUE.
         HAVE_UNCAL(QUW) = .TRUE.
         HAVE_UNCAL(QUTSonic) = .TRUE.
      ENDIF
      IF (EC_T_STRLEN(HygName) .GT. 0) THEN
         HAVE_UNCAL(QUHumidity) = .TRUE.
      ENDIF
      IF (EC_T_STRLEN(CoupName) .GT. 0) THEN
         HAVE_UNCAL(QUTCouple) = .TRUE.
      ENDIF
      IF (EC_T_STRLEN(CO2Name) .GT. 0) THEN
         HAVE_UNCAL(QUCO2) = .TRUE.
      ENDIF
C
C Open an output file for the fluxes and write headers to it
C
      OPEN(FluxFile,FILE=FluxName,FORM='FORMATTED')
C
C The call to FDATE may not be portable: works on gnu-fortran.
C just remove the call, if you have problems. Or find something
C better or more portable.
C
      IF (Outputs(OSVers)) THEN
         WRITE(FluxFile,'(A, A, A, A)') 'EC_NCDF Version: ', 
     &                     VERSION, ' run at ', FDATE()
      ENDIF
      CALL  EC_F_WFlux(FluxFile, .TRUE.,
     &                StartTime, StopTime,
     &                M, Mok, Mean, TolMean, Cov, TolCov,
     &                QPhys, dQPhys, Std, dStd, Struct,
     &                dStruct, R, dR, DiagFlag,
     &                OutMean, OutCov, OutPh, 
     &                OutStd, OutNum, OutStr, Outputs)
      
C
C Number of quantities involved in this experiment
C
      N = NNMax  ! Number of calibrated quantities you will follow
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
C NOTE: this mapping is NOT done in the calibration step!
C

      Delay(QUU)        = NINT(CalSonic(QQDelay)*
     &                         0.001D0*ExpVar(QEFreq))
      Delay(QUV)        = NINT(CalSonic(QQDelay)*
     &                         0.001D0*ExpVar(QEFreq))
      Delay(QUW)        = NINT(CalSonic(QQDelay)*
     &                         0.001D0*ExpVar(QEFreq))
      Delay(QUHumidity) = NINT(CalHyg(QQDelay)*
     &                         0.001D0*ExpVar(QEFreq))
      Delay(QUTcouple)  = NINT(CalTherm(QQDelay)*
     &                         0.001D0*ExpVar(QEFreq))
      Delay(QUTsonic)   = NINT(CalSonic(QQDelay)*
     &                         0.001D0*ExpVar(QEFreq))
      Delay(QUCO2)      = NINT(CalCO2(QQDelay)*
     &                         0.001D0*ExpVar(QEFreq))

      Gain(QUU)        = CalSonic(QQGain)
      Gain(QUV)        = CalSonic(QQGain)
      Gain(QUW)        = CalSonic(QQGain)
      Gain(QUHumidity) = CalHyg(QQGain)
      Gain(QUTcouple)  = CalTherm(QQGain)
      Gain(QUTsonic)   = CalSonic(QQGain)
      Gain(QUCO2)      = CalCO2(QQGain)

      Offset(QUU)        = CalSonic(QQOffset)
      Offset(QUV)        = CalSonic(QQOffset)
      Offset(QUW)        = CalSonic(QQOffset)
      Offset(QUHumidity) = CalHyg(QQOffset)
      Offset(QUTcouple)  = CalTherm(QQOffset)
      Offset(QUTsonic)   = CalSonic(QQOffset)
      Offset(QUCO2)      = CalCO2(QQOffset)

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
C If planar fit tilt correction , try to get rotation matrix from file
C
      IF (DoCorr(QCPF)) THEN
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
     &    PCal,PIndep, Psychro,CalSonic,CalTherm,CalHyg,
     &    CalCO2, P, Calibrat,
     &    Sample,Flag,Mok,Cok,MIndep,CIndep,Rc,BadTc,
     &    DoCorr, PCorr, ExpVar, CorrPar,
     &    DirYaw,
     &    DirPitch,
     &    DirRoll,
     &    Apf, 
     &    SonFactr,
     &    O2Factor,
     &    FrCor,
     &    Mean,TolMean,Cov,TolCov,
     &    QPhys, dQPhys,
     &    HAVE_UNCAL, HAVE_CAL, DiagFlag, StartTime(1))
C
C Calculate structure parameters
C
        IF (DoStruct) THEN
          DO I=1,NNMax
	     DO J=1,NNMax
	        IF ((HAVE_CAL(I) .AND. HAVE_CAL(J)) .AND.
     &              OutStr(I,J)) THEN
	           R = ExpVar(QEStructSep)
                   CALL EC_Ph_Struct(Sample,NNMax,MMMax,M,Flag,
     &                  I,J, R,dR,ExpVar(QEFreq),CIndep,
     &                  Struct(I,J), dStruct(I,J))
                   Struct(J,I) = Struct(I,J)
                   dStruct(J,I) = dStruct(I,J)
                ELSE
                   Struct(I,J) = DUMMY
                   dStruct(I,J) = DUMMY
                ENDIF
             ENDDO
          ENDDO
        ELSE
          DO I=1,NNMax
	     DO J=1,NNMax
                R = DUMMY
                dR = DUMMY
                Struct(I,J) = DUMMY
                dStruct(I,J) = DUMMY
             ENDDO
          ENDDO
        ENDIF
C
C Determine standard deviations
C
        DO I=1,NNMax
	   IF (HAVE_CAL(I) .AND. HAVE_CAL(J)) THEN
	     Std(I) = DUMMY
	     dStd(I) = DUMMY
	   ELSE
	     Std(I) = SQRT(Cov(I,I))
	     dStd(I) = 0.5D0*TolCov(I,I)/
     &              SQRT(Cov(I,I))
           ENDIF
        ENDDO
C
C Print fluxes
C     
      ENDIF 
      CALL  EC_F_WFlux(FluxFile, .FALSE.,
     &                StartTime, StopTime,
     &                M, Mok, Mean, TolMean, Cov, TolCov,
     &                QPhys, dQPhys, Std, dStd, Struct,
     &                dStruct, R, dR, DiagFlag,
     &                OutMean, OutCov, OutPh, 
     &                OutStd, OutNum, OutStr, Outputs)

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








