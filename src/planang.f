C $Id$
C  planang.f, part of ecpack
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
C Function  : planang
C Purpose   : To compute planar fit angles for periods that are read from
C             a file.
C Interface : a configuration file is read from standard input
C Author    : Arnold Moene (but based on ec_ncdf)
C Date      : December 5, 2002
C...........................................................................
      PROGRAM PLANANG

      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'
      INCLUDE 'version.inc'

      INTEGER MaxPF, NPF
      PARAMETER (MaxPF = 100000)
      CHARACTER*255 FNAME, DatDir, OutDir, ParmDir, FluxName,
     &              ParmName,InterName, PlfName, PlfIntName
      CHARACTER*255 OutName,DumName1,DumName2
      LOGICAL PRaw,PYaw,PPitch,PRoll,PFreq,PO2,PWebb,PCal,PIndep,
     &  DoYaw,DoPitch,DoRoll,DoFreq,DoO2,DoWebb,DoCrMean,PSonic,
     &  DoSonic,DoTilt,PTilt,DoPrint,Flag(NNMAx,MMMax),DoDetren,
     &  PDetrend,DoStruct,BadTc, DoPF, PPF
      INTEGER N,i,j,M,MIndep(NNMax),CIndep(NNMax,NNMax),FOO,
     &  Channels,Delay(NNNMax),Mok(NNMax),Cok(NNMax,NNMax), FirstDay,
     &  NDelta, NLock, NHigh, NLow,
     &  StartTime(3),StopTime(3),
     &  LStartTime(3),LStopTime(3)

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
     &  Alpha, Beta, Gamma, WBias,
     &  UMean(3, MaxPF), Apf(3,3), PFValid
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ), CalCO2(NQQ),
     &  R,dR,CTSon2,CTCop2,Cq2,CTSonq,CTCopq,
     &  dCTSon2,dCTCop2,dCq2,dCTSonq,dCTCopq
      REAL*8 StdTs, dStdTs,
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
      LOGICAL       HAVE_UNCAL(NNNMax), EndInter, PastStart, BeforeEnd,
     &              SingleRun, DoWBias, LastInter

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
      WRITE(*,*) 'PLANFIT Version: ', VERSION
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
C Check which calibration files we have
C
      IF (EC_T_STRLEN(SonName) .GT. 0) THEN
C For now assume it is always a 3-D sonic
         HAVE_UNCAL(U) = .TRUE.
         HAVE_UNCAL(V) = .TRUE.
         HAVE_UNCAL(W) = .TRUE.
         HAVE_UNCAL(TSonic) = .TRUE.
      ENDIF
C
C Number of quantities involved in this experiment
C
      N = 3  ! Number of calibrated quantities you will follow
C
C Read parameters from file with name ParmName
C
      CALL EC_F_Params(ParmName,
     &  Freq,PitchLim,RollLim,PreYaw,PrePitch,PreRoll,
     &  LLimit,ULimit,DoCrMean,DoDetren,DoSonic,
     &  DoTilt,DoYaw,DoPitch,DoRoll,DoFreq,DoO2,DoWebb,DoStruct,DoPF,
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
C
C Set delays, gains and offset of all channels in raw datafile
C At reading time all channels are corrected by mapping:
C                x --> (x/Gain) + Offset
C

      Delay(U)        = NINT(CalSonic(QQDelay)*0.001D0*Freq)
      Delay(V)        = NINT(CalSonic(QQDelay)*0.001D0*Freq)
      Delay(W)        = NINT(CalSonic(QQDelay)*0.001D0*Freq)

      Gain(U)        = CalSonic(QQGain)
      Gain(V)        = CalSonic(QQGain)
      Gain(W)        = CalSonic(QQGain)

      Offset(U)        = CalSonic(QQOffset)
      Offset(V)        = CalSonic(QQOffset)
      Offset(W)        = CalSonic(QQOffset)

C
C Open interval file
C
      OPEN(IntervalFile,FILE = InterName,IOSTAT=FOO,STATUS='OLD')
      IF (FOO.NE.0) THEN
        WRITE(*,*) 'Could not open interval-file ',InterName
        STOP'- Fatal: no time info could be read -'
      ENDIF
C
C Get planar fit interval file
C 
      OPEN(PlfIntFile,FILE = PlfIntName,IOSTAT=FOO,STATUS='OLD')
      IF (FOO.NE.0) THEN
        WRITE(*,*) 'Could not open planar fit interval-file ',PlfIntName
        STOP'- Fatal: planarfit intervals could be read -'
      ENDIF
C
C Open Planar fit angles file
C
      OPEN(PlfFile,FILE = PlfName,IOSTAT=FOO,STATUS='Unknown')
      IF (FOO.NE.0) THEN
        WRITE(*,*) 'Could not open planar fit angles-file ',PlfName
        STOP'- Fatal: no planar fit angles could be written -'
      ENDIF

C
C######################################################################
C######################################################################
C
C Here starts the loop over the different Planar fit intervals
C
C######################################################################
C######################################################################
C
C
C Read time-interval from file
C
 36   READ(PlfIntFile,*,END=37) (StartTime(i),i=1,3),
     &  (StopTime(i),i=1,3)
     
      FName   = DumName1
      OutName = DumName2

      EndInter = .FALSE.
      NPF = 0
C
C Now cycle through Interval file to find intervals that are within
C this Planar fit interval
C
      OPEN(IntervalFile,FILE = InterName,IOSTAT=FOO,STATUS='OLD')
      IF (FOO.NE.0) THEN
         WRITE(*,*) 'Could not open interval-file ',InterName
          STOP'- Fatal: no time info could be read -'
      ENDIF
C Read first line blindly (comment)
      READ(IntervalFile,*)
      
      DO WHILE (.NOT. EndInter)
        READ(IntervalFile,*,END=37) (LStartTime(i),i=1,3),
     &       (LStopTime(i),i=1,3),Psychro,p,DumName1, DumName2
C
C Check if we are past start
C
        IF ((LStartTime(1) .GT. StartTime(1)) .OR.
     &      ((LStartTime(1) .EQ. StartTime(1)) .AND.
     &       (LStartTime(2) .GT. StartTime(2))) .OR.
     &      ((LStartTime(1) .EQ. StartTime(1)) .AND.
     &       (LStartTime(2) .EQ. StartTime(2)) .AND.
     &       (LStartTime(3) .GE. StartTime(3)))) THEN
           PastStart = .TRUE.
        ELSE
           PastStart = .FALSE.
        ENDIF
C
C Check if we are before end
C
        IF ((LStopTime(1) .LT. StopTime(1)) .OR.
     &      ((LStopTime(1) .EQ. StopTime(1)) .AND.
     &       (LStopTime(2) .LT. StopTime(2))) .OR.
     &      ((LStopTime(1) .EQ. StopTime(1)) .AND.
     &       (LStopTime(2) .EQ. StopTime(2)) .AND.
     &       (LStopTime(3) .LE. StopTime(3)))) THEN
           BeforeEnd = .TRUE.
        ELSE
           BeforeEnd = .FALSE.
        ENDIF
C
C Check if we are at last interval (might be last interval of file,
C that's why we check it): compute angles either when this
C is the last interval of the total Planar fit interval,
C or when we are beyond the interval (somehow we missed the
C last interval
C
        IF ((LStopTime(1) .EQ. StopTime(1)) .AND.
     &      (LStopTime(2) .EQ. StopTime(2)) .AND.
     &      (LStopTime(3) .EQ. StopTime(3))) THEN
           LastInter = .TRUE.
        ELSE
           LastInter = .FALSE.
        ENDIF

        IF (PastStart .AND. BeforeEnd) THEN
           NPF = NPF + 1
           
           DUMSTRING=DatDir
           Fname = DumName1
           CALL EC_T_STRCAT(DUMSTRING, '/')
           CALL EC_T_STRCAT(DUMSTRING, Fname)
           Fname=DUMSTRING
C
C Show which file is currently being analysed
C
           WRITE(*,951) FName(:EC_T_STRLEN(Fname)),(LStartTime(i),i=1,3)
 951  FORMAT(a,1X,3(I5,1X))
C
C For subroutine calibrat below find the first day; transfer via interface 
C via EC_G_Main
C 
C
           FirstDay = StartTime(1)
C      
C Set Number of uncalibrated signals
C
           Channels = 3
C      
C Set Number of calibrated signals
C
           N = 3
C
C Reset Have_Uncal
C
           DO I=1,Channels
               HAVE_UNCAL(I) = .FALSE.
           ENDDO
C
C Read raw data from file
C
           CALL EC_F_ReadNCDF(FName,LStartTime,LStopTime,
     &         Delay,Gain,Offset,RawSampl,NNNMax,MMMax,M,
     &         NCVarName, HAVE_UNCAL)
C 
C Do calibration (per sample)
C 
           N = 3
           DO i=1,M
               CALL Calibrat(RawSampl(1,i),Channels,P,CorMean,
     &	          CalSonic,CalTherm,CalHyg,CalCO2, BadTc,Sample(1,i),N,Flag(1,i),
     &            Have_Uncal, FirstDay)
           ENDDO
C
C Determine mean velocities, check if there are enough valid samples
C
    	   CALL EC_M_Averag(Sample, NNMax, N, MMMax, M, Flag,
     &                      Mean, TolMean, Cov, TolCov, MIndep, CIndep,
     &                      Mok, Cok)
           IF ((DBLE(MOK(U))/DBLE(M) .LT. PFValid) .OR.
     &         (DBLE(MOK(U))/DBLE(M) .LT. PFValid) .OR.
     &         (DBLE(MOK(W))/DBLE(M) .LT. PFValid)) THEN
               NPF = NPF - 1
               write(*,*) 'Interval rejected: ', 
     &             Min(DBLE(MOK(U))/DBLE(M), DBLE(MOK(U))/DBLE(M), 
     &                 DBLE(MOK(W))/DBLE(M)), ' < ', PFValid
           ELSE
               UMean(U,NPF) = Mean(U)
               UMean(V,NPF) = Mean(V)
               UMean(W,NPF) = Mean(W)
           ENDIF
        ENDIF
        IF (PastStart .AND. 
     &      ((.NOT. BeforeEnd) .OR. (LastInter))) THEN
           EndInter = .TRUE.
           SingleRun = (.FALSE.)
	   DoWBias = (.FALSE.)
           CALL EC_C_T01(DoWBias, SingleRun, UMean, NPF,
     &                   Apf,Alpha,Beta,Gamma,WBias)
           Write(PlfFile,100) (StartTime(i),i=1,3), 
     &           (StopTime(i),i=1,3), Alpha, Beta, Gamma, WBias,
     &           ((Apf(i,j),i=1,3),j=1,3)
 100       FORMAT(6(I8,1X),1X,13(F20.15))
        ENDIF
      ENDDO
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
      CLOSE(PlfFile)
      CLOSE(PlfIntFile)

      END








