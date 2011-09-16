C $Id$
C  ec_file.f, part of ecpack
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


C     ****f* ec_file.f/EC_F_GetConf
C NAME
C     EC_F_GetConf
C SYNOPSIS
C     CALL EC_F_GetConf(ConfUnit,
C               DatDir, OutDir,  ParmDir,
C               FluxName, ParmName, InterName, PlfIntName,
C               PlfName, SonName, CoupName, HygName, CO2Name,
C               NCVarName, NCNameLen,
C               OutMean, OutCov, OutPh, OutTime, OutStd, OutNum, OutStr,
C               OutFrcor,
C               Outputs, DoCorr, CorrPars, ExpVar, DoStruct, DoPrint,
C               PCorr, PRaw, PCal, PIndep)
C INPUTS
C     ConfUnit : [INTEGER]
C               file unit
C     NCNameLen: [INTEGER]
C                Length of array NCVarName
C OUTPUT
C     DatDir   : [CHARACTER*(*)]
C                Directory with data files
C     OutDir   : [CHARACTER*(*)]
C                Directory for output
C     ParmDir  : [CHARACTER*(*)]
C                Directory for parameter/calibration files
C     FluxName : [CHARACTER*(*)]
C                file name for flux file
C     ParmName : [CHARACTER*(*)]
C                file name for parameter file (which corrections etc.)
C     InterName: [CHARACTER*(*)]
C                file name for interval file
C     PlfIntName: [CHARACTER*(*)]
C                file name for planar fit interval file
C     PlfName  : [CHARACTER*(*)]
C                file name for planar fit angles file
C     SonName  : [CHARACTER*(*)]
C                Sonic calibration file name
C     CoupName : [CHARACTER*(*)]
C                Thermocouple calibration file name
C     HygName  : [CHARACTER*(*)]
c                Hygrometer calibration file name
C     CO2Name  : [CHARACTER*(*)]
C                CO2 sensor calibration file
C     NCVarName: [CHARACTER*(*)](NCNameLen)
C                array with NetCDF variables names
C     OutMean  : [LOGICAL](NNMax)
C                which means to write to output
C     OutCov   : [LOGICAL](NNMax, NNMax)
C                which covariances to write to output
C     OutPh    : [LOGICAL](NMaxPhys)
C                which physical quantities to write to output 
C     OutTime  : [LOGICAL](NMaxOST)
C                which time variables to output 
C     OutStd   : [LOGICAL](NNMax)
C                which standard deviations to write to output
C     OutNum   : [LOGICAL](NNMax)
C                which number of samples to write to output
C     OutStr   : [LOGICAL](NNMax, NNMax)
C                which structure parameters to write to output,
C     OutFrCor : [LOGICAL](NNMax, NNMax)
C                which frequency correction factors to write to output,
C     OutPuts  : [LOGICAL](NMaxOS)
C                further info about output
C     DoCorr   : [LOGICAL](NMaxCorr)
C                Which corrections to do
C     CorrPars : [REAL*8](NMaxCorrPar)
C                Parameters of corrections
C     ExpVar   : [REAL*8](NMaxExp)
C                array with Experimental settings
C     DoStruct : [LOGICAL]
C                compute structure parameters
C     DoPrint  : [LOGICAL]
C                print any output
C     DoCorr   : [LOGICAL](NMaxCorr)
C                Switch for corrections
C     PCorr    : [LOGICAL](NMaxCorr)
C                Switch for printing intermediate results after
C                corrections
C     PRaw     : [LOGICAL]
C                write intermediate results with respect to
C                raw data
C     PCal     : [LOGICAL]
C                write intermediate results with respect to
C                calibrated data
C     PIndep     : [LOGICAL]
C                write intermediate results with respect to
C                number of independent samples
C FUNCTION
C     Routine to read info with respect to input and output files
C AUTHOR
C     Arnold Moene
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     EC_T_ClearSTR
C     EC_T_StripSTR
C     EC_T_UpCase
C     EC_T_StrCat
C     EC_T_GOut1
C     EC_T_GOut2
C     EC_T_GPhys
C     EC_T_GOutS
C     EC_T_GCorr
C     EC_T_GCorrPars
C     parcnst.inc
C     ***

      SUBROUTINE EC_F_GetConf(ConfUnit,
     &           DatDir, OutDir,  ParmDir,
     &           FluxName, ParmName, InterName, PlfIntName,
     &           PlfName, SonName, CoupName, HygName, CO2Name,
     &           NCVarName, NCNameLen,
     &           OutMean, OutCov, OutPh, OutTime, 
     &           OutStd, OutNum, OutStr,
     &           OutFrcor,
     &           Outputs, DoCorr, CorrPars, ExpVar, 
     &           DoStruct, DoPrint,
     &           PCorr, PRaw, PCal, PIndep)
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'
      
      
      LOGICAL         DoPrint,PRaw,PCal,PIndep,
     &                 DoStruct, PCorr(NMaxCorr)
c     CHARACTER*(*)   InName
      INTEGER         ConfUnit, NCNameLen, I, J
      CHARACTER*(*)   DatDir, OutDir, ParmDir, FluxName, ParmName,
     &                InterName, SonName, CoupName, HygName, 
     &                CO2Name, PlfName, PlfIntName
      CHARACTER*(*)   NCVarName(NCNameLen)
      LOGICAL         OutMean(NNMax), OutCov(NNMax, NNMax), 
     &                OutStd(NNMax), OutNum(NNMax),  
     &                OutStr(NNMax, NNMax), 
     &                OutFrcor(NNMax, NNMax), 
     &                OutPh(NMaxPhys), OutTime(NMaxOST),
     &                Outputs(NMaxOS),
     &                Docorr(NMaxCorr) 
      REAL*8          CorrPars(NMaxCorrPar), ExpVar(NMaxExp)

      INTEGER         IOCODE, KINDEX
      CHARACTER*255   LINE, TOKLINE, VALLINE, DUMSTRING
      INTEGER         EC_T_STRLEN
      EXTERNAL        EC_T_STRLEN

C
C Set defaults
C
      DO 1000, I=1,NCNameLen
         CALL EC_T_CLEARSTR(NCVarName(I))
 1000 CONTINUE
      CALL EC_T_CLEARSTR(DatDir)
      CALL EC_T_CLEARSTR(OutDir)
      CALL EC_T_CLEARSTR(ParmDir)
      CALL EC_T_CLEARSTR(FluxName)
      CALL EC_T_CLEARSTR(ParmName)
      CALL EC_T_CLEARSTR(InterName)
      CALL EC_T_CLEARSTR(PlfIntName)
      CALL EC_T_CLEARSTR(SonName)
      CALL EC_T_CLEARSTR(CoupName)
      CALL EC_T_CLEARSTR(HygName)
      CALL EC_T_CLEARSTR(CO2Name)
      CorrPars(QCPBlock) = 1



      IOCODE = 0
      DO 4000, WHILE (IOCODE .EQ. 0)
C     Read one line from configuration file
         READ(Confunit,'(A)', IOSTAT=IOCODE, END=9000) LINE
         IF (IOCODE .NE. 0) THEN
            WRITE(*,*) 'ERROR in reading of configuration file'
            STOP
         ENDIF
         CALL EC_T_STRIPSTR(LINE)
C     Check for comment (should start with //
         KINDEX = INDEX(LINE, '//')
         IF (KINDEX .EQ. 0) THEN
C     Find the equality sign
            KINDEX = INDEX(LINE, '=')
            IF (KINDEX .EQ. 0) THEN
                WRITE(*,*) 'WARNING: no equality sign found in ', LINE
                WRITE(*,*) 'assuming an empty line'
            ENDIF
c     Split into token and value
            WRITE(TOKLINE,*) LINE(:KINDEX-1)
            WRITE(VALLINE,*) LINE(KINDEX+1:)
            CALL EC_T_STRIPSTR(TOKLINE)
            CALL EC_T_STRIPSTR(VALLINE)
            CALL EC_T_UPCASE(TOKLINE)
C     See which token this is
            IF (INDEX(TOKLINE, 'DATDIR') .GT. 0) THEN
               DatDir = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'OUTDIR') .GT. 0) THEN
               OutDir = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'PARMDIR') .GT. 0) THEN
               ParmDir = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'FLUXNAME') .GT. 0) THEN
               FluxName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'PARMNAME') .GT. 0) THEN
               ParmName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'INTERNAME') .GT. 0) THEN
               InterName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'PLFINTNAME') .GT. 0) THEN
               PlfIntName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'SONNAME') .GT. 0) THEN
               SonName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'COUPNAME') .GT. 0) THEN
               CoupName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'HYGNAME') .GT. 0) THEN
               HygName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'CO2NAME') .GT. 0) THEN
               CO2Name = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'PLFNAME') .GT. 0) THEN
               PlfName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
C NetCdF variable names
            ELSE IF (INDEX(TOKLINE, 'U_VAR') .GT. 0) THEN
               NCVarName(QUU) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'V_VAR') .GT. 0) THEN
               NCVarName(QUV) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'W_VAR') .GT. 0) THEN
               NCVarName(QUW) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'TSONIC_VAR') .GT. 0) THEN
               NCVarName(QUTSonic) = 
     &            VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'TCOUPLE_VAR') .GT. 0) THEN
               NCVarName(QUTCouple) = 
     &            VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'HUMIDITY_VAR') .GT. 0) THEN
               NCVarName(QUHUMIDITY) = 
     &            VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'CO2_VAR') .GT. 0) THEN
               NCVarName(QUCO2) = 
     &            VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'DOY_VAR') .GT. 0) THEN
               NCVarName(QUDoy) = 
     &            VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'HOURMIN_VAR') .GT. 0) THEN
               NCVarName(QUHourMin) = 
     &            VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'SEC_VAR') .GT. 0) THEN
               NCVarName(QUSec) = 
     &            VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'DIAG_VAR') .GT. 0) THEN
               NCVarName(QUDiagnost) = 
     &            VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'TREF_VAR') .GT. 0) THEN
               NCVarName(QUTref) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'SSPEED_VAR') .GT. 0) THEN
               NCVarName(QUSSpeed) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
C Correction definitions
            ELSE IF (INDEX(TOKLINE, 'CORRECT') .GT. 0) THEN
               CALL EC_T_GCorr(VALLINE(:INDEX(VALLINE, CHAR(0))-2),
     &                            DoCorr)
C Correction parameters
            ELSE IF (INDEX(TOKLINE, 'CORRPAR') .GT. 0) THEN
               CALL EC_T_GCorrPar(VALLINE(:INDEX(VALLINE, CHAR(0))-2),
     &                            CorrPars)
C Output definitions
            ELSE IF (INDEX(TOKLINE, 'OUT_MEAN') .GT. 0) THEN
               CALL EC_T_GOut1(VALLINE(:INDEX(VALLINE, CHAR(0))-1),
     &                            OutMean)
            ELSE IF (INDEX(TOKLINE, 'OUT_COV') .GT. 0) THEN
               CALL EC_T_GOut2(VALLINE(:INDEX(VALLINE, CHAR(0))-1),
     &                            OutCov)
            ELSE IF (INDEX(TOKLINE, 'OUT_PHYS') .GT. 0) THEN
               CALL EC_T_GPhys(VALLINE(:INDEX(VALLINE, CHAR(0))-1),
     &                            OutPh)
            ELSE IF (INDEX(TOKLINE, 'OUT_TIME') .GT. 0) THEN
               CALL EC_T_GTime(VALLINE(:INDEX(VALLINE, CHAR(0))-1),
     &                            OutTime)
            ELSE IF (INDEX(TOKLINE, 'OUT_STRUCT') .GT. 0) THEN
               CALL EC_T_GOut2(VALLINE(:INDEX(VALLINE, CHAR(0))-1),
     &                            OutStr)
            ELSE IF (INDEX(TOKLINE, 'OUT_FRCOR') .GT. 0) THEN
               CALL EC_T_GOut2(VALLINE(:INDEX(VALLINE, CHAR(0))-1),
     &                            OutFrcor)
            ELSE IF (INDEX(TOKLINE, 'OUT_STD') .GT. 0) THEN
               CALL EC_T_GOut1(VALLINE(:INDEX(VALLINE, CHAR(0))-1),
     &                            OutStd)
            ELSE IF (INDEX(TOKLINE, 'OUT_NUM') .GT. 0) THEN
               CALL EC_T_GOut1(VALLINE(:INDEX(VALLINE, CHAR(0))-1),
     &                            OutNum)
            ELSE IF (INDEX(TOKLINE, 'OUTPUT') .GT. 0) THEN
               CALL EC_T_GOutS(VALLINE(:INDEX(VALLINE, CHAR(0))-1),
     &                            Outputs)
            ELSE 
               write(*,*) 'ERROR: unkown token: ', 
     &                    TOKLINE(:INDEX(TOKLINE, CHAR(0))-1)
               STOP
            ENDIF
         ENDIF
 4000 CONTINUE
 9000 CONTINUE

      CALL EC_T_CLEARSTR(DUMSTRING)
      IF (EC_T_STRLEN(PlfName) .GT. 0) THEN
         DUMSTRING = OutDir
         CALL EC_T_STRCAT(DUMSTRING,'/')
         CALL EC_T_STRCAT(DUMSTRING,PlfName)
      ENDIF
      PlfName = DUMSTRING

      CALL EC_T_CLEARSTR(DUMSTRING)
      IF (EC_T_STRLEN(FluxName) .GT. 0) THEN
         DUMSTRING = OutDir
         CALL EC_T_STRCAT(DUMSTRING,'/')
         CALL EC_T_STRCAT(DUMSTRING,FluxName)
      ENDIF
      FluxName = DUMSTRING

      CALL EC_T_CLEARSTR(DUMSTRING)
      IF (EC_T_STRLEN(ParmName) .GT. 0) THEN
         DUMSTRING = ParmDir
         CALL EC_T_STRCAT(DUMSTRING,'/')
         CALL EC_T_STRCAT(DUMSTRING,ParmName)
      ENDIF
      ParmName = DUMSTRING

      CALL EC_T_CLEARSTR(DUMSTRING)
      IF (EC_T_STRLEN(InterName) .GT. 0) THEN
         DUMSTRING = ParmDir
         CALL EC_T_STRCAT(DUMSTRING,'/')
         CALL EC_T_STRCAT(DUMSTRING,InterName)
      ENDIF
      InterName = DUMSTRING
   
      CALL EC_T_CLEARSTR(DUMSTRING)
      IF (EC_T_STRLEN(PlfIntName) .GT. 0) THEN
         DUMSTRING = ParmDir
         CALL EC_T_STRCAT(DUMSTRING,'/')
         CALL EC_T_STRCAT(DUMSTRING,PlfIntName)
      ENDIF
      PlfIntName = DUMSTRING

      CALL EC_T_CLEARSTR(DUMSTRING)
      IF (EC_T_STRLEN(SonName) .GT. 0) THEN
         DUMSTRING = ParmDir
         CALL EC_T_STRCAT(DUMSTRING,'/')
         CALL EC_T_STRCAT(DUMSTRING,SonName)
      ENDIF
      SonName = DUMSTRING
      
      CALL EC_T_CLEARSTR(DUMSTRING)
      IF (EC_T_STRLEN(CoupName) .GT. 0) THEN
         DUMSTRING = ParmDir
         CALL EC_T_STRCAT(DUMSTRING,'/')
         CALL EC_T_STRCAT(DUMSTRING,CoupName)
      ENDIF
      CoupName = DUMSTRING
      
      CALL EC_T_CLEARSTR(DUMSTRING)
      IF (EC_T_STRLEN(HygName) .GT. 0) THEN
         DUMSTRING = ParmDir
         CALL EC_T_STRCAT(DUMSTRING,'/')
         CALL EC_T_STRCAT(DUMSTRING,HygName)
      ENDIF
      HygName = DUMSTRING

      CALL EC_T_CLEARSTR(DUMSTRING)
      IF (EC_T_STRLEN(CO2Name) .GT. 0) THEN
         DUMSTRING = ParmDir
         CALL EC_T_STRCAT(DUMSTRING,'/')
         CALL EC_T_STRCAT(DUMSTRING,CO2Name)
      ENDIF
      CO2Name = DUMSTRING

C
C Set output flags for structure functions, if 'old' output is needed
C (structure functions are too expensive to just compute it without
C using it
C
      IF (Outputs(OSOld)) THEN
         DO I=1,NNMax
            DO J=1,NNMax
               OutStr(I,J) = .FALSE.
            ENDDO
         ENDDO
         OutStr(TSonic,TSonic) = .TRUE.
         OutStr(TCouple,TCouple) = .TRUE.
         OutStr(SpecHum,SpecHum) = .TRUE.
         OutStr(TCouple,Spechum) = .TRUE.
         OutStr(TSonic,Spechum) = .TRUE. 
      ENDIF

C 
C Now call EC_F_Params
C
      CALL EC_F_Params(ParmName, ExpVar,
     &    DoCorr, PCorr,
     &    DoStruct, DoPrint,
     &    PRaw,PCal,PIndep)

      END





C     ****f* ec_file.f/EC_F_Params
C NAME
C     EC_F_Params
C SYNOPSIS
C     CALL EC_F_Params(InName, ExpVar,
C             DoCorr, PCorr,
C             DoStruct, DoPrint,
C             PRaw,PCal,PIndep)
C INPUTS
C     InName   : [CHARACTER*(*)]
C                Name of parameter file
C OUTPUT
C     ExpVar   : [REAL*8](NMaxExp)
C                array with Experimental settings
C     DoCorr   : [LOGICAL](NMaxCorr)
C                Switch for corrections
C     PCorr    : [LOGICAL](NMaxCorr)
C                Switch for printing intermediate results after
C                corrections
C     DoStruct : [LOGICAL]
C                compute structure parameters
C     DoPrint  : [LOGICAL]
C                print any output
C     PRaw     : [LOGICAL]
C                write intermediate results with respect to
C                raw data
C     PCal     : [LOGICAL]
C                write intermediate results with respect to
C                calibrated data
C     PIndep     : [LOGICAL]
C                write intermediate results with respect to
C                number of independent samples
C FUNCTION
C     Routine to read info with respect to corrections and writing
C     of intermediate results
C AUTHOR
C     Arnold Moene
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     parcnst.inc
C     physcnst.inc
C     ***
      SUBROUTINE EC_F_Params(InName, ExpVar,
     &    DoCorr, PCorr,
     &    DoStruct, DoPrint,
     &    PRaw,PCal,PIndep)
C
C Purpose : Read settings for this particular session from file
C
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      REAL*8 ExpVar(NMaxExp)
      LOGICAL 
     &    DoPrint,PRaw,PCal,PIndep,
     &    DoStruct, DoCorr(NMaxCorr), PCorr(NMaxCorr)
      CHARACTER*(*) InName
      INTEGER      IOCODE

      OPEN(TempFile,FILE=InName, STATUS='old', IOSTAT=IOCODE)
      IF (IOCODE .NE. 0) THEN
         WRITE(*,*) 'ERROR: could not open parameter file ', InName
         STOP
      ENDIF
C
C Default values for planar fit settings
C
      DoCorr(QCPF) = .FALSE.
      PCorr(QCPF) = .FALSE.
      ExpVar(QEPFValid) = 1.0D0
      READ(TempFile,*) ExpVar(QEFreq)      ! [Hz] Sample frequency datalogger
      READ(TempFile,*)
      READ(TempFile,*) ExpVar(QEPitchLim)      ! [degree] Limit when Mean(W) is turned to zero
      READ(TempFile,*) ExpVar(QERollLim)      ! [degree] Limit when Cov(V,W) is turned to zero
      READ(TempFile,*)
      READ(TempFile,*) ExpVar(QEPreYaw)   ! Fixed yaw angle for known tilt-correction
      READ(TempFile,*) ExpVar(QEPrePitch)      ! Fixed pitch angle for known tilt-correction
      READ(TempFile,*) ExpVar(QEPreRoll)      ! Fixed roll angle for known tilt-correction
      READ(TempFile,*)
      READ(TempFile,*) ExpVar(QELLimit)      ! Smallest acceptable frequency-response
correction factor
      READ(TempFile,*) ExpVar(QEULimit)      ! Largest acceptable frequency-response
correction factor
      READ(TempFile,*)
      READ(TempFile,*) DoCorr(QCMean)     ! Replace mean quantities by better estimates
      READ(TempFile,*) DoCorr(QCDetrend)  ! Correct data for linear trend
      READ(TempFile,*) DoCorr(QCSonic)      ! Correct sonic temperature for humidity
      READ(TempFile,*) DoCorr(QCTilt)      ! Perform true tilt-correction with known angles
      READ(TempFile,*) DoCorr(QCYaw)      ! Turn system such that Mean(V) --> 0
      READ(TempFile,*) DoCorr(QCPitch)      ! Turn system such that Mean(W) --> 0
      READ(TempFile,*) DoCorr(QCRoll)      ! Turn System such that Cov(W,V) --> 0
      READ(TempFile,*) DoCorr(QCFreq)      ! Correct for poor frequency response
      READ(TempFile,*) DoCorr(QCO2)      ! Correct hygrometer for oxygen-sensitivity
      READ(TempFile,*) DoCorr(QCWebb)      ! Calculate mean velocity according to Webb
      READ(TempFile,*) DoStruct      ! Calculate structure parameters
      READ(TempFile,*)
      READ(TempFile,*) DoPrint      ! Skip printing intermediate results or not?
      READ(TempFile,*)
      READ(TempFile,*) PRaw      ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCal      ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCDetrend)   ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PIndep      ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCTilt)      ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCYaw)      ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCPitch)      ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCRoll)      ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCSonic)      ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCO2)      ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCFreq)      ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCWebb)      ! Indicator if these intermediate results are wanted
      READ(TempFile,*) ExpVar(QEStructSep)! Separation (meter) for which to calculate structure parameter
      READ(TempFile, *, END=5500) DoCorr(QCPF) ! Do Planar fit ?
      READ(TempFile, *, END=5500) PCorr(QCPF) ! Print planar fit intermediate results ?
      READ(TempFile, *, END=5500) ExpVar(QEPFValid) ! Minimum part (<= 1) of samples that should be valid to incorporate interval in angle PF calculation
      
5500  CONTINUE
      CLOSE(TempFile)

      RETURN
      END




C
C ########################################################################
C
C This routine reads specifications of a calibrated apparatus into an
C array. At the moment BOTH calibration constants AND position of
C the apparatus in the setup are read from the same file. (It may be useful
C to split this in a future version of this library into two files, such
C that setup-dependent constants don't require modification of the
C calibration file.)
C
C ########################################################################
C





      SUBROUTINE EC_F_ReadAp(InName,CalSpec)
C
C Read specifications of an apparatus from file
C
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      CHARACTER*(*) InName
      REAL*8 CalSpec(NQQ)
      INTEGER i, IOCODE
      INTEGER EC_T_STRLEN
      EXTERNAL EC_T_STRLEN

      INTEGER SPECS_READ

      OPEN(TempFile,FILE=InName,STATUS='OLD', IOSTAT=IOCODE)
      IF (IOCODE .NE. 0) THEN
         WRITE(*,5100) InName(:EC_T_STRLEN(InName))
         STOP
      ENDIF
 5100 FORMAT ('ERROR: can not open calibration file ', (A))

C First get the type of apparatus
      READ(TempFile,*) CalSpec(1)

C Now read the correct number of specs
      SPECS_READ=1
      DO i=2,ApNQQ(INT(CalSpec(1)))
           READ(TempFile,*,IOSTAT=IOCODE, END = 9000) CalSpec(i)
         IF (IOCODE .NE. 0) THEN
            WRITE(*,*) 'ERROR in reading of calibration file ', 
     &               InName(:EC_T_STRLEN(InName))
            STOP
         ENDIF
         SPECS_READ = SPECS_READ + 1
      ENDDO
 9000 CONTINUE
      IF (SPECS_READ .NE. ApNQQ(INT(CalSpec(1)))) THEN
         WRITE(*,*) 'ERROR: did not find correct number of specs in ',
     &               InName(:EC_T_STRLEN(InName))
      ENDIF

      CLOSE(TempFile)

      RETURN
      END






      SUBROUTINE EC_F_ReadDt(InName,Delay,Gain,Offset,x,NMax,N,MMax,M)
C
C Read raw data from file
C
C input : InName : CHARACTER*40 : Name of the datafile (ASCII)
C      NMax : INTEGER : Maximum number of quantities in this array
C      MMax : INTEGER : Maximum number of samples in this array
C      N : INTEGER : Number of quantities in file
C
C output : x(NMax,MMax) : REAL*8 Raw data array.
C              First index counts quantities; second counter
C              counts samples. Only the first N quantities
C              and the first M samples are used.
C       M : Number of samples in file
C
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      CHARACTER*40 InName
      INTEGER NMax,N,MMax,M,j,Delay(NMax),MaxDelay
      REAL*8 x(NMax,MMax),Dum(NNNMax),Gain(NMax),Offset(NMax)

      MaxDelay = 0
      DO j=1,N
         MaxDelay = MAX(Delay(j),MaxDelay)
      ENDDO

      OPEN(InFile,FILE=InName)
      M = 0

 10   READ(InFile,*,END=20) (Dum(j),j=1,N)

      M = M+1
      DO j=1,N
         IF ((M-Delay(j)).GT.0)
     &       x(j,M-Delay(j)) = Dum(j)/Gain(j) + Offset(j)
      ENDDO

      IF (M.LT.MMMax) GOTO 10

 20   CLOSE(InFile)
      M = M-MaxDelay

      RETURN
      END

      SUBROUTINE EC_F_ReadNCDF(InName,StartTime,StopTime,
     &  Delay,Gain,Offset,x,NMax,MMax,M, NCVarName,
     &  Have_Uncal)
C
C Read raw data from NETCDF-file
C
C input : InName : CHARACTER*40 : Name of the datafile (ASCII)
C      NMax : INTEGER : Maximum number of quantities in this array
C      MMax : INTEGER : Maximum number of samples in this array
C
C output : x(NMax,MMax) : REAL*8 Raw data array.
C              First index counts quantities; second counter
C              counts samples. Only the first N quantities
C              and the first M samples are used.
C       M : Number of samples in file
C Revision 28-05-2001: added info on whether uncalibrated data are
C                      available for a given variable (mainly important
C                      for sonic and/or Couple temperature since that
C                      is used for various corrections)
C Revision 19-9-2002:  removed N (number of quantities in file) from
C                      interface
C
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'
      INCLUDE 'netcdf.inc'

      CHARACTER*(*) InNAME
      CHARACTER*40 NAME ! Maximum length of names is given by NF_MAX_NAME
      INTEGER STATUS,
     &        NCID,       ! ID for file
     &        NDIMS,      ! # of dimensions
     &        NVARS,      ! # of variables
     &        NGATTS,     ! # of global attributes
     &        UNLIMITED,  ! ID of unlimited dimension
     &        XTYPE,      ! number type of variable
     &        VARDIMS,    ! # of dimensions for variable
     &        DIMIDS(NF_MAX_VAR_DIMS), ! IDs of dimensions
     &        NATTS,      ! # of attributes
     &        DIMLEN      ! length of a dimension   
      INTEGER I, J, K
      INTEGER NMax,MMax,M,Delay(NMax),HMStart,HMStop
      LOGICAL ok
      REAL*8 x(NMax,MMax),Gain(NMax),Offset(NMax)
      INTEGER StartTime(3), StopTime(3)
      CHARACTER*(*) NCVarName(NMax)
      LOGICAL       Have_Uncal(NMax)

      CHARACTER*255 ATT_NAME
      INTEGER       NCVarID(NNNMax), DoyStart, DoyStop,
     +              STARTIND, STOPIND, NSAMPLE
      LOGICAL       HAS_SCALE(NMax), HAS_OFFSET(NMax)
      REAL*8        NC_SCALE(NMax), NC_OFFSET(NMax)
      INTEGER       EC_T_STRLEN
      EXTERNAL      EC_T_STRLEN

      HMStart = 100*StartTime(2)+StartTime(3)
      HMStop  = 100*StopTime( 2)+StopTime( 3)
C
C Open file (NF_NOWRITE and NF_NOERR are defined in netcdf.inc)
C
      STATUS = NF_OPEN(InNAME, NF_NOWRITE, NCID)
      IF (STATUS .NE. NF_NOERR) CALL EC_NCDF_HANDLE_ERR(STATUS)
C
C Inquire about dimensions etc.
C
      STATUS = NF_INQ(NCID, NDIMS, NVARS, NGATTS, UNLIMITED)
      IF (STATUS .NE. NF_NOERR) CALL EC_NCDF_HANDLE_ERR(STATUS)

      M = 0
C
C Inquire about variables
C
      DO I=1,NMax
         NCVarID(I) = 0
         HAS_SCALE(I) = .FALSE.
         HAS_OFFSET(I) = .FALSE.
         Have_Uncal(I) = .FALSE.
      ENDDO
      DO I=1,NVARS
        STATUS = NF_INQ_VAR(NCID,I,NAME,XTYPE,VARDIMS,DIMIDS,NATTS)
        IF (STATUS .NE. NF_NOERR) CALL EC_NCDF_HANDLE_ERR(STATUS)
        DO J=1,NMax
           IF (Name .EQ. NCVarName(J)) THEN
              NCVarID(J) = I
              Have_Uncal(J) = .TRUE.
              DO K=1,NATTS
                 STATUS = NF_INQ_ATTNAME(NCID, I, K, ATT_NAME)
                 IF (STATUS .NE. NF_NOERR)
     &              CALL EC_NCDF_HANDLE_ERR(STATUS)
                 IF (INDEX(ATT_NAME, 'scale_factor') .GT. 0) THEN
                    STATUS = NF_GET_ATT_DOUBLE(NCID, I, 'scale_factor',
     &                       NC_SCALE(J))
                    IF (STATUS .NE. NF_NOERR)
     &                    CALL EC_NCDF_HANDLE_ERR(STATUS)
                    HAS_SCALE(J) = .TRUE.
                 ENDIF
                 IF (INDEX(ATT_NAME, 'add_offset') .GT. 0) THEN
                    STATUS = NF_GET_ATT_DOUBLE(NCID, I, 'add_offset',
     &                       NC_OFFSET(J))
                    IF (STATUS .NE. NF_NOERR)
     &                  CALL EC_NCDF_HANDLE_ERR(STATUS)
                    HAS_OFFSET(J) = .TRUE.
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
      ENDDO
C
C Check whether all variables found
C
      OK = (.TRUE.)
      DO I=1,NMax
        IF ((EC_T_STRLEN(NCVarName(I)) .GT. 0) .AND.
     &      (NCVarID(I) .EQ. 0)) THEN
            WRITE(*,5110) NCVarName(I)(:EC_T_STRLEN(NCVarName(I)))
            OK = (.FALSE.)
        ENDIF
      ENDDO
 5110 FORMAT('ERROR: variable ',A,' not found')
      IF (.NOT. OK) STOP

C
C Check whether data have samples
C
      DO I=1,NMax
        IF (NCVarID(I) .NE. 0) THEN
            STATUS = NF_INQ_VAR(NCID,NCVARID(I),
     &                          NAME,XTYPE,VARDIMS,DIMIDS,NATTS)
            IF (STATUS .NE. NF_NOERR) CALL EC_NCDF_HANDLE_ERR(STATUS)
            STATUS = NF_INQ_DIMLEN(NCID, DIMIDS(1), DIMLEN)
            IF (STATUS .NE. NF_NOERR) CALL EC_NCDF_HANDLE_ERR(STATUS)
        IF (DIMLEN .LT. 1) THEN
               WRITE(*,5120) NCVarName(I)(:EC_T_STRLEN(NCVarName(I)))
               STATUS = NF_CLOSE(NCID)
               IF (STATUS .NE. NF_NOERR) CALL EC_NCDF_HANDLE_ERR(STATUS)
               RETURN
        ENDIF
        ENDIF
      ENDDO
 5120 FORMAT('ERROR: no samples in variable ',A)
      

C
C Find start and stop index in file
C
      DOYStart = StartTime(1)
      DOYStop = StopTime(1)
C For security: always use time search routine without seconds:
C is more robust
      CALL EC_NCDF_FINDNOSEC(DoyStart, HMSTART, NCID, NCVarID(QUDoy),
     +                 NCVarID(QUHourMin),
     +                 STARTIND)
      CALL EC_NCDF_FINDNOSEC(DoyStop, HMSTOP, NCID, NCVarID(QUDoy),
     +                 NCVarID(QUHourMin),
     +                 STOPIND)
      NSAMPLE = STOPIND - STARTIND
      IF ((STARTIND .EQ. -1) .OR. (STOPIND .EQ. -1))  NSAMPLE = 0
      IF ((NSAMPLE .EQ. 0) .OR.
     +    (STOPIND .LT. STARTIND)) THEN
          WRITE(*,*) 'ERROR: No data found for given interval'
          STATUS = NF_CLOSE(NCID)
          IF (STATUS .NE. NF_NOERR) CALL EC_NCDF_HANDLE_ERR(STATUS)
          RETURN
      ENDIF
      IF (NSAMPLE .GT. MMMAX) THEN
         WRITE(*,5500) MMMAX, NSAMPLE
         STOP
      ENDIF
      Do I=1,NMax
         IF (STOPIND + Delay(i) .GT. DIMLEN) THEN
            NSAMPLE = STOPIND - Delay(i) - STARTIND
            write(*,*) 'WARNING: Shortened all series, since ', 
     +                 NCVarName(I)(:EC_T_STRLEN(NCVarName(I))), 
     +                ' overruns array by extra delay of ', Delay(i),
     +                ' samples'
         ENDIF
      ENDDO
C
C Read all data from file
C
      OK = (.TRUE.)
C
      M = NSAMPLE
      DO i=1,NMax
C
C Try to read one channel at a time
C
        IF (NCVarID(I) .GT. 0) THEN
           DO J=1,NSAMPLE
              STATUS = NF_GET_VAR1_DOUBLE(NCID,NCVarID(I),
     &                                  STARTIND+J-1+Delay(i),
     &                                  x(i,J))
              OK = (STATUS .EQ. NF_NOERR)
           ENDDO
C
C If something went wrong, then abort further execution of the program
C
           IF (.NOT. ok) THEN
             CALL EC_NCDF_HANDLE_ERR(STATUS)
             STOP
           ENDIF
C
C Apply NetCDF gain and offset first, if present
C
           IF (HAS_SCALE(I)) THEN
              DO J=1,NSAMPLE
                 X(I,J) = X(I,J) * NC_SCALE(I)
              ENDDO
           ENDIF

           IF (HAS_OFFSET(I)) THEN
              DO J=1,NSAMPLE
                 X(I,J) = X(I,J) +  NC_OFFSET(I)
              ENDDO
           ENDIF
C
C Apply gain and offset given in calibration file
C
           DO J=1,NSAMPLE
             x(i,J) = x(i,J)/Gain(i) + Offset(i)
           ENDDO
        ENDIF
      ENDDO


 5500 FORMAT('WARNING: stopped reading data because compiled ',
     &              ' size of storage (',I6,') is unsufficient ',
     &              ' for requested data (',I6,')')
C
C Close file again
C
      STATUS = NF_CLOSE(NCID)
      IF (STATUS .NE. NF_NOERR) CALL EC_NCDF_HANDLE_ERR(STATUS)

      RETURN
      END





      SUBROUTINE EC_F_ReadWou(InName,Delay,Gain,Offset,x,NMax,N,MMax,M)
C
C Read raw data from file
C
C input : InName : CHARACTER*40 : Name of the datafile (binary Wouter format)
C      NMax : INTEGER : Maximum number of quantities in this array
C      MMax : INTEGER : Maximum number of samples in this array
C      N : INTEGER : Number of quantities in file
C
C output : x(NMax,MMax) : REAL*8 Raw data array.
C              First index counts quantities; second counter
C              counts samples. Only the first N quantities
C              and the first M samples are used.
C       M : Number of samples in file
C
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      CHARACTER*40 InName
      INTEGER NMax,N,MMax,M,j,Delay(NMax),MaxDelay,CRT,NHead,NErr
      INTEGER*2 Buffer(NNNMax)
      BYTE ByteBuf(2*NNNMax),C
      CHARACTER Char
      REAL*8 x(NMax,MMax),Gain(NMax),Offset(NMax)
      LOGICAL Ready
      EQUIVALENCE(Buffer,ByteBuf)

      MaxDelay = 0
      DO j=1,N
         MaxDelay = MAX(Delay(j),MaxDelay)
      ENDDO

      OPEN(InFile,FILE=InName,STATUS='OLD',ACCESS='direct',recl=1)
C
C Skip header
C
      CRT = 0
      NHead = 0
      DO WHILE (CRT.LT.3)
        NHead = NHead + 1
        READ(InFile,REC=NHead) Char
        IF (ICHAR(Char).EQ.13) CRT = CRT+1
        ENDDO
      NHead = NHead + 1

      M = 0
      Ready = (.FALSE.)

 10   CONTINUE
      DO j=1,2*N
        READ(InFile,REC=NHead+M*2*N+j,IOSTAT=NErr) ByteBuf(j)
        Ready = (Ready.OR.(NErr.NE.0))
      ENDDO

      IF (.NOT.Ready) THEN
        M = M+1
        DO j=1,N
          C = ByteBuf(2*(j-1)+1)
          ByteBuf(2*(j-1)+1) = ByteBuf(2*(j-1)+2)
          ByteBuf(2*(j-1)+2) = C

      IF ((M-Delay(j)).GT.0)
     &      x(j,M-Delay(j)) = DBLE(Buffer(j))/Gain(j) + Offset(j)
        ENDDO
      ENDIF

      IF ((M.LT.MMMax).AND.(.NOT.READY)) GOTO 10

 20   CLOSE(InFile)
      M = M-MaxDelay

      RETURN
      END

      



      
C     ****f* ec_file.f/EC_F_GetPF
C NAME
C     EC_F_GetPF
C SYNOPSIS
C     CALL  EC_F_GetPF(PlfName, StartTime, StopTime, 
C                      Alpha, Beta, Gamma, 
C                      WBias, Apf)C     
C INPUTS
C     PlfName : [CHARACTER(*)]
C               name of file with planar fit matrix
C     StartTime: [INTEGER(3)]
C               Starttime of interval
C     StoptTime: [INTEGER(3)]
C               Ebdtime of interval
C OUTPUT
C     Alpha    : [REAL*8]
C               first tilt angle (degrees)
C     Beta     : [REAL*8]
C               second tilt angle (degrees)
C     Gamma    : [REAL*8]
C               third tilt angle (degrees)
C     Wbias    : [REAL*8]
C               bias in vertical windspeed (m/s)
C     Apf     : [REAL*8(3,3)]
C               untilt matrix
C     AnglesFound: [LOGICAL]
C               true if angles for this interval found
C FUNCTION
C     Routine to read planar fit until matrix from file.
C     File has probably been produced by program planang
C AUTHOR
C     Arnold Moene
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     parcnst.inc
C     ***
C  
      SUBROUTINE EC_F_GetPF(PlfName, StartTime, StopTime, 
     &                  Alpha, Beta, Gamma, 
     &                  WBias, Apf, AnglesFound)
      IMPLICIT NONE     
      INCLUDE 'parcnst.inc'
      
      CHARACTER*(*) PlfName
      REAL*8        Alpha, Beta, Gamma,
     &              Wbias, Apf(3,3)
      INTEGER       FOO, I, J, StartTime(3), StopTime(3), 
     &              LStartTime(3), LStopTime(3)
      LOGiCAL       PastStart, BeforeEnd, AnglesFound
     
C
C By default we didn't find anything
C
      AnglesFound = .FALSE.
C
C Open Planar fit angles file
C
      OPEN(PlfFile,FILE = PlfName,IOSTAT=FOO,STATUS='Unknown')
      IF (FOO.NE.0) THEN
        WRITE(*,*) 'Could not open planar fit angles-file ',PlfName
        STOP'- Fatal: no planar fit angles could be read -'
      ENDIF
      
C
C Loop over file to see if we have an untilt matrix valid
C for our interval
C
      DO While (.NOT. AnglesFound)
            READ(PlfFile,*,END=37) (LStartTime(i),i=1,3),
     &       (LStopTime(i),i=1,3),Alpha, Beta, Gamma, WBias,
     &       ((Apf(I,J),I=1,3), J=1,3)
C
C Check if we are past start
C
        IF ((StartTime(1) .GT. LStartTime(1)) .OR.
     &      ((StartTime(1) .EQ. LStartTime(1)) .AND.
     &       (StartTime(2) .GT. LStartTime(2))) .OR.
     &      ((StartTime(1) .EQ. LStartTime(1)) .AND.
     &       (StartTime(2) .EQ. LStartTime(2)) .AND.
     &       (StartTime(3) .GE. LStartTime(3)))) THEN
           PastStart = .TRUE.
        ELSE
           PastStart = .FALSE.
        ENDIF
C
C Check if we are before end
C
        IF ((StopTime(1) .LT. LStopTime(1)) .OR.
     &      ((StopTime(1) .EQ. LStopTime(1)) .AND.
     &       (StopTime(2) .LT. LStopTime(2))) .OR.
     &      ((StopTime(1) .EQ. LStopTime(1)) .AND.
     &       (StopTime(2) .EQ. LStopTime(2)) .AND.
     &       (StopTime(3) .LE. LStopTime(3)))) THEN
           BeforeEnd = .TRUE.
        ELSE
           BeforeEnd = .FALSE.
        ENDIF
        IF (PastStart .AND. BeforeEnd) THEN
           AnglesFound = .TRUE.
        ENDIF
      ENDDO
 37   CONTINUE
C     Check if these are valid angles (any of the angles will be invalid; assume that it is sufficient to check Alpha
      IF (AnglesFound) THEN
          IF (INT(Alpha) .EQ. INT(DUMMY)) THEN
             AnglesFound = .FALSE.
          ENDIF
      ENDIF
      IF (.NOT. AnglesFound) THEN
         WRITE(*,*) 'No planar fit angles found for requested interval'
      ENDIF      
C
C     Close the planar fit angles file again
C
      Close(PlfFile)
      END
      
C     ****f* ec_file.f/EC_F_GetPF
C NAME
C     EC_F_WFlux
C SYNOPSIS
C     CALL  EC_F_WFlux(Outunit, Header,
C                      StartTime, StopTime,
C                      M, Mok, Mean, dMean, Cov, dCov,
C                      Phys, dPhys, Std, dStd, Struct,
C                      dStruct, R, dR, DiagFlag, FrCor,
C                      OutMean, OutCov, OutPh, OutTime,
C                      OutStd, OutNum, OutStr, OutFrcor, Outputs)
C INPUTS
C     Outunit   : [INTEGER]
C                 output file unit
c     Header    : [LOGICAL]
C                 write header in this call?
C     StartTime : [REAL*8](3)
C                 starting time of interval
C     StopTime  : [REAL*8](3)
C                 end time of interval
C     M         : [INTEGER]
C                 total number of samples in interval
C     Mok       : [INTEGER](NNMax)
C                 number of valid samples for each calibrated channel
C     Mean      : [REAL*8](NNMax)
C                 mean for each calibrated channel
C     dMean     : [REAL*8](NNMax)
C                 tol. in mean for each calibrated channel
C     Cov       : [REAL*8](NNMax,NNMax)
C                 covariances
C     dCov      : [REAL*8](NNMax,NNMax)
C                 tol. in covariances
C     Phys      : [REAL*8](NMaxPhys)
C                 physical (derived) quantities
C     dPhys     : [REAL*8](NMaxPhys)
C                 tol. in physical (derived) quantities
C     Std       : [REAL*8](NNMax)
C                 standard deviations
C     dStd      : [REAL*8](NNMax)
C                 tol. in standard deviations
C     Struct    : [REAL*8](NNMax,NNMax)
C                 (cross-) structure parameters
C     dStruct   : [REAL*8](NNMax,NNMax)
C                 tol. in (cross-) structure parameters
C     R         : [REAL*8] 
C                 separation in meters at which one wants to estimate
C     dR        : [REAL*8]  
C                 separation in meters corresponding with a delay
C                 of one sample (i.e. tolerance in R)
C     DiagFlag  : [INTEGER](NMaxDiag)
C                 count of error flags of CSAT
C     FrCor     : [REAL*8](NNMax, NNMax)
C                 frequency response correction factors for variances/covariances
C     OutMean   : [LOGICAL](NNMax)
C                 flags for output of means
C     OutCov    : [LOGICAL](NNMax,NNMax)
C                 flags for output of covariances
C     OutPh,    : [LOGICAL](NMaxPhys)
C                 flags for output of physical quantities
C     OutTime,  : [LOGICAL](NMaxOST)
C                 flags for output of time variables
C     OutStd    : [LOGICAL](NNMax)
C                 flags for output of standard deviations
C     OutNum    : [LOGICAL](NNMax)
C                 flags for output of number of samples 
C     OutStr    : [LOGICAL](NNMax,NNMax)
C                 flags for output of structure parameters
C     OutFrcor  : [LOGICAL](NNMax,NNMax)
C                 flags for output of frequency response correction
C     Outputs   : [LOGICAL](NMAxOS)
C                 general output flags
C FUNCTION
C     Routine to write results
C AUTHOR
C     Arnold Moene
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     parcnst.inc
C     EC_T_STRLEN
C     EC_T_STRIPSTR
C     EC_M_DECDAY
C     EC_M_DECHOUR
C     EC_M_DECSEC
C     ***
C
      SUBROUTINE  EC_F_WFlux(Outunit, Header,
     &                 StartTime, StopTime,
     &                 M, Mok, Mean, dMean, Cov, dCov,
     &                 Phys, dPhys, Std, dStd, Struct,
     &                 dStruct, R, dR,DiagFlag, FrCor,
     &                 OutMean, OutCov, OutPh, OutTime,
     &                 OutStd, OutNum, OutStr, OutFrcor, Outputs)
      IMPLICIT NONE
      include 'parcnst.inc'

      INTEGER OutUnit, StartTime(3), StopTime(3), M, Mok(NNMax),
     &        I, J, DiagFlag(NMaxDiag)
      LOGICAL Header, OutMean(NNMax), OutCov(NNMax,NNMax),
     &        OutPh(NMaxPhys), OutTime(NMaxOST), 
     &        OutStd(NNMax), OutNum(NNMax),
     &        OutStr(NNMax,NNMax),
     &        Outputs(NMaxOS), AnyStruct,
     &        OutFrcor(NNmax, NNMax)
      REAL*8  Mean(NNMax), dMean(NNMax),
     &        Cov(NNMax,NNMax), dCov(NNMax,NNMax),
     &        Std(NNMax),dStd(NNMax),
     &        Struct(NNMax,NNMax), dStruct(NNMax,NNMax),
     &        Phys(NMaxPHys), dPhys(NMaxPhys), R, dR,
     &        FrCor(NNMax, NNMax)
      CHARACTER*30  FRM, NAME, NAME1, NAME2, OUTNAME

      INTEGER EC_T_STRLEN
      EXTERNAL EC_T_STRLEN

      REAL*8 EC_M_DECDAY, EC_M_DECHOUR, EC_M_DECSEC
      EXTERNAL EC_M_DECDAY, EC_M_DECHOUR, EC_M_DECSEC

C You really want the old format !!?
      IF (Outputs(OSOLD)) THEN
         IF (Header) THEN
           WRITE(FluxFile,56)
           WRITE(FluxFile,57)
         ELSE
           IF (M .GT. 0) THEN
             WRITE(OutUnit,55)
     &       (StartTime(i),i=1,3),
     &       (StopTime(i),i=1,3),
     &       M,(Mok(i),i=1,7), Mok(TTime),
     &       NINT(Phys(QPDirFrom)), NINT(dPhys(QPDirFrom)),
     &       Phys(QPVectWind), dPhys(QPVectWind),
     &       Mean(TSonic),dMean(TSonic),
     &       Mean(TCouple),dMean(TCouple),
     &       Mean(SpecHum),dMean(SpecHum),
     &       Std(Tsonic), dStd(TSonic),
     &       Std(Tcouple), dStd(Tcouple),
     &       Std(SpecHum), dStd(SpecHum),
     &       Std(U), dStd(U),
     &       Std(V), dStd(V),
     &       Std(W), dStd(W),
     &       Cov(TSonic,SpecHum),dCov(TSonic,SpecHum),
     &       Cov(TCouple,SpecHum),dCov(TCouple,SpecHum),
     &       Cov(TSonic,U),dCov(TSonic,U),
     &       Cov(TCouple,U),dCov(TCouple,U),
     &       Cov(SpecHum,U),dCov(SpecHum,U),
     &       Phys(QPHSonic),dPhys(QPHSonic),
     &       Phys(QPHTc),dPhys(QPHTc),
     &       Phys(QPLvE),dPhys(QPLvE),
     &       Phys(QPLvEWebb),dPhys(QPLvEWebb),
     &       Phys(QPSumLvE), dPhys(QPSumLvE),
     &       Phys(QPUstar),dPhys(QPUstar),
     &       Phys(QPTau),dPhys(QPTau),
     &       R,dR, Struct(TSonic,TSonic), dStruct(TSonic,TSonic),
     &       Struct(TCouple,TCouple), dStruct(TCouple,TCouple),
     &       Struct(SpecHum,SpecHum), dStruct(SpecHum,SpecHum),
     &       Struct(TSonic,SpecHum), dStruct(TSonic,SpecHum),
     &       Struct(TCouple,SpecHum), dStruct(TCouple,SpecHum),
     &       Phys(QPMeanW),dPhys(QPMeanW),
     &       Mok(CO2),
     &       Mean(CO2), dMean(CO2),
     &       Mean(specCO2), dMean(specCO2),
     &       Std(CO2), dStd(CO2),
     &       Std(specCO2), dStd(specCO2),
     &       Phys(QPFCO2),dPhys(QPFCO2),
     &       Phys(QPFCO2Webb),dPhys(QPFCO2Webb),
     &       Phys(QPSumFCO2), dPhys(QPSumFCO2),
     &       DiagFlag(QDDelta), DiagFlag(QDLock), 
     &       DiagFlag(QDHigh), DiagFlag(QDLow),
     &       (Mean(I), I=1,3),
     &       ((COV(I,J),I=1,3),J=1,3)  
           ELSE
             WRITE(OutUnit,55)
     &       (StartTime(i),i=1,3),
     &       (StopTime(i),i=1,3),
     &       M,(0,i=1,8), INT(DUMMY),INT(DUMMY),(DUMMY,i=1,58),
     &       0, (DUMMY,i=1,14), (INT(DUMMY), i=1,4),
     &       ((DUMMY,I=1,3),J=1,3), (DUMMY,I=1,3)
           ENDIF
        ENDIF
C Great, you want the NEW format
      ELSE
C Check if any structure parameter needs to be written
        AnyStruct = .FALSE.
        DO I=1,NNMax
           DO J=1,NNMax
              IF (OutStr(I,J)) AnyStruct = .TRUE.
           ENDDO
        ENDDO
        IF (HEADER) THEN
C First line of header
C Special format for first column
           IF (Outputs(OSCSV)) THEN
              FRM='(A,$)'
           ELSE
              FRM='(A,$)'
           ENDIF
           write(OutUnit, FRM) 'Doy'
C General format for other columns
           IF (Outputs(OSCSV)) THEN
              FRM='(",",A,$)'
           ELSE
              FRM='(" ",A20,$)'
           ENDIF
C Rest of time info
           write(OutUnit, FRM) 'Hour'
           write(OutUnit, FRM) 'Min'
           write(OutUnit, FRM) 'Doy'
           write(OutUnit, FRM) 'Hour'
           write(OutUnit, FRM) 'Min'
           DO I=1,NMaxOST
             IF (OutTime(I)) THEN
                NAME = QTName(I)
                CALL EC_T_STRIPSTR(NAME)
                WRITE(OUTNAME, "(A)") NAME(:EC_T_STRLEN(NAME))
                write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
             ENDIF
           ENDDO
C Number of samples
           write(OutUnit, FRM) '#Samples'
           DO I=1,NNMax
             IF (OutNum(I)) THEN
                NAME = QName(I)
                CALL EC_T_STRIPSTR(NAME)
                WRITE(OUTNAME, "('#',A)") NAME(:EC_T_STRLEN(NAME))
                write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
             ENDIF
           ENDDO
C Diagnostics
           IF (Outputs(OSDiag)) THEN
              DO I=1,NMaxDiag
                NAME = QDName(I)
                CALL EC_T_STRIPSTR(NAME)
                WRITE(OUTNAME, "(A)")
     &                NAME(:EC_T_STRLEN(NAME))
                write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
              ENDDO
           ENDIF
C Means
           DO I=1,NNMax
             IF (OutMean(I)) THEN
                NAME = QName(I)
                CALL EC_T_STRIPSTR(NAME)
                WRITE(OUTNAME, "('Mean(',A,')')")
     &                NAME(:EC_T_STRLEN(NAME))
                write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                IF (Outputs(OSTol)) THEN
                   NAME = QName(I)
                   CALL EC_T_STRIPSTR(NAME)
                   WRITE(OUTNAME, "('TolMean(',A,')')")
     &                NAME(:EC_T_STRLEN(NAME))
                   write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                ENDIF
             ENDIF
           ENDDO
C Standard deviations
           DO I=1,NNMax
             IF (OutStd(I)) THEN
                NAME = QName(I)
                CALL EC_T_STRIPSTR(NAME)
                WRITE(OUTNAME, "('Std(',A,')')")
     &                NAME(:EC_T_STRLEN(NAME))
                write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                IF (Outputs(OSTol)) THEN
                   WRITE(OUTNAME, "('TolStd(',A,')')")
     &                NAME(:EC_T_STRLEN(NAME))
                    write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                ENDIF
             ENDIF
           ENDDO
C Covariances
           DO I=1,NNMax
             DO J=1,NNMax
                IF (OutCov(I,J)) THEN
                   NAME1 = QName(I)
                   NAME2 = QName(J)
                   CALL EC_T_STRIPSTR(NAME1)
                   CALL EC_T_STRIPSTR(NAME2)
                   WRITE(OUTNAME, "('Cov(',A,'*',A,')')")
     &                NAME1(:EC_T_STRLEN(NAME1)),
     &                NAME2(:EC_T_STRLEN(NAME2))
                   write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                   IF (Outputs(OSTol)) THEN
                      WRITE(OUTNAME, "('TolCov(',A,'*',A,')')")
     &                   NAME1(:EC_T_STRLEN(NAME1)),
     &                   NAME2(:EC_T_STRLEN(NAME2))
                      write(OutUnit, FRM)
     &                  OUTNAME(:EC_T_STRLEN(OUTNAME))
                   ENDIF
                ENDIF
             ENDDO
           ENDDO
C Structure parameters
           IF (AnyStruct) THEN
             write(OutUnit, FRM) 'R'
             write(OutUnit, FRM) 'dR'
             DO I=1,NNMax
                DO J=1,NNMax
                  IF (OutStr(I,J)) THEN
                     NAME1 = QName(I)
                     NAME2 = QName(J)
                     CALL EC_T_STRIPSTR(NAME1)
                     CALL EC_T_STRIPSTR(NAME2)
                     WRITE(OUTNAME, "('C_(',A,'*',A,')')")
     &                  NAME1(:EC_T_STRLEN(NAME1)),
     &                  NAME2(:EC_T_STRLEN(NAME2))
                     write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                     IF (Outputs(OSTol)) THEN
                        WRITE(OUTNAME, "('TolC_(',A,'*',A,')')")
     &                     NAME1(:EC_T_STRLEN(NAME1)),
     &                     NAME2(:EC_T_STRLEN(NAME2))
                        write(OutUnit, FRM)
     &                    OUTNAME(:EC_T_STRLEN(OUTNAME))
                     ENDIF
                  ENDIF
                ENDDO
             ENDDO
           ENDIF
C Physical quantities
           DO I=1,NMaxPhys
             IF (OutPh(I)) THEN
                NAME = QPName(I)
                CALL EC_T_STRIPSTR(NAME)
                WRITE(OUTNAME, "(A)")
     &                NAME(:EC_T_STRLEN(NAME))
                write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                IF (Outputs(OSTol)) THEN
                   WRITE(OUTNAME, "('Tol(',A,')')")
     &                NAME(:EC_T_STRLEN(NAME))
                   write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                ENDIF
             ENDIF
           ENDDO
C Frequency response factors
           DO I=1,NNMax
             DO J=1,NNMax
                IF (OutFrcor(I,J)) THEN
                   NAME1 = QName(I)
                   NAME2 = QName(J)
                   CALL EC_T_STRIPSTR(NAME1)
                   CALL EC_T_STRIPSTR(NAME2)
                   WRITE(OUTNAME, "('FrC(',A,'*',A,')')")
     &                NAME1(:EC_T_STRLEN(NAME1)),
     &                NAME2(:EC_T_STRLEN(NAME2))
                   write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                ENDIF
             ENDDO
           ENDDO
C New line
           WRITE(OutUnit,*)
C Second line of header
C Special format for first column
           IF (Outputs(OSCSV)) THEN
              FRM='(A,$)'
           ELSE
              FRM='(A,$)'
           ENDIF
           write(OutUnit, FRM) '[-]'
C General format for other columns
           IF (Outputs(OSCSV)) THEN
              FRM='(",",A,$)'
           ELSE
              FRM='(" ",A20,$)'
           ENDIF
C Rest of time info
           write(OutUnit, FRM) '[-]'
           write(OutUnit, FRM) '[-]'
           write(OutUnit, FRM) '[-]'
           write(OutUnit, FRM) '[-]'
           write(OutUnit, FRM) '[-]'
C Extra time columns
           DO I=1,NMaxOST
             IF (OutTime(I)) THEN
                NAME = UTName(I)
                CALL EC_T_STRIPSTR(NAME)
                WRITE(OUTNAME, "(A)") NAME(:EC_T_STRLEN(NAME))
                write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
             ENDIF
           ENDDO
C # samples
           write(OutUnit, FRM) '[-]'
           DO I=1,NNMax
             if (OutNum(I)) THEN
                 write(OutUnit, FRM) '[-]'
             ENDIF
           ENDDO
C Diagnostics
           IF (Outputs(OSDiag)) THEN
              DO I=1,NMaxDiag
                NAME = UDName(I)
                CALL EC_T_STRIPSTR(NAME)
                WRITE(OUTNAME, "(A)")
     &                NAME(:EC_T_STRLEN(NAME))
                write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
              ENDDO
           ENDIF
C Means
           DO I=1,NNMax
             IF (OutMean(I)) THEN
                NAME = UName(I)
                CALL EC_T_STRIPSTR(NAME)
                WRITE(OUTNAME, "(A)")  NAME(:EC_T_STRLEN(NAME))
                write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                IF (Outputs(OSTol)) THEN
                   write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                ENDIF
             ENDIF
           ENDDO
C Standard deviations
           DO I=1,NNMax
             IF (OutStd(I)) THEN
                NAME = UName(I)
                CALL EC_T_STRIPSTR(NAME)
                WRITE(OUTNAME, "(A)")  NAME(:EC_T_STRLEN(NAME))
                write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                IF (Outputs(OSTol)) THEN
                    write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                ENDIF
             ENDIF
           ENDDO
C Covariances
           DO I=1,NNMax
             DO J=1,NNMax
                IF (OutCov(I,J)) THEN
                   NAME1 = UName(I)
                   NAME2 = UName(J)
                   CALL EC_T_STRIPSTR(NAME1)
                   CALL EC_T_STRIPSTR(NAME2)
                   WRITE(OUTNAME, "(A,A)")
     &                NAME1(:EC_T_STRLEN(NAME1)),
     &                NAME2(:EC_T_STRLEN(NAME2))
                   write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                   IF (Outputs(OSTol)) THEN
                      write(OutUnit, FRM)
     &                  OUTNAME(:EC_T_STRLEN(OUTNAME))
                   ENDIF
                ENDIF
             ENDDO
           ENDDO
C Structure parameters
           IF (AnyStruct) THEN
             write(OutUnit, FRM) '[m]'
             write(OutUnit, FRM) '[m]'
             DO I=1,NNMax
                DO J=1,NNMax
                  IF (OutStr(I,J)) THEN
                     NAME1 = UName(I)
                     NAME2 = UName(J)
                     CALL EC_T_STRIPSTR(NAME1)
                     CALL EC_T_STRIPSTR(NAME2)
                     WRITE(OUTNAME, "(A,A,A)")
     &                  NAME1(:EC_T_STRLEN(NAME1)),
     &                  NAME2(:EC_T_STRLEN(NAME2)),
     &                  '[m^-2/3]'
                     write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                     IF (Outputs(OSTol)) THEN
                        write(OutUnit, FRM)
     &                    OUTNAME(:EC_T_STRLEN(OUTNAME))
                     ENDIF
                  ENDIF
                ENDDO
             ENDDO
           ENDIF
C Physical quantities
           DO I=1,NMaxPhys
             IF (OutPh(I)) THEN
                NAME = UPName(I)
                CALL EC_T_STRIPSTR(NAME)
                WRITE(OUTNAME, "(A)")
     &                NAME(:EC_T_STRLEN(NAME))
                write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                IF (Outputs(OSTol)) THEN
                   write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                ENDIF
             ENDIF
           ENDDO
C Covariances
           DO I=1,NNMax
             DO J=1,NNMax
                IF (OutFrcor(I,J)) THEN
                   NAME1 = UName(I)
                   NAME2 = UName(J)
                   CALL EC_T_STRIPSTR(NAME1)
                   CALL EC_T_STRIPSTR(NAME2)
                   WRITE(OUTNAME, "(A,A)")
     &                NAME1(:EC_T_STRLEN(NAME1)),
     &                NAME2(:EC_T_STRLEN(NAME2))
                   write(OutUnit, FRM) OUTNAME(:EC_T_STRLEN(OUTNAME))
                ENDIF
             ENDDO
           ENDDO
C New line
           WRITE(OutUnit,*)
C The real data!
        ELSE
           write(OUTNAME, '(I3)') StartTime(1)
           CALL EC_T_STRIPSTR(OUTNAME)
           WRITE(OutUnit,'(A,$)') OUTNAME(:EC_T_STRLEN(OUTNAME))
C General format for other columns
           IF (Outputs(OSCSV)) THEN
              FRM='(",",I6,$)'
           ELSE
              FRM='(" ",I20,$)'
           ENDIF
C Time info
           write(OutUnit, FRM) StartTime(2)
           write(OutUnit, FRM) StartTime(3)
           write(OutUnit, FRM) StopTime(1)
           write(OutUnit, FRM) StopTime(2)
           write(OutUnit, FRM) StopTime(3)
C Extra time columns
           IF (Outputs(OSCSV)) THEN
              FRM='(",",G16.8,$)'
           ELSE
              FRM='(" ",G20.8,$)'
           ENDIF
           DO I=1,NMaxOST
             IF (OutTime(I)) THEN
                IF (I .EQ. OSTDecDayB)
     &             write(OutUnit, FRM) 
     &             EC_M_DecDay(StartTime(1), StartTime(2), StartTime(3))
                IF (I .EQ. OSTDecDayE)
     &             write(OutUnit, FRM) 
     &             EC_M_DecDay(StopTime(1), StopTime(2), StopTime(3))
                IF (I .EQ. OSTDecHrB)
     &             write(OutUnit, FRM) 
     &             EC_M_DecHour(StartTime(2), StartTime(3))
                IF (I .EQ. OSTDecHrE)
     &             write(OutUnit, FRM) 
     &             EC_M_DecHour(StopTime(2), StopTime(3))
                IF (I .EQ. OSTDecSecB)
     &             write(OutUnit, FRM) 
     &             EC_M_DecSec(StartTime(2), StartTime(3))
                IF (I .EQ. OSTDecSecE)
     &             write(OutUnit, FRM) 
     &             EC_M_DecSec(StopTime(2), StopTime(3))
             ENDIF
           ENDDO
C Number of samples
           IF (Outputs(OSCSV)) THEN
              FRM='(",",I6,$)'
           ELSE
              FRM='(" ",I20,$)'
           ENDIF
           write(OutUnit, FRM) M
           DO I=1,NNMax
             IF (OutNum(I)) THEN
                CALL EC_F_WINT(OutUnit, FRM, Mok(I), M)
             ENDIF
           ENDDO
           IF (Outputs(OSDiag)) THEN
              IF (Outputs(OSCSV)) THEN
                 FRM='(",",I6,$)'
              ELSE
                 FRM='(" ",I20,$)'
              ENDIF
              DO I=1,NMaxDiag
                CALL EC_F_WINT(OutUnit, FRM, DiagFlag(I), M)
              ENDDO
           ENDIF
           IF (Outputs(OSCSV)) THEN
              FRM='(",",G16.8,$)'
           ELSE
              FRM='(" ",G20.5,$)'
           ENDIF
C Means
           DO I=1,NNMax
             IF (OutMean(I)) THEN
                CALL EC_F_WREAL(OutUnit, FRM, Mean(I), M)
                IF (Outputs(OSTol)) THEN
                   CALL EC_F_WREAL(OutUnit, FRM, dMean(I), M)
               ENDIF
             ENDIF
           ENDDO
C Standard deviations
           DO I=1,NNMax
             IF (OutStd(I)) THEN
                CALL EC_F_WREAL(OutUnit, FRM, Std(I), M)
                IF (Outputs(OSTol)) THEN
                    CALL EC_F_WREAL(OutUnit, FRM, dStd(I), M)
                ENDIF
             ENDIF
           ENDDO
C Covariances
           DO I=1,NNMax
             DO J=1,NNMax
                IF (OutCov(I,J)) THEN
                   CALL EC_F_WREAL(OutUnit, FRM, Cov(I,J), M)
                   IF (Outputs(OSTol)) THEN
                      CALL EC_F_WREAL(OutUnit, FRM, dCov(I,J), M)
                   ENDIF
                ENDIF
             ENDDO
           ENDDO
C Structure parameters
           IF (AnyStruct) THEN
             CALL EC_F_WREAL(OutUnit, FRM, R, M)
             CALL EC_F_WREAL(OutUnit, FRM, dR, M)
             DO I=1,NNMax
                DO J=1,NNMax
                  IF (OutStr(I,J)) THEN
                     CALL EC_F_WREAL(OutUnit, FRM, Struct(I,J), M)
                     IF (Outputs(OSTol)) THEN
                        CALL EC_F_WREAL(OutUnit, FRM, dStruct(I,J), M)
                     ENDIF
                  ENDIF
                ENDDO
             ENDDO
           ENDIF
C Physical quantities
           DO I=1,NMaxPhys
             IF (OutPh(I)) THEN
                CALL EC_F_WREAL(OutUnit, FRM, Phys(I), M)
                IF (Outputs(OSTol)) THEN
                   CALL EC_F_WREAL(OutUnit, FRM, dPhys(I), M)
                ENDIF
             ENDIF
           ENDDO
C Frequency response correction factors
           DO I=1,NNMax
             DO J=1,NNMax
                IF (OutFrcor(I,J)) THEN
                   CALL EC_F_WREAL(OutUnit, FRM, FrCor(I,J), M)
                ENDIF
             ENDDO
           ENDDO
C New line
           WRITE(OutUnit,*)
        ENDIF
      ENDIF
 55   FORMAT(2(I3,1X,2(I2,1X)),9(I6,1X),2(I5,1X),
     &        29(2(G15.5:,1X)),
     &        I13,  1X,
     &        7(2(G15.5:,1X)), 4(I15,1X), 12(G15.5:,1X))
 56   FORMAT(
     &  'DOY Hr Mn  ',
     &  'DOY Hr Mn ',
     &  '#samples',
     &  '  #U     #V     #W  #TSon  #TCop  #RhoV     #q  #time   ',
     &  'Dir d(dir)   ',
     &  'Mean(vectorU)  dMean(vectorU)   ',
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
     &  '#AmplHigh      #AmplLow         ',
     &  'Mean(U)        Mean(V)          ',
     &  'Mean(W)        Cov(U,U)         ',
     &  'Cov(U,V)       Cov(U,W)         ',
     &  'Cov(V,U)       Cov(V,V)         ',
     &  'Cov(V,W)       Cov(W,U)         ',
     &  'Cov(W,V)       Cov(W,W)         ')
  57   FORMAT(     
     &  '-  -  -   ',
     &  '-  -  -  ',
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
     &  '[-]            [-]              ',
     &  '[m/s]          [m/s]            ',
     &  '[m/s]          [m/s]**2         ',
     &  '[m/s]**2       [m/s]**2         ',
     &  '[m/s]**2       [m/s]**2         ',     
     &  '[m/s]**2       [m/s]**2         ',
     &  '[m/s]**2       [m/s]**2         ')     
      END


      SUBROUTINE EC_F_WREAL(OUTUNIT, FRM, NUMBER, NSAMP)
      IMPLICIT NONE
      INCLUDE 'parcnst.inc'
      REAL*8 NUMBER
      CHARACTER*(*) FRM
      INTEGER NSAMP, OUTUNIT

      IF (NSAMP .GT. 0) THEN
          write(OutUnit, FRM) NUMBER
      ELSE
          write(OutUnit, FRM) DUMMY
      ENDIF
      END

      SUBROUTINE EC_F_WINT(OUTUNIT, FRM, NUMBER, NSAMP)
      IMPLICIT NONE
      INCLUDE 'parcnst.inc'

      CHARACTER*(*) FRM
      INTEGER NSAMP, OUTUNIT, NUMBER

      IF (NSAMP .GT. 0) THEN
          write(OutUnit, FRM) NUMBER
      ELSE
          write(OutUnit, FRM) INT(DUMMY)
      ENDIF
      END
