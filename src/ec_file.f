C $Id$
C  ec_file.f, part of ecpack
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

      SUBROUTINE EC_F_GetConf(ConfUnit,
     &           DatDir, OutDir,  ParmDir,
     &           FluxName, ParmName, InterName, PlfIntName,
     &           PlfName,
     &           SonName, CoupName, HygName, CO2Name,
     &           NCVarName,
     &           NCNameLen)
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'
      
      
      INTEGER         ConfUnit, NCNameLen, I
      CHARACTER*(*)   DatDir, OutDir, ParmDir, FluxName, ParmName,
     &                InterName, SonName, CoupName, HygName, 
     &                CO2Name, PlfName, PlfIntName
      CHARACTER*(*)   NCVarName(NCNameLen)

      INTEGER         IOCODE, KINDEX
      CHARACTER*255   LINE, TOKLINE, VALLINE, DUMSTRING
      INTEGER         EC_T_STRLEN
      EXTERNAL        EC_T_STRLEN


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
      CALL EC_T_CLEARSTR(PlfName)
      
      IOCODE = 0
      DO 4000, WHILE (IOCODE .EQ. 0)
C     Read one line from configuration file
         READ(Confunit, '(A)', IOSTAT=IOCODE, END=9000) LINE
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
                WRITE(*,*) 'ERROR no equality sign found'
                STOP
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


      END
C###########################################################################
C
C
C
C      EEEEE   CCCC        PPPPP        A         CCCC   K    K
C      E      C    C       P    P      A A       C    C  K   K
C      E      C	           P     P    A   A     C	 K  K
C      EEEE   C     ------ PPPPPP    A     A    C	 K K
C      E      C		   P	    AAAAAAAAA   C	 KK K
C      E      C    C	   P	   A         A   C    C  K   K
C      EEEEE   CCCC 	   P      A	      A   CCCC   K    K
C
C	     Library for processing Eddy-Correlation data
C     EC Special Interest Group of Wag-UR-METAIR Wageningen and KNMI
C
C
C
C Version of release    : 1.10
C (note that version numbers of subroutines are not maintained)
C Date	   : September 26 2000
C Author		: Arjan van Dijk
C             Arnold Moene
C For : EC Special Interest Group of Wag-UR-METAIR Wageningen
C	and KNMI (Henk de Bruin, Arjan van Dijk, Wim Kohsiek,
C	Fred Bosveld, Cor Jacobs and Bart van de Hurk)
C Contact address	: Duivendaal 2
C			  6701 AP  Wageningen
C			  The Netherlands
C			  WWW.MetAir.WAU.nl
C			  WWW.MetAir.Wag-UR.nl (in due time...)
C			  Tel. +31 317 483981 (secretary)
C			  Fax. +31 317 482811
C
C Purpose :
C
C The aim of the special interest group is to develop a standard processing
C method for eddycorrelation data. With such a standard method results will
C become better comparable. All important corrections can be carried out.
C These are:
C
C  - Correction of sonic-temperature for lateral velocity (deferred to
C    calibration routine).
C  - Correction of sonic-temperature for presence of humidity.
C  - Correction of sonic path length to make mean sonic temperature match
C    mean temperature according to a different device (e.g. thermocouple).
C  - Time-delay in the datafile between the different columns
C  - Yaw-correction to make Mean(V) --> 0
C  - Pitch-correction to make Mean(W) --> 0
C  - Roll-correction to make Cov(W,V) --> 0
C  - Frequency-reponse corrections for slow apparatus and path length
C    integration.
C  - Oxygen-correction for hygrometers which are sensitive for O2.
C  - Inclusion of the mean vertical velocity according to Webb.
C
C If you process your data using EC-Pack, please make a reference
C of the kind: "This data has been processed using EC-PACK version XX".
C This will make your data better exchangeable and comparable.
C
C If you have suggestions for additional functionalities, comments
C or questions concerning EC-Pack, please contact us at the above
C address. To prevent the formation of EC-Pack dialects, and to bring
C bright insights to the attention and disposal of fellow researchers:
C		  Share your knowledge via us!
C
C DISCLAIMER/COPYRIGHT: The routines in EC-Pack are provided for free
C and may freely be copied. We would appreciate acknowledgement in your
C publications. Furthermore: EC-Pack is provided as is. No responsibility
C can be taken for consequences of the use of EC-Pack. Any use of EC-Pack
C is done at your own risk.
C
C
C Developments :
C - Better accessibility of frequency-response corrections
C - Separation of calibration routine from library
C - Detection of different trends and subsequent handling
C - Better documentation of all input/output variables and functionality
C   of subroutines
C
C
C###########################################################################

C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C
C General routines
C
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################

      SUBROUTINE EC_F_Params(InName, ExpVar,
     &	DoCorr, PCorr,
     &	DoStruct, DoPrint,
     &	PRaw,PCal,PIndep)
C
C Purpose : Read settings for this particular session from file
C
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      REAL*8 ExpVar(NMaxExp)
      LOGICAL 
     &	DoPrint,PRaw,PCal,PIndep,
     &	DoStruct, DoCorr(NMaxCorr), PCorr(NMaxCorr)
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
      ExpVar(QEPFValid) = 0.0D0
      READ(TempFile,*) ExpVar(QEFreq)	  ! [Hz] Sample frequency datalogger
      READ(TempFile,*)
      READ(TempFile,*) ExpVar(QEPitchLim)	  ! [degree] Limit when Mean(W) is turned to zero
      READ(TempFile,*) ExpVar(QERollLim)	  ! [degree] Limit when Cov(V,W) is turned to zero
      READ(TempFile,*)
      READ(TempFile,*) ExpVar(QEPreYaw)   ! Fixed yaw angle for known tilt-correction
      READ(TempFile,*) ExpVar(QEPrePitch)	  ! Fixed pitch angle for known tilt-correction
      READ(TempFile,*) ExpVar(QEPreRoll)	  ! Fixed roll angle for known tilt-correction
      READ(TempFile,*)
      READ(TempFile,*) ExpVar(QELLimit)	  ! Smallest acceptable frequency-response
correction factor
      READ(TempFile,*) ExpVar(QEULimit)	  ! Largest acceptable frequency-response
correction factor
      READ(TempFile,*)
      READ(TempFile,*) DoCorr(QCMean)     ! Replace mean quantities by better estimates
      READ(TempFile,*) DoCorr(QCDetrend)  ! Correct data for linear trend
      READ(TempFile,*) DoCorr(QCSonic)	  ! Correct sonic temperature for humidity
      READ(TempFile,*) DoCorr(QCTilt)	  ! Perform true tilt-correction with known angles
      READ(TempFile,*) DoCorr(QCYaw)	  ! Turn system such that Mean(V) --> 0
      READ(TempFile,*) DoCorr(QCPitch)	  ! Turn system such that Mean(W) --> 0
      READ(TempFile,*) DoCorr(QCRoll)	  ! Turn System such that Cov(W,V) --> 0
      READ(TempFile,*) DoCorr(QCFreq)	  ! Correct for poor frequency response
      READ(TempFile,*) DoCorr(QCO2)	  ! Correct hygrometer for oxygen-sensitivity
      READ(TempFile,*) DoCorr(QCWebb)	  ! Calculate mean velocity according to Webb
      READ(TempFile,*) DoStruct	  ! Calculate structure parameters
      READ(TempFile,*)
      READ(TempFile,*) DoPrint	  ! Skip printing intermediate results or not?
      READ(TempFile,*)
      READ(TempFile,*) PRaw	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCal	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCDetrend)   ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PIndep	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCTilt)	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCYaw)	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCPitch)	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCRoll)	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCSonic)	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCO2)	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCFreq)	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCorr(QCWebb)	  ! Indicator if these intermediate results are wanted
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

      OPEN(TempFile,FILE=InName,STATUS='OLD', IOSTAT=IOCODE)
      IF (IOCODE .NE. 0) THEN
         WRITE(*,5100) InName(:EC_T_STRLEN(InName))
         STOP
      ENDIF
 5100 FORMAT ('ERROR: can not open calibration file ', (A))

C First get the type of apparatus
      READ(TempFile,*) CalSpec(1)

C Now read the correct number of specs
      DO i=2,ApNQQ(INT(CalSpec(1)))
       	READ(TempFile,*,IOSTAT=IOCODE, END = 9000) CalSpec(i)
         IF (IOCODE .NE. 0) THEN
            WRITE(*,*) 'ERROR in reading of configuration file'
            STOP
         ENDIF
      ENDDO
 9000 CONTINUE
      IF ((I-1) .NE. ApNQQ(INT(CalSpec(1)))) THEN
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
C	  NMax : INTEGER : Maximum number of quantities in this array
C	  MMax : INTEGER : Maximum number of samples in this array
C	  N : INTEGER : Number of quantities in file
C
C output : x(NMax,MMax) : REAL*8 Raw data array.
C			  First index counts quantities; second counter
C			  counts samples. Only the first N quantities
C			  and the first M samples are used.
C	   M : Number of samples in file
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
     &    x(j,M-Delay(j)) = Dum(j)/Gain(j) + Offset(j)
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
C	  NMax : INTEGER : Maximum number of quantities in this array
C	  MMax : INTEGER : Maximum number of samples in this array
C
C output : x(NMax,MMax) : REAL*8 Raw data array.
C			  First index counts quantities; second counter
C			  counts samples. Only the first N quantities
C			  and the first M samples are used.
C	   M : Number of samples in file
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
C	  NMax : INTEGER : Maximum number of quantities in this array
C	  MMax : INTEGER : Maximum number of samples in this array
C	  N : INTEGER : Number of quantities in file
C
C output : x(NMax,MMax) : REAL*8 Raw data array.
C			  First index counts quantities; second counter
C			  counts samples. Only the first N quantities
C			  and the first M samples are used.
C	   M : Number of samples in file
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
      IF (.NOT. AnglesFound) THEN
         WRITE(*,*) 'No planar fit angles found for requested interval'
      ENDIF      
      END
