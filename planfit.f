      PROGRAM PLANFIT

      IMPLICIT None
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'
      INCLUDE 'calcomm.inc'

C
C Test of the planar fit method for tilt-correction by Wilczak
C
      CHARACTER*255 FNAME, DatDir, OutDir, ParmDir, FluxName,
     &              ParmName,InterName, PlfIntName, PlfName,
     &              OutName
      CHARACTER*255 SonName,CoupName,HygName, DUMSTRING
      CHARACTER*255 NCvarname(NNNMax)
      INTEGER N,M,MIndep(NNMax),CIndep(NNMax,NNMax),FOO,
     &  Channels,Delay(NNNMax),Mok(NNMax),Cok(NNMax,NNMax)

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
     &  StructSep
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ),
     &  R,dR,CTSon2,CTCop2,Cq2,CTSonq,CTCopq,
     &  dCTSon2,dCTCop2,dCq2,dCTSonq,dCTCopq



      INTEGER RunLength,NRuns,  MaxRuns
      PARAMETER (MaxRuns=1000)
      INTEGER i,j,k,l
      REAL*8 UMean(3, MaxRuns),Apf(3,3),Alpha,Beta,Gamma,WBias,Dum(4),
     &  UUSum(3,3),RApf(3,3),RAlpha,RBeta,RGamma,RWBias,USum(3),
     &  Swing,SinSwing,CosSwing,UHor, RunStart(3, MaxRuns),
     &  RunStop(3, MaxRuns)

      LOGICAL       HAVE_UNCAL(NNNMax),Flag(NNMAx,MMMax)
      LOGICAL ENDINTERVAL, BadTc, DoPrint, HAVE_SAMP
      INTEGER STRLEN
C
C The name of the calibration routine specified at the end of this file
C
      EXTERNAL Calibrat, STRLEN


C...........................................................................
C Start of executable statements
C...........................................................................
C Set DoPrint to true
      DoPrint = .TRUE.
C
C Open and read configuration file
C

      IF (ECConfFile .NE. 5) THEN
         OPEN(ECConfFile, FILE='ecconf')
      ENDIF

C Report version number (VERSION is defined in version.inc)
C Give some RCS info (do not edit this!!, RCS does it for us) 
      WRITE(*,*) '$Name$'
      WRITE(*,*) '$Date$'
      WRITE(*,*) '$Revision$'
      
      CALL GetConf(ECConfFile,
     &             DatDir, OutDir, ParmDir,
     &             FluxName, ParmName, InterName,
     &             PlfName,
     &             SonName, CoupName, HygName,
     &             NCVarname, NNNMax)
C
C Check which calibration files we have
C
      IF (STRLEN(SonName) .GT. 0) THEN
         HAVE_SONCAL = .TRUE.
      ELSE
         HAVE_SONCAL = .FALSE.
      ENDIF
      IF (STRLEN(HygName) .GT. 0) THEN
         HAVE_HYGCAL = .TRUE.
      ELSE
         HAVE_HYGCAL = .FALSE.
      ENDIF
      IF (STRLEN(CoupName) .GT. 0) THEN
         HAVE_TCCAL = .TRUE.
      ELSE
         HAVE_TCCAL = .FALSE.
      ENDIF

C
C Read calibration data from files :
C
      DO i=1,NNNMax
        Delay( i) = 0
        Gain(  i) = 1.D0
        Offset(i) = 0.D0
      ENDDO
      IF (HAVE_SONCAL)
     &     CALL ECReadAp(SonName,CalSonic) ! The sonic
      IF (HAVE_TCCAL)
     &     CALL ECReadAp(CoupName,CalTherm) ! A thermo-couple
      IF (HAVE_HYGCAL)
     &     CALL ECReadAp(HygName,CalHyg) ! A Krypton hygrometer
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

      Gain(U)        = CalSonic(QQGain)
      Gain(V)        = CalSonic(QQGain)
      Gain(W)        = CalSonic(QQGain)
      Gain(Humidity) = CalHyg(QQGain)
      Gain(Tcouple)  = CalTherm(QQGain)
      Gain(Tsonic)   = CalSonic(QQGain)
      
      Offset(U)        = CalSonic(QQOffset)
      Offset(V)        = CalSonic(QQOffset)
      Offset(W)        = CalSonic(QQOffset)
      Offset(Humidity) = CalHyg(QQOffset)
      Offset(Tcouple)  = CalTherm(QQOffset)
      Offset(Tsonic)   = CalSonic(QQOffset)
      
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
C Now cycle the interval file to read names of intermediate files
C and clear those files
      ENDINTERVAL = .FALSE.
      DO WHILE (.NOT. ENDINTERVAL)
         READ(IntervalFile,*,END=900) (StartTime(i), i=1,3),
     &        (StopTime(i),i=1,3), Psychro, p, Fname, OutName,
     &        PlfIntName
         DUMSTRING =  OutDir
         CALL STRCAT(DUMSTRING, "/")
         CALL STRCAT(DUMSTRING, PlfIntName)
         PlfIntName = DUMSTRING
         OPEN(PlfIntFile,FILE=PlfIntName,STATUS='REPLACE')
         CLOSE(PlfIntFile)
      ENDDO
 900  CONTINUE
      CLOSE(IntervalFile)

C Open file where results of Planar fit method are written
      OPEN(PlfFile, FILE=PlfName, status='Replace',
     &     FORM='FORMATTED')
      WRITE(PlfFile,*) 'Doy  Hourmin  Sec Doy Hourmin  Sec  ',
     &                 '   Alpha   Beta    Gamma WBias   Swing'

C Open interval file again
      OPEN(IntervalFile,FILE = InterName,IOSTAT=FOO,STATUS='OLD')
      IF (FOO.NE.0) THEN
        WRITE(*,*) 'Could not open interval-file ',InterName
        STOP'- Fatal: no time info could be read -'
      ENDIF
C Read one line of comments (blindly, it is supposed to be a comment line)
      READ(IntervalFile,*)
      ENDINTERVAL = .FALSE.
      DO WHILE (.NOT. ENDINTERVAL)
         READ(IntervalFile,*,END=1000) (StartTime(i), i=1,3),
     &        (StopTime(i),i=1,3), Psychro, p, Fname, OutName,
     &        PlfIntName
         DUMSTRING =  DatDir
         CALL STRCAT(DUMSTRING, "/")
         CALL STRCAT(DUMSTRING, Fname)
         Fname = DUMSTRING

         DUMSTRING =  OutDir
         CALL STRCAT(DUMSTRING, "/")
         CALL STRCAT(DUMSTRING, OutName)
         OutName = DUMSTRING
         OPEN(OutFile, FILE=OutName, Status='UNKNOWN')

         DUMSTRING =  OutDir
         CALL STRCAT(DUMSTRING, "/")
         CALL STRCAT(DUMSTRING, PlfIntName)
         PlfIntName = DUMSTRING
C
C Show which file is currently being analysed
C
         WRITE(*,951) FName(:STRLEN(Fname)),(NINT(StartTime(i)),i=1,3)
 951     FORMAT(a,1X,3(I5,1X))
C
C Read raw data from file
C
         CALL ECReadNCDF(FName,StartTime,StopTime,
     &     Delay,Gain,Offset,RawSampl,NNNMax,Channels,MMMax,M,
     &     NCVarName, HAVE_UNCAL)
C 
C Do calibration (per sample)
C
         DO i=1,M
            CALL Calibrat(RawSampl(1,i),Channels,P,CorMean,
     &	      CalSonic,CalTherm,CalHyg,BadTc,Sample(1,i),N,Flag(1,i),
     &        Have_Uncal)
         ENDDO

C Determine mean velocity and write them to an intermediate file
         IF (M. GT. 0) THEN
           DO k=1,3
	     UMean(k,1) = 0.D0
             Mok(k) = 0
	   ENDDO
           DO j=1,M
	     DO k=1,3
               IF (.NOT. Flag(k,j)) THEN
	         UMean(k,1) = UMean(k,1) + Sample(k,j)
                 Mok(k) = Mok(k) + 1
               ENDIF
	     ENDDO
	   ENDDO
	   DO k=1,3
	     UMean(k,1) = UMean(k,1)/Mok(k)
	   ENDDO
           OPEN(PlfIntFile,FILE=PlfIntName,ACCESS='append',
     &          FORM='formatted')
           WRITE(PlfIntFile,100) (StartTime(i),i=1,3), 
     &           (StopTime(i),i=1,3),(Umean(k,1), k=1,3)
 100  FORMAT(2(F8.0,1X,F8.0,1X,F8.2),1X,3(F20.15))
           CLOSE(PlfIntFile)
	 ENDIF
      ENDDO
 1000 CONTINUE
      CLOSE(IntervalFile)

C
C Now start over again, and get the mean velocities from the intermediate files
C
      OPEN(IntervalFile,FILE = InterName,IOSTAT=FOO,STATUS='OLD')
      IF (FOO.NE.0) THEN
        WRITE(*,*) 'Could not open interval-file ',InterName
        STOP'- Fatal: no time info could be read -'
      ENDIF
C Read one line of comments (blindly, it is supposed to be a comment line)
      READ(IntervalFile,*)
      ENDINTERVAL = .FALSE.
      DO WHILE (.NOT. ENDINTERVAL)
         READ(IntervalFile,*,END=1100) (StartTime(i), i=1,3),
     &        (StopTime(i),i=1,3), Psychro, p, Fname, OutName,
     &        PlfIntName
         DUMSTRING =  DatDir
         CALL STRCAT(DUMSTRING, "/")
         CALL STRCAT(DUMSTRING, Fname)
         Fname = DUMSTRING

         DUMSTRING =  OutDir
         CALL STRCAT(DUMSTRING, "/")
         CALL STRCAT(DUMSTRING, OutName)
         OutName = DUMSTRING

         DUMSTRING =  OutDir
         CALL STRCAT(DUMSTRING, "/")
         CALL STRCAT(DUMSTRING, PlfIntName)
         PlfIntName = DUMSTRING
C
         OPEN(PlfIntFile,FILE=PlfIntName,STATUS='OLD')
         i=0
         DO WHILE (.TRUE.)
            IF (i+1 .GT. MaxRuns) THEN
               WRITE(*,*) 'Maximum number of runs exceeded: ',
     &                    MAXRUNS
               STOP
            ENDIF
            READ(PlfIntFile,*,END=1200) 
     &                  (RunStart(k,i+1), k=1,3),
     &                  (RunStop(k,i+1), k=1,3),
     &                  (Umean(k,i+1), k=1,3)
            i=i+1
         ENDDO
 1200    CONTINUE
         CLOSE(PlfIntFile)
         Nruns = i
         CALL PlanarFit(uMean,Nruns,Apf,Alpha,Beta,Gamma,WBias)

         IF (DoPrint) THEN
           WRITE(OutFile,*)
           WRITE(OutFile, 10) Alpha
 10        FORMAT('Alpha = ',F7.2,' degrees')
           WRITE(OutFile,20) Beta
 20        FORMAT('Beta  = ',F7.2,' degrees')
           WRITE(OutFile,30) Gamma
 30        FORMAT('Gamma = ',F7.2,' degrees')
           WRITE(OutFile, 40) WBias
 40        FORMAT('WBias = ',F7.4,' m/s')
           WRITE(OutFile, *) 'Matrix = '
           WRITE(OutFile,50) Apf(1,1),Apf(1,2),Apf(1,3)
           WRITE(OutFile,50) Apf(2,1),Apf(2,2),Apf(2,3)
           WRITE(OutFile,50) Apf(3,1),Apf(3,2),Apf(3,3)
           WRITE(OutFile,*)
 50        FORMAT(3(F7.3,1X))
         ENDIF
CC
CC Repeat the whole planar fit calculation as a check if it worked
CC
CC
CC Apply the planar fit matrix plus WBias-correction to all mean vectors
CC
C      WRITE(*,*)
C      WRITE(*,*) 'After planar fit correction : '
C      WRITE(*,*)
C      DO i=1,NRuns
C        DO k=1,3
C	  Dum(k) = UMean(i,k)
C	ENDDO
C	Dum(3) = Dum(3) - WBias
C	CALL ECMapVec(Apf,Dum,Dum)
C        DO k=1,3
C	  UMean(i,k) = Dum(k)
C	ENDDO
C      ENDDO
C      CALL PlanarFit(uMean,NRuns,Apf,Alpha,Beta,Gamma,WBias)
C
C      WRITE(*,5)
C      WRITE(*,10) Alpha
C      WRITE(*,20) Beta
C      WRITE(*,30) Gamma
C      WRITE(*,40) WBias
C      WRITE(*,*) 'Matrix = '
C      WRITE(*,50) Apf(1,1),Apf(1,2),Apf(1,3)
C      WRITE(*,50) Apf(2,1),Apf(2,2),Apf(2,3)
C      WRITE(*,50) Apf(3,1),Apf(3,2),Apf(3,3)
C
C Now perform individual tilt-corrections per run
C

C
C Read raw data from file
C
         CALL ECReadNCDF(FName,StartTime,StopTime,
     &     Delay,Gain,Offset,RawSampl,NNNMax,Channels,MMMax,M,
     &     NCVarName, HAVE_UNCAL)
C 
C Do calibration (per sample)
C
         DO i=1,M
            CALL Calibrat(RawSampl(1,i),Channels,P,CorMean,
     &	      CalSonic,CalTherm,CalHyg,BadTc,Sample(1,i),N,Flag(1,i),
     &        Have_Uncal)
         ENDDO
         Runlength = M
         DO k=1,3
            USum(k) = 0.D0
            Mok(k) = 0
            DO l=1,3
               UUSum(k,l) = 0.D0
               Cok(k,l) = 0
            ENDDO
         ENDDO
         DO j=1,RunLength
            DO k=1,3
C Since U=1, V=2, and W=3 !!
              Dum(k) = Sample(k,j)
            ENDDO
            Dum(3) = Dum(3) - WBias

	    CALL ECMapVec(Apf,Dum,Dum)
	    DO k=1,3
               IF (.NOT. Flag(k,j)) THEN
                  Mok(k) = Mok(k) + 1
	          USum(k) = USum(k) + Dum(k)
               ENDIF
	       DO l=1,3
                  IF ((.NOT. Flag(k,j)) .AND.
     &                (.NOT. Flag(l,j))) THEN
	             UUSum(k,l) = UUSum(k,l) + Dum(k)*Dum(l)
                     Cok(k,l) = Cok(k,l) + 1
                  ENDIF
	       ENDDO
	    ENDDO
	 ENDDO
	 DO k=1,3
c	    USum(k) = USum(k)/Mok(k)
	    USum(k) = USum(k)
	    DO l=1,3
c	      UUSum(k,l) = UUSum(k,l)/Cok(k,l)
 	      UUSum(k,l) = UUSum(k,l)
	    ENDDO
	 ENDDO


         CALL PFit(uSum,UUSum,RApf,RAlpha,RBeta,RGamma,RWBias)


         UHor = (RAlpha**2+RBeta**2)**0.5
         CosSwing = RBeta/UHor
         SinSwing = RAlpha/UHor
         Swing = 180.D0*ACOS(CosSwing)/PI
         IF (SinSwing.LT.0.D0) Swing = 360.D0-Swing
         IF (Swing.GT.180.D0) Swing = Swing-360.D0
       
         WRITE(OutFile,*)
         WRITE(OutFile,*) 'This run: ', (StartTime(i),i=1,3),
     &                    (StopTime(i),i=1,3)
         WRITE(OutFile,*) 'For this run the remaining',
     *                    ' tilts are:'
         WRITE(OutFile,*) 'Alpha   Beta    Gamma WBias Swing'
         WRITE(OutFile, 60)  RAlpha,RBeta,RGamma,RWBias,Swing
         WRITE(PlfFile, 61) (StartTime(i),i=1,3), (StopTime(i),i=1,3),
     &                       RAlpha,RBeta,RGamma,RWBias,Swing
 60      FORMAT(3(F7.2,1X),F7.4,1X,F7.1)
 61      FORMAT(6(F5.0,1X),1X,3(F7.2,1X),F7.4,1X,F7.1)
      ENDDO
 1100 CONTINUE
      CLOSE(IntervalFile)
      CLOSE(PlfFile)
      CLOSE(OutFile)

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
C ########################################################################
C
      SUBROUTINE Calibrat(RawSampl,Channels,P,CorMean,
     &  CalSonic,CalTherm,CalHyg,BadTc,Sample,N,Error, Have_Uncal)

      IMPLICIT NONE
      
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'
      INCLUDE 'calcomm.inc'

      INTEGER N,i,ColU,ColV,ColW,ColTSonic,ColTCple,ColHum,
     &  ColDay,Channels,ColHrMin,ColScnds,ColDiagnostic, ColTref
      LOGICAL Error(N),BadTc, Have_Uncal(N)
      REAL*8 RawSampl(Channels),Sample(N),P,UDum,VDum,Hook,Dum,
     &  CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ),CorMean(N),
     &  Hours,Minutes,Days,Secnds, TsCorr, WDum
      REAL*8 ECQ,ECBaseF ! External function calls
      
      REAL*8 C(0:NMaxOrder)
      INTEGER IFLG, DUMORDER, DUMTYP
      REAL*8 ABSERROR, RELERROR, ESTIM, DNLIMIT, UPLIMIT

      REAL*8 DUMFUNC, ECRawSchot
      EXTERNAL DUMFUNC, ECRawSchot
C
C To communicate with function of which zero is searched
C
      COMMON /DUMCOM/ C, DUM, DUMORDER, DUMTYP

C
C Which quantity can be found in which column of the raw data?
C
      PARAMETER(ColDay        = Doy)
      PARAMETER(ColHrMin      = HourMin)
      PARAMETER(ColScnds      = Sec)
      PARAMETER(ColU          = U)
      PARAMETER(ColV          = V)
      PARAMETER(ColW          = W)
      PARAMETER(ColTSonic     = TSonic)
      PARAMETER(ColDiagnostic = Diagnost)
      PARAMETER(ColTCple      = Tcouple)
      PARAMETER(ColHum        = Humidity)
      PARAMETER(ColTref        = Tref)
C
C In principle no apparatus has shown to give a valid signal.
C
      DO i=1,N
        Error(i) = (.TRUE.)
      ENDDO
C
C Calculate time from the columns in the raw data
C
      Days = RawSampl(ColDay) - DBLE(FirstDay)
      Hours = INT(0.01D0*RawSampl(ColHrMin))
      Minutes = RawSampl(ColHrMin) - 100.D0*Hours
      Secnds = RawSampl(ColScnds)
      Sample(TTime) =
     &  Secnds + 60.D0*(Minutes +60.D0*(Hours + 24.D0*Days))

      Error(TTime) = (.FALSE.)
C
C Take the velocity from the raw data
C
      IF (HAVE_SONCAL) THEN
         UDum      =  RawSampl(ColU) ! U [m/s]
         VDum      =  RawSampl(ColV) ! V [m/s]
         Sample(W) =  RawSampl(ColW) ! W [m/s]
         
         Error(U) = (ABS(UDum).GT.40.D0)
         Error(V) = (ABS(VDum).GT.40.D0)
         Error(W) = (ABS(Sample(W)).GT.20.D0)
C
C If this is a wind tunnel calibrated sonic,
C we can apply the KNMI calibrations
C
         IF (CalSonic(1) .EQ. ApSon3Dcal) THEN
            WDum = Sample(W)
            CALL ECSCal(CalSonic,
     &                  UDUM, VDUM, WDum,
     &                  ERROR(U), ERROR(V), ERROR(W))
            Sample(W) = WDum
         ENDIF
   
C This is the construction when we have a diagnostic variable
C      IF (((Error(U).OR.Error(V)).OR.Error(W)).OR.
C     &  ((RawSampl(ColDiagnostic).LT.0.D0).OR.
C     &   (RawSampl(ColDiagnostic).GT.100.D0))) THEN
         IF ((Error(U).OR.Error(V)).OR.Error(W)) THEN
           Error(U) = (.TRUE.)
           Error(V) = (.TRUE.)
           Error(W) = (.TRUE.)
         ENDIF
   
         IF (.NOT.((Error(U).OR.Error(V)).OR.Error(W))) THEN

C
C Rotate velocity according to angle of setup's north relative real north
C
           Hook = PI*CalSonic(QQYaw)/180.D0
           Sample(U) =  COS(Hook)*UDum + SIN(Hook)*VDum
           Sample(V) = -SIN(Hook)*UDum + COS(Hook)*VDum

           Sample(TSonic) = RawSampl(ColTSonic) + Kelvin
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
           Sample(TSonic) = Sample(TSonic)
           Error(Tsonic) = ((Sample(TSonic).GT.(Kelvin+MaxT))
     &       .OR. (Sample(TSonic).LT.(Kelvin+MinT)))
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
      IF (HAVE_TCCAL) THEN
         IF (HAVE_TREF) THEN

C
C This is calibration according to Campbells P14 instruction
C      T =  Calibration(voltage +
C                       inverse calibration(reference temperature))
C Reference temperature is supposed to be in CELCIUS !!
C
            DO i=0,CalTherm(QQOrder)
              c(i) = CalTherm(QQC0+i)
            ENDDO
            Dum = RawSampl(ColTref)
            Error(TCouple) = 
     &           (((Kelvin + DUM) .GT.(Kelvin+MaxT))
     &             .OR. ((Kelvin + DUM) .LT.(Kelvin+MinT)))
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
     &             ECBaseF(DNLIMIT + RawSampl(Tcouple),
     &              NINT(CalTherm(QQFunc)),
     &              NINT(CalTherm(QQOrder)),c)
            ENDIF
         ELSE
C
C Suppose that sample is already temperature
C
             Sample(TCouple) = RawSampl(ColTCple)+Kelvin
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
      IF (HAVE_HYGCAL) THEN
         Dum = RawSampl(ColHum)
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

            Sample(Humidity) = CorMean(Humidity) +
     &         ECBaseF(Dum,NINT(CalHyg(QQFunc)),
     &           NINT(CalHyg(QQOrder)),c)/1000.D0

            Error(Humidity) = ((Sample(Humidity).GT.MaxRhoV)
     &          .OR. (Sample(Humidity).LT.MinRhoV))
            Error(SpecHum) = Error(Humidity)
C
C Calculate specific humidity associated with this absolute humidity
C
            IF (.NOT.Error(Humidity)) THEN
              IF (.NOT.BadTc) THEN
                IF (.NOT.Error(TCouple)) THEN
                  Sample(SpecHum)=ECQ(Sample(Humidity),
     &                                Sample(TCouple),P)
                ELSE IF (.NOT.Error(TSonic)) THEN
                  TsCorr = ECRawSchot(Sample(TSonic),
     &                                Sample(Humidity), P)
                  Sample(SpecHum)=ECQ(Sample(Humidity),
     &                                TsCorr,P)
                ELSE
                  Error(SpecHum) = (.TRUE.)
                ENDIF
              ELSE
                IF (.NOT.Error(TSonic)) THEN
                  TsCorr = ECRawSchot(Sample(TSonic),
     &                                Sample(Humidity), P)
                  Sample(SpecHum)=ECQ(Sample(Humidity),
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
C Time is used for detrending. Therefore we take out all samples with
C a defect in one of the velocities.
C I wouldnt know why so test has been removed AM
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

      REAL*8 ECBaseF
      EXTERNAL ECBaseF

C
C To communicate with function of which zero is searched
C
      COMMON /DUMCOM/ C, DUM, DUMORDER, DUMTYP

      DUMFUNC = ECBaseF(X, DUMTYP, DUMORDER, C) - DUM
      END
      



