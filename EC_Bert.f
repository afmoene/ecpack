C Contains changes due to Arnold Moene:
C - Use sonic temp. rather than thermocouple for calculation of
C   specific humidity
C - Check for existence of time file
      PROGRAM EC

      INCLUDE 'physcnst.for'

      CHARACTER*255 FNAME
      CHARACTER*40 OutName,FluxName,DumName1,DumName2
      LOGICAL PRaw,PPitch,PYaw,PRoll,PFreq,PO2,PWebb,PCal,PIndep,
     &  DoPitch,DoYaw,DoRoll,DoFreq,DoO2,DoWebb,DoCrMean,PSonic,
     &  DoSonic,DoTilt,PTilt,DoPrint,Flag(NNMAx,MMMax),DoDetren,
     &  PDetrend,DoStruct
      INTEGER N,i,M,MIndep(NNMax),CIndep(NNMax,NNMax),
     &  Channels,Delay(NNNMax),FirstDay,Mok(NNMax),Cok(NNMax,NNMax)
      REAL*8 RawSampl(NNNMax,MMMax),Sample(NNMax,MMMax),P,Psychro,
     &  Mean(NNMax),TolMean(NNMax),Cov(NNMax,NNMax),
     &  TolCov(NNMax,NNMax),LLimit,ULimit,FrCor(NNMax,NNMax),
     &  YawLim,RollLim,DirPitch,RC(NNMax),O2Factor(NNMax),
     &  Freq,DirYaw,DirRoll,
     &  SonFactr(NNMax),PrePitch,PreYaw,PreRoll,
     &  HSonic,dHSonic,HTc,dHTc,
     &  LvE,dLvE,LvEWebb,dLvEWebb,
     &  UStar,dUStar,Tau,dTau,CorMean(NNMax),
     &  Gain(NNNMax),Offset(NNNMax),StartTime(3),StopTime(3)
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ),
     &  R,dR,CTSon2,CTCop2,Cq2,CTSonq,CTCopq,
     &  dCTSon2,dCTCop2,dCq2,dCTSonq,dCTCopq
      CHARACTER*6 QName(NNMax)
      CHARACTER*9 UName(NNMax)
      CHARACTER*12 SonName,CoupName,HygName,ParmName,InterName

      INTEGER FOO
C
C The name of the calibration routine specified at the end of this file
C
      EXTERNAL Calibrat
C
C Common block to pass specific data to calibration routine below
C
      COMMON/CalOffset/ FirstDay
C
C File to which the final fluxes are to be written
C
      FluxName = 'c:/Rapid Analysis/EC_bert/EC_Bert.flx'
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
     &  'Mean(U)        dMean(U)          ',
     &  'Mean(TSon)     dMean(TSon)       ',
     &  'Mean(TCop)     dMean(TCop)       ',
     &  'Mean(q)        dMean(q)          ',
     &  'sigma(TSon)    dsigma(TSon)      ',
     &  'sigma(TCop)    dsigma(TCop)      ',
     &  'sigma(q)       dsigma(q)         ',
     &  'sigma(u)       dsigma(u)         ',
     &  'sigma(v)       dsigma(v)         ',
     &  'sigma(w)       dsigma(w)         ',
     &  'cov(TSon,q)    dcov(TSon,q)      ',
     &  'cov(TCop,q)    dcov(TCop,q)      ',
     &  'cov(TSon,U)    dcov(TSon,U)      ',
     &  'cov(TCop,U)    dcov(TCop,U)      ',
     &  'cov(q,U)       dcov(q,U)         ',
     &  'H(Sonic)       Tol(HSonic)       ',
     &  'H(TCouple)     Tol(HTCouple)     ',
     &  'LvECov         dLvECov           ',
     &  'LvEWebb        dLvEWebb          ',
     &  'LvE            dLvE              ',
     &  'UStar          dUStar            ',
     &  'Tau            dTau              ',
     &  'R(delay)       dR                ',
     &  'CTSon2         dCTSon2           ',
     &  'CTCop2         dCTCop2           ',
     &  'Cq2            dCq2              ',
     &  'CTSonq         dCTSonq           ',
     &  'CTCopq         dCTCopq           '
     &  )
C
C Number of quantities involved in this experiment
C
      N = 8  ! Number of calibrated quantities you will follow
C
C Read specifications about setup from file
C

       ParmName = 'specbert.dat'

      CALL ECParams(ParmName,
     &  Freq,YawLim,RollLim,PrePitch,PreYaw,PreRoll,
     &  LLimit,ULimit,DoCrMean,DoDetren,DoSonic,
     &  DoTilt,DoPitch,DoYaw,DoRoll,DoFreq,DoO2,DoWebb,DoStruct,DoPrint,
     &  PRaw,PCal,PDetrend,PIndep,PTilt,PPitch,PYaw,
     &  PRoll,PSonic,PO2,PFreq,PWebb,SonName,CoupName,HygName,
     &  InterName)
C
C Read calibration data from files :
C
      DO i=1,NNNMax
        Delay( i) = 0
        Gain(  i) = 1.D0
        Offset(i) = 0.D0
      ENDDO

      CALL ECReadAp(SonName,CalSonic) ! The Solent sonic
      CALL ECReadAp(CoupName,CalTherm) ! A thermo-couple
      CALL ECReadAp(HygName,CalHyg) ! A Krypton hygrometer
C
C Set delays, gains and offset of all channels in raw datafile
C At reading time all channels are corrected by mapping:
C                x --> (x/Gain) + Offset
C
      Delay( 4) = 2 ! For lower sonic : Always two samples.
      Delay( 5) = 2 ! Index of "delay" refers to columns in raw datafile
      Delay( 6) = 2
      Delay( 7) = 2
      Delay( 8) = 2
      Delay( 9) = NINT(CalTherm(QQDelay)*0.001D0*Freq)
      Delay(10) = NINT(CalHyg(QQDelay)*0.001D0*Freq)

      Gain( 4) = CalSonic(QQGain)
      Gain( 5) = CalSonic(QQGain)
      Gain( 6) = CalSonic(QQGain)
      Gain( 9) = CalTherm(QQGain)
      Gain(10) = CalHyg(QQGain)

      Offset( 4) = CalSonic(QQOffset)
      Offset( 5) = CalSonic(QQOffset)
      Offset( 6) = CalSonic(QQOffset)
      Offset( 9) = CalTherm(QQOffset)
      Offset(10) = CalHyg(QQOffset)
C
C Get names of files
C
      OPEN(IntervalFile, IOSTAT = FOO, FILE = InterName,
     +              STATUS = 'old')
      IF (FOO .NE. 0) THEN
         WRITE(*,*) 'Could not open file ', Intername
         STOP '-- Fatal: no time info could be read --'
      ENDIF

c      OPEN(IntervalFile,FILE = InterName)
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
C Read time-interval from file, plus value of psychrometer RhoV and atm. pressure
C
 36   READ(IntervalFile,*,END=37) (StartTime(i),i=1,3),
     &  (StopTime(i),i=1,3),Psychro,p,DumName1,DumName2
      FName   = DumName1
      OutName = DumName2
C
C Here you can insert the name of the directory of the raw data
C and the name of the directory for intermediate results
C Do mind! Use forward-slashes in the directory names!!!!!!!
C
C The following lines are examples:
C
C      FName = 'e:/suikerbieten/ec_split/'//DumName1
C      OutName = 'd:/wouter/'//DumName2

        FName = 'c:/Rapid Data/EC_bert/'//DumName1
        OutName = 'c:/Rapid Analysis/EC_bert/'//DumName2


C
C Show which file is currently being analysed
C
      WRITE(*,951) DumName1,(NINT(StartTime(i)),i=1,2)
 951  FORMAT(a,1X,2(I5,1X))
C
C If wanted, make a file for printing intermediate results
C
      IF (DoPrint) THEN
        OPEN(OutFile,FILE=OutName,FORM='FORMATTED')
        WRITE(OutFile,54) FName
 54     FORMAT('File = ',a)
        WRITE(OutFile,*) 'From (day,hour,minute)',
     &    (StartTime(i),' ',i=1,3),' to ',(StopTime(i),i=1,3)
        WRITE(OutFile,*) 'At RhoV(slow) = ',Psychro,' [kg m^{-3}]',
     &    ' and atm. pressure = ',p,' [Pascal]'
      ENDIF
C
C For subroutine calibrat below find the first day; transfer via common-block
C
      FirstDay = StartTime(1)
C
C Read raw data from file
C
      CALL ECReadNCDF(FName,StartTime,StopTime,
     &  Delay,Gain,Offset,RawSampl,NNNMax,Channels,MMMax,M)

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
        CALL ECMain(OutFile,DoPrint,
     &    RawSampl,NNNMax,Channels,NNMax,N,MMMax,M,
     &    DoCrMean,PCal,PIndep,
     &    Psychro,Freq,CalSonic,CalTherm,CalHyg,P,
     &    Calibrat,
     &    Sample,Flag,Mok,Cok,MIndep,CIndep,Rc,
     &    QName,UName,
     &    DoDetren,PDetrend,
     &    DoTilt,PTilt,PrePitch,PreYaw,PreRoll,
     &    DoPitch,PPitch,DirPitch,
     &    DoYaw,PYaw,YawLim,DirYaw,
     &    DoRoll,PRoll,RollLim,DirRoll,
     &    DoSonic,PSonic,SonFactr,
     &    DoO2,PO2,O2Factor,
     &    DoFreq,PFreq,LLimit,ULimit,FrCor,
     &    DoWebb,PWebb,
     &    Mean,TolMean,Cov,TolCov,
     &    HSonic,dHSonic,HTc,dHTc,
     &    LvE,dLvE,LvEWebb,dLvEWebb,
     &    UStar,dUStar,Tau,dTau)
C
C Calculate structure parameters
C
        IF (DoStruct) THEN
          R = 1.00D0
          CALL ECStruct(Sample,NNMax,N,MMMax,M,Flag,TSonic ,TSonic ,
     &      R,dR,Freq,CIndep,CTSon2,dCTSon2)
          CALL ECStruct(Sample,NNMax,N,MMMax,M,Flag,TCouple,TCouple,
     &      R,dR,Freq,CIndep,CTCop2,dCTCop2)
          CALL ECStruct(Sample,NNMax,N,MMMax,M,Flag,SpecHum,SpecHum,
     &      R,dR,Freq,CIndep,Cq2   ,dCq2   )
          CALL ECStruct(Sample,NNMax,N,MMMax,M,Flag,TSonic ,SpecHum,
     &      R,dR,Freq,CIndep,CTSonq,dCTSonq)
          CALL ECStruct(Sample,NNMax,N,MMMax,M,Flag,TCouple,SpecHum,
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
        WRITE(FluxFile,55)
     &    (NINT(StartTime(i)),i=1,3),
     &    (NINT(StopTime(i)),i=1,3),
     &    M,(Mok(i),i=1,8),
     &    NINT(DirPitch),
     &    NINT(2.D0*180.D0*ATAN(SQRT(Cov(V,V))/Mean(U))/Pi),
     &    Mean(U),TolMean(U),
     &    Mean(TSonic),TolMean(TSonic),
     &    Mean(TCouple),TolMean(TCouple),
     &    Mean(SpecHum),TolMean(SpecHum),
     &    SQRT(Cov(TSonic,TSonic)),
     &    0.5D0*TolCov(TSonic,TSonic)/SQRT(Cov(TSonic,TSonic)),
     &    SQRT(Cov(TCouple,TCouple)),
     &    0.5D0*TolCov(TCouple,TCouple)/SQRT(Cov(TCouple,TCouple)),
     &    SQRT(Cov(SpecHum,SpecHum)),
     &    0.5D0*TolCov(SpecHum,SpecHum)/SQRT(Cov(SpecHum,SpecHum)),
     &    SQRT(Cov(U,U)),
     &    0.5D0*TolCov(U,U)/SQRT(Cov(U,U)),
     &    SQRT(Cov(V,V)),
     &    0.5D0*TolCov(V,V)/SQRT(Cov(V,V)),
     &    SQRT(Cov(W,W)),
     &    0.5D0*TolCov(W,W)/SQRT(Cov(W,W)),
     &    Cov(TSonic,SpecHum),TolCov(TSonic,SpecHum),
     &    Cov(TCouple,SpecHum),TolCov(TCouple,SpecHum),
     &    Cov(TSonic,U),TolCov(TSonic,U),
     &    Cov(TCouple,U),TolCov(TCouple,U),
     &    Cov(SpecHum,U),TolCov(SpecHum,U),
     &    HSonic,dHSonic,HTc,dHTc,
     &    LvE,dLvE,LvEWebb,dLvEWebb,
     &    LvE+LvEWebb,
     &    SQRT(dLvE**2.D0 + dLvEWebb**2.D0),
     &    UStar,dUStar,Tau,dTau,
     &    R,dR,
     &    CTSon2,dCTSon2,CTCop2,dCTCop2,Cq2,dCq2,CTSonq,dCTSonq,
     &    CTCopq,dCTCopq
 55     FORMAT(2(I3,1X,2(I2,1X)),9(I6,1X),2(I5,1X),62(2(G14.5:,1X),3X))
      ELSE
        WRITE(FluxFile,55)
     &    (NINT(StartTime(i)),i=1,3),
     &    (NINT(StopTime(i)),i=1,3),
     &    M,(0,i=1,8),-9999,-9999,(-9999.0,i=1,56)
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
C OutPut: Sample   : REAL*8(N) : Calibrated sample
C         Error    : LOGICAL(N) : Indicator if quantities are valid
C                    after calibration
C
C ########################################################################
C
      SUBROUTINE Calibrat(RawSampl,Channels,P,CorMean,
     &  CalSonic,CalTherm,CalHyg,Sample,N,Error)

      INCLUDE 'physcnst.for'

      INTEGER N,i,ColU,ColV,ColW,ColTSonic,ColTCple,ColHum,
     &  ColDay,Channels,ColHrMin,ColScnds,FirstDay,ColDiagnostic
      LOGICAL Error(N)
      REAL*8 RawSampl(Channels),Sample(N),P,UDum,VDum,Hook,Dum,
     &  CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ),CorMean(N),C(0:5),
     &  Hours,Minutes,Days,Secnds
      REAL*8 ECQ,ECBaseF ! External function calls
      COMMON/CalOffset/ FirstDay
C
C Which quantity can be found in which column of the raw data?
C
      PARAMETER(ColDay        =  1)
      PARAMETER(ColHrMin      =  2)
      PARAMETER(ColScnds      =  3)
      PARAMETER(ColU          =  4)
      PARAMETER(ColV          =  5)
      PARAMETER(ColW          =  6)
      PARAMETER(ColTSonic     =  7)
      PARAMETER(ColDiagnostic =  8)
      PARAMETER(ColTCple      =  9)
      PARAMETER(ColHum        = 10)
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
      UDum      =  RawSampl(ColU) ! U [m/s]
      VDum      =  RawSampl(ColV) ! V [m/s]
      Sample(W) =  RawSampl(ColW) ! W [m/s]

      Error(U) = (ABS(UDum).GT.40.D0)
      Error(V) = (ABS(VDum).GT.40.D0)
      Error(W) = (ABS(Sample(W)).GT.20.D0)
      IF (((Error(U).OR.Error(V)).OR.Error(W)).OR.
     &  ((RawSampl(ColDiagnostic).LT.0.D0).OR.
     &   (RawSampl(ColDiagnostic).GT.100.D0))) THEN
        Error(U) = (.TRUE.)
        Error(V) = (.TRUE.)
        Error(W) = (.TRUE.)
      ENDIF

      IF (.NOT.((Error(U).OR.Error(V)).OR.Error(W))) THEN
C
C Rotate velocity according to angle of setup's north relative real north
C
        Hook = PI*CalSonic(QQPitch)/180.D0
        Sample(U) =  COS(Hook)*UDum + SIN(Hook)*VDum
        Sample(V) = -SIN(Hook)*UDum + COS(Hook)*VDum

        Sample(TSonic) = RawSampl(ColTSonic) + Kelvin
C
C Here the Schotanus et al. correction for sensitivity of the sonic for
C lateral velocity is applied directly to the raw data. This is only
C half of the Schotanus-correction. Humidity-correction comes later.
C
        Sample(TSonic) = Sample(TSonic)
     &    + ((Sample(U))**2 + (Sample(V))**2)/GammaR
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
     &    .OR. (Sample(TSonic).LT.(Kelvin+MinT)))
C
C Calibrate thermocouple
C
        Sample(TCouple) = RawSampl(ColTCple)+Kelvin
        Error(TCouple) = ((Sample(TCouple).GT.(Kelvin+MaxT))
     &    .OR. (Sample(TCouple).LT.(Kelvin+MinT)))
C
C Calibrate hygrometer
C
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
     &      ECBaseF(Dum,NINT(CalHyg(QQFunc)),
     &        NINT(CalHyg(QQOrder)),c)/1000.D0

          Error(Humidity) = ((Sample(Humidity).GT.MaxRhoV)
     &        .OR. (Sample(Humidity).LT.MinRhoV))
C
C Calculate specific humidity associated with this absolute humidity
C
C          Sample(SpecHum) = ECQ(Sample(Humidity),Sample(TCouple),P)
          Sample(SpecHum) = ECQ(Sample(Humidity),Sample(TSonic),P)
C          Error(SpecHum) = Error(Humidity)
          Error(SpecHum) = (Error(Humidity) .OR. Error(Tsonic))

        ENDIF
      ENDIF
C
C Time is used for detrending. Therefore we take out all samples with
C a defect in one of the velocities.
C
      Error(TTime) = (.FALSE.)
      DO i=U,W
        Error(TTime) = (Error(TTime).OR.Error(i))
      ENDDO

      RETURN
      END





      INCLUDE 'ecpack.f'


