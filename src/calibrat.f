
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
C         CalCO2   : Calibration array REAL*8(NQQ) of the CO2-sensor
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
      INTEGER IFLG, DUMORDER, DUMTYP, DIAG_WORD
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
      Days = RawSampl(QUDoy) - DBLE(FirstDay)
      Hours = NINT(0.01D0*RawSampl(QUHourMin))
      Minutes = RawSampl(QUHourMin) - 100.D0*Hours
      Secnds = RawSampl(QUSec)
      Sample(TTime) =
     &  Secnds + 60.D0*(Minutes +60.D0*(Hours + 24.D0*Days))

      Error(TTime) = (.FALSE.)
C
C Take the velocity from the raw data (again assume it's a 3D sonic)
C
      IF (HAVE_UNCAL(QUU) .AND. 
     &    HAVE_UNCAL(QUV) .AND. 
     &    HAVE_UNCAL(QUW)) THEN
         UDum  =  RawSampl(QUU) ! U [m/s]
         VDum  =  RawSampl(QUV) ! V [m/s]
         Sample(W) = RawSampl(QUW) ! W [m/s]

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
         IF (Have_Uncal(QUDiagnost)) THEN
            DIAG_WORD = INT(RawSampl(QUDiagnost)/4096)
            IF ((Error(U).OR.Error(V).OR.Error(W)).OR.
     &          (DIAG_WORD.NE.0.D0)) THEN
               Error(U) = (.TRUE.)
               Error(V) = (.TRUE.)
               Error(W) = (.TRUE.)
               Error(TSonic) = (.TRUE.)               
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

           Sample(TSonic) = RawSampl(QUTSonic) + Kelvin
C
C Here the Schotanus et al. correction for sensitivity of the sonic for
C lateral velocity is applied directly to the raw data. This is only
C half of the Schotanus-correction. Humidity-correction comes later.
C
           Sample(TSonic) = Sample(TSonic)
     &         + (Sample(U)**2 + Sample(V)**2)/GammaR
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
	 Error(U) = (Error(U) .OR. (.NOT. Have_Uncal(QUU)))
	 Error(V) = (Error(V) .OR. (.NOT. Have_Uncal(QUV)))
	 Error(W) = (Error(W) .OR. (.NOT. Have_Uncal(QUW)))
	 Error(TSonic) = (Error(TSonic) .OR. 
     &                    (.NOT. Have_Uncal(QUTSonic)))
      ELSE
         Error(U) = .TRUE.
         Error(V) = .TRUE.
         Error(W) = .TRUE.
         Error(Tsonic) = .TRUE.
      ENDIF
C
C Calibrate thermocouple
C
      IF (HAVE_UNCAL(QUTCouple)) THEN
         IF (HAVE_UNCAL(QUTref)) THEN

C
C This is calibration according to Campbells P14 instruction
C      T =  Calibration(voltage +
C                       inverse calibration(reference temperature))
C Reference temperature is supposed to be in Celcius !!
C
            DO i=0,NINT(CalTherm(QQOrder))
              c(i) = CalTherm(QQC0+i)
            ENDDO
            Dum = RawSampl(QUTref) 
            Error(TCouple) =
     &           ((DUM .GT.(MaxT))
     &             .OR. (DUM .LT.(MinT)))
            IF (.NOT. Error(TCouple)) THEN
               UPLIMIT = (DUM + 3.0D0 - CalTherm(QQC0))/CalTherm(QQC1)
               DNLIMIT = (DUM - 3.0D0 - CalTherm(QQC0))/CalTherm(QQC1)
               ESTIM = (Dum - CalTherm(QQC0))/CalTherm(QQC1)
               RELERROR = 1D-4
               ABSERROR = 0.001D0/CalTherm(QQC1)
               DUMORDER = NINT(CalTherm(QQOrder))
               DUMTYP = NINT(CalTherm(QQFunc))
               CALL DFZERO(DUMFUNC, DNLIMIT, UPLIMIT, ESTIM, RELERROR,
     &                    ABSERROR, IFLG)
               IF (IFLG .GT. 3) THEN
	                ERROR(TCouple) = .TRUE.
	            ENDIF
               Sample(Tcouple) = Kelvin +
     &             EC_M_BaseF(DNLIMIT + RawSampl(QUTcouple),
     &              NINT(CalTherm(QQFunc)),
     &              NINT(CalTherm(QQOrder)),c)
            ENDIF
         ELSE
C
C Suppose that sample is already temperature (now in Kelvin, Sept. 18, 2002!)
C
             Sample(TCouple) = RawSampl(QUTCouple) + Kelvin
         ENDIF

         Error(TCouple) = ((Sample(TCouple).GT.(MaxT))
     &     .OR. (Sample(TCouple).LT.(MinT)))
         Error(TCouple) = (Error(TCouple) .OR.
     &                        (.NOT. Have_Uncal(QUTCouple)))
      ELSE
         Error(Tcouple) = .TRUE.
      ENDIF
C
C Calibrate hygrometer
C
      IF (HAVE_UNCAL(QUHumidity)) THEN
         Dum = RawSampl(QUHumidity)
         IF (Dum .GE. Epsilon) THEN
C
C Add an optional correction to the krypton's humidity to compensate for
C drift in the hygrometer (e.g. dirt). Compensation may be provided by
C the difference between the mean humidity by a (slow) psychrometer's
C output and the mean of the humidity estimated by the hygrometer.
C
            DO i=0,NINT(CalHyg(QQOrder))
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
     &                        (.NOT. Have_Uncal(QUHumidity)))
         Error(SpecHum) = (Error(SpecHum) .OR. Error(Humidity))
      ELSE
         Error(Humidity) = .TRUE.
         Error(SpecHum) = .TRUE.
      ENDIF
C
C Get CO2 sample
C
      IF (HAVE_UNCAL(QUCO2)) THEN
         DO i=0,INT(CalCO2(QQOrder))
              c(i) = CalCO2(QQC0+i)
         ENDDO
	 Sample(CO2) = CorMean(CO2) + 
     &         EC_M_BaseF(RawSampl(QUCO2),NINT(CalCO2(QQFunc)),
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
     &                        (.NOT. Have_Uncal(QUCO2)))
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
