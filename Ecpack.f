C $Id$
C###########################################################################
C
C
C
C      EEEEE	   CCCC        PPPPP	       A	        CCCC   K    K
C      E	       C    C	      P    P      A A       C    C  K   K
C      E	      C		         P     P    A   A     C	      K  K
C      EEEE	   C    	  ------ PPPPPP    A	    A    C	      K K
C      E        C		         P	      AAAAAAAAA   C	      KK K
C      E	        C    C	      P	     A	      A   C    C  K   K
C      EEEEE	   CCCC 	      P      A	          A   CCCC   K    K
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
C  - Pitch-correction to make Mean(V) --> 0
C  - Yaw-correction to make Mean(W) --> 0
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
C - Flow distortion correction
C - Better standardisation of different types of functions
C
C
C###########################################################################

      SUBROUTINE ECMain(OutF,DoPrint,
     &	RawSampl,MaxChan,Channels,NMax,N,MMax,M,DoCrMean,PCal,PIndep,
     &	Psychro,Freq,CalSonic,CalTherm,CalHyg,P,
     &	Calibr,
     &	Sample,Flag,Mok,Cok,MIndep,CIndep,Rc,BadTc,
     &	QName,UName,
     &	DoDetren,PDetrend,
     &	DoTilt,PTilt,PrePitch,PreYaw,PreRoll,
     &	DoPitch,PPitch,DirPitch,
     &	DoYaw,PYaw,YawLim,DirYaw,
     &	DoRoll,PRoll,RollLim,DirRoll,
     &	DoSonic,PSonic,SonFactr,
     &	DoO2,PO2,O2Factor,
     &	DoFreq,PFreq,LLimit,ULimit,FrCor,
     &	DoWebb,PWebb,
     &	Mean,TolMean,Cov,TolCov,
     &	HSonic,dHSonic,HTc,dHTc,
     &	LvE,dLvE,LvEWebb,dLvEWebb,
     &	UStar,dUStar,Tau,dTau,
     &  MEANW, TOLMEANW, HAVE_UNCAL)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECMain
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Purpose : Integrated routine which:
C  - Calibrates raw samples
C  - Estimates mean values and covariances with respective tolerances
C  - Corrects the resulting mean values and covariances for all
C    effects selected by the user
C  - Estimates, from the final mean values and covariances, the surface-
C    fluxes with tolerances.
C
C Revision: 03-04-2001: get mean W and its tolerance before tilt
C                       correction; export via interface (AM)
C Revision: 28-05-2001: added passing of info on whether uncalibrated
C                       data are available (AM)

      INCLUDE 'physcnst.for'
      INCLUDE 'calcomm.inc'

      INTEGER NMax,N,i,M,j,MIndep(NMax),CIndep(NMax,NMax),
     &	OutF,MMax,MaxChan,Channels,Mok(NMax),Cok(NMax,NMax),NTcOut
      LOGICAL PPitch,PYaw,PRoll,PFreq,PO2,PWebb,PCal,PIndep,
     &	DoPitch,DoYaw,DoRoll,DoFreq,DoO2,DoWebb,DoCrMean,PSonic,
     &	DoSonic,DoTilt,PTilt,DoPrint,Flag(NMAx,MMax),DoDetren,
     &	PDetrend,AnyTilt,DumPitch,DumYaw,DumRoll,DumTilt,BadTc
      REAL*8 RawSampl(MaxChan,MMax)
      REAL*8 P,Psychro,Sample(NMax,MMax),
     &	Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),LLimit,ULimit,FrCor(NMax,NMax),
     &	YawLim,RollLim,DirPitch,RC(NMax),O2Factor(NMax),
     &	Freq,DirYaw,DirRoll,
     &	SonFactr(NMax),CorMean(NNMax),PrePitch,PreYaw,PreRoll,
     &	HSonic,dHSonic,HTc,dHTc,
     &	LvE,dLvE,LvEWebb,dLvEWebb,
     &	UStar,dUStar,Tau,dTau,Speed(3),DumCov(3,3)
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ)
      REAL*8 MEANW, TOLMEANW
      CHARACTER*6 QName(NMax)
      CHARACTER*9 UName(NMax)
      LOGICAL HAVE_CAL(NMax)
      LOGICAL HAVE_UNCAL(NMax)
C
C Check whether we have the needed calibration info
C
      DO I=1,NMax
         HAVE_CAL(i) = .FALSE.
      ENDDO
      IF (HAVE_SONCAL) THEN
         HAVE_CAL(U) = .TRUE.
         HAVE_CAL(V) = .TRUE.
         HAVE_CAL(W) = .TRUE.
         HAVE_CAL(TSonic) = .TRUE.
      ENDIF
      IF (HAVE_HYGCAL) THEN
         HAVE_CAL(Humidity) = .TRUE.
         HAVE_CAL(SpecHum) = .TRUE.
      ENDIF
      IF (HAVE_TCCAL) THEN
         HAVE_CAL(TCouple) = .TRUE.
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
C Names of quantities and units
C
      QName(U	    ) = '     U'
      QName(V	    ) = '     V'
      QName(W	    ) = '     W'
      QName(TSonic  ) = 'T(son)'
      QName(TCouple ) = 'T(cpl)'
      QName(Humidity) = '  RhoV'
      QName(SpecHum ) = '     q'
      QName(TTime   ) = '  Time'

      UName(U	    ) = '[m/s]'
      UName(V	    ) = '[m/s]'
      UName(W	    ) = '[m/s]'
      UName(TSonic  ) = '[K]'
      UName(TCouple ) = '[K]'
      UName(Humidity) = '[kg m^-3]'
      UName(SpecHum ) = '[kg/kg]'
      UName(TTime   ) = '[s]'
C
C
C Calibrate the raw samples for the first time ignoring drift in apparatus.
C
C
      DO i=1,N
       	CorMean(i) = 0.D0
      ENDDO

      NTcOut = 0
      BadTc = (.FALSE.)
      DO i=1,M
         CALL Calibr(RawSampl(1,i),Channels,P,CorMean,
     &	  CalSonic,CalTherm,CalHyg,BadTc,Sample(1,i),N,Flag(1,i),
     &    Have_Uncal)
        IF (Flag(TCouple,i)) NTcOut = NTcOut + 1
      ENDDO
      BadTc = (NTcOut.GE.(M/2))
      IF (BadTc) THEN
        DO i=1,M
          CALL Calibr(RawSampl(1,i),Channels,P,CorMean,
     &      CalSonic,CalTherm,CalHyg,BadTc,Sample(1,i),N,Flag(1,i),
     &      Have_Uncal)
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
	CALL ECAverag(Sample,NMax,N,MMax,M,Flag,
     &	  Mean,TolMean,Cov,TolCov,MIndep,CIndep,Mok,Cok)

	CorMean(Humidity) = Psychro - Mean(Humidity)
	IF (DoPrint.AND.PCal) THEN
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'Added to humidity : ',
     &	  (Psychro - Mean(Humidity)),' [kg m^{-3}]'
	  WRITE(OutF,*)
        ENDIF

	DO i=1,M
	  CALL Calibr(RawSampl(1,i),Channels,P,CorMean,
     &	    CalSonic,CalTherm,CalHyg,BadTc,Sample(1,i),N,Flag(1,i),
     &      Have_Uncal)
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
      CALL ECAverag(Sample,NMax,N,MMax,M,Flag,
     &	Mean,TolMean,Cov,TolCov,MIndep,CIndep,Mok,Cok)

C
C Print number of independent contributions
C
      IF (DoPrint.AND.PIndep) THEN
	CALL ECShwInd(OutF,QName,MIndep,CIndep,NMax,N,M,Freq)
      ENDIF
C
C Print averages of raw, calibrated data
C
      IF (DoPrint.AND.PCal) THEN
	WRITE(OutF,*)
	WRITE(OutF,*) 'For raw calibrated data : '
	WRITE(OutF,*)
	CALL ECShow(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
      ENDIF
C
C
C Subtract a linear trend from all quantities
C
C
      IF (DoDetren) THEN

	CALL ECDetren(Sample,NMax,N,MMAx,M,Mean,Cov,Sample,RC)
C
C
C Estimate mean values, covariances and tolerances of both
C for the detrended dataset
C
C
	CALL ECAverag(Sample,NMax,N,MMax,M,Flag,
     &	  Mean,TolMean,Cov,TolCov,MIndep,CIndep,Mok,Cok)

	IF (DoPrint.AND.PDetrend) THEN
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'After detrending : '
	  WRITE(OutF,*)

	  IF (PIndep) THEN
	    CALL ECShwInd(OutF,QName,MIndep,CIndep,NMax,N,M,Freq)
	  ENDIF

	  CALL ECShow(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
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
      CALL ECCorrec(OutF,
     &	DoPrint,
     &	Mean,NMax,N,TolMean,
     &	Cov,TolCov,
     &	QName,UName,
     &  BadTc,
     &	DoTilt,PTilt,PrePitch,PreYaw,PreRoll,
     &	DoPitch,PPitch,DirPitch,
     &	DoYaw,PYaw,YawLim,DirYaw,
     &	DoRoll,PRoll,RollLim,DirRoll,
     &	DoSonic,PSonic,SonFactr,
     &	DoO2,PO2,O2Factor,
     &	DoFreq,PFreq,LLimit,ULimit,Freq,CalSonic,CalTherm,CalHyg,FrCor,
     &	DoWebb,PWebb,P,Have_Uncal)
C
C If any transformation of coordinates was required (one of the options
C DoTilt, DoPitch, DoYaw or DoRoll was selected), then the numbers
C of independent samples of the respective quantities has to be
C re-estimated. It is not possible to "tilt the error-bars" on basis of
C the quantities which have been calculated until here.
C Therefore we return to the calibrated time series and
C make the transformations BEFORE averaging. Then averaging and corrections
C are repeated all over again.
C
      AnyTilt = ((DoTilt.OR.DoPitch).OR.(DoYaw.OR.DoRoll))

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
            CALL ECKPitch(Speed,NMax,N,DumCov,PrePitch)
            CALL ECKYaw(  Speed,NMax,N,DumCov,PreYaw  )
            CALL ECKRoll( Speed,NMax,N,DumCov,PreRoll )
          ENDIF

          IF (DoPitch) CALL ECKPitch(Speed,3,3,DumCov,DirPitch)
          IF (DoYaw  ) CALL ECKYaw(  Speed,3,3,DumCov,DirYaw  )
          IF (DoRoll ) CALL ECKRoll( Speed,3,3,DumCov,DirRoll )

          Sample(U,i) = Speed(1)
          Sample(V,i) = Speed(2)
          Sample(W,i) = Speed(3)

        ENDDO
C
C Reestablish the averages and covariances in the correct frame of reference.
C
	CALL ECAverag(Sample,NMax,N,MMax,M,Flag,
     &	  Mean,TolMean,Cov,TolCov,MIndep,CIndep,Mok,Cok)

	IF (DoPrint) THEN
	  WRITE(OutF,*)
	  WRITE(OutF,*) 'After untilting raw data : '
	  WRITE(OutF,*)

	  IF (PIndep) THEN
	    CALL ECShwInd(OutF,QName,MIndep,CIndep,NMax,N,M,Freq)
	  ENDIF

	  CALL ECShow(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
	ENDIF

        DumTilt  = (.FALSE.)
        DumPitch = (.FALSE.)
        DumYaw   = (.FALSE.)
        DumRoll  = (.FALSE.)
C
C Perform all necessary corrections on the mean values and (co-)variances.
C
        CALL ECCorrec(OutF,
     &	  DoPrint,
     &	  Mean,NMax,N,TolMean,
     &	  Cov,TolCov,
     &	  QName,UName,
     &    BadTc,
     &	  DumTilt,PTilt,PrePitch,PreYaw,PreRoll,
     &	  DumPitch,PPitch,DirPitch,
     &	  DumYaw,PYaw,YawLim,DirYaw,
     &	  DumRoll,PRoll,RollLim,DirRoll,
     &	  DoSonic,PSonic,SonFactr,
     &	  DoO2,PO2,O2Factor,
     &	  DoFreq,PFreq,LLimit,ULimit,Freq,
     &    CalSonic,CalTherm,CalHyg,FrCor,
     &	  DoWebb,PWebb,P, Have_Uncal)

      ENDIF
C
C
C Calculate fluxes from processed mean values and covariances.
C
C
      CALL ECFlux(Mean,NMax,Cov,TolMean,TolCov,p,BadTc,
     &	HSonic,dHSonic,HTc,dHTc,
     &	LvE,dLvE,LvEWebb,dLvEWebb,
     &	UStar,dUStar,Tau,dTau)

      RETURN
      END





      SUBROUTINE ECParams(InName,
     &	Freq,YawLim,RollLim,PrePitch,PreYaw,PreRoll,
     &	LLimit,ULimit,DoCrMean,DoDetren,DoSonic,DoTilt,DoPitch,DoYaw,
     &	DoRoll,DoFreq,DoO2,DoWebb,DoStruct,DoPrint,
     &	PRaw,PCal,PDetrend,PIndep,PTilt,PPitch,PYaw,
     &	PRoll,PSonic,PO2,PFreq,PWebb, StructSep)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECParams
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Purpose : Read settings for this particular session from file
C

      INCLUDE 'physcnst.for'

      REAL*8 Freq,YawLim,RollLim,PrePitch,PreYaw,PreRoll,LLimit,ULimit,
     &       StructSep
      LOGICAL DoCrMean,DoSonic,DoTilt,DoPitch,DoYaw,DoRoll,DoFreq,DoO2,
     &	DoWebb,DoPrint,PRaw,PCal,PIndep,PTilt,PPitch,PYaw,PRoll,PSonic,
     &	PO2,PFreq,PWebb,DoDetren,PDetrend,DoStruct
      CHARACTER*(*) InName
      INTEGER      IOCODE

      OPEN(TempFile,FILE=InName, STATUS='old', IOSTAT=IOCODE)
      IF (IOCODE .NE. 0) THEN
         WRITE(*,*) 'ERROR: could not open parameter file ', InName
         STOP
      ENDIF

      READ(TempFile,*) Freq	  ! [Hz] Sample frequency datalogger
      READ(TempFile,*)
      READ(TempFile,*) YawLim	  ! [degree] Limit when Mean(W) is turned to zero
      READ(TempFile,*) RollLim	  ! [degree] Limit when Cov(V,W) is turned to zero
      READ(TempFile,*)
      READ(TempFile,*) PrePitch   ! Fixed pitch angle for known tilt-correction
      READ(TempFile,*) PreYaw	  ! Fixed yaw angle for known tilt-correction
      READ(TempFile,*) PreRoll	  ! Fixed roll angle for known tilt-correction
      READ(TempFile,*)
      READ(TempFile,*) LLimit	  ! Smallest acceptable frequency-response correction factor
      READ(TempFile,*) ULimit	  ! Largest acceptable frequency-response correction factor
      READ(TempFile,*)
      READ(TempFile,*) DoCrMean   ! Replace mean quantities by better estimates
      READ(TempFile,*) DoDetren   ! Correct data for linear trend
      READ(TempFile,*) DoSonic	  ! Correct sonic temperature for humidity
      READ(TempFile,*) DoTilt	  ! Perform true tilt-correction with known angles
      READ(TempFile,*) DoPitch	  ! Turn system such that Mean(V) --> 0
      READ(TempFile,*) DoYaw	  ! Turn system such that Mean(W) --> 0
      READ(TempFile,*) DoRoll	  ! Turn System such that Cov(W,V) --> 0
      READ(TempFile,*) DoFreq	  ! Correct for poor frequency response
      READ(TempFile,*) DoO2	  ! Correct hygrometer for oxygen-sensitivity
      READ(TempFile,*) DoWebb	  ! Calculate mean velocity according to Webb
      READ(TempFile,*) DoStruct	  ! Calculate structure parameters
      READ(TempFile,*)
      READ(TempFile,*) DoPrint	  ! Skip printing intermediate results or not?
      READ(TempFile,*)
      READ(TempFile,*) PRaw	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PCal	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PDetrend   ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PIndep	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PTilt	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PPitch	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PYaw	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PRoll	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PSonic	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PO2	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PFreq	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) PWebb	  ! Indicator if these intermediate results are wanted
      READ(TempFile,*) StructSep! Separation (meter) for which to calculate structure parameter
      

      CLOSE(TempFile)

      RETURN
      END





      SUBROUTINE ECCorrec(OutF,
     &	DoPrint,
     &	Mean,NMax,N,TolMean,
     &	Cov,TolCov,
     &	QName,UName,
     &  BadTc,
     &	DoTilt,PTilt,PrePitch,PreYaw,PreRoll,
     &	DoPitch,PPitch,DirPitch,
     &	DoYaw,PYaw,YawLim,DirYaw,
     &	DoRoll,PRoll,RollLim,DirRoll,
     &	DoSonic,PSonic,SonFactr,
     &	DoO2,PO2,O2Factor,
     &	DoFreq,PFreq,LLimit,ULimit,Freq,CalSonic,CalTherm,CalHyg,FrCor,
     &	DoWebb,PWebb,P, Have_Uncal)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECCorrec
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Purpose : Integrated correction routine applying ALL (user-selected)
C	    corrections in this library on mean values and covariances.
C	    All intermediate results can be output to a file.
C	    Moreover they are returned to the calling routine in
C	    respective variables.
C
C Revision 28-05-2001: added info on whether uncalibrated data are
C                      available for a given variable (mainly important
C                      for sonic and/or Couple temperature since that
C                      is used for various corrections)
      INCLUDE 'physcnst.for'

      INTEGER N,NInt,NMAx,OutF, WhichTemp
      LOGICAL PPitch,PYaw,PRoll,PFreq,PO2,PWebb,
     &	DoPitch,DoYaw,DoRoll,DoFreq,DoO2,DoWebb,PSonic,
     &	DoSonic,DoTilt,PTilt,DoPrint,BadTc,
     &	QPitch,QYaw,QRoll,QFreq,QO2,QWebb,QSonic,QTilt,
     &  Have_Uncal(NMax)
      REAL*8 P,Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),LLimit,ULimit,FrCor(NMax,NMax),
     &	YawLim,RollLim,DirPitch,O2Factor(NMax),
     &	NSTA,NEND,TAUV,TauD,Freq,DirYaw,DirRoll,
     &	SonFactr(NMax),PrePitch,PreYaw,PreRoll
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ)
      CHARACTER*6 QName(NMax)
      CHARACTER*9 UName(NMax)
C
C
C Only print intermediate results if desired
C
C
      QTilt  = (PTilt  .AND. DoPrint)
      QPitch = (PPitch .AND. DoPrint)
      QYaw   = (PYaw   .AND. DoPrint)
      QRoll  = (PRoll  .AND. DoPrint)
      QSonic = (PSonic .AND. DoPrint)
      QO2    = (PO2    .AND. DoPrint)
      QFreq  = (PFreq  .AND. DoPrint)
      QWebb  = (PWebb  .AND. DoPrint)
C
C
C Perform a tilt-correction of the whole system using KNOWN (!!!!!!) angles.
C Used to compensate for KNOWN miss-alignment!!!!!
C The classic "tilt-corrections", which turn Mean(V), Mean(W) and Cov(W,V)
C to zero, to are carried out later!!!!!!
C If no fixed tilt-angles are known, this subroutine may be skipped.
C
C
      IF (DoTilt) CALL ECPTilt(QTilt,OutF,Mean,TolMean,
     &	NMax,N,Cov,TolCov,QName,UName,PrePitch,PreYaw,PreRoll)
C
C
C Perform classic pitch-correction : Mean(V) --> 0
C
C
      IF (DoPitch) CALL ECPPitch(QPitch,OutF,Mean,TolMean,
     &	NMax,N,Cov,TolCov,QName,UName,DirPitch)
C
C
C Perform classic yaw-correction : Mean(W) --> 0
C
C
      IF (DoYaw) CALL ECPYaw(QYaw,OutF,Mean,TolMean,NMax,
     &	N,Cov,TolCov,YawLim,QName,UName,DirYaw)
C
C
C Perform classic roll-correction : Cov(W,V) --> 0
C
C
      IF (DoRoll) CALL ECPRoll(QRoll,OutF,Mean,TolMean,
     &	NMax,N,Cov,TolCov,RollLim,QName,UName,DirRoll)
C
C
C Correct sonic temperature and all covariances with sonic temperature
C for humidity. This is half the Schotanus-correction. Side-wind
C correction is done at calibration time directly after reading raw data.
C
C Now we need to know if and which temperature we have. Default to
C Sonic temperature
C
      IF (HAVE_UNCAL(TSonic)) THEN
         WhichTemp = Tsonic
      ELSE IF (Have_UNCAL(TCouple)) THEN
         WhichTemp = TCouple
      ELSE
         WhichTemp = -1
      ENDIF

      IF (DoSonic) CALL ECPSchot(QSonic,OutF,Mean,TolMean,
     &	   NMax,N,Cov,TolCov,QName,UName,SonFactr)
C
C
C Perform correction for oxygen-sensitivity of hygrometer
C
C
      IF (DoO2) THEN
        IF (WhichTemp .GT. 0) THEN
            CALL ECPO2(QO2,OutF,Mean,TolMean,NMax,N,
     &	        Cov,TolCov,QName,UName,P,O2Factor, WhichTemp)
        ELSE
	    WRITE(*,*) 'ERROR: can not perform O2 correction without ',
     &                 'a temperature'
        ENDIF
      ENDIF
C
C
C Perform correction for poor frequency response and large paths
C
C

C
C Constants for integration routine in correction frequency response
C
      NSTA   = -5.D0	 ! [?] start frequency numerical integration
      NEND   = 0.69897D0 ! [?] end frequency numerical integration (LOG(5))
      NINT   = 19	 ! [1] number of intervals
      TAUV   = 0.D0	 ! [?] Low pass filter time constant
      TauD   = 0.D0	 ! [?] interval length for running mean

      IF (DoFreq) THEN
        IF (WhichTemp .GT. 0) THEN
          CALL ECPFreq(QFreq,OutF,Mean,TolMean,
     &    NMax,N,Cov,TolCov,QName,UName,LLimit,ULimit,NSta,NEnd,
     &    NInt,Freq,TauD,TauV,CalSonic,CalTherm,CalHyg,FrCor, WhichTemp)
        ELSE
	    WRITE(*,*) 'ERROR: can not perform freq. response correction',
     &                 ' without a temperature'
	ENDIF
      ENDIF
C
C
C Calculate mean vertical velocity according to Webb
C
C
      IF (DoWebb) THEN
         IF (WhichTemp .GT. 0) THEN
	    CALL ECPWebb(QWebb,OutF,Mean,TolMean,
     &	        NMax,N,Cov,TolCov,P,QName,UName, WhichTemp)
         ELSE
	    WRITE(*,*) 'ERROR: can not perform Webb correction',
     &                 ' without a temperature'
	 ENDIF
      ENDIF

      RETURN
      END





      SUBROUTINE ECReadDt(InName,Delay,Gain,Offset,x,NMax,N,MMax,M)
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
      INCLUDE 'physcnst.for'

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





      SUBROUTINE ECReadWou(InName,Delay,Gain,Offset,x,NMax,N,MMax,M)
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
      INCLUDE 'physcnst.for'

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


C...........................................................................
C FUNCTION:  CSI2HOUR
C PURPOSE:   To convert Campbell hour/min to decimal hours
C INTERFACE: HOURMIN    IN     Campbell hour/minutes (REAL)
C            return value      decimal hours (REAL)
C AUTHOR:    Arnold Moene
C DATE:      May 28, 2001
C...........................................................................
      REAL FUNCTION CSI2HOUR(HOURMIN)

      REAL HOURMIN
      INTEGER HOUR, MIN

      HOUR = INT(HOURMIN/100.0)
      MIN = INT(HOURMIN-HOUR*100)
      
      CSI2HOUR = HOUR + MIN/60.0
      END
      
C...........................................................................
C SUBROUTINE:  FINDTIME
C PURPOSE   :  To find the index in a NetCDF File for a given DOY and time
C...........................................................................
      SUBROUTINE FINDTIME(FDOY, FHOURMIN, NCID,
     +                    NCdoyID, NChmID, NCsecID,IND)
      INTEGER FDOY, FHOURMIN, NCID, NCdoyID, NChmID, NCsecID, IND,
     +        DIMID, MAXIND, ATTEMPTS, PREVCH
      REAL SEC1, SEC2, DOY1, DOY2, HM1, HM2, SAMPF, TIME1, TIME2,
     +     DOYD, HOURD, SECD, TIMED

      STATUS = NF_INQ_VARDIMID(NCID,NCsecID,DIMID)
      STATUS = NF_INQ_DIMLEN(NCID, DIMID,MAXIND )
      STATUS = NF_GET_VAR1_REAL(NCID,NCsecID, 1, SEC1)
      STATUS = NF_GET_VAR1_REAL(NCID,NCsecID, 2, SEC2)
      STATUS = NF_GET_VAR1_REAL(NCID,NCdoyID, 1, DOY1)
      STATUS = NF_GET_VAR1_REAL(NCID,NCdoyID, 2, DOY2)
      STATUS = NF_GET_VAR1_REAL(NCID,NChmID, 1, HM1)
      STATUS = NF_GET_VAR1_REAL(NCID,NChmID, 2, HM2)
      IF ((DOY1 .EQ. DOY2) .AND. (HM1 .EQ. HM2)) THEN
          SAMPF = SEC2-SEC1
      ELSE IF (DOY1 .EQ. DOY2) THEN
          TIME1 = CSI2HOUR(HM1)
          TIME2 = CSI2HOUR(HM2)
          SAMPF = (TIME2*3600+SEC2) - (TIME1*3600 + SEC1)
      ELSE
          TIME1 = CSI2HOUR(HM1)
          TIME2 = CSI2HOUR(HM2)
          SAMPF = (DOY2-DOY1)*3600*24 +
     +         (TIME2*3600+SEC2) - (TIME1*3600 + SEC1)
      ENDIF


      DOYD = FDOY - DOY1
      HOURD = CSI2HOUR(REAL(FHOURMIN)) - CSI2HOUR(HM1)
      SECD = 0 - SEC1
      TIMED = 3600.0*24*DOYD + HOURD*3600.0 + SECD
      ATTEMPTS = 0
      PREVCH = TIMED/SAMPF+1
      CURCH = TIMED/SAMPF
      IND = 0
      DO WHILE ((ABS(CURCH) .GT. 0) .AND. (ATTEMPTS .LT. 20))
C Determine time difference between current sample and required time
         IF (ABS(CURCH) .GT. ABS(PREVCH)) THEN
             CURCH = CURCH*0.8
         ENDIF
         IND = IND + CURCH
         PREVCH = CURCH
         IF ((IND .LT. 1) .OR. (IND .GT. MAXIND)) THEN
             WRITE(*,*)
     +         'ERROR: requested data appear to be not in this file'
             IND = -1
             RETURN
         ENDIF
         STATUS = NF_GET_VAR1_REAL(NCID,NCsecID, IND, SEC1)
         STATUS = NF_GET_VAR1_REAL(NCID,NCdoyID, IND, DOY1)
         STATUS = NF_GET_VAR1_REAL(NCID,NChmID, IND, HM1)
         STATUS = NF_GET_VAR1_REAL(NCID,NCsecID, IND+1, SEC2)
         STATUS = NF_GET_VAR1_REAL(NCID,NCdoyID, IND+1, DOY2)
         STATUS = NF_GET_VAR1_REAL(NCID,NChmID, IND+1, HM2)
         IF ((DOY1 .EQ. DOY2) .AND. (HM1 .EQ. HM2)) THEN
             SAMPF = SEC2-SEC1
         ELSE IF (DOY1 .EQ. DOY2) THEN
             TIME1 = CSI2HOUR(HM1)
             TIME2 = CSI2HOUR(HM2)
             SAMPF = (TIME2*3600+SEC2) - (TIME1*3600 + SEC1)
         ELSE
             TIME1 = CSI2HOUR(HM1)
             TIME2 = CSI2HOUR(HM2)
             SAMPF = (DOY2-DOY1)*3600*24+
     +         (TIME2*3600+SEC2) - (TIME1*3600 + SEC1)
         ENDIF
         DOYD = FDOY - DOY1
         HOURD = CSI2HOUR(REAL(FHOURMIN)) - CSI2HOUR(HM1)
         SECD = 0 - SEC1
         TIMED = 3600.0*24*DOYD + HOURD*3600.0 + SECD
         ATTEMPTS = ATTEMPTS + 1
         CURCH = TIMED/SAMPF
C
C A fix to get the right sample (needed in some cases)
         IF ((SECD .LT. 0) .AND. 
     +       ((INT(FHOURMIN) - INT(HM1)) .EQ. 1)) THEN
            CURCH = CURCH + 1
	 ENDIF
      ENDDO
      IF (ATTEMPTS .EQ. 20) THEN
         WRITE(*,*) 'ERROR: Gave up on searching time in file'
         IND = -1
      ENDIF
      END
      
      SUBROUTINE ECReadNCDF(InName,StartTime,StopTime,
     &  Delay,Gain,Offset,x,NMax,N,MMax,M, NCVarName,
     &  Have_Uncal)
C
C Read raw data from NETCDF-file
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
C Revision 28-05-2001: added info on whether uncalibrated data are
C                      available for a given variable (mainly important
C                      for sonic and/or Couple temperature since that
C                      is used for various corrections)
C
      INCLUDE 'physcnst.for'
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
     &        NATTS       ! # of attributes
      INTEGER I, J, K
      INTEGER NMax,N,MMax,M,Delay(NMax),HMStart,HMStop,Counter
      REAL Dum
      LOGICAL ok,Ready,Started,NotStopped
      REAL*8 x(NMax,MMax),Gain(NMax),Offset(NMax),StartTime(3),
     &  StopTime(3)
      CHARACTER*(*) NCVarName(NMax)
      LOGICAL       Have_Uncal(NMax)

      CHARACTER*255 ATT_NAME
      INTEGER       NCVarID(NNNMax), DoyStart, DoyStop,
     +              STARTIND, STOPIND, NSAMPLE
      LOGICAL       HAS_SCALE(NMax), HAS_OFFSET(NMax)
      REAL*8        NC_SCALE(NMax), NC_OFFSET(NMax)
      INTEGER       STRLEN
      EXTERNAL      STRLEN

      HMStart = 100*ANINT(StartTime(2))+ANINT(StartTime(3))
      HMStop  = 100*ANINT(StopTime( 2))+ANINT(StopTime( 3))
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

      N = NVars
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
        IF ((STRLEN(NCVarName(I)) .GT. 0) .AND.
     &      (NCVarID(I) .EQ. 0)) THEN
            WRITE(*,5110) NCVarName(I)(:STRLEN(NCVarName(I)))
            OK = (.FALSE.)
        ENDIF
      ENDDO
 5110 FORMAT('ERROR: variable ',A,' not found')
      IF (.NOT. OK) STOP

C
C Find start and stop index in file
C
      DOYStart = ANINT(StartTime(1))
      DOYStop = ANINT(StopTime(1))
      CALL FINDTIME(DoyStart, HMSTART, NCID, NCVarID(Doy),
     +             NCVarID(HourMin),
     +             NCVarID(Sec),
     +             STARTIND)
      CALL FINDTIME(DoyStop, HMSTOP, NCID, NCVarID(Doy),
     +             NCVarID(HourMin),
     +             NCVarID(Sec),
     +             STOPIND)
      IF (((STARTIND .LE. 0) .OR. (STOPIND .LE. 0)) .OR.
     +    (STOPIND .LT. STARTIND)) THEN
          WRITE(*,*) 'ERROR: No data found for given interval'
          STATUS = NF_CLOSE(NCID)
          IF (STATUS .NE. NF_NOERR) CALL EC_NCDF_HANDLE_ERR(STATUS)
          RETURN
      ENDIF
      NSAMPLE = STOPIND - STARTIND 
      IF ((STOPIND - STARTIND) .GT. MMMAX) THEN
         WRITE(*,5500) MMMAX, NSAMPLE
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
              Ready = (STATUS .EQ. nf_einvalcoords)
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
     &              ' size of storage (',I4,') is unsufficient ',
     &              ' for requested data (',I4,')')
C
C Close file again
C
      STATUS = NF_CLOSE(NCID)
      IF (STATUS .NE. NF_NOERR) CALL EC_NCDF_HANDLE_ERR(STATUS)

      RETURN
      END





      SUBROUTINE EC_NCDF_HANDLE_ERR(STATUS)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: EC_NCDF_HANDLE_ERR
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C To display error message of NETCDF-functions and stop program
C
      IMPLICIT NONE

      INCLUDE 'netcdf.inc'

      INTEGER STATUS
      IF (STATUS .NE. NF_NOERR) THEN
         WRITE(*,*) NF_STRERROR(STATUS)
         STOP 'Stopped'
      ENDIF
      END





C
C ########################################################################
C
C Small routine to estimate surface fluxes (with their tolerances!) from
C processed mean values, covariances and tolerances.
C
C ########################################################################
C





      SUBROUTINE ECFlux(Mean,NMax,Cov,TolMean,TolCov,p,BadTc,
     &	HSonic,dHSonic,HTc,dHTc,LvE,dLvE,LvEWebb,dLvEWebb,
     &	UStar,dUStar,Tau,dTau)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECFlux
C Version of subroutine : 1.03
C Date			: 21 October 1999
C Author		: Arjan van Dijk
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
C Construct estimates for surface fluxes from mean values and covariances
C
      INCLUDE 'physcnst.for'

      LOGICAL BadTc
      INTEGER NMax
      REAL*8 ECRhoWet,Mean(NMax),Cov(NMax,NMax),HSonic,dHSonic,HTc,dHTc,
     &	LvE,dLvE,LvEWebb,dLvEWebb,Tau,dTau,RhoSon,RhoTc,Frac1,Frac2,Sgn,
     &	p,TolMean(NMax),TolCov(NMax,NMAx),UStar,dUStar
C
C Sensible heat flux [W m^{-2}]
C
      RhoSon = ECRhoWet(Mean(Humidity),Mean(TSonic),p)
      HSonic = Cp*RhoSon*Cov(W,TSonic)
      dHSonic = Cp*RhoSon*TolCov(W,TSonic)

      IF (.NOT.BadTc) THEN
        RhoTc = ECRhoWet(Mean(Humidity),Mean(TCouple),p)
        HTc = Cp*RhoTc*Cov(W,TCouple)
        dHTc = Cp*RhoTc*TolCov(W,TCouple)
      ELSE
        HTc = -9999.D0
        dHTc = -9999.D0
      ENDIF
C
C Latent heat flux [W m^{-2}]
C
      LvE = Lv*Cov(W,Humidity)
      dLvE = Lv*TolCov(W,Humidity)
      LvEWebb = Lv*Mean(W)*Mean(Humidity)
C
C These few statements are eliminated to make sure that the
C error in the mean velocity is NOT YET taken into the error
C of the sensible heat. By uncommenting the following four commented
C statements, this is restored.
C
C      IF (ABS(Mean(W)).GT.1.D-10) THEN    ! statement 1
C	Frac1 = ABS(TolMean(W)/Mean(W))    ! statement 2
C      ELSE                                ! statement 3
	Frac1 = 0.D0
C      ENDIF                               ! statement 4

      Frac2 = ABS(TolMean(Humidity)/Mean(Humidity))

      dLvEWebb = LvEWebb*(Frac1**2.D0+Frac2**2.D0)**0.5D0
C
C Friction velocity
C
      IF (Cov(W,U).GT.0.D0) THEN
	Sgn = -1.D0
	UStar = -SQRT(Cov(W,U))
      ELSE
	Sgn = 1.D0
	UStar = SQRT(-Cov(W,U))
      ENDIF

      dUStar = 0.5*ABS(TolCov(W,U)/Cov(W,U))*UStar
C
C Friction force [N m^{-2}]
C
      IF (.NOT.BadTc) THEN
        Tau = Sgn*0.5*(RhoSon+RhoTc)*(ABS(UStar))**2.D0
      ELSE
        Tau = Sgn*(RhoSon)*(ABS(UStar))**2.D0
      ENDIF
      dTau = ABS(TolCov(W,U)/Cov(W,U))*Tau

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





      SUBROUTINE ECReadAp(InName,CalSpec)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECReadAp
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Read specifications of an apparatus from file
C
      INCLUDE 'physcnst.for'

      CHARACTER*(*) InName
      REAL*8 CalSpec(NQQ)
      INTEGER i, IOCODE
      INTEGER STRLEN
      EXTERNAL STRLEN

      OPEN(TempFile,FILE=InName,STATUS='OLD', IOSTAT=IOCODE)
      IF (IOCODE .NE. 0) THEN
         WRITE(*,5100) InName(:STRLEN(InName))
         STOP
      ENDIF
 5100 FORMAT ('ERROR: can not open calibration file ', (A))

C First get the type of apparatus
      READ(TempFile,*) CalSpec(1)

C Now read the correct number of specs
      DO i=2,ApNQQ(CalSpec(1))
       	READ(TempFile,*,IOSTAT=IOCODE, END = 9000) CalSpec(i)
         IF (IOCODE .NE. 0) THEN
            WRITE(*,*) 'ERROR in reading of configuration file'
            STOP
         ENDIF
      ENDDO
 9000 CONTINUE
      IF ((I-1) .NE. ApNQQ(CalSpec(1))) THEN
         WRITE(*,*) 'ERROR: did not find correct number of specs in ',
     &               InName(:STRLEN(InName))
      ENDIF
      
      CLOSE(TempFile)

      RETURN
      END





C
C ########################################################################
C
C Some silly routines to print intermediate results to file. Handy
C for checking the origin of unexpected results!
C
C ########################################################################
C





      SUBROUTINE ECShow(OutF,QName,UName,Mean,TolMean,Cov,
     &	TolCov,NMax,N)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECShow
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Prints mean values and (co-)variances plus respective tolerances
C
      INCLUDE 'physcnst.for'

      INTEGER i,j,N,NMax,OutF
      REAL*8 Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),Dum(NNMax)
      CHARACTER*6 QName(NMax)
      CHARACTER*9 UName(NMax)

      WRITE(OutF,*)
      WRITE(OutF,*) 'Mean quantities : '
      WRITE(OutF,*)

      DO i=1,N
	WRITE(OutF,10) QName(i),Mean(i),TolMean(i),UName(i)
      ENDDO
 10   FORMAT(a6,' = ',F20.10,' +/- ',F20.10,1X,a9)

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
 40   FORMAT(a6,20(a8:,1X))

      DO i=1,N
	DO j=1,N
	  IF (ABS(Cov(i,i)*Cov(j,j)).GT.Epsilon) THEN
	    Dum(j) = Cov(i,j)/SQRT(Cov(i,i)*Cov(j,j))
	  ELSE
	    Dum(j) = 0.D0
	  ENDIF
	ENDDO
	WRITE(OutF,50) QName(i),(Dum(j),j=1,N)
 50	FORMAT(a6,20(F8.4:,1X))
      ENDDO

      WRITE(OutF,*)
      WRITE(OutF,*) '###############################################'
      WRITE(OutF,*)

      RETURN
      END





      SUBROUTINE ECShwInd(OutF,QName,MIndep,
     &	CIndep,NMax,N,M,Freq)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECShwInd
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Prints number of independent observations
C
      INCLUDE 'physcnst.for'

      INTEGER i,j,N,NMax,OutF,MIndep(NMax),CIndep(NMax,NMax),M
      REAL*8 Freq
      CHARACTER*6 QName(NMax)

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





      SUBROUTINE ECShwFrq(OutF,QName,FrCor,NMax,N)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECShwFrq
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Prints correction factors associated with frequency response
C
      INCLUDE 'physcnst.for'

      INTEGER i,j,N,NMax,OutF
      REAL*8 FrCor(NMax,NMax)
      CHARACTER*6 QName(NMax)

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





C
C ########################################################################
C
C The following (large) routine is the main body of the data-reduction.
C From large time-series we extract mean quantities and co-variances.
C For both mean values and for covariances the tolerances are estimated.
C
C ########################################################################
C





      SUBROUTINE ECAverag(x,NMax,N,MMax,M,Flag,
     &	Mean,TolMean,Cov,TolCov,MIndep,CIndep,Mok,Cok)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECAverag
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C From the given set of calibrated samples, calculate the averages,
C variances, covariances and their tolerances.
C
C input : x : array with physical dimensions (NMax,MMax)
C	      First index counts quantities; second counter
C	      counts samples. Only the first N quantities
C	      and the first M samples are used.
C	  flag : LOGICAL array with physical dimensions NMAx by MMAx,
C	      out of which only N by M are used. If flag(j,i) is true, then
C	      quantity j in sample i is not ok.
C
C output : Mean : array with physical dimension NMax out of
C		   which only N are used. Mean contains the average value
C		   of array x. Only samples with Flag = 0 are used.
C	   TolMean : tolerance of Mean, defined as
C		       2 * sigma / sqrt(NIndep)
C		   where sigma is the standarddeviation of the quantities
C		   and where the number of independent samples is estimated
C		   as twice the number of sign-changes of the fluctuations
C		   of the respective quantities around their means.
C	   Cov : REAL*8 (NMax,NMax) : covariances
C	   TolCov : tolerances of Cov, estimated as tolerances of mean.
C	   MIndep : INTEGER(NMAx) : Number of independent samples
C		   in time series from which means are calculated
C	   CIndep : INTEGER(NMAx,NMAx) : Number of independent samples
C		   in time series from which covariances are calculated
C
      INCLUDE 'physcnst.for'

      INTEGER NMax,N,MMax,M,i,j,k,NChange(NNMax,NNMax),
     &	MIndep(NMax),CIndep(NMax,NMax),NMin(NNMax,NNMax),
     &  Mok(NMax),Cok(NMax,NMax),NPlus(NNMax,NNMax)
      LOGICAL Flag(NMax,MMax)
      REAL*8 x(NMax,MMax),RawMean(NNMax),Mean(NMax),TolMean(NMax),
     &	Cov(NMax,NMax),TolCov(NMax,NMax),xPrime(NNMax),
     &	PrevPrime(NNMax),dTolCov(NNMax,NNMax),PSwap
C
C Initialise the arrays
C
      DO i=1,NMax
	RawMean(i) = 0.D0
	Mean(i) = 0.D0
	TolMean(i) = 0.D0
	Mok(i) = 0
	DO j=1,NMax
	  Cov(i,j) = 0.D0
	  TolCov(i,j) = 0.D0
	  Cok(i,j) = 0
	ENDDO
      ENDDO
C
C Find a rough estimate for the mean
C
      DO i=1,M
	DO j=1,N
	  IF (.NOT.Flag(j,i)) THEN
	    Mok(j) = Mok(j) + 1
	    RawMean(j) = RawMean(j) + x(j,i)
	  ENDIF
	ENDDO
      ENDDO

      DO j=1,N
	IF (Mok(j).GT.0) RawMean(j) = RawMean(j)/DBLE(Mok(j))
      ENDDO
C
C Find final estimates for mean
C
      DO i=1,M
	DO j=1,N
	  IF (.NOT.Flag(j,i)) THEN
	    Mean(j) = Mean(j) + (x(j,i) - RawMean(j))
	  ENDIF
	ENDDO
      ENDDO

      DO j=1,N
	IF (Mok(j).GT.0) Mean(j) = Mean(j)/DBLE(Mok(j))
	Mean(j) = Mean(j) + RawMean(j)
      ENDDO
C
C Find (co-)variances and from them the tolerances of the mean
C
      DO j=1,N
	NChange(j,j) = 0
        NMin(j,j) = 0
        NPlus(j,j) = 0
      ENDDO

      DO i=1,M
	DO j=1,N
	  IF (.NOT.Flag(j,i)) THEN
	    xPrime(j) = x(j,i) - Mean(j)
            IF (xPrime(j).LT.0.D0) NMin(j,j) = NMin(j,j)+1
            IF (xPrime(j).GT.0.D0) NPlus(j,j) = NPlus(j,j)+1
	    IF (i.GT.1) THEN
	      IF (PrevPrime(j)*xPrime(j).LE.0.D0)
     &		NChange(j,j) = NChange(j,j) + 1
	    ENDIF
	    PrevPrime(j) = xPrime(j)
	  ENDIF
	ENDDO


	DO j=1,N
	  DO k=j,N
	    IF (.NOT.(Flag(j,i).OR.Flag(k,i))) THEN
	      Cok(j,k) = Cok(j,k) + 1
	      Cov(j,k) = Cov(j,k) + xPrime(j)*xPrime(k)
	    ENDIF
	  ENDDO
	ENDDO
      ENDDO

      DO j=1,N
	DO k=j,N
	  IF (Cok(j,k).GT.0) Cov(j,k) = Cov(j,k)/DBLE(Cok(j,k))
	ENDDO
      ENDDO

      DO j=1,N
	DO k=1,N
	  IF (k.LT.j) Cov(j,k) = Cov(k,j)
	ENDDO
      ENDDO

      DO j=1,N
        PSwap = 2.D0*DBLE(NMin(j,j))*DBLE(NPlus(j,j))/
     &    DBLE(NMin(j,j) + NPlus(j,j))**2.D0
	MIndep(j) = ANINT(DBLE(NChange(j,j))/PSwap) - 1
	MIndep(j) = MAX(MIndep(j),1)
	TolMean(j) = 2.D0*(Cov(j,j)/DBLE(MIndep(j)))**0.5D0
      ENDDO
C
C Find tolerances for (co-)variances
C
      DO i=1,M
        DO j=1,N
	  IF (.NOT.Flag(j,i)) THEN
	    xPrime(j) = x(j,i) - Mean(j)
	  ENDIF
	ENDDO
	DO j=1,N
	  DO k=j,N
	    IF (.NOT.(Flag(j,i).OR.Flag(k,i))) THEN
	      dTolCov(j,k)=xPrime(j)*xPrime(k)-Cov(j,k)
	      TolCov(j,k)=TolCov(j,k)+(dTolCov(j,k))**2
	    ENDIF
	  ENDDO
	ENDDO
      ENDDO

      DO j=1,N
	DO k=j,N
C
C Calculate the standard deviation of the instantaneous contributions
C to the covariances.
C
	  IF (Cok(j,k).GT.0) TolCov(j,k)=TolCov(j,k)/DBLE(Cok(j,k))
C
C Here we estimate the number of independent contributions to the
C covariance by counting the least number of independent contributions
C to either of the factors.
C
	  CIndep(j,k) = MIN(MIndep(j),MIndep(k))
	  TolCov(j,k) = TolCov(j,k)/DBLE(CIndep(j,k))
C
C Tolerance is defined as 2*sigma, where sigma is standard deviation
C
	  TolCov(j,k) = 2.D0*(TolCov(j,k))**0.5
	ENDDO
      ENDDO

      DO j=1,N
	DO k=1,(j-1)
	  TolCov(j,k) = TolCov(k,j)
	  CIndep(j,k) = CIndep(k,j)
	ENDDO
      ENDDO

      RETURN
      END





C
C ########################################################################
C
C This routine constructs a trend-corrected time-series from a given
C time-series: The model is an additive model. The straight line from
C least squares regression (linear trend) is subtracted from the dataset.
C
C ########################################################################
C





      SUBROUTINE ECDetren(x,NMax,N,MMAx,M,Mean,Cov,y,RC)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECDetren
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Construct a linearly detrended dataset from a given dataset
C
C input : x : REAL*8(NMax,MMax) : x(i,j) = quantity i in sample j
C		only the first N quantities and the first M samples
C		are used.
C	  Cov : REAL*8 (NMax,NMAx) : covariances of quantities.
C		used to find trend.
C output : y : REAL*8(NMAx,MMax) : detrended timeseries.
C	   RC : REAL*8(NMAx) : Directional coefficients of linear
C		regression trend-lines.
C
      INCLUDE 'physcnst.for'

      INTEGER NMax,N,i,MMax,M,j
      REAL*8 x(NMax,MMax),Mean(NMax),Cov(NMax,NMax),RC(NMax),
     &	y(NMax,MMax),Trend

      DO j = 1,N
	RC(j) = Cov(TTime,j)/Cov(TTime,TTime)
      ENDDO

      DO i= 1,M
	DO j= 1,N
	  IF (j.NE.TTime) THEN
	    Trend = RC(j)*(x(TTime,i)-Mean(TTime))
	    y(j,i) = x(j,i) - Trend
	  ENDIF
	ENDDO
      ENDDO

      RETURN
      END





C
C ########################################################################
C
C A routine to calculate the value of certain simple basefunctions
C
C ########################################################################
C





      REAL*8 FUNCTION ECBaseF(x,FType,Order,C)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECBaseF
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Purpose : Calculate simple functions. Currently implemented:
C  - Ordinary polynomials
C  - Polynomials in the natural logarithm of x
C
      INCLUDE 'physcnst.for'

      INTEGER FType,Order,i
      REAL*8 x,C(0:Order),Dum,LogX

      Dum = 0.D0
      IF (FType.EQ.NormPoly) THEN
	DO i=0,Order
	  Dum = Dum + C(i)*x**i
	ENDDO
      ELSE IF (FType.EQ.LogPoly) THEN
	LogX = LOG(x)
	DO i=0,Order
	  Dum = Dum + C(i)*LogX**i
	ENDDO
      ENDIF

      ECBaseF = Dum

      RETURN
      END





C
C ########################################################################
C
C The following set of (small) routines use the ideal gas law to
C calculate densities and humidities for wet air.
C
C ########################################################################
C





      REAL*8 FUNCTION ECRhoDry(RhoV,T,P)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECRhoDry
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Calculate the density of dry air component in wet air
C
C input : RhoV : Density of water vapour [kg m^{-3}]
C	  T    : Temperature		 [K]
C	  P    : Pressure		 [Pa]
C
C output : RhoDry : Density of dry air component [kg m^{-3}]
C
C Via Dalton's law : Pressure is sum of partial pressures :
C	      P = RhoV*Rv*T + RhoD*Rd*T
C
      INCLUDE 'physcnst.for'

      REAL*8 RhoV,T,P

      ECRhoDry = P/(Rd*T) - RhoV*Rv/Rd	  ! [kg m^{-3}]

      RETURN
      END





      REAL*8 FUNCTION ECRhoWet(RhoV,T,P)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECRhoWet
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Calculate the density of wet air
C
C input : RhoV : Density of water vapour [kg m^{-3}]
C	  T    : Temperature		 [K]
C	  P    : Pressure		 [Pa]
C
C output : RhoWet : Density of wet air [kg m^{-3}]
C
      INCLUDE 'physcnst.for'

      REAL*8 RhoV,T,P,ECRhoDry

      ECRhoWet = ECRhoDry(RhoV,T,P) + RhoV ! [kg m^{-3}]

      RETURN
      END





      REAL*8 FUNCTION ECQ(RhoV,T,P)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECQ
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Calculate the specific humidity of wet air
C
C input : RhoV : Density of water vapour [kg m^{-3}]
C	  T    : Temperature		 [K]
C	  P    : Pressure		 [Pa]
C
C output : ECQ : specific humidity	 [kg kg^{-1}]
C
      INCLUDE 'physcnst.for'

      REAL*8 RhoV,T,P,ECRhoWet,RhoWet

      RhoWet = ECRhoWet(RhoV,T,P)
      ECQ = RhoV/RhoWet 		 ! [kg/kg]

      RETURN
      END





C
C ########################################################################
C
C Here comes a set routines that actually implement all corrections
C necessary to process eddy-correlation measurements.
C They do contain interesting models and physics/mathematics!
C All following routines affect mean quantities and co-variances only.
C No raw data-series necessary. Therefore one can use these routines
C on pre-processed data (raw mean values and raw co-variances).
C
C ########################################################################
C





      SUBROUTINE ECMapVec(a,x,y)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECMapVec
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Calculates the image of "x" under the map "a"; y(i) = a(ji)x(j)
C
C input : a : REAL*8 (3,3) : the mapping matrix
C	  x : REAL*8 (3)   : the vector to be mapped
C output : y : REAL*8 (3)  : the image of the map
C
      IMPLICIT NONE

      REAL*8 a(3,3),x(3),y(3),Dum(3)
      INTEGER I,J

      DO I=1,3
	Dum(I) = 0.D0
      ENDDO

      DO I=1,3
	DO J=1,3
	  Dum(I) = Dum(I) + a(J,I)*x(J)
	ENDDO
      ENDDO

      DO I=1,3
	y(I) = Dum(I)
      ENDDO

      RETURN
      END





      SUBROUTINE ECMapMtx(a,x,y)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECMapMtx
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Calculates the image of "x" under the map "a";
C y(ji) = a(ki)a(lj)x(lk)
C
C input : a : REAL*8 (3,3) : the mapping matrix
C	  x : REAL*8 (3,3) : the tensor to be mapped
C output : y : REAL*8 (3,3) : the image of the map
C
      IMPLICIT NONE

      REAL*8 a(3,3),x(3,3),y(3,3),Dum(3,3)
      INTEGER I,J,K,L

      DO I=1,3
	 DO J=1,3
	    Dum(I,J) = 0.D0
	 ENDDO
      ENDDO

      DO I=1,3
	 DO J=1,3
	    DO K=1,3
	       DO L=1,3
		  Dum(J,I) = Dum(J,I) + a(K,I)*a(L,J)*x(L,K)
	       ENDDO
	    ENDDO
	 ENDDO
      ENDDO

      DO I=1,3
	 DO J=1,3
	    y(J,I) = Dum(J,I)
	 ENDDO
      ENDDO

      RETURN
      END





      SUBROUTINE ECShake(Mean,NMax,N,Cov,Speed,Stress,DumVecs,NN)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECShake
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Help routine for routines for correction of coordinate system
C
      INCLUDE 'physcnst.for'

      INTEGER i,j,NMax,N,NN
      REAL*8 Stress(3,3),Speed(3),DumVecs(3,4:NN),Mean(NMax),
     &	Cov(NMax,NMax)
C
C Loop : All velocity components
C
      DO i = U,W
	Speed(i) = Mean(i)
	DO j = U,W
	  Stress(i,j) = Cov(i,j)
	ENDDO
	DO j = 4,N
	  DumVecs(i,j) = Cov(i,j)
	ENDDO
      ENDDO

      RETURN
      END





      SUBROUTINE ECUnShak(Speed,Stress,DumVecs,NN,Mean,NMax,N,Cov)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECUnShak
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Help routine for routines for correction of coordinate system
C
      INCLUDE 'physcnst.for'

      INTEGER i,j,NMax,N,NN
      REAL*8 Stress(3,3),Speed(3),DumVecs(3,4:NN),Mean(NMax),
     &	Cov(NMax,NMax)
C
C Loop : All velocity components
C
      DO i = U,W
	Mean(i) = Speed(i)
	DO j = U,W
	  Cov(i,j) = Stress(i,j)
	ENDDO
	DO j = 4,N
	  Cov(i,j) = DumVecs(i,j)
	  Cov(j,i) = Cov(i,j)
	ENDDO
      ENDDO

      RETURN
      END





      SUBROUTINE ECRotate(Mean,NMax,N,Cov,Map)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECRotate
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Routine to change coordinate system according to tensor "Map".
C This routine is called by all tilt-correction procedures.
C Both the mean velocity and the Reynoldsstresses and the
C covariances of all velocity components with other quantities
C are rotated.
C
      INCLUDE 'physcnst.for'

      INTEGER NMax,N,j
      REAL*8 Map(3,3),Stress(3,3),Speed(3),
     &  DumVecs(3,4:NNmax),Mean(NMax),Cov(NMax,NMax)

      CALL ECShake(Mean,NMax,N,Cov,Speed,Stress,DumVecs,NNMax)
      CALL ECMapVec(Map,Speed,Speed)
      CALL ECMapMtx(Map,Stress,Stress)
      DO j = 4,N
        CALL ECMapVec(Map,DumVecs(1,j),DumVecs(1,j))
      ENDDO
      CALL ECUnShak(Speed,Stress,DumVecs,NNMax,Mean,NMax,N,Cov)

      RETURN
      END





      SUBROUTINE ECKPitch(Mean,NMax,N,Cov,Direction)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECKPitch
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Rotate coordinate system about a KNOWN pitch-angle around the vertical
C
      INCLUDE 'physcnst.for'

      INTEGER NMax,N
      REAL*8 SinPhi,CosPhi,Pitch(3,3),Mean(NMax),Cov(NMax,NMax),
     &	Direction

      SinPhi = SIN(PI*Direction/180.D0)
      CosPhi = COS(PI*Direction/180.D0)

      Pitch(1,1) = CosPhi
      Pitch(1,2) = -SinPhi
      Pitch(1,3) = 0.D0
      Pitch(2,1) = SinPhi
      Pitch(2,2) = CosPhi
      Pitch(2,3) = 0.D0
      Pitch(3,1) = 0.D0
      Pitch(3,2) = 0.D0
      Pitch(3,3) = 1.D0

      CALL ECRotate(Mean,NMax,N,Cov,Pitch)

      RETURN
      END





      SUBROUTINE ECPitch(Mean,NMax,N,Cov,Direction)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECPitch
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Transform coordinate system such that v_mean = 0
C (no mean lateral horizontal velocity component)
C
      INCLUDE 'physcnst.for'

      INTEGER NMax,N
      REAL*8 UHor,SinPhi,CosPhi,Mean(NMax),Cov(NMax,NMax),Direction

      UHor = ((Mean(U))**2+(Mean(V))**2)**0.5
      SinPhi = Mean(V)/UHor
      CosPhi = Mean(U)/UHor
      Direction = 180.D0*ACOS(CosPhi)/PI
      IF (SinPhi.LT.0.D0) Direction = 360.D0-Direction

      CALL ECKPitch(Mean,NMax,N,Cov,Direction)

      RETURN
      END





      SUBROUTINE ECKYaw(Mean,NMax,N,Cov,Direction)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECKYaw
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Rotate coordinate system about a KNOWN yaw-angle around vector (0,1,0).
C
      INCLUDE 'physcnst.for'

      INTEGER NMax,N
      REAL*8 SinTheta,CosTheta,Yaw(3,3),Mean(NMax),Cov(NMax,NMax),
     &	Direction

      SinTheta = SIN(PI*Direction/180.D0)
      CosTheta = COS(PI*Direction/180.D0)

      Yaw(1,1) = CosTheta
      Yaw(1,2) = 0.D0
      Yaw(1,3) = -SinTheta
      Yaw(2,1) = 0.D0
      Yaw(2,2) = 1.D0
      Yaw(2,3) = 0.D0
      Yaw(3,1) = SinTheta
      Yaw(3,2) = 0.D0
      Yaw(3,3) = CosTheta

      CALL ECRotate(Mean,NMax,N,Cov,Yaw)

      RETURN
      END





      SUBROUTINE ECYaw(Mean,NMax,N,Cov,YawLim,Direction)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECYaw
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Transform coordinate system such that w_mean = 0
C (no mean vertical velocity component)
C
      INCLUDE 'physcnst.for'

      INTEGER NMax,N
      REAL*8 SinTheta,Mean(NMax),Cov(NMax,NMax),Direction,UTot,YawLim

      UTot = ((Mean(U))**2+(Mean(W))**2)**0.5
      SinTheta = Mean(W)/UTot
      Direction = 180.D0*ASIN(SinTheta)/PI

      IF (ABS(Direction) .LE. YawLim)
     &	CALL ECKYaw(Mean,NMax,N,Cov,Direction)

      RETURN
      END





      SUBROUTINE ECKRoll(Mean,NMax,N,Cov,Direction)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECKRoll
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Rotate coordinate system about a KNOWN roll-angle around vector (1,0,0).
C
      INCLUDE 'physcnst.for'

      INTEGER NMax,N
      REAL*8 SinRoll,CosRoll,Roll(3,3),Mean(NMax),Cov(NMax,NMax),
     &	Direction

      SinRoll = SIN(PI*Direction/180.D0)
      CosRoll = COS(PI*Direction/180.D0)

      Roll(1,1) = 1.D0
      Roll(1,2) = 0.D0
      Roll(1,3) = 0.D0
      Roll(2,1) = 0.D0
      Roll(2,2) = CosRoll
      Roll(2,3) = -SinRoll
      Roll(3,1) = 0.D0
      Roll(3,2) = SinRoll
      Roll(3,3) = CosRoll

      CALL ECRotate(Mean,NMax,N,Cov,Roll)

      RETURN
      END





      SUBROUTINE ECRoll(Mean,NMax,N,Cov,RollLim,Direction)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECRoll
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Transform coordinate system such that Cov(V,W) = 0
C (vertical velocity fluctuations are independent from horizontal
C fluctuations)
C
      INCLUDE 'physcnst.for'

      INTEGER NMax,N
      REAL*8 Mean(NMax),Cov(NMax,NMax),Direction,RollLim,RollAngl,
     &	Arg1,Arg2

      Arg1 = 2*Cov(V,W)
      Arg2 = -Cov(V,V)+Cov(W,W)
      IF (Arg2.LT.0.D0) THEN
        RollAngl =  0.5D0*ATAN2(Arg1,-Arg2)
      ELSE
        RollAngl = -0.5D0*ATAN2(Arg1, Arg2)
      ENDIF
      Direction = 180.D0*RollAngl/Pi

      IF (ABS(Direction) .LE. RollLim) THEN
	CALL ECKRoll(Mean,NMax,N,Cov,Direction)
      ENDIF

      RETURN
      END





      SUBROUTINE ECFrResp(Mean,Cov,TolCov,NMax,NSize,Lower,Upper,
     &	NSta,NEnd,NInt,NS,TauD,TauV,CalSonic,CalTherm,CalHyg,WXT,
     & WhichTemp)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECFrResp
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
C Based on old routine by Van de Hurk, Elbers and Nieveen.
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
C Calculate frequency response corrections for sensor response, path
C length averaging, sensor separation, signal processing and dampening
C of fluctuations in a tube. Based on publications;
C
C	Moore, C.J. (1986): 'Frequency Response Corrections for Eddy
C	Correlation Systems'. Boundary Layer Met. 37: 17-35.
C
C	Philip, J.R. (1963): 'The Damping of Fluctuating Concentration
C	by Continuous Sampling Through a tube' Aust. J. Phys. 16: 454-463.
C
C	Leuning, R. and K.M. King (1991): 'Comparison of Eddy-Covariance
C	Measurements of CO2 Fluxes by open- and closed-path CO2 analysers'
C	(unpublished)
C
C Revision 28-05-2001: added info on which temperature should be used
C                      in corrections (Sonic or thermocouple)
      INCLUDE 'physcnst.for'

      INTEGER I1,I2,I3,I4,I,J,NINT,NMax,NSize, WhichTemp
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ),
     &	GAIND,GAINV,GAINT1,GAINT2,GAINUV,GAINW,GAINQ1,
     &	TRANA,TRANQ1,TRANUV,TRANW,TRANT1,TRANT2,TRANWT1,TRANWT2,TRANWQ1,
     &	TRANWUV,ATT,AWW,AUU,AWT,AUW,BTT,BWW,BUU,BWT,BUW,TAUD,TAUV,TAUT,
     &	XI,ZETA,CW,CU,WT,UW,TT,UU,WW,F,NS,C,NSTA,NEND,N,LF,X,ZL,
     &	INTTS(15),G(15),WXT(NMax,NMax),Mean(NMax),Cov(NMax,NMax),UStar,
     &	Lower,Upper,TolCov(NMax,NMax)
C
C-- Declarations of constants ----------------------------------------
C
      TAUT = CalTherm(QQTime)/
     &	(1.D0+4.9D0*(CalTherm(QQTime)**0.5*Mean(U))**0.45D0)
      LF   = (NEND-NSTA)/DBLE(NINT-1)	    ! Interval width num. integr.
      C    = 1.D0/(LOG(10.D0)*LF/3.D0)	    ! Constant of integration, eq.1
      LF   = 10.D0**LF
      N    = 10.D0**NSTA
      DO I=1,15
	INTTS(I)=0.D0
      ENDDO

      Ustar = ((Cov(U,W))**2.D0)**0.25D0
      ZL = (CalSonic(QQZ)*Karman*GG*Cov(W,WhichTemp))/
     &	(-Mean(WhichTemp)*Ustar**3.D0)

      IF(ZL.GE.0.D0) THEN
C
C-- Stable spectral and co-spectral function factors, eq. 20b, 21b
C
	ATT=0.0961D0+0.644D0*ZL**0.6D0
	AWW=0.838D0+1.172D0*ZL
	AUU=0.2D0*AWW
	AWT=0.284D0*(1.D0+6.4D0*ZL)**0.75D0
	AUW=0.124D0*(1.D0+7.9D0*ZL)**0.75D0
	BTT=3.124D0/ATT**0.667D0
	BWW=3.124D0/AWW**0.667D0
	BUU=3.124D0/AUU**0.667D0
	BWT=2.34D0/AWT**1.1D0
	BUW=2.34D0/AUW**1.1D0
      ELSE
C
C-- Unstable spectral and co-spectral function factors, eq. 23 -----
C
	XI   = (-ZL)**0.667D0
	ZETA = (0.001D0*CalSonic(QQZ))**1.667D0
	CW   = 0.7285D0+1.4115D0*XI
	CU   = 9.546D0+1.235D0*XI/(ZETA**0.4D0)
      END IF
C
C-- Start Simpson's rule numerical integration -----------------------
C     From -5 to log(5) in NInt steps
C
      I3=0
      DO I1=1,10
	DO I2=2,4,2
	  I3=I3+1
	  IF (I3.GT.NInt) GOTO 200
	  I4=I2
	  IF (I3.EQ.1)	I4=I4-1
	  IF (I3.EQ.NInt) I4=I4-1
C
C-- Skip integration if frequency exceeds Nyquist frequency ------
C
	  IF (N.GT.(0.5D0*NS)) GOTO 999

C
C-- ELECTRONIC / DIGITAL FILTERS !!
C-- Auto Regressive Moving Average filter response gain, eq. 19 --
C
	  GAIND=1.D0
	  X=2.D0*PI*N*TAUD
	  IF (X.LT.6.D0 .AND. X.GT.0.D0) GAIND=X*X/(1.D0+X*X)
C
C-- Butterworth filter frequency response gain, eq. 5 ------------
C
	  GAINV = 1.D0
	  GAINV=(2.D0*PI*N*TAUV)**2.D0
	  GAINV=1.D0+GAINV*GAINV
C
C-- Data aquisition co-spectral transfer function, eq. 15 (b=3) --
C
	  TRANA=1.D0+(N/(NS-N))**3.D0
C
C-- LUW THERMOCOUPLE TEMPERATURE !!
C-- Thermocouple frequency response gain, eq. 2 ------------------
C
	  GAINT1=1.D0+(2.D0*PI*N*TAUT)**2.D0
C
C-- Thermocouple spatial averaging transfer function, eq. 7 ------
C
	  TRANT1=1.D0	       !spatial averaging negligible
C
C-- W-T lateral separation transfer function, eq. 11 -------------
C
	  TRANWT1=1.D0
	  X=N/Mean(U)*(CalTherm(QQX)-CalSonic(QQX))
	  IF (X.GT.0.01D0) TRANWT1=EXP(-9.9D0*X**1.5D0)
C
C-- SOLENT SONIC TEMPERATURE !!
C-- Sonic temperature frequency response gain, eq. 2 -------------
C
	  GAINT2=1.D0	       !sonic temperature gain negligible
C
C-- Sonic temperature spatial averaging transfer function, eq.7 --
C
	  TRANT2=1.D0
	  X=2.D0*PI*N/Mean(U)*CalSonic(QQExt4)
	  IF (X.GT.0.02D0) TRANT2=(3.D0+EXP(-X)-4.D0*(1.D0-EXP(-X))/X)/X
C
C-- W-T lateral separation transfer function, eq. 11 -------------
C
	  TRANWT2=1.D0	       !there is no lateral separation

C
C-- SOLENT SONIC HORIZONTAL WINDSPEED !!
C-- UV-sensor frequency response gain, eq. 2 ---------------------
C
	  GAINUV=1.D0	       !sonic UV windspeed gain negligible
C
C-- UV-sensor spatial averaging transfer function, eq. 7 ---------
C
	  TRANUV=1.D0
	  X=2.D0*PI*N/Mean(U)*CalSonic(QQPath)
	  IF (X.GT.0.02D0) TRANUV=(3.D0+EXP(-X)-4.D0*(1.D0-EXP(-X))/X)/X
C
C-- W-UV lateral separation transfer function, eq. 11 ------------
C
	  TRANWUV=1.D0
	  X = N/Mean(U)*CalSonic(QQExt3)
	  IF (X .GT.0.01D0) TRANWUV = EXP(-9.9D0*X**1.5D0)

C
C-- SOLENT SONIC VERTICAL WINDSPEED !!
C-- W-sensor frequency response gain, eq. 2 ----------------------
C
	  GAINW=1.D0	       !sonic W windspeed gain neglectible
C
C-- W-sensor spatial averaging transfer function, eq. 9 ----------
C
	  TRANW=1.D0
	  X=2.D0*PI*N/Mean(U)*CalSonic(QQPath)
	  IF (X.GT.0.04D0)
     &	    TRANW=(1.D0+(EXP(-X)-3.D0*(1.D0-EXP(-X))/X)/2.D0)*4.D0/X

C
C-- LYMANN-ALPHA OPEN PATH HYGROMETER !!
C-- hygrometer frequency response gain, eq. 2 ------------
C
	  GAINQ1=1.D0	       !hygrometer gain neglectible
C
C-- lymann-alpha spatial averaging transfer function, eq.7 ------------
C
	  TRANQ1=1.D0
	  X=2.D0*PI*N/Mean(U)*CalHyg(QQPath)
	  IF (X.GT.0.02D0) TRANQ1=(3.D0+EXP(-X)-4.D0*(1.D0-EXP(-X))/X)/X
C
C-- W-Q lateral separation transfer function, eq. 11 -------------
C
	  TRANWQ1=1.D0
	  X=N/Mean(U)*CalHyg(QQX)
	  IF (X.GT.0.01D0) TRANWQ1=EXP(-9.9D0*X**1.5D0)
C
C-- Composite transfer functions, eq. 28 -------------------------
C
	  DO I=1,15
	    G(I) = 0.D0
	  ENDDO

	  G(1)	= TRANUV / GAINUV		    !UU
	  G(2)	= TRANUV / GAINUV		    !VV
	  G(3)	= TRANW  / GAINW		    !WW
	  G(4)	= TRANT1 / GAINT1		    !TTth
	  G(5)	= TRANT2 / GAINT2		    !TTso
	  G(6)	= TRANQ1 / GAINQ1		    !QQkr
	  G(9)	= TRANWUV * SQRT (G(3) * G(1))	    !WU
	  G(10) = TRANWUV * SQRT (G(3) * G(2))	    !WV
	  G(11) = TRANWT1 * SQRT (G(3) * G(4))	    !WTth
	  G(12) = TRANWT2 * SQRT (G(3) * G(5))	    !WTso
	  G(13) = TRANWQ1 * SQRT (G(3) * G(6))	    !WQkr
	  DO I=1,15
	    G(I) = G(I)*TRANA/GAINV*GAIND
	  ENDDO

	  F=N/Mean(U)*CalSonic(QQZ)
	  IF (ZL.LT.0.D0) THEN
C
C-- Unstable normalised spectral and co-spectral forms ---------
C
	    UU=210.0D0*F/(1.D0+33.D0*F)**1.667D0		       !eq. 23
	    UU=UU+F*XI/(ZETA+2.2D0*F**1.667D0)
	    UU=UU/CU
	    WW=16.D0*F*XI/(1.D0+17.D0*F)**1.667D0		       !eq. 22
	    WW=WW+F/(1.D0+5.3D0*F**1.667D0)
	    WW=WW/CW
	    TT=14.94D0*F/(1.+24.D0*F)**1.667D0			     !eq. 24
	    IF (F.GE.0.15D0) TT=6.827D0*F/(1.D0+12.5D0*F)**1.667D0
	    WT=12.92D0*F/(1.D0+26.7D0*F)**1.375D0		       !eq. 25
	    IF (F.GE.0.54D0) WT=4.378D0*F/(1.D0+3.8D0*F)**2.4D0
	    UW=20.78D0*F/(1.D0+31.D0*F)**1.575D0		       !eq. 26
	    IF (F.GE.0.24D0) UW=12.66D0*F/(1.D0+9.6D0*F)**2.4D0
	  ELSE
C
C-- Stable normalised spectral and co-spectral forms, eq. 20a --
C
	    UU=F/(AUU+BUU*F**1.667D0)
	    WW=F/(AWW+BWW*F**1.667D0)
	    TT=F/(ATT+BTT*F**1.667D0)
	    WT=F/(AWT+BWT*F**2.1D0)
	    UW=F/(AUW+BUW*F**2.1D0)
	  ENDIF
C
C-- Integral of co-spectral transfer function * atm. co-spectrum -
C
	  INTTS(1)  = INTTS(1)	+ I4*G(1) *UU	   !UU
	  INTTS(2)  = INTTS(2)	+ I4*G(2) *UU	   !VV
	  INTTS(3)  = INTTS(3)	+ I4*G(3) *WW	   !WW
	  INTTS(4)  = INTTS(4)	+ I4*G(4) *TT	   !TTth
	  INTTS(5)  = INTTS(5)	+ I4*G(5) *TT	   !TTso
	  INTTS(6)  = INTTS(6)	+ I4*G(6) *TT	   !QQkr
	  INTTS(9)  = INTTS(9)	+ I4*G(9) *UW	   !WU
	  INTTS(10) = INTTS(10) + I4*G(10)*UW	   !WV
	  INTTS(11) = INTTS(11) + I4*G(11)*WT	   !WTth
	  INTTS(12) = INTTS(12) + I4*G(12)*WT	   !WTso
	  INTTS(13) = INTTS(13) + I4*G(13)*WT	   !WQkr

C
C-- Increase frequency by interval width -------------------------
C
	  N=N*LF
	ENDDO
      ENDDO

200   CONTINUE

C
C-- Final flux loss correction factors !! ----------------------------
C

      DO I=1,NSize
	DO J=1,NSize
	  WXT(I,J) = 1.0D0
	ENDDO
      ENDDO

999   CONTINUE
      WXT(U,U) = C/INTTS(1)		   !UU
      WXT(V,V) = C/INTTS(2)		   !VV
      WXT(U,V) = WXT(U,U)		   !UV
      WXT(V,U) = WXT(U,V)		   !VU
      WXT(W,W) = C/INTTS(3)		   !WW
      WXT(W,U) = C/INTTS(9)		   !WU
      WXT(U,W) = WXT(W,U)
      WXT(W,V) = C/INTTS(10)		   !WV
      WXT(V,W) = WXT(W,V)

      WXT(W,TSonic) = C/INTTS(12)	 !WTso
      WXT(TSonic,W) = WXT(W,TSonic)
      WXT(TSonic,TSonic) = C/INTTS(5)	 !TTso

      WXT(W,TCouple) = C/INTTS(11)	  !WTth
      WXT(TCouple,W) = WXT(W,TCouple)
      WXT(TCouple,TCouple) = C/INTTS(4)    !TTth

      WXT(W,Humidity) = C/INTTS(13)	  !WQly
      WXT(Humidity,W) = WXT(W,Humidity)
      WXT(Humidity,Humidity) = C/INTTS(6)  !QQly

      WXT(W,SpecHum)	   = WXT(W,Humidity)
      WXT(SpecHum,W)	   = WXT(Humidity,W)
      WXT(SpecHum,SpecHum) = WXT(Humidity,Humidity)

      DO i=1,NSize
	DO j=1,NSize
	  IF ((WXT(i,j).GE.Lower).AND.(WXT(i,j).LE.Upper)) THEN
	    Cov(i,j) = Cov(i,j)*WXT(i,j)
C
C The tolerance is augmented by half the amount added by frequency-
C response correction.
C
	    TolCov(i,j) = SQRT(TolCov(i,j)**2.D0+
     &	      (0.5D0*Cov(i,j)*(WXT(i,j)-1.D0))**2.D0)
	  ENDIF
	ENDDO
      ENDDO

      RETURN
      END





      SUBROUTINE ECOxygen(Mean,NMax,N,Cov,P,Factor, WhichTemp)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECOxygen
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Contribution of other gases especially oxygen absorb
C some of the radiation at the wavelengths at which the
C Krypton works.
C
C Revision 28-05-2001: added info on which temperature should be used
C                      in corrections (Sonic or thermocouple)
      INCLUDE 'physcnst.for'

      INTEGER NMax,N,i, WhichTemp
      REAL*8 Mean(NMax),Cov(NMax,NMax),Factor(NMax),OXC,P

      OXC = FracO2*MO2*P*KoK/(RGas*(Mean(WhichTemp))**2.D0*KwK)

      DO i = 1,N
	Factor(i) = 1.D0 + OXC*Cov(i,WhichTemp)/Cov(i,Humidity)
	IF ((i.EQ.Humidity).OR.(i.EQ.SpecHum))
     &	  Factor(i) = Factor(i)**2.D0
	Cov(i,Humidity) = Factor(i)*Cov(i,Humidity)
	Cov(i,SpecHum ) = Factor(i)*Cov(i,SpecHum )
	Cov(Humidity,i) = Cov(i,Humidity)
	Cov(SpecHum ,i) = Cov(i,SpecHum )
      ENDDO

      RETURN
      END





      SUBROUTINE ECSchot(Mean,NMax,N,Cov,Factor)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECSchot
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Partial Schotanus et al. correction: correction for humidity
C of sonic temperature, and of all covariances with sonic temperature.
C Sidewind-correction has already been applied in the
C routine where the sonic signal is calibrated.
C
      INCLUDE 'physcnst.for'

      INTEGER NMax,N,i
      REAL*8 Mean(NMax),Cov(NMax,NMax),Factor(NMax)

      DO i = 1,N
	Factor(i) = 1.D0 - 0.51D0*Mean(SpecHum)
     &	  -0.51D0*Mean(TSonic)*Cov(i,SpecHum)/Cov(i,TSonic)
	IF (i.EQ.TSonic) Factor(i) = Factor(i)**2.D0
	Cov(i,TSonic) = Factor(i)*Cov(i,TSonic)
	Cov(TSonic,i) = Cov(i,TSonic)
      ENDDO

      Mean(TSonic) = Mean(TSonic)/(1.D0+0.51D0*Mean(SpecHum))

      RETURN
      END

      REAL*8 FUNCTION ECRawSchot(Temp, Rhov, Press)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECRawSchot
C Version of subroutine : 1.0
C Date			: July 27 2000
C Author		: Arnold Moene
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
C To do humidity part of Schotanus et al. correction on a
C raw sample of sonic temperature, with Rhov that was not yet
C corrected
C
C Sidewind-correction has already been applied in the
C routine where the sonic signal is calibrated.
C
      INCLUDE 'physcnst.for'

      INTEGER NMax,N,i
      REAL*8   Temp, Rhov, Press, SPHUM, SPHUMNEW, NEWTEMP
      REAL*8 ECQ,ECBaseF ! External function calls

      SPHUM = ECQ(Temp, Rhov, Press)
      SPHUMNEW = 1.0
      NEWTEMP = Temp
C Dit convergeert dus helemaal niet !!!
      DO WHILE (ABS((SPHUM - SPHUMNEW)/SPHUM) .GT. 0.0001)

         NEWTEMP = TEMP/(1+0.51*SPHUM)
         SPHUMNEW = SPHUM
         SPHUM = ECQ(NEWTEMP, Rhov, Press)
      ENDDO
      ECRAWSCHOT = NEWTEMP
      RETURN
      END




      SUBROUTINE ECWebb(Mean,NMax,Cov,P, WhichTemp)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECWebb
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C Mean vertical velocity according to Webb, Pearman and Leuning
C
C Revision 28-05-2001: added info on which temperature should be used
C                      in corrections (Sonic or thermocouple)
      INCLUDE 'physcnst.for'

      INTEGER NMax, WhichTemp
      REAL*8 Mean(NMax),Cov(NMax,NMax),Sigma,ECRhoDry,P
C
C Ratio of mean densities of water vapour and dry air
C
      Sigma = Mean(Humidity)/ECRhoDry(Mean(Humidity),Mean(WhichTemp),P)
      Mean(W) = (1.D0+Mu*Sigma)*Cov(W,WhichTemp)/Mean(WhichTemp) +
     &	Mu*Sigma*Cov(W,Humidity)/Mean(Humidity)

      RETURN
      END





C
C ########################################################################
C
C Here comes a set of ready-to-use routines, performing corrections
C and printing results on demand. They are entirely based on routines
C described earlier in this library. No important physics or models...
C
C ########################################################################
C





       SUBROUTINE ECPTilt(DoPrint,OutF,Mean,TolMean,
     &	NMax,N,Cov,TolCov,QName,UName,PrePitch,PreYaw,PreRoll)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECPTilt
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C All-in procedure to tilt the system using KNOWN pitch, yaw and roll angles
C and print results
C
      INCLUDE 'physcnst.for'

      LOGICAL DoPrint
      INTEGER NMax,OutF,N
      REAL*8 Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),PrePitch,PreYaw,PreRoll
      CHARACTER*6 QName(NMax)
      CHARACTER*9 UName(NMax)

      CALL ECKPitch(Mean,NMax,N,Cov,PrePitch)
      CALL ECKYaw(  Mean,NMax,N,Cov,PreYaw  )
      CALL ECKRoll( Mean,NMax,N,Cov,PreRoll )

      IF (DoPrint) THEN
	WRITE(OutF,*)
	WRITE(OutF,*) 'Averages of data after fixed ',
     &	  'tilt-correction: '
	WRITE(OutF,*)
	CALL ECShow(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
      ENDIF

      RETURN
      END





      SUBROUTINE ECPPitch(DoPrint,OutF,Mean,TolMean,NMax,N,Cov,
     &	TolCov,QName,UName,Dirs)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECPPitch
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C All-in procedure to tilt the system such that Mean(V) --> 0
C and print results
C
      INCLUDE 'physcnst.for'

      LOGICAL DoPrint
      INTEGER NMax,OutF,N
      REAL*8 Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),Dirs
      CHARACTER*6 QName(NMax)
      CHARACTER*9 UName(NMax)

      CALL ECPitch(Mean,NMax,N,Cov,Dirs)

      IF (DoPrint) THEN
	WRITE(OutF,*)
	WRITE(OutF,*) 'Pitch angle = ',Dirs,' degrees'
	WRITE(OutF,*)
	WRITE(OutF,*) 'After pitch-correction (Mean(V) -> 0) : '
	WRITE(OutF,*)
	CALL ECShow(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
      ENDIF

      RETURN
      END





      SUBROUTINE ECPYaw(DoPrint,OutF,Mean,TolMean,NMax,N,Cov,
     &	TolCov,YawLim,QName,UName,Dirs)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECPYaw
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C All-in procedure to tilt the system such that Mean(W) --> 0
C and print results
C
      INCLUDE 'physcnst.for'

      LOGICAL DoPrint
      INTEGER NMax,OutF,N
      REAL*8 Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),Dirs,YawLim
      CHARACTER*6 QName(NMax)
      CHARACTER*9 UName(NMax)

      CALL ECYaw(Mean,NMax,N,Cov,YawLim,Dirs)

      IF (DoPrint) THEN
	WRITE(OutF,*)
	WRITE(OutF,*) 'Yaw angle = ',Dirs,' degrees'
	WRITE(OutF,*)
	WRITE(OutF,*) 'After yaw-correction (Mean(W) -> 0) : '
	WRITE(OutF,*)
	CALL ECShow(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
      ENDIF

      RETURN
      END





      SUBROUTINE ECPRoll(DoPrint,OutF,Mean,TolMean,NMax,N,Cov,
     &	TolCov,RollLim,QName,UName,Dirs)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECPRoll
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C All-in procedure to tilt the system such that Cov(W,V) --> 0
C and print results
C
      INCLUDE 'physcnst.for'

      LOGICAL DoPrint
      INTEGER NMax,OutF,N
      REAL*8 Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),Dirs,RollLim
      CHARACTER*6 QName(NMax)
      CHARACTER*9 UName(NMax)

      CALL ECRoll(Mean,NMax,N,Cov,RollLim,Dirs)

      IF (DoPrint) THEN
	WRITE(OutF,*)
	WRITE(OutF,*) 'Roll angle = ',Dirs,' degrees'
	WRITE(OutF,*)
	WRITE(OutF,*) 'After roll-correction (Cov(V,W) -> 0) : '
	WRITE(OutF,*)
	CALL ECShow(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
      ENDIF

      RETURN
      END





      SUBROUTINE ECPSchot(DoPrint,OutF,Mean,TolMean,NMax,N,Cov,
     &	TolCov,QName,UName,SonFactr)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECPSchot
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C All-in procedure to correct for water sensitivity of sonic temperature
C and print results
C
      INCLUDE 'physcnst.for'

      LOGICAL DoPrint
      INTEGER NMax,OutF,N,i
      REAL*8 Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),SonFactr(NMax)
      CHARACTER*6 QName(NMax)
      CHARACTER*9 UName(NMax)

      CALL ECSchot(Mean,NMax,N,Cov,SonFactr)

      IF (DoPrint) THEN
	WRITE(OutF,*)
	WRITE(OutF,*) 'Correction factors for H2O-sensitivity of ',
     &	  'sonic'
	DO i=1,N
	  WRITE(OutF,45) QName(i),SonFactr(i)
	ENDDO
 45	FORMAT('For (',a6,',T(son)) : ',F10.5)
	WRITE(OutF,*)
	WRITE(OutF,*) 'After H2O-correction for sonic : '
	WRITE(OutF,*)
	CALL ECShow(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
      ENDIF

      RETURN
      END





      SUBROUTINE ECPO2(DoPrint,OutF,Mean,TolMean,NMax,N,Cov,
     &	TolCov,QName,UName,P,O2Factor, WhichTemp)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECPO2
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C All-in procedure to correct for oxygen sensitivity of hygrometer
C and print results
C
C Revision 28-05-2001: added info on which temperature should be used
C                      in corrections (Sonic or thermocouple)
      INCLUDE 'physcnst.for'

      LOGICAL DoPrint
      INTEGER NMax,OutF,N,i, WhichTemp
      REAL*8 Mean(NMax),TolMean(NMax),Cov(NMax,NMax),
     &	TolCov(NMax,NMax),O2Factor(NMax),P
      CHARACTER*6 QName(NMax)
      CHARACTER*9 UName(NMax)

      CALL ECOxygen(Mean,NMax,N,Cov,P,O2Factor, WhichTemp)

      IF (DoPrint) THEN
	WRITE(OutF,*)
	WRITE(OutF,*) 'Correction factors for O2-sensitivity of ',
     &	  'hygrometer'
	DO i=1,N
	  WRITE(OutF,46) QName(i),O2Factor(i)
	ENDDO
 46	FORMAT('For (',a6,',RhoV or q) : ',F10.5)
	WRITE(OutF,*)
	WRITE(OutF,*) 'After oxygen-correction for hygrometer : '
	WRITE(OutF,*)
	CALL ECShow(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
      ENDIF

      RETURN
      END




      SUBROUTINE ECPFreq(DoPrint,OutF,Mean,TolMean,NMax,N,
     &	Cov,TolCov,QName,UName,LLimit,ULimit,NSta,NEnd,NInt,Freq,
     &	TauD,TauV,CalSonic,CalTherm,CalHyg,FrCor, WhichTemp)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECPFreq
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C All-in procedure to perform frequency response correction and print results
C
C input : DoPrint [LOGICAL] : Indicator if intermediate results must be
C           printed to file OutF.
C         OutF [INTEGER] : Unit number of output-file for intermediate
C           results. OutF must be an open file.
C         Mean [REAL*8(NMax)] : Array of mean values of the quantities in
C           this experiment (only the first N quantities are used).
C         TolMean [REAL*8(NMax)] : Tolerances of Mean.
C         NMax [INTEGER] : Physical dimension of array Mean
C         N [INTEGER] : Number of quantities actually involved in this
C           experiment.
C         Cov [REAL*8(NMax,NMax)] : covariances of the fluctuations.
C         TolCov [REAL*8(NMax,NMax)] : Tolerances of covariances (2 sigma).
C         QName [CHARACTER*9(NMax)] : Names of the quantities.
C         UName [CHARACTER*9(NMax)] : Names of the units of the quantities.
C         LLimit [REAL*8] : Lower acceptance limit for frequency-response
C           factors. Correction factors smaller than LLimit are set to 1.
C         ULimit [REAL*8] : Upper acceptance limit for frequency-response
C           factors. Correction factors larger than ULimit are set to 1.
C         NSta [REAL*8] : Start frequency numerical integration.
C           Popular value: -5.D0 [unit?].
C         NEND [REAL*8] : End frequency numerical integration.
C           Popular value: LOG(5) = 0.69897D0 [unit?].
C         NINT [INTEGER] : Number of intervals in integration.
C           Popular value: 19.
C         Freq [REAL*8] : Sampling frequency [Hz].
C         TauD [REAL*8] : Interval length for running mean.
C           Popular value: 0.D0 [unit?].
C         TAUV [REAL*8] : Low pass filter time constant.
C           Popular value: 0.D0 [unit?]
C         CalSonic [REAL*8(NQQ)] : Calibration specification array of
C           sonic anemometer.
C         CalTherm [REAL*8(NQQ)] : Calibration specification array of
C           thermometer.
C         CalHyg [REAL*8(NQQ)] : Calibration specification array of
C           hygrometer.
C
C output : FrCor [REAL*8(NMax,NMax)] : Correction factors for covariances.
C
C Revision 28-05-2001: added info on which temperature should be used
C                      in corrections (Sonic or thermocouple)
      INCLUDE 'physcnst.for'

      LOGICAL DoPrint
      INTEGER NMax,OutF,N,NInt, WhichTemp
      REAL*8 Mean(NMax),TolMean(NMax),Cov(NMax,NMax),TolCov(NMax,NMax)
      CHARACTER*6 QName(NMax)
      CHARACTER*9 UName(NMax)
      REAL*8 LLimit,ULimit,FrCor(NMax,NMax),NSTA,NEND,TAUV,TauD,Freq
      REAL*8 CalSonic(NQQ),CalTherm(NQQ),CalHyg(NQQ)

      CALL ECFrResp(Mean,Cov,TolCov,NMax,N,
     &	LLimit,ULimit,NSta,NEnd,NInt,Freq,TauD,TauV,
     &	CalSonic,CalTherm,CalHyg,FrCor, WhichTemp)

      IF (DoPrint) THEN
	CALL ECShwFrq(OutF,QName,FrCor,NMax,N)
	WRITE(OutF,*)
	WRITE(OutF,*) 'After frequency-correction : '
	WRITE(OutF,*) 'Factors accepted between ',LLimit,
     &	  ' and ',ULimit
	WRITE(OutF,*)
	CALL ECShow(OutF,QName,UName,Mean,TolMean,Cov,TolCov,NMax,N)
      ENDIF

      RETURN
      END





      SUBROUTINE ECPWebb(DoPrint,OutF,Mean,TolMean,NMax,N,Cov,TolCov,
     &	P,QName,UName, WhichTemp)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECPWebb
C Version of subroutine : 1.02
C Date			: 27 August 1999
C Author		: Arjan van Dijk
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
C All-in procedure to calculate Webb-velocity and print results
C
C input : DoPrint [LOGICAL] : Indicator if intermediate results must be
C           printed to file OutF.
C         OutF [INTEGER] : Unit number of output-file for intermediate
C           results. OutF must be an open file.
C         Mean [REAL*8(NMax)] : Array of mean values of the quantities in
C           this experiment (only the first N quantities are used).
C         TolMean [REAL*8(NMax)] : Tolerances of Mean.
C         NMax [INTEGER] : Physical dimension of array Mean
C         N [INTEGER] : Number of quantities actually involved in this
C           experiment.
C         Cov [REAL*8(NMax,NMax)] : covariances of the fluctuations.
C         TolCov [REAL*8(NMax,NMax)] : Tolerances of covariances (2 sigma).
C         P [REAL*8] : Ambient air pressure in Pascal.
C         QName [CHARACTER*9(NMax)] : Names of the quantities.
C         UName [CHARACTER*9(NMax)] : Names of the units of the quantities.
C output : This routine sets the mean vertical velocity "Mean(w)" to the
C         Webb-velocity.
C
C Revision 28-05-2001: added info on which temperature should be used
C                      in corrections (Sonic or thermocouple)
      INCLUDE 'physcnst.for'

      LOGICAL DoPrint
      INTEGER NMax,OutF,N, WhichTemp
      REAL*8 Mean(NMax),TolMean(NMax),Cov(NMax,NMax),P,
     &	TolCov(NMax,NMax)
      CHARACTER*6 QName(NMax)
      CHARACTER*9 UName(NMax)

      CALL ECWebb(Mean,NMax,Cov,P, WhichTemp)

      IF (DoPrint) THEN
	WRITE(OutF,*)
	WRITE(OutF,*) 'Webb-velocity (vertical) = ',Mean(W),' m/s'
	WRITE(OutF,*)
	WRITE(OutF,*) 'After addition of Webb-term : '
	WRITE(OutF,*)
	CALL ECShow(OutF,QName,UName,Mean,TolMean,Cov,TolCov,
     &              NMax,N)
      ENDIF

      RETURN
      END





C
C
C Routines to calculate structure parameters
C
C





      SUBROUTINE ECStruct(Sample,NMax,N,MMax,M,Flag,XIndex,YIndex,
     &  R,dR,Freq,CIndep,Cxy,dCxy)
C
C EC-Pack Library for processing of Eddy-Correlation data
C Subroutine		: ECStruct
C Version of subroutine : 1.02
C Date			: 26 August 1999
C Author		: Arjan van Dijk
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
C Calculate structure parameters <(x(r)-x(r+R))*(y(r)-y(r+R))>/R^2/3
C
C input : Sample [REAL*8(NMax,MMax)] : The first N out of NMax quantities
C           and the first M samples in Sample are used.
C         NMax [INTEGER] : physical first dimension of array Sample.
C         N [INTEGER] : actual number of meaningful quantities in array
C           Sample.
C         MMax [INTEGER] : physical second dimension of array Sample.
C         M [INTEGER] : actual number of meaningful samples in array Sample.
C         Flag [LOGICAL(NMax,MMax)] : if Flag(i,j) is true, then something
C           is wrong with quantity i in sample j.
C         XIndex [INTEGER] : indicator of first quantity involved in
C           structure function.
C         YIndex [INTEGER] : indicator of second quantity involved in
C           structure function.
C         R [REAL*8] : separation in meters at which one wants to estimate
C           the structure function.
C         Freq [REAL*8] : Sampling frequency in s^-1.
C         CIndep [INTEGER[NMax,NMax]) : Number of independent contributions
C           by array Sample to covariance between quantities selected with
C           XIndex and YIndex.
C output: dR [REAL*8] : separation in meters corresponding with a delay
C           of one sample.
C         cxy [REAL*8] : Structure parameter.
C         dcxy [REAL*8] : Tolerance of cxy.
C
      INCLUDE 'physcnst.for'
      INTEGER NMax,N,i,M,MMax,XIndex,YIndex,NOk,NSeparate,
     &  CIndep(NMax,NMax)
      LOGICAL Flag(NMax,MMax),ok
      REAL*8 R,Sample(NMax,MMax),UMean,Freq,TwoThird,Cxy,dCxy,Dum,
     &  Increment,dR

      TwoThird = 2.D0/3.D0
C
C Calculate the average of the length of the velocity (and not
C the length of the average velocity!!!)
C
      UMean = 0.D0
      NOk = 0
      DO i=1,M
        IF ( (.NOT.(Flag(U,i))).AND.
     &      ((.NOT.(Flag(V,i))).AND.
     &       (.NOT.(Flag(W,i))))) THEN
          NOk = NOk + 1
          UMean = UMean + Sample(U,i)**2+Sample(V,i)**2+Sample(W,i)**2
        ENDIF
      ENDDO
      IF (NOk.GT.0) UMean = UMean/DBLE(NOk)
      UMean = SQRT(UMean)
C
C Estimate how many samples delay one must go to let Taylor's hypothesis
C of frozen turbulence give the correct spatial separation R.
C The discrete nature of sampling may call for strong rounding off of the
C delay distance. The rounded off value for R is returned to the calling
C rourine.
C
      dR = UMean/Freq
      NSeparate = NINT(R/dR)
      R = R*NSeparate/(R/dR)
C
C Calculate structure parameter
C
      Cxy = 0.D0
      NOk = 0
      DO i=1,(M-NSeparate)
        ok = (((.NOT.Flag(XIndex, i           )).AND.
     &         (.NOT.Flag(XIndex,(i+NSeparate))))
     &        .AND.
     &        ((.NOT.Flag(YIndex, i           )).AND.
     &         (.NOT.Flag(YIndex,(i+NSeparate)))))
        IF (ok) THEN
          NOk = NOk + 1
          Cxy = Cxy +
     &      (Sample(XIndex,i)-Sample(XIndex,(i+NSeparate)))*
     &      (Sample(YIndex,i)-Sample(YIndex,(i+NSeparate)))
        ENDIF
      ENDDO
      IF (NOk.GT.0) Cxy = Cxy/DBLE(Nok)
      Dum = Cxy ! For use in tolerance estimation loop
      IF (NSeparate.GT.0) Cxy = Cxy/R**TwoThird
C
C Estimate tolerance of structure parameter
C
      dCxy = 0.D0
      IF (NSeparate.GT.0) THEN
        DO i=1,(M-NSeparate)
          ok = (((.NOT.Flag(XIndex, i           )).AND.
     &           (.NOT.Flag(XIndex,(i+NSeparate))))
     &          .AND.
     &          ((.NOT.Flag(YIndex, i           )).AND.
     &           (.NOT.Flag(YIndex,(i+NSeparate)))))
          IF (ok) THEN
            Increment =
     &        (Sample(XIndex,i)-Sample(XIndex,(i+NSeparate)))*
     &        (Sample(YIndex,i)-Sample(YIndex,(i+NSeparate)))
            dCxy = dCxy + (Increment - Dum)**2.D0
          ENDIF
        ENDDO
        IF (NOk.GT.0) dCxy = dCxy/DBLE(Nok)
        dCxy = SQRT(dCxy)/R**TwoThird ! Standard deviation
        dCxy = 2.D0*dCxy/SQRT(DBLE(CIndep(XIndex,YIndex))) ! Tolerance
      ENDIF
C
C We use the number of independent samples found earlier in the estimation
C of the covariances
C
      RETURN
      END

C...........................................................................
C Function  : ECSCal
C Purpose   : to calibrate a sonic signal according to wind tunnel
C             calibration (for th moment this works for apparatus 6,
C             i.e. a wind tunnel calibrated sonic)
C Interface : Cal     IN        array of length NQQ with calibation info
C             UDum    IN/OUT    one horizontal component (on exit: calibrated)
C             VDum    IN/OUT    another horizontal component (on exit: calibrated)
C             WDum    IN/OUT    vertical component (on exit: calibrated)
C             UERROR  IN/OUT    error flag for U (.TRUE. if wrong data)
C             VERROR  IN/OUT    error flag for V (.TRUE. if wrong data)
C             WERROR  IN/OUT    error flag for W (.TRUE. if wrong data)
C Author    : Arnold Moene
C Date      : September 26, 2000
C Remarks   : The method is based on a wind tunnel calibration of the sonic
C             The real velocity components can be derived from the
C             measured components and the real azimuth and elevation angle.
C             But the latter are not known and have to be determined
C             iteratively from the measured components. The relationship
C             between the real components and the measured components is:
C
C               Ureal =  Umeas/(UC1*(1 - 0.5*
C                                    ((Azi + (Elev/0.5236)*UC2)*
C                                     (1 - UC3*Abs(Elev/0.5236)))**2 ))
C               Vreal =  Vmeas*(1 - VC1*Abs(Elev/0.5236))
C               Wreal =  Wmeas/(WC1*(1 - 0.5*(Azi*WC2)**2))
C
C             and
C               Azi = arctan(V/U)
C               Elev = arctan(W/sqrt(U**2 + V**2))
C
C             where UC1, UC2, UC3, VC1, WC1, WC2 are fitting coefficients.
C             An azimuth angle of zero is supposed to refer to a wind
C             direction from the most optimal direction (i.e. the 'open'
C             side of a sonic). Samples with an absolute azimuth angle of
C             more than 40 degrees are rejected.
C...........................................................................
      SUBROUTINE ECScal(Cal, UDum, VDum, WDum,
     &                  UError, VError, WError)
     
      INCLUDE 'parcnst.inc'

      REAL*8   Cal(NQQ), UDum, VDum, WDum
      LOGICAL  UError, VError, WError

      REAL*8   UCorr, VCorr, WCorr, AziNew, AziOld, ElevNew, ElevOld,
     &         UC1, UC2, UC3, VC1, WC1, WC2
      INTEGER  I, NITER, ITMAX


      UC1 = Cal(QQExt6)
      UC2 = Cal(QQExt7)
      UC3 = Cal(QQExt8)
      VC1 = Cal(QQExt9)
      WC1 = Cal(QQExt10)
      WC2 = Cal(QQExt11)

      ITMAX = 20
      NITER = 0
      AziOld = 9999D0
      ElevOld = 9999D0
      AziNew = ATAN(VDum/UDum)
      ElevNew = ATAN(WDum/SQRT(UDum**2 + VDum**2))

      UCorr = UDum
      VCorr = VDum
      WCorr = WDum

      IF (ABS(AziNew) .LE. 0.698) THEN
         DO WHILE (((ABS(AziNew-AziOld) .GT. 1./60) .OR.
     &              (ABS(ElevNew-ElevOld) .GT. 1./60)) .AND.
     &           (NITER .LT. ITMAX))
            UCorr =  UDum/(UC1*(1 - 0.5*
     &              ((AziNew + (ElevNew/0.5236D0)*UC2)*
     &              (1 - UC3*Abs(ElevNew/0.5236D0))
     &              )**2        )
     &                     )
            VCorr =  VDum*(1 - VC1*Abs(ElevNew/0.5236D0))
            WCorr =  WDum/(WC1*(1 - 0.5*(ElevNew*WC2)**2))
   
            AziOld = AziNew
            ElevOld = ElevNew
   
            AziNew = ATAN(Vcorr/UCorr)
            ElevNew = ATAN(Wcorr/SQRT(UCorr**2 + VCorr**2))
   
            NITER = NITER + 1
         ENDDO
      ENDIF
      
      IF ((NITER .EQ. ITMAX) .OR. (ABS(AziNew) .GT. 0.698)) THEN
         UError = .TRUE.
         VError = .TRUE.
         WError = .TRUE.
      ELSE
         UDum = UCorr
         VDum = VCorr
         WDum = WCorr
      ENDIF

      END
