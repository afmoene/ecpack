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
     &  Channels,Delay(NNNMax),Mok(NNMax),Cok(NNMax,NNMax),
     &  ii, jj

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
      REAL*8 UMean(3,MaxRuns),Apf(3,3),Alpha,Beta,Gamma,WBias,Dum(3),
     &  Dumout(3),
     &  UUSum(3,3),RApf(3,3),RAlpha,RBeta,RGamma,RWBias,USum(3),
     &  Swing,SinSwing,CosSwing,UHor, RunStart(3, MaxRuns),
     &  RunStop(3, MaxRuns), 
     &  MEAN1(3), STRESS1(3,3),
     &  MEAN2(3), STRESS2(3,3), Yaw1, Roll1, Pitch1,
     &  dYaw1, dRoll1, dPitch1, Speed(3), DumCov(3,3),
     &  GMEAN(NNMAX), GCOV(NNMax,NNMax),
     &  GTolMEAN(NNMAX), GTolCOV(NNMax,NNMax),
     &  Yaw(3,3), InvYaw(3,3), SinGamma, CosGamma, MAP(3,3)


      LOGICAL       HAVE_UNCAL(NNNMax),Flag(NNMAx,MMMax),
     &        DoWBias, SingleRun
      LOGICAL ENDINTERVAL, BadTc, DoPrint, HAVE_SAMP, DoPf, DoClassic
      INTEGER EC_T_STRLEN
C
C The name of the calibration routine specified at the end of this file
C
      EXTERNAL Calibrat, EC_T_STRLEN


C...........................................................................
C Start of executable statements
C...........................................................................
C Set DoPrint to true
      DoPrint = .TRUE.
      DoPF = .TRUE.
      DoClassic = .FALSE.
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
      
      CALL EC_F_GetConf(ECConfFile,
     &             DatDir, OutDir, ParmDir,
     &             FluxName, ParmName, InterName,
     &             PlfName,
     &             SonName, CoupName, HygName,
     &             NCVarname, NNNMax)
C
C Check which calibration files we have
C
      IF (EC_T_STRLEN(SonName) .GT. 0) THEN
         HAVE_SONCAL = .TRUE.
      ELSE
         HAVE_SONCAL = .FALSE.
      ENDIF
      IF (EC_T_STRLEN(HygName) .GT. 0) THEN
         HAVE_HYGCAL = .TRUE.
      ELSE
         HAVE_HYGCAL = .FALSE.
      ENDIF
      IF (EC_T_STRLEN(CoupName) .GT. 0) THEN
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
     &     CALL EC_F_ReadAp(SonName,CalSonic) ! The sonic
      IF (HAVE_TCCAL)
     &     CALL EC_F_ReadAp(CoupName,CalTherm) ! A thermo-couple
      IF (HAVE_HYGCAL)
     &     CALL EC_F_ReadAp(HygName,CalHyg) ! A Krypton hygrometer
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
         CALL EC_T_STRCAT(DUMSTRING, "/")
         CALL EC_T_STRCAT(DUMSTRING, PlfIntName)
         PlfIntName = DUMSTRING
         OPEN(PlfIntFile,FILE=PlfIntName,STATUS='REPLACE')
         CLOSE(PlfIntFile)
      ENDDO
 900  CONTINUE
      CLOSE(IntervalFile)

C Open file where results of Planar fit method are written
      OPEN(PlfFile, FILE=PlfName, status='Replace',
     &     FORM='FORMATTED')
      WRITE(PlfFile,62) 'Doy', 'Hourmin', 'Sec', 'Doy', 'Hourmin',
     &                  'Sec', '#sample', '#u', '#v', '#w',
     &                  'Alpha',  'Beta', 'Gamma', 'WBias',
     &                  'Alpha1',  'Beta1', 'Gamma1', 'WBias1',
     &                  'Swing', 'meanu1',  'meanv1', 'meanw1',
     &                  'uu1', 'uv1', 'uw1', 'vu1', 'vv1', 'vw1',
     &                  'wu1','wv1','ww1',
     &                  'yaw','pitch','roll',
     &                  'meanu2',  'meanv2', 'meanw2',
     &                  'uu2', 'uv2', 'uw2', 'vu2', 'vv2', 'vw2',
     &                  'wu2', 'wv2', 'ww2', 'dyaw', 'dpitch', 'droll'


      if (.not. DOPF) GOTO 4000
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
         CALL EC_T_STRCAT(DUMSTRING, "/")
         CALL EC_T_STRCAT(DUMSTRING, Fname)
         Fname = DUMSTRING

         DUMSTRING =  OutDir
         CALL EC_T_STRCAT(DUMSTRING, "/")
         CALL EC_T_STRCAT(DUMSTRING, OutName)
         OutName = DUMSTRING
         OPEN(OutFile, FILE=OutName, Status='UNKNOWN')

         DUMSTRING =  OutDir
         CALL EC_T_STRCAT(DUMSTRING, "/")
         CALL EC_T_STRCAT(DUMSTRING, PlfIntName)
         PlfIntName = DUMSTRING
C
C Show which file is currently being analysed
C
         WRITE(*,951) FName(:EC_T_STRLEN(Fname)),(NINT(StartTime(i)),i=1,3)
 951     FORMAT(a,1X,3(I5,1X))
C
C Read raw data from file
C
         CALL EC_F_ReadNCDF(FName,StartTime,StopTime,
     &     Delay,Gain,Offset,RawSampl,NNNMax,Channels,MMMax,M,
     &     NCVarName, HAVE_UNCAL)
C 
C Do calibration (per sample)
C 
         N = 3
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
 100       FORMAT(2(F8.0,1X,F8.0,1X,F8.2),1X,3(F20.15))
           CLOSE(PlfIntFile)
	 ENDIF
      ENDDO
 1000 CONTINUE
      CLOSE(IntervalFile)

C
C Now start over again, and get the mean velocities from the intermediate files
C
 4000 CONTINUE
      OPEN(IntervalFile,FILE = InterName,IOSTAT=FOO,STATUS='OLD')
      IF (FOO.NE.0) THEN
        WRITE(*,*) 'Could not open interval-file ',InterName
        STOP'- Fatal: no time info could be read -'
      ENDIF
C Read one line of comments (blindly, it is supposed to be a comment line)
      READ(IntervalFile,*)
C End of skip
      ENDINTERVAL = .FALSE.

C
      DO WHILE (.NOT. ENDINTERVAL)
         READ(IntervalFile,*,END=1100) (StartTime(i), i=1,3),
     &        (StopTime(i),i=1,3), Psychro, p, Fname, OutName,
     &        PlfIntName
         write(*,*) (starttime(i),i=1,3)
         DUMSTRING =  DatDir
         CALL EC_T_STRCAT(DUMSTRING, "/")
         CALL EC_T_STRCAT(DUMSTRING, Fname)
         Fname = DUMSTRING

         DUMSTRING =  OutDir
         CALL EC_T_STRCAT(DUMSTRING, "/")
         CALL EC_T_STRCAT(DUMSTRING, OutName)
         OutName = DUMSTRING

         DUMSTRING =  OutDir
         CALL EC_T_STRCAT(DUMSTRING, "/")
         CALL EC_T_STRCAT(DUMSTRING, PlfIntName)
         PlfIntName = DUMSTRING
C
         IF (DOPF) THEN
            OPEN(PlfIntFile,FILE=PlfIntName,STATUS='OLD')
            i=0
            DO WHILE (.TRUE.)
               IF (i+1 .GT. MaxRuns) THEN
                  WRITE(*,*) 'Maximum number of runs exceeded: ',
     &                       MAXRUNS
                  STOP
               ENDIF
               READ(PlfIntFile,*,END=1200) 
     &                     (RunStart(k,i+1), k=1,3),
     &                     (RunStop(k,i+1), k=1,3),
     &                     (Umean(k,i+1), k=1,3)
               i=i+1
            ENDDO
 1200       CONTINUE
            CLOSE(PlfIntFile)
            Nruns = i
	    SingleRun = (.FALSE.)
	    DoWBias = (.TRUE.)
	    DoWBias = (.FALSE.)
            CALL EC_C_T01(DoWBias, SingleRun, uMean,Nruns,
     &                  Apf,Alpha,Beta,Gamma,WBias)
C Turn back the Gamma rotation
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

            IF (DoPrint) THEN
              WRITE(OutFile,*)
              WRITE(OutFile, 10) Alpha
 10           FORMAT('Alpha = ',F7.2,' degrees')
              WRITE(OutFile,20) Beta
 20           FORMAT('Beta  = ',F7.2,' degrees')
              WRITE(OutFile,30) Gamma
 30           FORMAT('Gamma = ',F7.2,' degrees')
              WRITE(OutFile, 40) WBias
 40           FORMAT('WBias = ',F7.4,' m/s')
              WRITE(OutFile, *) 'Matrix = '
              WRITE(OutFile,50) Apf(1,1),Apf(1,2),Apf(1,3)
              WRITE(OutFile,50) Apf(2,1),Apf(2,2),Apf(2,3)
              WRITE(OutFile,50) Apf(3,1),Apf(3,2),Apf(3,3)
              WRITE(OutFile,*)
 50           FORMAT(3(F7.3,1X))
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


C ENDIF of IF (DOPF)
         ENDIF
         SingleRun = (.TRUE.)

C
C Read raw data from file
C
         CALL EC_F_ReadNCDF(FName,StartTime,StopTime,
     &     Delay,Gain,Offset,RawSampl,NNNMax,Channels,MMMax,M,
     &     NCVarName, HAVE_UNCAL)
C 
C Do calibration (per sample)
C
         N = 3
         DO i=1,M
            CALL Calibrat(RawSampl(1,i),Channels,P,CorMean,
     &	      CalSonic,CalTherm,CalHyg,BadTc,Sample(1,i),N,Flag(1,i),
     &        Have_Uncal)
         ENDDO
	 CALL EC_M_Averag(Sample, NNMax, N, MMMax, M, Flag,
     &                 GMean, GTolMean, GCov, GTolCov, MIndep, CIndep,
     &                 Mok, Cok)
         Runlength = M
         IF (DOPF) THEN
C
C To which numbers should N and M be set ??? M has been set already above
            N = 3
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

	       CALL EC_M_MapVec(Apf,Dum,Dumout)

C Feed planarfit rotated sample back into Sample array
               do k=1,3
	          Sample(k,j) = Dumout(k)
	       enddo

	       DO k=1,3
                  IF (.NOT. Flag(k,j)) THEN
                     Mok(k) = Mok(k) + 1
	             USum(k) = USum(k) + Dumout(k)
                  ENDIF
	          DO l=1,3
                     IF ((.NOT. Flag(k,j)) .AND.
     &                   (.NOT. Flag(l,j))) THEN
	                UUSum(k,l) = UUSum(k,l) + Dumout(k)*Dumout(l)
                        Cok(k,l) = Cok(k,l) + 1
                     ENDIF
	          ENDDO
	       ENDDO
	    ENDDO
	    DO k=1,3
 	       USum(k) = USum(k)/Mok(k)
	       DO l=1,3
 	         UUSum(k,l) = UUSum(k,l)/Cok(k,l)
	       ENDDO
	    ENDDO
	    CALL EC_M_Averag(Sample, NNMax, N, MMMax, M, Flag,
     &                 GMean, GTolMean, GCov, GTolCov, MIndep, CIndep,
     &                 Mok, Cok)
C
C Compute the residual angles per run
C
            IF (RUNLENGTH .GT. 0) THEN
               CALL EC_C_T02(DoWBias, SingleRun,uSum,UUSum,RApf,
     &                RAlpha,
     &                RBeta,RGamma,RWBias)
            ENDIF


            UHor = (RAlpha**2+RBeta**2)**0.5
            CosSwing = RBeta/UHor
            SinSwing = RAlpha/UHor
            Swing = 180.D0*ACOS(CosSwing)/PI
            IF (SinSwing.LT.0.D0) Swing = 360.D0-Swing
            IF (Swing.GT.180.D0) Swing = Swing-360.D0
	    if (Runlength .EQ. 0) THEN
	        RALPHA = DUMMY
	        RBETA = DUMMY
	        RGAMMA = DUMMY
	        RWBIAS = DUMMY
	        SWING = DUMMY
	        DO I=1,3
	   	   MEAN1(I) = DUMMY
	           DO J=1,3
		      STRESS1(I,J) = DUMMY
		   ENDDO
                ENDDO
	    ELSE
	        DO I=1,3
C	   	   MEAN1(I) = USUM(I)
	   	   MEAN1(I) = GMEAN(I)
	           DO J=1,3
C		      STRESS1(I,J) = UUSUM(I,J)-USUM(I)*USUM(J)
		      STRESS1(I,J) = GCOV(I,J)
		   ENDDO
                ENDDO
	    ENDIF
        ELSE
C END OF DOPF
	    RALPHA = DUMMY
	    RBETA = DUMMY
	    RGAMMA = DUMMY
	    RWBIAS = DUMMY
	    SWING = DUMMY
	    DO I=1,3
	       MEAN1(I) = DUMMY
	       DO J=1,3
	          STRESS1(I,J) = DUMMY
	       ENDDO
            ENDDO
        ENDIF
C
C Compute the rotations using classic rotations
C
	 IF (RUNLENGTH .GT. 0) THEN
C Determine means and covariances after planarfit,
C tolerance in Yaw angle, compute yaw angle and
C apply it to raw data
	    CALL EC_M_Averag(Sample, NNMax, N, MMMax, M, Flag,
     &                 GMean, GTolMean, GCov, GTolCov, MIndep, CIndep,
     &                 Mok, Cok)
	    dYaw1 = (1/(1+(GMean(1)/GMean(2))**2))*
     &              abs(GMean(1)/GMean(2))*
     &              (abs(GTolMean(1)/GMean(1)) + 
     &               abs(GTolMean(2)/GMean(2)))
            DO i=1,3
	       Mean2(I) = Gmean(I)
	       DO j=1,3
	          Stress2(i,j) = GCov(i,j)
	       ENDDO
	    ENDDO
	    CALL EC_C_T07(Mean2(U), MEAN2(V), Yaw1)
	    CALL EC_C_T06(Yaw1, MAP)
	    CALL EC_C_T05(Mean2, 3, 3, STRESS2, MAP)
c	    CALL ECYaw(MEAN2, 3, 3, STRESS2, Yaw1)
	    DO I=1,RunLength
	       DO J=1,3
	         Speed(J)=Sample(j,i)
	       ENDDO
	       CALL EC_C_T05(Speed,3,3, DumCov, MAP)
	       DO J=1,3
	         Sample(j,i)=Speed(J)
	       ENDDO
	    ENDDO
C Recompute means, covariances and tolerances
C Determine tolerance in Pitch angle, compute Pitch angle and
C apply it to raw data
	    CALL EC_M_Averag(Sample, NNMax, N, MMMax, M, Flag,
     &                 Gmean, GTolMean, GCov, GTolCov, MIndep, CIndep,
     &                 Mok, Cok)
            IF (DoClassic) THEN
	    dPitch1 = (1/(1+(GMean(3)/GMean(1))**2))*
     &              abs(GMean(3)/GMean(1))*
     &              (abs(GTolMean(1)/GMean(1)) + 
     &               abs(GTolMean(3)/GMean(3)))
            DO i=1,3
	       Mean2(I) = Gmean(I)
	       DO j=1,3
	          Stress2(i,j) = GCov(i,j)
	       ENDDO
	    ENDDO
	    CALL EC_C_T09(Mean2(U), MEAN2(W), Pitch1)
	    CALL EC_C_T08(Pitch1, MAP)
	    CALL EC_C_T05(Mean2, 3, 3, STRESS2, MAP)
C	    CALL ECPitch(MEAN2, 3, 3, STRESS2, 30.0D0, Pitch1)
	    DO I=1,RunLength
	       DO J=1,3
	         Speed(J)=Sample(j,i)
	       ENDDO
	       CALL EC_C_T05(Speed,3,3, DumCov, MAP)
	       DO J=1,3
	         Sample(j,i)=Speed(J)
	       ENDDO
	    ENDDO
C Recompute means and coveriances and their tolerances
C Determine tolerance in Roll angle, compute Roll angle and
C apply it to raw data
	    CALL EC_M_Averag(Sample, NNMax, N, MMMax, M, Flag,
     &                 GMean, GTolMean, GCov, GTolCov, MIndep, CIndep,
     &                 Mok, Cok)
            dRoll1 = (1/(1+(2*GCov(2,3)/
     &                (GCov(2,2)-GCov(3,3)))**2)) *
     &               abs(2*GCov(2,3)/(GCov(2,2)-GCov(3,3)))*
     &               (abs(GTolCov(2,3)/GCov(2,3)) +
     &                abs((GTolCov(2,2)+GTolCov(3,3))/
     &                    (GCov(2,2)-GCov(3,3))))
            DO i=1,3
	       Mean2(I) = Gmean(I)
	       DO j=1,3
	          Stress2(i,j) = GCov(i,j)
	       ENDDO
	    ENDDO
	    CALL EC_C_T11(STRESS2(V,V),STRESS2(V,W), STRESS2(W,W), Roll1)
	    CALL EC_C_T10(Roll1, MAP)
	    CALL EC_C_T05(Mean2, 3, 3, STRESS2, MAP)
C	    CALL ECRoll(MEAN2, 3, 3, STRESS2, 30.0D0, Roll1)
	    DO I=1,RunLength
	       DO J=1,3
	         Speed(J)=Sample(j,i)
	       ENDDO
	       CALL EC_C_T05(Speed, 3, 3, DUMCOV, MAP)
	       DO J=1,3
	         Sample(j,i)=Speed(J)
	       ENDDO
	    ENDDO
	    CALL EC_M_Averag(Sample, NNMax, N, MMMax, M, Flag,
     &                 GMean, GTolMean, GCov, GTolCov, MIndep, CIndep,
     &                 Mok, Cok)
            write(*,*) (GMean(i),i=1,3)
            dYaw1 = dYaw1*(180.0/Pi)
            dPitch1 = dPitch1*(180.0/Pi)
            dRoll1 = dRoll1*(180.0/Pi)
C End of DoClassic
         ENDIF 
	 ELSE
	    DO I=1,3
	       MEAN2(I) = DUMMY
	       DO J=1,3
	          STRESS2(I,J) = DUMMY
               ENDDO
	    ENDDO
	    Yaw1 = DUMMY
	    Roll1 = DUMMY
	    Pitch1 = DUMMY
         ENDIF      

         WRITE(OutFile,*)
         WRITE(OutFile,*) 'This run: ', (StartTime(i),i=1,3),
     &                    (StopTime(i),i=1,3)
         WRITE(OutFile,*) 'For this run the remaining',
     *                    ' tilts are:'
         WRITE(OutFile,*) 'Alpha   Beta    Gamma WBias Swing'
         WRITE(OutFile, 60)  RAlpha,RBeta,RGamma,RWBias,Swing
         WRITE(PlfFile, 61) (StartTime(i),i=1,3), (StopTime(i),i=1,3),
     &                       M, (Mok(i),i=1,3),
     &                       Alpha,Beta,Gamma,WBias,
     &                       RAlpha,RBeta,RGamma,RWBias,Swing,
     &                       (MEAN1(ii), ii=1,3),
     &                       ((STRESS1(ii,jj), jj=1,3), ii=1,3),
     &                       Yaw1, Pitch1, Roll1,
     &                       (MEAN2(ii), ii=1,3),
     &                       ((STRESS2(ii,jj), jj=1,3), ii=1,3),
     &                       dYaw1, dPitch1, dRoll1
 60      FORMAT(3(F10.2,1X),F10.4,1X,F10.1)
 61      FORMAT(6(F10.0,1X),4(I10,1x),39(G14.6:,1X))
 62      FORMAT('#', 6(A10,1X), 43(A14,1X))
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
      REAL*8 EC_PH_Q,EC_M_BaseF ! External function calls
      
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
            CALL EC_C_SCal(CalSonic,
     &                  UDUM, VDUM, WDum,
     &                  ERROR(U), ERROR(V), ERROR(W))
            Sample(W) = WDum
         ENDIF
C
C If this is one of Kaijo Denki's, turn it an 90 degrees
C
         IF ((CalSonic(QQType) .EQ. ApKaijoTR90). .OR.
     &       (CalSonic(QQType) .EQ. ApKaijoTR61)) THEN
           Hook = PI*90.D0/180.D0
           Sample(U) =  COS(Hook)*UDum + SIN(Hook)*VDum
           Sample(V) = -SIN(Hook)*UDum + COS(Hook)*VDum
         ENDIF
       
C This is the construction when we have a diagnostic variable
         IF (Have_Uncal(Diagnost)) THEN
            IF ((Error(U).OR.Error(V).OR.Error(W)).OR.
     &          ((RawSampl(ColDiagnostic).LT.0.D0).OR.
     &           (RawSampl(ColDiagnostic).GT.100.D0))) THEN
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
     &             EC_M_BaseF(DNLIMIT + RawSampl(Tcouple),
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
     &         EC_M_BaseF(Dum,NINT(CalHyg(QQFunc)),
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

      REAL*8 EC_M_BaseF
      EXTERNAL EC_M_BaseF

C
C To communicate with function of which zero is searched
C
      COMMON /DUMCOM/ C, DUM, DUMORDER, DUMTYP

      DUMFUNC = EC_M_BaseF(X, DUMTYP, DUMORDER, C) - DUM
      END
      



