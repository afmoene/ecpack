C $Id$
C
C Constants to indicate files
C
      INTEGER InFile,OutFile,TempFile,DumpFile,FluxFile,DataFile,
     &	TrendFile,IntervalFile, ECConfFile
      PARAMETER(InFile	  = 13)
      PARAMETER(OutFile   = 14)
      PARAMETER(TempFile  = 15)
      PARAMETER(DumpFile  = 16)
      PARAMETER(FluxFile  = 17)
      PARAMETER(DataFile  = 18)
      PARAMETER(TrendFile = 19)
      PARAMETER(IntervalFile = 20)
C If Ecconffile is set to 5, input is read from standard input
C else a file named 'ecconf' is opened
      PARAMETER(ECConfFile = 5)
C
C These constants are used to indicate quantities
C
      INTEGER U,V,W,TSonic,Tcouple,Humidity,SpecHum,TTime,
     &        Doy, HourMin, Sec, Diagnost, Tref

      PARAMETER(U	 = 1)
      PARAMETER(V	 = 2)
      PARAMETER(W	 = 3)
      PARAMETER(TSonic	 = 4)
      PARAMETER(Tcouple  = 5)
      PARAMETER(Humidity = 6)
      PARAMETER(SpecHum  = 7)
      PARAMETER(TTime	 = 8)
      PARAMETER(Doy   	 = 9)
      PARAMETER(HourMin	 = 10)
      PARAMETER(Sec   	 = 11)
      PARAMETER(Diagnost = 12)
      PARAMETER(Tref     = 13)
C
C Constants used for maximum array dimensions
C
      INTEGER MMMax,NNMax,MMaxDelay,NNNMax, NMaxOrder

      PARAMETER(MMMax  = 50000) ! Maximum number of samples allowed
      PARAMETER(NNNMax	  = 13) ! Maximum number of uncalibrated quantities in 1 sample
      PARAMETER(NNMax	  =  8) ! Maximum number of calibrated quantities in 1 sample
      PARAMETER(MMaxDelay =  8) ! Largest allowed delay between channels [number of samples]
      PARAMETER(NMaxOrder =  6) ! Maximum order of calibration polynomial
C
C
C Constants used for function types
C
C
      INTEGER NormPoly,LogPoly
      PARAMETER(NormPoly = 1)
      PARAMETER(LogPoly  = 2)
C
C Constants to characterize types of measurement apparatus, all starting
C with ''Ap''
C
      INTEGER ApGillSonic,ApTCouple,ApCampKrypton,ApPt100,ApPsychro,
     &        ApSon3Dcal

      PARAMETER(ApGillSonic   = 1)
      PARAMETER(ApTCouple     = 2)
      PARAMETER(ApCampKrypton = 3)
      PARAMETER(ApPt100       = 4)
      PARAMETER(ApPsychro     = 5)
      PARAMETER(ApSon3Dcal    = 6)

C
C Number of calibration data for each apparatus type
C
      INTEGER  ApNQQ(6) /27,27,27,27,27,33/
C
C Constants to indicate calibration data, all starting with QQ
C
      INTEGER QQType,QQIdent,QQGain,QQOffset,QQX,QQY,QQZ,QQC0,QQC1,
     &	QQC2,QQC3,QQC4,QQC5,QQC6,QQPitch,QQYaw,QQRoll,QQOrder,QQFunc,
     &	QQPath,QQTime,QQExt1,QQExt2,QQExt3,QQExt4,QQExt5,NQQ,QQDelay,
     &   QQExt6, QQExt7, QQExt8, QQExt9, QQExt10, QQExt11

      PARAMETER(NQQ	 = 33) ! Maximum number of calibration coefficients

      PARAMETER(QQType	 =  1) ! Either of the ApXXXXXX numbers
      PARAMETER(QQIdent  =  2) ! True series-number of the apparatus
      PARAMETER(QQGain	 =  3) ! Gain in this experiment
      PARAMETER(QQOffset =  4) ! Offset in this experiment
      PARAMETER(QQX	 =  5) ! Absolute x-position of the focal point
      PARAMETER(QQY	 =  6) ! Absolute y-position of the focal point
      PARAMETER(QQZ	 =  7) ! Absolute z-position of the focal point
      PARAMETER(QQC0	 =  8) ! Zeroth order coefficient of calibration
      PARAMETER(QQC1	 =  9) ! First order coefficient of calibration
      PARAMETER(QQC2	 = 10) ! Second order coefficient of calibration
      PARAMETER(QQC3	 = 11) ! Third order coefficient of calibration
      PARAMETER(QQC4	 = 12) ! Fourth order coefficient of calibration
      PARAMETER(QQC5	 = 13) ! Fifth order coefficient of calibration
      PARAMETER(QQC6	 = 14) ! Sixth order coefficient of calibration
      PARAMETER(QQPitch  = 15) ! Pitch angle of apparatus relative to north
      PARAMETER(QQYaw	 = 16) ! Yaw angle
      PARAMETER(QQRoll	 = 17) ! Roll angle
      PARAMETER(QQOrder  = 18) ! Order of the calibration model
      PARAMETER(QQFunc	 = 19) ! Type of fitfunction
      PARAMETER(QQPath	 = 20) ! Path length over which sensor integrates
      PARAMETER(QQTime	 = 21) ! Time constant
      PARAMETER(QQDelay  = 22) ! [ms] Delay between sampling and output
      PARAMETER(QQExt1	 = 23) ! Additional calibration coefficient
      PARAMETER(QQExt2	 = 24) ! Additional calibration coefficient
      PARAMETER(QQExt3	 = 25) ! Additional calibration coefficient
      PARAMETER(QQExt4	 = 26) ! Additional calibration coefficient
      PARAMETER(QQExt5	 = 27) ! Additional calibration coefficient
      PARAMETER(QQExt6	 = 28) ! Additional calibration coefficient
      PARAMETER(QQExt7	 = 29) ! Additional calibration coefficient
      PARAMETER(QQExt8	 = 30) ! Additional calibration coefficient
      PARAMETER(QQExt9	 = 31) ! Additional calibration coefficient
      PARAMETER(QQExt10	 = 32) ! Additional calibration coefficient
      PARAMETER(QQExt11	 = 33) ! Additional calibration coefficient
