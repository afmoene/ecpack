      IMPLICIT NONE
C
C Physics constants
C
      REAL*8 Rd 			    ! Gass constant for dry air
      PARAMETER(Rd = 287.04D0)		    ! [J/Kg.K]
      REAL*8 Rv 			    ! Gass constant for water vapour
      PARAMETER(Rv = 461.5D0)		    ! [J/Kg.K]
      REAL*8 RGas			    ! Universal gas constant
      PARAMETER(RGAs = 8314.D0) 	    ! [?]
      REAL*8 Epsilon			    ! Infinitesimal number
      PARAMETER(Epsilon = 1.D-20)	    ! [1]
      REAL*8 Pi 			    ! Pi
      PARAMETER(Pi=3.1415926535897932385D0) ! [1]
      REAL*8 Kelvin			    ! 0 degrees Celsius
      PARAMETER(Kelvin = 273.15D0)	    ! [K]
      REAL*8 GammaR			    ! ?????????????
      PARAMETER(GammaR = 403.D0)	    ! [m^{2}s^{-2}K^{-1}]
      REAL*8 DifCO2		   ! Molecular diffusivity CO2 at 30 degr
      PARAMETER(DifCO2 = 15.6D-6)  ! [m^{2}s^{-1}]
      REAL*8 DifH2O		   ! Molecular diffusivity H20 at 30 degr
      PARAMETER(DifH2O = 25.7D-6)  ! [m^{2}s^{-1}]
      REAL*8 Karman		   ! Von Karman constant
      PARAMETER(Karman = 0.4)	   ! [1]
      REAL*8 GG 		   ! Gravitation constant
      PARAMETER(GG = 9.81)	   ! [m s^{-2}]
      REAL*8 MO2		   ! molecular weight of oxygen
      PARAMETER(MO2 = 32.D0)	   ! [g mol^{-1}]
      REAL*8 MAir		   ! Molecular weight of dry air
      PARAMETER(MAir = 28.966D0)   ! [g/mol]
      REAL*8 MVapour		   ! Molecular weight of water vapour
      PARAMETER(MVapour = 18.016D0)! [g/mol]
      REAL*8 Mu 		   ! Ratio: Mu := Mair/Mv
      PARAMETER(Mu = 1.6078D0)	   ! [1]
      REAL*8 KoK		   ! extinction for oxygen for Krypton
      PARAMETER(KoK = 0.0085D0)    ! [m^{3} g^{-1} cm^{-1}]
      REAL*8 KwK		   ! extinction for water vapour Krypton
      PARAMETER(KwK = 0.143D0)	   ! [m^{3} g^{-1} cm^{-1}]
      REAL*8 FracO2		   ! fraction of O2 molecules in air
      PARAMETER(FracO2 = 0.21D0)   ! [1]
      REAL*8 Cp 		   ! Specific heat of air
      PARAMETER(Cp = 1004.67D0)    ! [J Kg^{-1} K^{-1}]
      REAL*8 Lv 		   ! Vaporisation heat of water
      PARAMETER(Lv = 2.45D6)	   ! [J Kg^{-1}]
C
C Limits for acceptance of humidity and temperature
C
      REAL*8 MinT,MaxT,MinRhoV,MaxRhoV
      PARAMETER(MinT	= -20.D0) ! [Celsius]
      PARAMETER(MaxT	=  60.D0) ! [Celsius]
      PARAMETER(MinRhoV =   0.D0) ! [kg m^{-3}]
      PARAMETER(MaxRhoV = 100.D0) ! [kg m^{-3}]
C
C Constants to indicate files
C
      INTEGER InFile,OutFile,TempFile,DumpFile,FluxFile,DataFile,
     &	TrendFile,IntervalFile
      PARAMETER(InFile	  = 13)
      PARAMETER(OutFile   = 14)
      PARAMETER(TempFile  = 15)
      PARAMETER(DumpFile  = 16)
      PARAMETER(FluxFile  = 17)
      PARAMETER(DataFile  = 18)
      PARAMETER(TrendFile = 19)
      PARAMETER(IntervalFile = 20)
C
C These constants are used to indicate quantities
C
      INTEGER U,V,W,TSonic,Tcouple,Humidity,SpecHum,TTime

      PARAMETER(U	 = 1)
      PARAMETER(V	 = 2)
      PARAMETER(W	 = 3)
      PARAMETER(TSonic	 = 4)
      PARAMETER(Tcouple  = 5)
      PARAMETER(Humidity = 6)
      PARAMETER(SpecHum  = 7)
      PARAMETER(TTime	 = 8)
C
C Constants used for maximum array dimensions
C
      INTEGER MMMax,NNMax,MMaxDelay,NNNMax

      PARAMETER(MMMax  = 50000) ! Maximum number of samples allowed
      PARAMETER(NNNMax	  = 20) ! Maximum number of uncalibrated quantities in 1 sample
      PARAMETER(NNMax	  =  8) ! Maximum number of calibrated quantities in 1 sample
      PARAMETER(MMaxDelay =  8) ! Largest allowed delay between channels [number of samples]
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
      INTEGER ApGillSonic,ApTCouple,ApCampKrypton,ApPt100,ApPsychro

      PARAMETER(ApGillSonic   = 1)
      PARAMETER(ApTCouple     = 2)
      PARAMETER(ApCampKrypton = 3)
      PARAMETER(ApPt100       = 4)
      PARAMETER(ApPsychro     = 5)
C
C Constants to indicate calibration data, all starting with QQ
C
      INTEGER QQType,QQIdent,QQGain,QQOffset,QQX,QQY,QQZ,QQC0,QQC1,
     &	QQC2,QQC3,QQC4,QQC5,QQC6,QQPitch,QQYaw,QQRoll,QQOrder,QQFunc,
     &	QQPath,QQTime,QQExt1,QQExt2,QQExt3,QQExt4,QQExt5,NQQ,QQDelay

      PARAMETER(NQQ	 = 27) ! Maximum number of calibration coefficients

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
