C $Id$
C  parcnst.inc, part of ecpack
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

C     ****f* parcnst.inc/parcnst.inc
C NAME
C     parcnst.inc
C FUNCTION
C     include file that defines various parameters (indicating
C     instruments etc.)
C     ***
C
C Constants to indicate files
C
      INTEGER InFile,OutFile,TempFile,DumpFile,FluxFile,DataFile,
     &	TrendFile,IntervalFile, ECConfFile, PlfIntFile, PlfFile
      PARAMETER(InFile	  = 13)
      PARAMETER(OutFile   = 14)
      PARAMETER(TempFile  = 15)
      PARAMETER(DumpFile  = 16)
      PARAMETER(FluxFile  = 17)
      PARAMETER(DataFile  = 18)
      PARAMETER(TrendFile = 19)
      PARAMETER(IntervalFile = 20)
      PARAMETER(PlfIntFile = 21)  ! File in which intermediate mean velocities are written
      PARAMETER(PlfFile = 22)     ! File to which planerfit results are written
C If Ecconffile is set to 5, input is read from standard input
C else a file named 'ecconf' is opened
      PARAMETER(ECConfFile = 5)
C
C Constants used for maximum array dimensions
C
      INTEGER MMMax,NNMax,MMaxDelay,NNNMax, NMaxOrder

      PARAMETER(MMMax = 100000) ! Maximum number of samples allowed
      PARAMETER(NNNMax	  = 15) ! Maximum number of uncalibrated quantities in 1 sample
      PARAMETER(NNMax	  =  10) ! Maximum number of calibrated quantities in 1 sample
      PARAMETER(MMaxDelay =  8) ! Largest allowed delay between channels [number of samples]
      PARAMETER(NMaxOrder =  6) ! Maximum order of calibration polynomial
C
C These constants are used to indicate quantities
C
      INTEGER U,V,W,TSonic,Tcouple,Humidity,SpecHum,TTime,
     &        Doy, HourMin, Sec, Diagnost, Tref, CO2, SpecCO2

      PARAMETER(U	 = 1)
      PARAMETER(V	 = 2)
      PARAMETER(W	 = 3)
      PARAMETER(TSonic	 = 4)
      PARAMETER(Tcouple  = 5)
      PARAMETER(Humidity = 6)
      PARAMETER(SpecHum  = 7)
      PARAMETER(CO2      = 8)
      PARAMETER(SpecCO2  = 9)
      PARAMETER(TTime	 = 10)
      PARAMETER(Doy   	 = 11)
      PARAMETER(HourMin	 = 12)
      PARAMETER(Sec   	 = 13)
      PARAMETER(Diagnost = 14)
      PARAMETER(Tref     = 15)


      CHARACTER*6 QName(NNNMax) 
     &                  /'     u',
     &                   '     v',
     &                   '     w',
     &                   'T(Son)',
     &                   'T(Cpl)',
     &                   '  rhoV',
     &                   '     q',
     &                   'RhoCO2',
     &                   '  qCO2',
     &                   '  time',
     &                   '   doy',
     &                   'hr min',
     &                   '   sec',
     &                   ' diagn',
     &                   '  Tref'/

      CHARACTER*9 UName(NNNMax) 
     &                  /'    [m/s]',
     &                   '    [m/s]',
     &                   '    [m/s]',
     &                   '      [K]',
     &                   '      [K]',
     &                   '[kg ^m-3]',
     &                   '  [kg/kg]',
     &                   '[kg m^-3]',
     &                   '  [kg/kg]',
     &                   '      [-]',
     &                   '      [-]',
     &                   '      [-]',
     &                   '      [-]',
     &                   '      [-]',
     &                   '      [K]'/
      
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
     &        ApSon3Dcal, ApMierijLyma, ApLiCor7500, ApKaijoTR90,
     &        ApKaijoTR61

      PARAMETER(ApGillSonic   = 1)
      PARAMETER(ApTCouple     = 2)
      PARAMETER(ApCampKrypton = 3)
      PARAMETER(ApPt100       = 4)
      PARAMETER(ApPsychro     = 5)
      PARAMETER(ApSon3Dcal    = 6)
      PARAMETER(ApMierijLyma  = 7)
      PARAMETER(ApLiCor7500   = 8)
      PARAMETER(ApKaijoTR90   = 9)
      PARAMETER(ApKaijoTR61   = 10)

C
C Number of calibration data for each apparatus type
C
      INTEGER  ApNQQ(8) /27,27,27,27,27,33,27,27/
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
      PARAMETER(QQYaw    = 15) ! Yaw angle of apparatus relative to north
      PARAMETER(QQPitch	 = 16) ! Pitch angle
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

      REAL DUMMY
      PARAMETER(DUMMY = -9999)


