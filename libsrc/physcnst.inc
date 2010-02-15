C $Id$
C  physcnst.inc, part of ecpack
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


C
C Physics constants
C
C     ****d* physcnst.inc/Rd
C NAME
C     Rd
C FUNCTION
C     Gas constant for dry air
C SOURCE
      REAL*8 Rd                 ! Gass constant for dry air
      PARAMETER(Rd = 287.04D0)            ! [J/Kg.K]
C     ***

C     ****d* physcnst.inc/Rv
C NAME
C     Rv
C FUNCTION
C     Gas constant for water vapour
C SOURCE
      REAL*8 Rv                 ! Gass constant for water vapour
      PARAMETER(Rv = 461.5D0)            ! [J/Kg.K]
C     ***

C     ****d* physcnst.inc/RGas
C NAME
C     RGas
C FUNCTION
C     Universal gas constant 
C SOURCE
      REAL*8 RGas                ! Universal gas constant
      PARAMETER(RGAs = 8314.D0)         ! [J/kmol.K]
C     ***

C     ****d* physcnst.inc/Epsilon
C NAME
C     Epsilon
C FUNCTION
C     Infinitisimal number 
C SOURCE
      REAL*8 Epsilon                ! Infinitesimal number
      PARAMETER(Epsilon = 1.D-30)        ! [1]
C     ***

C     ****d* physcnst.inc/Pi
C NAME
C     Pi
C FUNCTION
C     Pi
C SOURCE
      REAL*8 Pi                 ! Pi
      PARAMETER(Pi=3.1415926535897932385D0) ! [1]
C     ***

C     ****d* physcnst.inc/Kelvin
C NAME
C     Kelvin
C FUNCTION
C     Temperature of 0 degree Celsius in Kelvin
C SOURCE
      REAL*8 Kelvin                ! 0 degrees Celsius
      PARAMETER(Kelvin = 273.15D0)        ! [K]
C     ***

C     ****d* physcnst.inc/GammaR
C NAME
C     GammaR
C FUNCTION
C     Constant in correction of sonic temperature for cross-wind
C     (part of the Schotanus correction)
C SOURCE
      REAL*8 GammaR                ! gamma*R (in relationship between 
                                        ! sonic temperature and speed of sound
      PARAMETER(GammaR = 403.D0)        ! [m^{2}s^{-2}K^{-1}]
C     ***

C     ****d* physcnst.inc/DifCo2
C NAME
C     DiffCo2
C FUNCTION
C     Molecular diffusivity CO2 at 30 degree Celcius (?) 
C SOURCE
      REAL*8 DifCO2           ! Molecular diffusivity CO2 at 30 degr
      PARAMETER(DifCO2 = 15.6D-6)  ! [m^{2}s^{-1}]
C     ***

C     ****d* physcnst.inc/DifH20
C NAME
C     DiffH20
C FUNCTION
C     Molecular diffusivity H20 at 30 degree Celcius (?) 
C SOURCE
      REAL*8 DifH2O           ! Molecular diffusivity H20 at 30 degr
      PARAMETER(DifH2O = 25.7D-6)  ! [m^{2}s^{-1}]
C     ***

C     ****d* physcnst.inc/Karman
C NAME
C     Karman
C FUNCTION
C     Von Karman constant
C SOURCE
      REAL*8 Karman           ! Von Karman constant
      PARAMETER(Karman = 0.4D0)       ! [1]
C     ***

C     ****d* physcnst.inc/GG
C NAME
C     GG
C FUNCTION
C     Gravitational acceleration
C SOURCE
      REAL*8 GG            ! Gravitation constant
      PARAMETER(GG = 9.81D0)       ! [m s^{-2}]
C     ***

C     ****d* physcnst.inc/MO2
C NAME
C     MO2
C FUNCTION
C     molecular weight of oxygen
C SOURCE
      REAL*8 MO2           ! molecular weight of oxygen
      PARAMETER(MO2 = 32.D0)       ! [g mol^{-1}]
C     ***

C     ****d* physcnst.inc/MAir
C NAME
C     MAir
C FUNCTION
C     molecular weight of dry air
C SOURCE
      REAL*8 MAir           ! Molecular weight of dry air
      PARAMETER(MAir = 28.966D0)   ! [g/mol]
C     ***

C     ****d* physcnst.inc/MVapour
C NAME
C     MVapour
C FUNCTION
C     molecular weight of water vapour
C SOURCE
      REAL*8 MVapour           ! Molecular weight of water vapour
      PARAMETER(MVapour = 18.016D0)! [g/mol]
C     ***

C     ****d* physcnst.inc/Mu
C NAME
C     Mu
C FUNCTION
C     ratio: Mu = Mair/Mv
C SOURCE
      REAL*8 Mu            ! Ratio: Mu := Mair/Mv
      PARAMETER(Mu = 1.6078D0)       ! [1]
C     ***

C     ****d* physcnst.inc/Kok
C NAME
C     Kok
C FUNCTION
C     extinction for oxygen for Krypton hygrometer
C NOTES
C     First the value from work of Webb (1980) was used 
C     (KoK = 0.0085 [m^{3} g^{-1} cm^{-1}]). Now
C     a value due to van Dijk (1999) is used.
C     This constant may disappear from the library and may become
C     part of the user-supplied instrument specifications
C SOURCE
      REAL*8 KoK           ! extinction for oxygen for Krypton
      PARAMETER(KoK = 0.0038D0)      ! [m^{3} g^{-1} cm^{-1}]
C     ***

C     ****d* physcnst.inc/Kwk
C NAME
C     Kwk
C FUNCTION
C     extinction for water vapour for Krypton hygrometer
C SOURCE
      REAL*8 KwK           ! extinction for water vapour Krypton
      PARAMETER(KwK = 0.143D0)       ! [m^{3} g^{-1} cm^{-1}]
C     ***

C     ****d* physcnst.inc/Kola
C NAME
C     Kola
C FUNCTION
C     Extinction for oxygen for Lymann-alpha
C SOURCE
      REAL*8 KoLa           ! extinction for oxygen for Lymann-alpha
      PARAMETER(KoLa = 0.001085D0)      ! [m^{3} g^{-1} cm^{-1}]
C     ***

C     ****d* physcnst.inc/Kwla
C NAME
C     Kwla
C FUNCTION
C     Extinction for water vapour for Lymann-alpha
C SOURCE
      REAL*8 KwLa           ! extinction for water vapour Lymann-alpha
      PARAMETER(KwLa = 0.09125D0)       ! [m^{3} g^{-1} cm^{-1}]
C     ***

C     ****d* physcnst.inc/FracO2
C NAME
C     FracO2
C FUNCTION
C     Fraction of O2 molecules in air
C SOURCE
      REAL*8 FracO2           ! fraction of O2 molecules in air
      PARAMETER(FracO2 = 0.21D0)   ! [1]
C     ***

c     ****d* physcnst.inc/MinT
C NAME
C     MinT
C FUNCTION
C     Lower limit for acceptance of temperature sample
C SOURCE
      REAL*8 MinT
      PARAMETER(MinT    = 243.D0) ! [Kelvin]
C     ***

c     ****d* physcnst.inc/MaxT
C NAME
C     MaxT
C FUNCTION
C     Upper limit for acceptance of temperature sample
C SOURCE
      REAL*8 MaxT
      PARAMETER(MaxT    = 333.D0) ! [Kelvin]
C     ***

c     ****d* physcnst.inc/MinRhov
C NAME
C     MinRhov
C FUNCTION
C     Lower limit for acceptance of water vapour sample
C SOURCE
      REAL*8 MinRhoV
      PARAMETER(MinRhoV =   0.D0) ! [kg m^{-3}]
C     ***

c     ****d* physcnst.inc/MaxRhov
C NAME
C     MaxRhov
C FUNCTION
C     Upper limit for acceptance of water vapour sample
C SOURCE
      REAL*8 MaxRhoV
      PARAMETER(MaxRhoV =   1.D0) ! [kg m^{-3}]
C     *** 

c     ****d* physcnst.inc/MinRhoCO2
C NAME
C     MinRhoCO2
C FUNCTION
C     Lower limit for acceptance of CO2 sample
C SOURCE
      REAL*8 MinRhoCO2
      PARAMETER(MinRhoCO2 =   0.D0) ! [kg m^{-3}]
C     ***

c     ****d* physcnst.inc/MaxRhoCO2
C NAME
C     MaxRhoCO2
C FUNCTION
C     Upper limit for acceptance of CO2 sample
C SOURCE
      REAL*8 MaxRhoCO2
      PARAMETER(MaxRhoCO2 =   1.D0) ! [kg m^{-3}]
C     *** 
