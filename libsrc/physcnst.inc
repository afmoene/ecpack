C $Id$
C  physcnst.inc, part of ecpack
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
C First value is from work of Webb, 1980
c     PARAMETER(KoK = 0.0085D0)    ! [m^{3} g^{-1} cm^{-1}]
C This value is due to Arjan van Dijk, 1999
      PARAMETER(KoK = 0.0038)      ! [m^{3} g^{-1} cm^{-1}]

      REAL*8 KwK		   ! extinction for water vapour Krypton
      PARAMETER(KwK = 0.143D0)	   ! [m^{3} g^{-1} cm^{-1}]

      REAL*8 KoLa		   ! extinction for oxygen for Lymann-alpha
      PARAMETER(KoLa = 0.001085D0)      ! [m^{3} g^{-1} cm^{-1}]
      REAL*8 KwLa		   ! extinction for water vapour Lymann-alpha
      PARAMETER(KwLa = 0.09125D0)	   ! [m^{3} g^{-1} cm^{-1}]

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
      PARAMETER(MaxRhoV =   1.D0) ! [kg m^{-3}]
