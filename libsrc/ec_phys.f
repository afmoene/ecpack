C  ec_phys.f, part of ecpack
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




C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C
C Flux/physics routines
C
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################
C ########################################################################




C
C ########################################################################
C
C Small routine to estimate surface fluxes (with their tolerances!) from
C processed mean values, covariances and tolerances.
C
C ########################################################################
C





      SUBROUTINE EC_Ph_Flux(Mean,NMax,Cov,TolMean,TolCov,p,BadTc,
     &	HSonic,dHSonic,HTc,dHTc,LvE,dLvE,LvEWebb,dLvEWebb,
     &	UStar,dUStar,Tau,dTau, FCO2, dFCO2, FCO2Webb, dFCO2Webb, WebVel,
     &  VectWind, dVectWind, DirFrom, dDirFrom, dirYaw)
C     ****f* ec_phys.f/EC_Ph_Flux
C NAME
C     EC_Ph_Flux
C SYNOPSIS
C     CALL EC_Ph_Flux(Mean,NMax,Cov,TolMean,TolCov,p,BadTc,
C                     HSonic,dHSonic,HTc,dHTc,LvE,dLvE,LvEWebb,dLvEWebb,
C                     UStar,dUStar,Tau,dTau, 
C                     FCO2, dFCO2, FCO2Webb, dFCO2Webb,
C                     WebVel, VectWind, dVectWind, DirFrom, dDirFrom, dirYaw)
C FUNCTION
C     Construct estimates for surface fluxes from mean values and covariances
C INPUTS
C     Mean   : [REAL*8(NMax)]
C              Means of all variables
C     NMax   : [INTEGER]
C              Maximum number of variables
C     Cov    : [REAL*8(NMax,NMax)]
C              Covariances of all variables
C     TolMean: [REAL*8(NMax)]
C              Tolerances in means of all variables
C     TolCov : [REAL*8(NMax,NMax)]
C              Tolerances in covariances of all variables
C     p      : [REAL*8]
C              atmosperic pressure (Pa)
C     BadTc  : [LOGICAL]
C              indicator whether thermocouple temperature is corrupt
C     WebVel : [REAL*8]
C              Webb velocity (m/s)
C     DirYaw : [Real*8]
C              Yaw rotation angle (degrees)
C OUTPUT
C     HSonic : [REAL*8]
C              sensible heat flux with sonic temperature (W/m^2)
C     dHSonic: [REAL*8]
C              tolerance in sensible heat flux with sonic temperature (W/m^2)
C     HTc    : [REAL*8]
C              sensible heat flux with thermocouple temperature (W/m^2)
C     dHTc   : [REAL*8]
C              tolerance in sensible heat flux with thermocouple temperature (W/m^2)
C     LvE    : [REAL*8]
C              wq-covariance latent heat flux (W/m^2)
C     dLvE   : [REAL*8]
C              tolerance in wq-covariance latent heat flux (W/m^2)
C     LvEWebb: [REAL*8]
C              Webb term for latent heat  flux (W/m^2)
C     dLvEWebb: [REAL*8]
C              tolerance in Webb term for latent heat  flux (W/m^2)
C     UStar  : [REAL*8]
C              friction velocity (m/s)
C     dUStar : [REAL*8]
C              tolerance in friction velocity (m/s)
C     Tau    : [REAL*8]
C              surface shear stress (N/m^2)
C     dTau   : [REAL*8]
C              tolerance in surface shear stress (N/m^2)
C     FCO2   : [REAL*8]
C              w CO2-covariance CO2 flux (kg/m^2 s)
C     dFCO2  : [REAL*8]
C              tolerance in w CO2-covariance CO2 flux (kg/m^2 s)
C     FCO2Webb: [REAL*8]
C              Webb term for CO2  flux (kg/m^2 s)
C     dLvEWebb: [REAL*8]
C              tolerance in Webb term for CO2 flux (kg/m^2 s)
C     VectWind: [Real*8]
C              vector wind velocity (m/s)
C     dVectWind: [Real*8]
C              tolerance in vector wind velocity (m/s)
C     DirFrom : [Real*8]
C              wind direction (degrees)
C     dDirFrom : [Real*8]
C              tolerance in wind direction (degrees)
C AUTHOR
C     Arjan van Dijk, Arnold Moene
C HISTORY
C     07-10-2002: added CO2 fluxes and WebVel to interface. Webb-term
C                 is now computed with Webvel, rather than Mean(W)
C     $Name$ 
C     $Id$
C USES
C     EC_Ph_RhoWet
C     parcnst.inc
C     Cp
C     Lv
C     ***
      IMPLICIT NONE
      REAL*8 x(3),b(3),a(3,3),aInv(3,3)
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      LOGICAL BadTc
      INTEGER NMax
      REAL*8 EC_Ph_RhoWet,Mean(NMax),Cov(NMax,NMax),HSonic,dHSonic,HTc,
     &  dHTc,
     &	LvE,dLvE,LvEWebb,dLvEWebb,Tau,dTau,RhoSon,RhoTc,Frac1,Frac2,Sgn,
     &	p,TolMean(NMax),TolCov(NMax,NMAx),UStar,dUStar,
     &  FCO2, dFCO2, FCO2Webb, dFCO2Webb, WebVel, VectWind, dVectWind,
     &  DirFrom, dDirFrom, DirYaw, LocalDir
C
C Sensible heat flux [W m^{-2}]
C
      RhoSon = EC_Ph_RhoWet(Mean(Humidity),Mean(TSonic),p)
      HSonic = Cp*RhoSon*Cov(W,TSonic)
      dHSonic = Cp*RhoSon*TolCov(W,TSonic)

      IF (.NOT.BadTc) THEN
        RhoTc = EC_Ph_RhoWet(Mean(Humidity),Mean(TCouple),p)
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
      LvEWebb = Lv*WebVel*Mean(Humidity)
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
C Friction velocity (vector sum of two components
C
      IF (.NOT. ((COV(W,U) .EQ. DUMMY) .AND.
     &            (COV(W,V) .EQ. DUMMY))) THEN
	  UStar = SQRT(SQRT(Cov(W,U)**2+Cov(W,U)**2))
      ENDIF

      dUStar = (0.25*(2*TolCov(W,U)*ABS(COV(W,U)) + 
     &                 2*TolCov(W,V)*ABS(COV(W,V)))/UStar**2)*
     &           UStar
C
C Friction force [N m^{-2}]
C
      IF (.NOT.BadTc) THEN
        Tau = 0.5*(RhoSon+RhoTc)*(ABS(UStar))**2.D0
      ELSE
        Tau = (RhoSon)*(ABS(UStar))**2.D0
      ENDIF
      dTau = (0.5*(2*ABS(COV(W,U))*TolCov(W,U) + 
     &             2*ABS(COV(W,V))*TolCov(W,V))/UStar**2)*
     &         Ustar**2
C
C CO2 flux [kg m^{-2} s^{-1}]
C
      FCO2 = Cov(W,CO2)
      dFCO2 = TolCov(W,CO2)
      FCO2Webb = WebVel*Mean(CO2)
C Frac1 is set to zero, assuming no error in W (else error in 
C Webb term would be enormous
      Frac1 = 0.D0
      Frac2 = ABS(TolMean(CO2)/Mean(CO2))

      dFCO2Webb = FCO2Webb*(Frac1**2.D0+Frac2**2.D0)**0.5D0   
C 
C Vector wind
C
      VectWind = SQRT(Mean(U)**2+Mean(V)**2)
      dVectWind = 0.5 * (TolMean(U)*ABS(Mean(U)) + 
     &                   TolMean(V)*ABS(Mean(V)))
C      
C Wind direction
C
      IF (DirYaw.LT.180.D0) THEN
          DirFrom = DirYaw - 180.D0
      ELSE
          DirFrom = DirYaw + 180.D0
      ENDIF
      DirFrom = DirFrom - ATAN2(Mean(V), -Mean(U))*180/PI
      IF (DirFrom .LT. 0) THEN
         DirFrom = DirFrom + 360.0D0
      ENDIF
      IF (DirFrom .GT. 360.D0) THEN
         DirFrom = DirFrom - 360.0D0
      ENDIF
      
      LocalDir = ATAN2(Mean(U), Mean(V))
      dDirFrom = 2.D0*180.D0*ATAN(
     &              SQRT(Cos(LocalDir)*Cov(V,V)+
     &                   SIN(LocalDir)*Cov(U,U))/VectWind)/Pi
      RETURN
      END





      REAL*8 FUNCTION EC_Ph_Q(RhoV,T,P)
C     ****f* ec_phys.f/EC_Ph_Q
C NAME
C     EC_Ph_Q
C SYNOPSIS
C     Spec_hum = EC_Ph_Q(RhoV,T,P)
C FUNCTION
C     Calculate the specific humidity of wet air 
C INPUTS
C      Rhov : [REAL*8]
C             Density of air (kg/m^3)
C      T    : [REAL*8]
C             Temperature (K)
C      P    : [REAL*8]
C             Pressure (Pa)
C RETURN VALUE
C      return value : [REAL*8]
C             Specific humidity (kg/kg)
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     EC_Ph_RhoWet
C     ***
      IMPLICIT NONE
      REAL*8 RhoV,T,P,EC_Ph_RhoWet,RhoWet

      RhoWet = EC_Ph_RhoWet(RhoV,T,P)
      EC_Ph_Q = RhoV/RhoWet 		 ! [kg/kg]

      RETURN
      END


      REAL*8 FUNCTION EC_Ph_QCO2(RHOCO2,RHOV,T,P)
C     ****f* ec_phys.f/EC_Ph_QCO2
C NAME
C     EC_Ph_QCO2
C SYNOPSIS
C     Spec_CO2 = EC_Ph_QCO2(RHOCO2,RHOV,T,P)
C FUNCTION
C     Calculate the specific CO2 concentration of wet air 
C INPUTS
C      RhoCO2 : [REAL*8]
C             Density of CO2 (kg/m^3)
C      Rhov : [REAL*8]
C             Density of air (kg/m^3)
C      T    : [REAL*8]
C             Temperature (K)
C      P    : [REAL*8]
C             Pressure (Pa)
C RETURN VALUE
C      return value : [REAL*8]
C             Specific CO2 (kg/kg)
C AUTHOR
C     Arnold Moene
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     EC_Ph_RhoWet
C     ***
      IMPLICIT NONE
      REAL*8 RhoCO2, RhoV,T,P,EC_Ph_RhoWet,RhoWet

      RhoWet = EC_Ph_RhoWet(RhoV,T,P)
      EC_Ph_QCO2 = RhoCO2/RhoWet ! [kg/kg]

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





      REAL*8 FUNCTION EC_Ph_RhoDry(RhoV,T,P)
C     ****f* ec_phys.f/EC_Ph_RhoDry
C NAME
C     EC_Ph_RhoDry
C SYNOPSIS
C     Rho_dry = EC_Ph_RhoDry(RhoV,T,P)
C FUNCTION
C     Calculate the density of dry air component in wet air
C     Via Dalton's law : Pressure is sum of partial pressures :
C     P = RhoV*Rv*T + RhoD*Rd*T
C INPUTS
C      Rhov : [REAL*8]
C             Density of air (kg/m^3)
C      T    : [REAL*8]
C             Temperature (K)
C      P    : [REAL*8]
C             Pressure (Pa)
C RETURN VALUE
C      return value : [REAL*8]
C             Density of dry part of air (kg/m^3)
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     Rd
C     Rv
C     ***
      IMPLICIT NONE
      INCLUDE 'physcnst.inc'

      REAL*8 RhoV,T,P

      EC_Ph_RhoDry = P/(Rd*T) - RhoV*Rv/Rd	  ! [kg m^{-3}]

      RETURN
      END





      REAL*8 FUNCTION EC_Ph_RhoWet(RhoV,T,P)
C     ****f* ec_phys.f/EC_Ph_RhoWet
C NAME
C     EC_Ph_RhoWet
C SYNOPSIS
C     Rho_dry = EC_Ph_RhoWet(RhoV,T,P)
C FUNCTION
C     Calculate the density of wet air
C INPUTS
C      Rhov : [REAL*8]
C             Density of air (kg/m^3)
C      T    : [REAL*8]
C             Temperature (K)
C      P    : [REAL*8]
C             Pressure (Pa)
C RETURN VALUE
C      return value : [REAL*8]
C             Density of wet air (kg/m^3)
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$ 
C     $Id$
C USES
C     EC_Ph_RhoDry
C     ***

      REAL*8 RhoV,T,P,EC_Ph_RhoDry

      EC_Ph_RhoWet = EC_Ph_RhoDry(RhoV,T,P) + RhoV ! [kg m^{-3}]

      RETURN
      END
C
C
C Routines to calculate structure parameters
C
C


      SUBROUTINE EC_Ph_Struct(Sample,NMax,N,MMax,M,Flag,XIndex,YIndex,
     &  R,dR,Freq,CIndep,Cxy,dCxy)
C     ****f* ec_phys.f/EC_Ph_Struct
C NAME
C     EC_Ph_Struct
C SYNOPSIS
C     CALL EC_Ph_Struct(Sample,NMax,N,MMax,M,Flag,
C                       XIndex,YIndex,
C                       R,dR,Freq,CIndep,Cxy,dCxy)
C FUNCTION
C     Calculate structure parameters <(x(r)-x(r+R))*(y(r)-y(r+R))>/R^2/3
C INPUTS
C     Sample   : [REAL*8(NMax,MMax)]  
C                Samples (quantities in first dimension, samples in 
C                second dimension)
C     NMax     : [INTEGER]  
C                physical first dimension of array Sample.
C     N        : [INTEGER] 
C                actual number of meaningful quantities in array Sample.
C     MMax     : [INTEGER]  
C                physical second dimension of array Sample.
C     M        : [INTEGER]  
C                actual number of meaningful samples in array Sample.
C     Flag     : [LOGICAL(NMax,MMax)] 
C                if Flag(i,j) is true, then something is wrong with 
C                quantity i in sample j.
C     XIndex   : [INTEGER]  
C                indicator of first quantity involved in 
C                structure function.
C     YIndex   : [INTEGER]  
C                indicator of second quantity involved in
C                structure function.
C     R        : [REAL*8] : 
C                separation in meters at which one wants to estimate
C                the structure function.
C     Freq     : [REAL*8] : 
C                Sampling frequency in s^-1.
C     CIndep   : [INTEGER[NMax,NMax])  
C                Number of independent contributions by array Sample 
C                to covariance between quantities selected with
C                XIndex and YIndex.
C     Rhov     : [REAL*8]
C                Density of air (kg/m^3)
C     T        : [REAL*8]
C                Temperature (K)
C     P        : [REAL*8]
C                Pressure (Pa)
C OUTPUT 
C     dR       : [REAL*8]  
C                separation in meters corresponding with a delay
C                of one sample (i.e. tolerance in R)
C     cxy      : [REAL*8] 
C                Structure parameter.
C     dcxy     : [REAL*8] 
C                Tolerance of cxy.
C AUTHOR
C     Arjan van Dijk
C HISTORY
C     $Name$ 
C     $Id$
C     Revision:  20-9-2002: added implicit none and inclusion of
C                           parcnst.inc (needed for constants U,V, and W)
C USES
C     parcnst.inc
C     ***
      IMPLICIT NONE
      INCLUDE 'parcnst.inc'

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
      NSeparate = MAX(1,NINT(R/dR))
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
      IF (Nok .EQ. 0) THEN
         Cxy = DUMMY
	 dCxy = DUMMY
      ENDIF
      RETURN
      END
