


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
     &	UStar,dUStar,Tau,dTau)
C
C Construct estimates for surface fluxes from mean values and covariances
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      LOGICAL BadTc
      INTEGER NMax
      REAL*8 EC_Ph_RhoWet,Mean(NMax),Cov(NMax,NMax),HSonic,dHSonic,HTc,dHTc,
     &	LvE,dLvE,LvEWebb,dLvEWebb,Tau,dTau,RhoSon,RhoTc,Frac1,Frac2,Sgn,
     &	p,TolMean(NMax),TolCov(NMax,NMAx),UStar,dUStar
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





      REAL*8 FUNCTION EC_Ph_Q(RhoV,T,P)
C
C Calculate the specific humidity of wet air
C
C input : RhoV : Density of water vapour [kg m^{-3}]
C	  T    : Temperature		 [K]
C	  P    : Pressure		 [Pa]
C
C output : EC_Ph_Q : specific humidity	 [kg kg^{-1}]
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      REAL*8 RhoV,T,P,EC_Ph_RhoWet,RhoWet

      RhoWet = EC_Ph_RhoWet(RhoV,T,P)
      EC_Ph_Q = RhoV/RhoWet 		 ! [kg/kg]

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
C
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
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

      REAL*8 RhoV,T,P

      EC_Ph_RhoDry = P/(Rd*T) - RhoV*Rv/Rd	  ! [kg m^{-3}]

      RETURN
      END





      REAL*8 FUNCTION EC_Ph_RhoWet(RhoV,T,P)
C
C Calculate the density of wet air
C
C input : RhoV : Density of water vapour [kg m^{-3}]
C	  T    : Temperature		 [K]
C	  P    : Pressure		 [Pa]
C
C output : RhoWet : Density of wet air [kg m^{-3}]
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'

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
C
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
      INCLUDE 'physcnst.inc'
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
      IF (Nok .EQ. 0) THEN
         Cxy = DUMMY
	 dCxy = DUMMY
      ENDIF
      RETURN
      END
