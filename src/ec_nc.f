C  ec_nc.f, part of ecpack
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



C...........................................................................
C FUNCTION:  EC_NCDF_CSI2HOUR
C PURPOSE:   To convert Campbell hour/min to decimal hours
C INTERFACE: HOURMIN    IN     Campbell hour/minutes (REAL)
C            return value      decimal hours (REAL)
C AUTHOR:    Arnold Moene
C DATE:      May 28, 2001
C...........................................................................
      REAL*8 FUNCTION EC_NCDF_CSI2HOUR(HOURMIN)

      REAL*8 HOURMIN
      INTEGER HOUR, MIN

      HOUR = INT(HOURMIN/100.0D0)
      MIN = INT(HOURMIN-HOUR*100D0)

      EC_NCDF_CSI2HOUR = DBLE(HOUR + MIN/60.0D0)
      END
C...........................................................................
C SUBROUTINE:  EC_NCDF_FINDNOSEC
C PURPOSE   :  To find the index in a NetCDF File for a given DOY and
C               time without having seconds
C...........................................................................
      SUBROUTINE EC_NCDF_FINDNOSEC(FDOY, FHOURMIN, NCID,
     +                    NCdoyID, NChmID, IND)
      INTEGER FDOY, FHOURMIN, NCID, NCdoyID, NChmID, IND,
     +        DIMID, MAXIND, ATTEMPTS, PREVCH, IND1, IND2, IND3,
     +        CURCH
      REAL*8 DOY1, DOY2, DOY3, HM1, HM2, HM3,
     +     SAMPF, TIME1, TIME2, TIME3, HM4, DOY4,
     +     DOYD, HOURD, TIMED
      REAL*8 EC_NCDF_CSI2HOUR
      LOGICAL OK

      STATUS = NF_INQ_VARDIMID(NCID,NChmID,DIMID)
      STATUS = NF_INQ_DIMLEN(NCID, DIMID,MAXIND )
      IF (MAXIND .GT. 0) THEN
         IND1 = 1
         IND3 = MAXIND
         OK = .FALSE.
         DO WHILE (.NOT. OK)
           STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND1, DOY1)
           STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND1, HM1)
	   IF ((DOY1 .GT. 0) .AND. (HM1 .GE. 0) .OR.
     +         (IND1 .EQ. MAXIND)) THEN
	      OK = .TRUE.
	   ELSE
	      IND1 = IND1 + 1
	   ENDIF
         ENDDO
         OK = .FALSE.
         DO WHILE (.NOT. OK)
           STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND3, DOY3)
           STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND3, HM3)
	   IF ((DOY3 .GT. 0) .AND. (HM3 .GE. 0) .OR.
     +         (IND3 .EQ. 1)) THEN
	      OK = .TRUE.
	   ELSE
	      IND3 = IND3 - 1
	   ENDIF
         ENDDO
         IND2 = ANINT((IND3+IND1)/2.0)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND1, DOY1)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND2, DOY2)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND3, DOY3)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND1, HM1)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND2, HM2)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND3, HM3)

C Check whether point is before start or beyond end of file
         DOYD = FDOY - DOY1
         HOURD = EC_NCDF_CSI2HOUR(DBLE(FHOURMIN)) - EC_NCDF_CSI2HOUR(HM1)
         TIMED = 3600.0*24*DOYD + HOURD*3600.0
         IF (TIMED .LT. 0) THEN
            FDOY = DOY1
	    FHOURMIN = HM1
         ENDIF
         DOYD = DOY3 - FDOY
         HOURD = EC_NCDF_CSI2HOUR(HM3) - EC_NCDF_CSI2HOUR(DBLE(FHOURMIN))
         TIMED = 3600.0*24*DOYD + HOURD*3600.0
         IF (TIMED .LT. 0) THEN
            FDOY = DOY3
	    FHOURMIN = HM3
         ENDIF

         CURCH = IND2
         ATTEMPTS = 0

         DO WHILE ((ABS(CURCH) .GT. 0) .AND. (ATTEMPTS .LT. 50))
C Determine time difference between current sample and required time
            PREVCH = CURCH
            DOYD = FDOY - DOY2
            HOURD = EC_NCDF_CSI2HOUR(DBLE(FHOURMIN)) - EC_NCDF_CSI2HOUR(HM2)
            TIMED = 3600.0*24*DOYD + HOURD*3600.0
            IF (ANINT(TIMED) .LT. 0) THEN
	       CURCH = IND3 - IND2
               IND3 = IND2
	       HM3 = HM2
	       DOY3 = DOY2
	    ELSE IF (ANINT(TIMED) .GT. 0) THEN
	       CURCH = IND1 - IND2
               IND1 = IND2
	       HM1 = HM2
	       DOY1= DOY2
	    ELSE
	       IND4=IND2-1
               STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND4, DOY4)
               STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND4, HM4)
	       IF (ABS(ANINT(HM2 - HM4)) .GT. 0) THEN
	         IND = IND2
	         RETURN
	       ELSE
	         CURCH = IND2 - IND3
                 IND3 = IND2
	         HM3 = HM2
	         DOY3 = DOY2
	       ENDIF
            ENDIF
            IND2 = ANINT((IND3+IND1)/2.0)
            STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND2, DOY2)
            STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND2, HM2)


	    IF (CURCH .EQ. 0) THEN
	       IND = IND2
	       RETURN
	    ENDIF
            ATTEMPTS = ATTEMPTS + 1
         ENDDO
      ELSE
         IND = -1
      ENDIF
      IF (ATTEMPTS .EQ. 50) THEN
         IND = -1
      ENDIF
      IF (IND2 .LT. 1) IND2 = 1
      IF (IND2 .GT. MAXIND) IND2 = MAXIND
      IF (ABS(TIMED) .GT. 1.5/60D0) IND2 = -1
      IND = IND2
      END

C...........................................................................
C SUBROUTINE:  EC_NCDF_FINDTIME
C PURPOSE   :  To find the index in a NetCDF File for a given DOY and time
C...........................................................................
      SUBROUTINE EC_NCDF_FINDTIME(FDOY, FHOURMIN, NCID,
     +                    NCdoyID, NChmID, NCsecID,IND)
      INTEGER FDOY, FHOURMIN, NCID, NCdoyID, NChmID, NCsecID, IND,
     +        DIMID, MAXIND, ATTEMPTS, PREVCH, IND1, IND2, IND3,
     +        CURCH
      REAL*8 SEC1, SEC2, SEC3, DOY1, DOY2, DOY3, HM1, HM2, HM3,
     +     SAMPF, TIME1, TIME2, TIME3, SEC4, HM4, DOY4,
     +     DOYD, HOURD, SECD, TIMED
      LOGICAL OK

      REAL*8 EC_NCDF_CSI2HOUR

      STATUS = NF_INQ_VARDIMID(NCID,NCsecID,DIMID)
      STATUS = NF_INQ_DIMLEN(NCID, DIMID,MAXIND )
      IF (MAXIND .GT. 0) THEN
         IND1 = 1
         IND3 = MAXIND
         OK = .FALSE.
         DO WHILE (.NOT. OK)
           STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND1, DOY1)
           STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND1, HM1)
	   IF (((DOY1 .GT. 0) .AND. (DOY1 .LE. 366) .AND.
     +          (HM1 .GE. 0) .AND. (HM1 .LE. 2400)) .OR.
     +         (IND1 .EQ. MAXIND)) THEN
	      OK = .TRUE.
	   ELSE
	      IND1 = IND1 + 1
	   ENDIF
         ENDDO
         OK = .FALSE.
         DO WHILE (.NOT. OK)
           STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND3, DOY3)
           STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND3, HM3)
	   IF (((DOY3 .GT. 0) .AND. (DOY3 .LE. 366) .AND.
     +          (HM3 .GE. 0) .AND. (HM3 .LE. 2400)) .OR.
     +         (IND3 .EQ. 1)) THEN
	      OK = .TRUE.
	   ELSE
	      IND3 = IND3 - 1
	   ENDIF
         ENDDO

         STATUS = NF_GET_VAR1_DOUBLE(NCID,NCsecID, IND1, SEC1)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NCsecID, IND3, SEC3)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND1, DOY1)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND3, DOY3)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND1, HM1)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND3, HM3)
C Get mean sample rate
         DOYD = DOY3 - DOY1
         HOURD = EC_NCDF_CSI2HOUR(HM3) - EC_NCDF_CSI2HOUR(HM1)
         SECD = SEC3 - SEC1
         TIMED = 3600.0*24*DOYD + HOURD*3600.0 + SECD
         SAMPF = MAXIND/ABS(TIMED)
C Check whether point is before start or beyond end of file
         DOYD = FDOY - DOY1
         HOURD = EC_NCDF_CSI2HOUR(DBLE(FHOURMIN)) - EC_NCDF_CSI2HOUR(HM1)
         SECD = 0 - SEC1
         TIMED = 3600.0*24*DOYD + HOURD*3600.0 + SECD
         IF (TIMED .LT. 0) THEN
            FDOY = DOY1
	    FHOURMIN = HM1
         ENDIF
         DOYD = DOY3 - FDOY
         HOURD = EC_NCDF_CSI2HOUR(HM3) - EC_NCDF_CSI2HOUR(DBLE(FHOURMIN))
         SECD = SEC3 - 0
         TIMED = 3600.0*24*DOYD + HOURD*3600.0 + SECD
         IF (TIMED .LT. 0) THEN
            FDOY = DOY3
	    FHOURMIN = HM3
         ENDIF

C Get time difference between start and requested point
         DOYD = FDOY - DOY1
         HOURD = EC_NCDF_CSI2HOUR(DBLE(FHOURMIN)) - EC_NCDF_CSI2HOUR(HM1)
         SECD = 0 - SEC1
         TIMED = 3600.0*24*DOYD + HOURD*3600.0 + SECD
         ATTEMPTS = 0
         IND2 = TIMED*SAMPF+1
         IF (IND2 .LT. 1) THEN
            IND2 = 1
         ELSE IF (IND2 .GT. MAXIND) THEN
            IND2 = MAXIND
         ENDIF
         CURCH = IND2

         DO WHILE ((ABS(CURCH) .GT. 0) .AND. (ATTEMPTS .LT. 20))
C Determine time difference between current sample and required time
            PREVCH = CURCH
            STATUS = NF_GET_VAR1_DOUBLE(NCID,NCsecID, IND2, SEC2)
            STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND2, DOY2)
            STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND2, HM2)
            STATUS = NF_GET_VAR1_DOUBLE(NCID,NCsecID, IND2+1, SEC4)
            STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND2+1, DOY4)
            STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND2+1, HM4)
C Get local sample rate
            DOYD = DOY4 - DOY2
            HOURD = EC_NCDF_CSI2HOUR(HM4) - EC_NCDF_CSI2HOUR(HM2)
            SECD = SEC4 - SEC2
            TIMED = 3600.0*24*DOYD + HOURD*3600.0 + SECD
	    SAMPF = 1D0/TIMED

C Get distance between current point and requested point
            DOYD = FDOY - DOY2
            HOURD = EC_NCDF_CSI2HOUR(DBLE(FHOURMIN)) - EC_NCDF_CSI2HOUR(HM2)
            SECD = 0D0 - SEC2
            TIMED = 3600.0*24*DOYD + HOURD*3600.0 + SECD
            CURCH = ANINT(TIMED*SAMPF)
C We're there
	    IF (CURCH .EQ. 0) THEN
	       IND = IND2
	       RETURN
C Shift forward, but not beyond end of file
	    ELSE IF (TIMED .GT. 0) THEN
	        IF (IND2 .EQ. MAXIND) THEN
	           IND = IND2
	   	   RETURN
	        ENDIF
	        IND1 = IND2
	        HM1 = HM2
	        SEC1 = SEC2
	        DOY1 = DOY2
                IND2 = IND2 + ANINT(TIMED*SAMPF)
C Shift backward, but not beyond start of file
	    ELSE IF  (TIMED .LT. 0) THEN
	        IF (IND2 .EQ. 1) THEN
	           IND = IND2
	   	   RETURN
	        ENDIF
	        IND3 = IND2
	        HM3 = HM2
	        SEC3 = SEC2
	        DOY3 = DOY2
                IND2 = IND2 + ANINT(TIMED*SAMPF)
	   ENDIF
C
C Get time difference between start of interval and requested point
           ATTEMPTS = ATTEMPTS + 1
         ENDDO
      ELSE
         IND = -1
      ENDIF

      IF (ATTEMPTS .EQ. 20) THEN
         IND = -1
      ENDIF
      IF (IND2 .LT. 1) IND2 = 1
      IF (IND2 .GT. MAXIND) IND2 = MAXIND
      IF (ABS(TIMED) .GT. 1./60D0) IND2 = -1
      IND = IND2
      END





      SUBROUTINE EC_NCDF_HANDLE_ERR(STATUS)
C
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
