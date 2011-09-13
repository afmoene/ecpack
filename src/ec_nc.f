C  ec_nc.f, part of ecpack
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



C...........................................................................
C FUNCTION:  EC_NCDF_CSI2HOUR
C PURPOSE:   To convert Campbell hour/min to decimal hours
C INTERFACE: HOURMIN    IN     Campbell hour/minutes (REAL)
C            return value      decimal hours (REAL)
C AUTHOR:    Arnold Moene
C DATE:      May 28, 2001
C...........................................................................
      REAL*8 FUNCTION EC_NCDF_CSI2HOUR(HOURMIN)
      IMPLICIT NONE
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
      IMPLICIT NONE     
      INCLUDE 'netcdf.inc'
      
      INTEGER FDOY, FHOURMIN, NCID, NCdoyID, NChmID, IND,
     +        DIMID, MAXIND, ATTEMPTS, IND1, IND2, IND3, IND4,
     +        CURCH, STATUS
      REAL*8 DOY1, DOY2, DOY3, HM1, HM2, HM3, HM4, DOY4,
     +       DOYD, HOURD, TIMED
      REAL*8 EC_NCDF_CSI2HOUR
      LOGICAL OK
      


C Initialize 
      STATUS = NF_INQ_VARDIMID(NCID,NChmID,DIMID)
      STATUS = NF_INQ_DIMLEN(NCID, DIMID,MAXIND )
      ATTEMPTS = 0

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
         IND2 = INT((IND3+IND1)/2.0)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND1, DOY1)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND2, DOY2)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NCdoyID, IND3, DOY3)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND1, HM1)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND2, HM2)
         STATUS = NF_GET_VAR1_DOUBLE(NCID,NChmID, IND3, HM3)

C Check whether point is before start or beyond end of file
         DOYD = FDOY - DOY1
         HOURD = EC_NCDF_CSI2HOUR(DBLE(FHOURMIN)) - 
     +           EC_NCDF_CSI2HOUR(HM1)
         TIMED = 3600.0D0*24.0D0*DOYD + HOURD*3600.0D0
         IF (TIMED .LT. 0) THEN
c           Data is before start of file so in this file the interval might start at first sample
c           FDOY = INT(DOY1)
c	    FHOURMIN = INT(HM1)
            IND = IND1
            RETURN
         ENDIF
         DOYD = DOY3 - FDOY
         HOURD = EC_NCDF_CSI2HOUR(HM3) - 
     +           EC_NCDF_CSI2HOUR(DBLE(FHOURMIN))
         TIMED = 3600.0D0*24.0D0*DOYD + HOURD*3600.0D0
         IF (TIMED .LT. 0) THEN
c           Data is after end of file so in this file the interval might end at last sample
c           FDOY = INT(DOY3)
c           FHOURMIN = INT(HM3)
            IND = IND3
            RETURN
         ENDIF

         CURCH = IND2

         DO WHILE ((ABS(CURCH) .GT. 0) .AND. (ATTEMPTS .LT. 50))
C Determine time difference between current sample and required time
            
            DOYD = FDOY - DOY2
            HOURD = EC_NCDF_CSI2HOUR(DBLE(FHOURMIN)) - 
     +              EC_NCDF_CSI2HOUR(HM2)
            TIMED = 3600.0D0*24.0D0*DOYD + HOURD*3600.0D0
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
            IND2 = INT((IND3+IND1)/2.0)
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
      IF (ABS(TIMED) .GT. 1.5D0/60D0) IND2 = -1
      IND = IND2
      END

C...........................................................................
C SUBROUTINE:  EC_NCDF_FINDTIME
C PURPOSE   :  To find the index in a NetCDF File for a given DOY and time
C...........................................................................
      SUBROUTINE EC_NCDF_FINDTIME(FDOY, FHOURMIN, NCID,
     +                    NCdoyID, NChmID, NCsecID,IND)
      IMPLICIT None
      INCLUDE 'netcdf.inc'
      
      INTEGER FDOY, FHOURMIN, NCID, NCdoyID, NChmID, NCsecID, IND,
     +        DIMID, MAXIND, ATTEMPTS, IND1, IND2, IND3,
     +        CURCH, STATUS
      REAL*8 SEC1, SEC2, SEC3, DOY1, DOY2, DOY3, HM1, HM2, HM3,
     +     SAMPF, SEC4, HM4, DOY4,
     +     DOYD, HOURD, SECD, TIMED
      LOGICAL OK

      REAL*8 EC_NCDF_CSI2HOUR

C Initialize
      STATUS = NF_INQ_VARDIMID(NCID,NCsecID,DIMID)
      STATUS = NF_INQ_DIMLEN(NCID, DIMID,MAXIND )
      ATTEMPTS = 0

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
         TIMED = 3600.0D0*24.0D0*DOYD + HOURD*3600.0D0 + SECD
         SAMPF = MAXIND/ABS(TIMED)
C Check whether point is before start or beyond end of file
         DOYD = FDOY - DOY1
         HOURD = EC_NCDF_CSI2HOUR(DBLE(FHOURMIN)) - 
     +           EC_NCDF_CSI2HOUR(HM1)
         SECD = 0 - SEC1
         TIMED = 3600.0D0*24.0D0*DOYD + HOURD*3600.0D0 + SECD
         IF (TIMED .LT. 0) THEN
            FDOY = INT(DOY1)
	    FHOURMIN = INT(HM1)
         ENDIF
         DOYD = DOY3 - FDOY
         HOURD = EC_NCDF_CSI2HOUR(HM3) - 
     +           EC_NCDF_CSI2HOUR(DBLE(FHOURMIN))
         SECD = SEC3 - 0
         TIMED = 3600.0D0*24.0D0*DOYD + HOURD*3600.0D0 + SECD
         IF (TIMED .LT. 0) THEN
            FDOY = INT(DOY3)
	    FHOURMIN = INT(HM3)
         ENDIF

C Get time difference between start and requested point
         DOYD = FDOY - DOY1
         HOURD = EC_NCDF_CSI2HOUR(DBLE(FHOURMIN)) - 
     +          EC_NCDF_CSI2HOUR(HM1)
         SECD = 0 - SEC1
         TIMED = 3600.0D0*24.0D0*DOYD + HOURD*3600.0D0 + SECD
         ATTEMPTS = 0
         IND2 = INT(TIMED*SAMPF+1)
         IF (IND2 .LT. 1) THEN
            IND2 = 1
         ELSE IF (IND2 .GT. MAXIND) THEN
            IND2 = MAXIND
         ENDIF
         CURCH = IND2

         DO WHILE ((ABS(CURCH) .GT. 0) .AND. (ATTEMPTS .LT. 20))
C Determine time difference between current sample and required time
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
            TIMED = 3600.0D0*24.0D0*DOYD + HOURD*3600.0D0 + SECD
	    SAMPF = 1D0/TIMED

C Get distance between current point and requested point
            DOYD = FDOY - DOY2
            HOURD = EC_NCDF_CSI2HOUR(DBLE(FHOURMIN)) - 
     +              EC_NCDF_CSI2HOUR(HM2)
            SECD = 0D0 - SEC2
            TIMED = 3600.0D0*24.0D0*DOYD + HOURD*3600.0D0 + SECD
            CURCH = INT(TIMED*SAMPF)
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
                IND2 = IND2 + INT(TIMED*SAMPF)
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
                IND2 = IND2 + INT(TIMED*SAMPF)
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
      IF (ABS(TIMED) .GT. 1.0D0/60D0) IND2 = -1
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
         WRITE(*,*) 'NetCDF lib error: ', NF_STRERROR(STATUS)
         STOP 'Stopped'
      ENDIF
      END
