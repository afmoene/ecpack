C  ec_text.f, part of ecpack
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
C Routine   : EC_T_CLEARSTR
C Purpose   : To set all characters in a string to CHAR(0)
C Interface : STRING      intent(INOUT)   string to be changed
C Author    : Arnold Moene
C Date      : May 25, 2000
C...........................................................................
      SUBROUTINE EC_T_CLEARSTR(STRING)

      IMPLICIT NONE
      CHARACTER*(*) STRING
      INTEGER       I

      DO 15, I=1,LEN(STRING)
         STRING(I:I) = CHAR(0)
   15 CONTINUE
      END

C...........................................................................
C Routine   : stripstr
C Purpose   : to replace all trailing spaces in a string by NULL characters
C             and strip all leading spaces. Furthermore, remove
C             opening and closing quotes
C Interface : STRING      intent(INOUT)   string to be checked
C Author    : Arnold Moene
C Date      : May 24, 2000
C...........................................................................
      SUBROUTINE EC_T_STRIPSTR(STRING)

      IMPLICIT NONE
      CHARACTER*(*) STRING
      INTEGER       I, J, K, CHANGED, EC_T_STRLEN
      EXTERNAL EC_T_STRLEN

      CHARACTER     QUOTE(2)
C Remove trailing spaces
      DO 15, I = LEN(STRING), 1, -1
        IF (STRING(I:I) .NE. ' ') GO TO 20
        IF (STRING(I:I) .EQ. ' ') STRING(I:I) = CHAR(0)
   15 CONTINUE
   20 CONTINUE

C Remove leading spaces
      CHANGED = 1      
      I = 1
      DO 30, WHILE (CHANGED .EQ. 1)
        CHANGED = 0
        IF ((STRING(I:I) .EQ. ' ') .OR.
     &      (STRING(I:I) .EQ. CHAR(0))) THEN
           DO 35, J = I+1,LEN(STRING)-1
              STRING(J-1:J-1) = STRING(J:J)
   35      CONTINUE
           CHANGED = 1
        ELSE
           I = I + 1
        ENDIF
   30 CONTINUE
C Remove outer quotes
      QUOTE(1) = "'"
      QUOTE(2) = '"'

      DO 75, K = 1,2
         IF (STRING(1:1) .EQ. QUOTE(K)) THEN
              DO 55, I = LEN(STRING),1,-1
                 IF (STRING(I:I) .EQ. QUOTE(K)) THEN
                     STRING(I:I) = CHAR(0)
                    DO 65, J = 1,LEN(STRING)
                        STRING(J:J) = STRING(J+1:J+1)
   65               CONTINUE
                 ENDIF
   55         CONTINUE
         ENDIF
   75 CONTINUE
      END


C...........................................................................
C Function  : EC_T_STRLEN
C Purpose   : to return the length of a string
C Interface : STRING      intent(IN)      string to be checked
C             return value                length of string (integer)
C Author    : after function in book of Clive Page (Professional programmer's
C             guide to Fortran 77)
C Date      : May 24, 2000
C...........................................................................
      INTEGER FUNCTION EC_T_STRLEN(STRING)

      IMPLICIT NONE
      CHARACTER*(*) STRING
      INTEGER       I

      DO 15, I = LEN(STRING), 1, -1
        IF ((STRING(I:I) .NE. ' ') .AND.
     &      (STRING(I:I) .NE. CHAR(0)))
     &     GO TO 20
   15 CONTINUE
   20 EC_T_STRLEN = I
      END

C...........................................................................
C routine   : strcat
C Purpose   : to concatenate a string to another
C Interface : STRING1     intent(INOUT)   string to be added to (supposed
C                                         to be of sufficient length)
C             STRING2     intent(IN)      string to add
C Author    : Arnold Moene
C Date      : May 24, 2000
C...........................................................................
      SUBROUTINE EC_T_STRCAT(STRING1, STRING2)

      IMPLICIT NONE
      CHARACTER*(*) STRING1, STRING2
      INTEGER       I, LEN1
      INTEGER       EC_T_STRLEN
      EXTERNAL      EC_T_STRLEN


      LEN1 = EC_T_STRLEN(STRING1)
      DO 15, I = 1, EC_T_STRLEN(STRING2)
        STRING1(LEN1+I:LEN1+I) = STRING2(I:I)
   15 CONTINUE
      END
      
C...........................................................................
C FUNCTION   : EC_T_EQSTRING
C Purpose   : To test whether strings are equal
C Interface : STRING1      intent(IN)   string to be checked
C             STRING2      intent(IN)   string to be checked
C             return value : LOGICAL
C Author    : Arnold Moene
C Date      : January 30, 2003
C...........................................................................
      LOGICAL FUNCTION EC_T_EQSTRING(STRING1, STRING2)
      IMPLICIT NONE
      
      CHARACTER*(*) STRING1, STRING2
      INTEGER       EC_T_STRLEN
      EXTERNAL      EC_T_STRLEN
      
      EC_T_EQSTRING = 
     &    ((INDEX(STRING1(:Ec_T_STRLEN(STRING1)),
     &            STRING2(:EC_T_STRLEN(STRING2))) .GT. 0) .AND.
     &       (EC_T_STRLEN(STRING1) .EQ. EC_T_STRLEN(STRING2))) 
      END
C...........................................................................
C Routine   : EC_T_UPCASE
C Purpose   : To convert a string to upper case
C Interface : STRING      intent(INOUT)   string to be changed
C Author    : Arnold Moene
C Date      : May 24, 2000
C...........................................................................
      SUBROUTINE EC_T_UPCASE(STRING)

      IMPLICIT NONE
      CHARACTER*(*) STRING
      INTEGER       I

      DO 15, I=1,LEN(STRING)
         IF ((ICHAR(STRING(I:I)) .GT. 96) .AND.
     &       (ICHAR(STRING(I:I)) .LT. 124))
     &       STRING(I:I) = CHAR(ICHAR(STRING(I:I))-32)
   15 CONTINUE
      END

C...........................................................................
C Routine   : EC_T_GOut1
C Purpose   : To set flag to give output of statistic of one variable
C Interface : STRING      intent(IN)      string to be parsed
C             FLAGS       intent(INOUT)   array of flags to be set
C Author    : Arnold Moene
C Date      : January 29, 2003
C...........................................................................
      SUBROUTINE EC_T_GOut1(String, Flags)
      IMPLICIT NONE
      include 'parcnst.inc'
      
      CHARACTER*(*) STRING
      CHARACTER*100 DUMSTRING
      LOGICAL       FLAGS(NNMax), OneFound
      INTEGER       I
      
      INTEGER EC_T_STRLEN
      LOGICAL EC_T_EQSTRING
      EXTERNAL EC_T_STRLEN, EC_T_EQSTRING
      

      CALL EC_T_STRIPSTR(STRING)
      CALL EC_T_UPCASE(STRING)
      OneFound = .FALSE.
      DO I=1,NNMAX
	 DUMSTRING = QName(I)
	 CALL EC_T_STRIPSTR(DUMSTRING)
         CALL EC_T_UPCASE(DUMSTRING)
         IF (EC_T_EQSTRING(STRING, DUMSTRING)) THEN
	     Flags(I) = .TRUE.
	     OneFound = .TRUE.
         ENDIF
         IF (EC_T_EQSTRING(STRING, 'ALL')) THEN
	     Flags(I) = .TRUE.
	     OneFound = .TRUE.
         ENDIF
      ENDDO
            
      IF (.NOT. OneFound) THEN
         WRITE(*,*) 'Cannot find variable ', STRING
         STOP
      ENDIF
      END
      
C...........................................................................
C Routine   : EC_T_GOut2
C Purpose   : To set flag to give output of statistic of two variables
C Interface : STRING      intent(IN)      string to be parsed
C             FLAGS       intent(INOUT)   array of flags to be set
C Author    : Arnold Moene
C Date      : January 29, 2003
C...........................................................................
      SUBROUTINE EC_T_GOut2(String, Flags)
      IMPLICIT NONE
      include 'parcnst.inc'
      
      CHARACTER*(*) STRING
      CHARACTER*255 DUMSTRING1, DUMSTRING2,  PART1, PART2
      LOGICAL       FLAGS(NNMax, NNMax), OneFound
      INTEGER       I, J, KINDEX
      
      INTEGER EC_T_STRLEN
      LOGICAL EC_T_EQSTRING
      EXTERNAL EC_T_STRLEN, EC_T_EQSTRING
      
      CALL EC_T_STRIPSTR(STRING)
      CALL EC_T_UPCASE(STRING)
C     Find the separating comma
      KINDEX = INDEX(STRING, ',')
      IF (KINDEX .EQ. 0) THEN
          WRITE(*,*) 'ERROR no comma sign found in ', STRING
          STOP
      ENDIF
c     Split into token and value
      WRITE(PART1,*) STRING(:KINDEX-1)
      WRITE(PART2,*) STRING(KINDEX+1:)
      CALL EC_T_STRIPSTR(PART1)
      CALL EC_T_STRIPSTR(PART2)
      CALL EC_T_UPCASE(PART1)
      CALL EC_T_UPCASE(PART2)

      OneFound = .FALSE.
      DO I=1,NNMax
         DUMSTRING1=QName(I)
         CALL EC_T_STRIPSTR(DUMSTRING1)
         CALL EC_T_UPCASE(DUMSTRING1)
         IF (EC_T_EQSTRING(PART1, DUMSTRING1) .OR.
     &       EC_T_EQSTRING(PART1, 'ALL')) THEN
             DO J=1,NNMax
	        DUMSTRING2=QName(J)
                CALL EC_T_STRIPSTR(DUMSTRING2)
                CALL EC_T_UPCASE(DUMSTRING2)
                IF (EC_T_EQSTRING(PART2, DUMSTRING2) .OR.
     &              EC_T_EQSTRING(PART2, 'ALL')) THEN
	           Flags(I,J) = .TRUE.
	           OneFound = .TRUE.
                ENDIF
	     ENDDO
         ENDIF
      ENDDO
      
      IF (.NOT. OneFound) THEN
         WRITE(*,*) 'Cannot find variable combination', STRING
      ENDIF
     
      END      
      
C...........................................................................
C Routine   : EC_T_GPhys
C Purpose   : To set flag to give output of a physical quantity
C Interface : STRING      intent(IN)      string to be parsed
C             FLAGS       intent(INOUT)   array of flags to be set
C Author    : Arnold Moene
C Date      : January 29, 2003
C...........................................................................
      SUBROUTINE EC_T_GPhys(String, Flags)
      IMPLICIT NONE
      include 'parcnst.inc'
      
      CHARACTER*(*) STRING
      CHARACTER*100 DUMSTRING, ALTDUMSTR
      LOGICAL       FLAGS(NMaxPhys), OneFound
      INTEGER       I
      INTEGER EC_T_STRLEN
      LOGICAL EC_T_EQSTRING
      EXTERNAL EC_T_STRLEN, EC_T_EQSTRING      

      CALL EC_T_STRIPSTR(STRING)
      CALL EC_T_UPCASE(STRING)
      OneFound = .FALSE.
      DO I=1,NMaxPhys
	 DUMSTRING = QPName(I)
	 ALTDUMSTR = QPAltN(I)
	 CALL EC_T_STRIPSTR(DUMSTRING)
	 CALL EC_T_STRIPSTR(ALTDUMSTR)
         CALL EC_T_UPCASE(DUMSTRING)
         CALL EC_T_UPCASE(ALTDUMSTR)
         IF ((EC_T_EQSTRING(DUMSTRING, STRING)) .OR.
     &       (EC_T_EQSTRING(ALTDUMSTR, STRING))) THEN
	     Flags(I) = .TRUE.
	     OneFound = .TRUE.
         ENDIF
         IF (EC_T_EQSTRING('ALL', STRING)) THEN
	     Flags(I) = .TRUE.
	     OneFound = .TRUE.
         ENDIF
      ENDDO
            
      IF (.NOT. OneFound) THEN
         WRITE(*,*) 'Cannot find variable ', STRING
         STOP
      ENDIF
      END
 
C...........................................................................
C Routine   : EC_T_GOutS
C Purpose   : To set flag to give info about output type
C Interface : STRING      intent(IN)      string to be parsed
C             FLAGS       intent(INOUT)   array of flags to be set
C Author    : Arnold Moene
C Date      : January 29, 2003
C...........................................................................
      SUBROUTINE EC_T_GOuts(String, Flags)
      IMPLICIT NONE
      include 'parcnst.inc'
      
      CHARACTER*(*) STRING
      LOGICAL       FLAGS(NMaxOS), OneFound
      
      INTEGER EC_T_STRLEN
      EXTERNAL EC_T_STRLEN

      CALL EC_T_STRIPSTR(STRING)
      CALL EC_T_UPCASE(STRING)
      OneFound = .FALSE.
C
C Old/new output
C
      IF (INDEX('OLD',
     &          STRING(:EC_T_STRLEN(STRING))) .GT. 0) THEN
         FLAGS(OSOld) = .TRUE.
         OneFound = .TRUE.
      ELSE IF (INDEX('NEW',
     &          STRING(:EC_T_STRLEN(STRING))) .GT. 0) THEN
         FLAGS(OSOld) = .FALSE.
         OneFound = .TRUE.
      ENDIF
C
C Tolerances yes/no
C
      IF (INDEX('TOL',
     &          STRING(:EC_T_STRLEN(STRING))) .GT. 0) THEN
         FLAGS(OSTol) = .TRUE.
         OneFound = .TRUE.
      ELSE IF (INDEX('NOTOL',
     &          STRING(:EC_T_STRLEN(STRING))) .GT. 0) THEN
         FLAGS(OSTOL) = .FALSE.
         OneFound = .TRUE.
      ENDIF
C
C CSV or SSV output yes/no
C
      IF (INDEX('CSV',
     &          STRING(:EC_T_STRLEN(STRING))) .GT. 0) THEN
         FLAGS(OScsv) = .TRUE.
         OneFound = .TRUE.
      ELSE IF (INDEX('SSV',
     &          STRING(:EC_T_STRLEN(STRING))) .GT. 0) THEN
         FLAGS(OScsv) = .FALSE.
         OneFound = .TRUE.
      ENDIF
C
C Diagnoaitics output yes/no
C
      IF (INDEX('DIAG',
     &          STRING(:EC_T_STRLEN(STRING))) .GT. 0) THEN
         FLAGS(OSDiag) = .TRUE.
         OneFound = .TRUE.
      ELSE IF (INDEX('NODIAG',
     &          STRING(:EC_T_STRLEN(STRING))) .GT. 0) THEN
         FLAGS(OSDiag) = .FALSE.
         OneFound = .TRUE.
      ENDIF
C
C Version number in output file?
C
      IF (INDEX('VERSION',
     &          STRING(:EC_T_STRLEN(STRING))) .GT. 0) THEN
         FLAGS(OSVers) = .TRUE.
         OneFound = .TRUE.
      ELSE IF (INDEX('NOVERSION',
     &          STRING(:EC_T_STRLEN(STRING))) .GT. 0) THEN
         FLAGS(OSVers) = .FALSE.
         OneFound = .TRUE.
      ENDIF
      
            
      IF (.NOT. OneFound) THEN
         WRITE(*,*) 'Cannot find variable ', STRING
         STOP
      ENDIF
      END
 
