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
                    DO 65, J = I+1,LEN(STRING)-1
                       STRING(J-1:J-1) = STRING(J:J)
   65               CONTINUE
                 ENDIF
   55         CONTINUE
         ENDIF
   75 CONTINUE
      END


C...........................................................................
C Function  : strlen
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
