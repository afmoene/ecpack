C $Id$
      SUBROUTINE EC_F_GetConf(ConfUnit,
     &           DatDir, OutDir,  ParmDir,
     &           FluxName, ParmName, InterName,
     &           PlfName,
     &           SonName, CoupName, HygName, NCVarName,
     &           NCNameLen)

      INCLUDE 'physcnst.inc'
      INCLUDE 'parcnst.inc'
      
      
      INTEGER         ConfUnit, NCNameLen
      CHARACTER*(*)   DatDir, OutDir, ParmDir, FluxName, ParmName,
     &                InterName, SonName, CoupName, HygName, PlfName
      CHARACTER*(*)   NCVarName(NCNameLen)

      INTEGER         IOCODE, KINDEX
      CHARACTER*255   LINE, TOKLINE, VALLINE, DUMSTRING
      CHARACTER       ONECHAR
      INTEGER         STRLEN
      EXTERNAL        STRLEN


      DO 1000, I=1,NCNameLen
         CALL CLEARSTR(NCVarName(I))
 1000 CONTINUE
      CALL CLEARSTR(DatDir)
      CALL CLEARSTR(OutDir)
      CALL CLEARSTR(ParmDir)
      CALL CLEARSTR(FluxName)
      CALL CLEARSTR(ParmName)
      CALL CLEARSTR(InterName)
      CALL CLEARSTR(SonName)
      CALL CLEARSTR(CoupName)
      CALL CLEARSTR(HygName)
      CALL CLEARSTR(PlfName)
      
      IOCODE = 0
      DO 4000, WHILE (IOCODE .EQ. 0)
C     Read one line from configuration file
         READ(Confunit, '(A)', IOSTAT=IOCODE, END=9000) LINE
         IF (IOCODE .NE. 0) THEN
            WRITE(*,*) 'ERROR in reading of configuration file'
            STOP
         ENDIF
         CALL STRIPSTR(LINE)
C     Check for comment (should start with //
         KINDEX = INDEX(LINE, '//')
         IF (KINDEX .EQ. 0) THEN
C     Find the equality sign
            KINDEX = INDEX(LINE, '=')
            IF (KINDEX .EQ. 0) THEN
                WRITE(*,*) 'ERROR no equality sign found'
                STOP
            ENDIF
c     Split into token and value
            WRITE(TOKLINE,*) LINE(:KINDEX-1)
            WRITE(VALLINE,*) LINE(KINDEX+1:)
            CALL STRIPSTR(TOKLINE)
            CALL STRIPSTR(VALLINE)
            CALL UPCASE(TOKLINE)
C     See which token this is
            IF (INDEX(TOKLINE, 'DATDIR') .GT. 0) THEN
               DatDir = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'OUTDIR') .GT. 0) THEN
               OutDir = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'PARMDIR') .GT. 0) THEN
               ParmDir = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'FLUXNAME') .GT. 0) THEN
               FluxName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'PARMNAME') .GT. 0) THEN
               ParmName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'INTERNAME') .GT. 0) THEN
               InterName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'SONNAME') .GT. 0) THEN
               SonName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'COUPNAME') .GT. 0) THEN
               CoupName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'HYGNAME') .GT. 0) THEN
               HygName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'PLFNAME') .GT. 0) THEN
               PlfName = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
C NetCdF variable names
            ELSE IF (INDEX(TOKLINE, 'U_VAR') .GT. 0) THEN
               NCVarName(U) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'V_VAR') .GT. 0) THEN
               NCVarName(V) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'W_VAR') .GT. 0) THEN
               NCVarName(W) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'TSONIC_VAR') .GT. 0) THEN
               NCVarName(TSonic) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'TCOUPLE_VAR') .GT. 0) THEN
               NCVarName(TCouple) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'HUMIDITY_VAR') .GT. 0) THEN
               NCVarName(HUMIDITY) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'DOY_VAR') .GT. 0) THEN
               NCVarName(Doy) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'HOURMIN_VAR') .GT. 0) THEN
               NCVarName(HourMin) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'SEC_VAR') .GT. 0) THEN
               NCVarName(Sec) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'DIAG_VAR') .GT. 0) THEN
               NCVarName(Diagnost) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ELSE IF (INDEX(TOKLINE, 'TREF_VAR') .GT. 0) THEN
               NCVarName(Tref) = VALLINE(:INDEX(VALLINE, CHAR(0))-1)
            ENDIF
         ENDIF
 4000 CONTINUE
 9000 CONTINUE

      CALL CLEARSTR(DUMSTRING)
      IF (STRLEN(PlfName) .GT. 0) THEN
         DUMSTRING = OutDir
         CALL STRCAT(DUMSTRING,'/')
         CALL STRCAT(DUMSTRING,PlfName)
      ENDIF
      PlfName = DUMSTRING

      CALL CLEARSTR(DUMSTRING)
      IF (STRLEN(FluxName) .GT. 0) THEN
         DUMSTRING = OutDir
         CALL STRCAT(DUMSTRING,'/')
         CALL STRCAT(DUMSTRING,FluxName)
      ENDIF
      FluxName = DUMSTRING

      CALL CLEARSTR(DUMSTRING)
      IF (STRLEN(ParmName) .GT. 0) THEN
         DUMSTRING = ParmDir
         CALL STRCAT(DUMSTRING,'/')
         CALL STRCAT(DUMSTRING,ParmName)
      ENDIF
      ParmName = DUMSTRING

      CALL CLEARSTR(DUMSTRING)
      IF (STRLEN(InterName) .GT. 0) THEN
         DUMSTRING = ParmDir
         CALL STRCAT(DUMSTRING,'/')
         CALL STRCAT(DUMSTRING,InterName)
      ENDIF
      InterName = DUMSTRING
      
      CALL CLEARSTR(DUMSTRING)
      IF (STRLEN(SonName) .GT. 0) THEN
         DUMSTRING = ParmDir
         CALL STRCAT(DUMSTRING,'/')
         CALL STRCAT(DUMSTRING,SonName)
      ENDIF
      SonName = DUMSTRING
      
      CALL CLEARSTR(DUMSTRING)
      IF (STRLEN(CoupName) .GT. 0) THEN
         DUMSTRING = ParmDir
         CALL STRCAT(DUMSTRING,'/')
         CALL STRCAT(DUMSTRING,CoupName)
      ENDIF
      CoupName = DUMSTRING
      
      CALL CLEARSTR(DUMSTRING)
      IF (STRLEN(HygName) .GT. 0) THEN
         DUMSTRING = ParmDir
         CALL STRCAT(DUMSTRING,'/')
         CALL STRCAT(DUMSTRING,HygName)
      ENDIF
      HygName = DUMSTRING


 5000 FORMAT(A)

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
      INTEGER FUNCTION STRLEN(STRING)

      IMPLICIT NONE
      CHARACTER*(*) STRING
      INTEGER       I

      DO 15, I = LEN(STRING), 1, -1
        IF ((STRING(I:I) .NE. ' ') .AND.
     &      (STRING(I:I) .NE. CHAR(0)))
     &     GO TO 20
   15 CONTINUE
   20 STRLEN = I
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
      SUBROUTINE STRCAT(STRING1, STRING2)

      IMPLICIT NONE
      CHARACTER*(*) STRING1, STRING2
      INTEGER       I, LEN1
      INTEGER       STRLEN
      EXTERNAL      STRLEN


      LEN1 = STRLEN(STRING1)
      DO 15, I = 1, STRLEN(STRING2)
        STRING1(LEN1+I:LEN1+I) = STRING2(I:I)
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
      SUBROUTINE STRIPSTR(STRING)

      IMPLICIT NONE
      CHARACTER*(*) STRING
      INTEGER       I, J, K, CHANGED, STRLEN
      EXTERNAL STRLEN

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
C Routine   : UPCASE
C Purpose   : To convert a string to upper case
C Interface : STRING      intent(INOUT)   string to be changed
C Author    : Arnold Moene
C Date      : May 24, 2000
C...........................................................................
      SUBROUTINE UPCASE(STRING)

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
C Routine   : CLEARSTR
C Purpose   : To set all characters in a string to CHAR(0)
C Interface : STRING      intent(INOUT)   string to be changed
C Author    : Arnold Moene
C Date      : May 25, 2000
C...........................................................................
      SUBROUTINE CLEARSTR(STRING)

      IMPLICIT NONE
      CHARACTER*(*) STRING
      INTEGER       I

      DO 15, I=1,LEN(STRING)
         STRING(I:I) = CHAR(0)
   15 CONTINUE
      END
      
