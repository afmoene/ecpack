C $Id$
C File: calcomm.inc
C Common to pass data to calibration routine
C
      INTEGER FirstDay
      LOGICAL HAVE_TREF
      LOGICAL HAVE_TCCAL, HAVE_SONCAL, HAVE_HYGCAL
      
      COMMON/CalOffset/ FirstDay
      SAVE/CalOffset/
      
      COMMON/TC_INFO/ HAVE_TREF
      SAVE/TC_INFO/

      COMMON/CALFILE/ HAVE_TCCAL, HAVE_SONCAL, HAVE_HYGCAL
      SAVE/CALFILE/
