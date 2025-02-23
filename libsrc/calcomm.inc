C $Id$
C  calcomm.inc, part of ecpack
C  
C  Copyright (C) 
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
