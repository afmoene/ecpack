#  Makefile, part of ecpack
#  
#  Copyright (C) 
#    1998-2000   Arjan van Dijk, Wageningen University, The Netherlands
#    2000-2002   Arjan van Dijk, IMAU, The Netherlands
#    1999-2002   Arnold Moene, Wageningen University, The Netherlands
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
# 


FC=g77
AR=ar
#
# For some Unix
#INCDIR=-I/usr/local/include
#LIBDIR=-L/usr/local/lib
#EXT=
#
# For DJGPP
# INCDIR=-Ic:/djgpp/include
# LIBDIR=-Lc:/djgpp/lib
# EXT=.exe
#
# For MingW
INCDIR=-Ic:/gcc-2.95.2/include
LIBDIR=-Lc:/gcc-2.95.2/lib
EXT=.exe

FFLAGS=$(INCDIR)  -ff2c -O3 -Wall -Wno-unused -fexpensive-optimizations -fomit-frame-pointer -ffixed-line-length-none
LDFLAGS=$(LIBDIR) -lnetcdf -L. -lecpack

all: ec_ncdf$(EXT) planfit$(EXT) libecpack.a

libecpack.a: slatec.o ec_corr.o ec_math.o ec_phys.o ec_gene.o
	$(AR) r $@ $?

ec_ncdf$(EXT): ec_ncdf.o  ec_file.o ec_nc.o libecpack.a
	$(FC) -o ec_ncdf$(EXT) ec_ncdf.o ec_file.o ec_nc.o  $(LDFLAGS)

planfit$(EXT): planfit.o ec_file.o ec_nc.o
	$(FC) -o planfit$(EXT) planfit.o ec_file.o ec_nc.o  $(LDFLAGS)

ec_ncdf.o: ec_ncdf.f physcnst.inc parcnst.inc calcomm.inc version.inc
	$(FC) $(FFLAGS) -c ec_ncdf.f

slatec.o: slatec.f
	$(FC) $(FFLAGS) -c slatec.f

ec_corr.o: ec_corr.f physcnst.inc parcnst.inc
	$(FC) $(FFLAGS) -c ec_corr.f

ec_file.o: ec_file.f physcnst.inc parcnst.inc 
	$(FC) $(FFLAGS) -c ec_file.f

ec_math.o: ec_math.f physcnst.inc parcnst.inc 
	$(FC) $(FFLAGS) -c ec_math.f

ec_phys.o: ec_phys.f physcnst.inc parcnst.inc 
	$(FC) $(FFLAGS) -c ec_phys.f

ec_gene.o: ec_gene.f physcnst.inc parcnst.inc 
	$(FC) $(FFLAGS) -c ec_gene.f

ec_nc.o: ec_nc.f physcnst.inc parcnst.inc 
	$(FC) $(FFLAGS) -c ec_nc.f

	
