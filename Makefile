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
#

# For some Unix
BUILD-LINUX=build-linux
LEXT=
LFC=g77
LAR=ar
LRANLIB=ranlib
LINCDIR='-I/home/arnold/include -I$(LIBSRC)/'
LLIBDIR=-L/home/arnold/lib

# For MingW
BUILD-WIN32=build-win32
WFC=/usr/local/bin/i386-mingw32-g77
WAR=/usr/local/bin/i386-mingw32-ar
WRANLIB=/usr/local/bin/i386-mingw32-ranlib
WEXT=.exe
WINCDIR='-I/usr/local/i386-mingw32/include -I$(LIBSRC)/'
WLIBDIR=-L/usr/local/i386-mingw32/lib
WEXT=.exe

FFLAGS=$(INCDIR)  -ff2c -O3 -Wall -Wno-unused -fexpensive-optimizations -fomit-frame-pointer -ffixed-line-length-none -g $(EXTRA_FFLAG)
LDFLAGS=$(LIBDIR) -lnetcdf -L. -lecpack
LIBSRC=../libsrc
SRC=../src

# For Robodoc
ROBODOC=robodoc
TMPDIR=`pwd`/tmp
DOCDIR=`pwd`/libdoc
MKDIR=mkdir -p
PROJECT=ecpack
LS=ls
SED=sed
RM=rm -f
RMDIR=rm -rf
MV=mv
DOCFILES=$(LIBSRC)/ec_corr.f $(LIBSRC)/ec_phys.f $(LIBSRC)/ec_math.f $(LIBSRC)/ec_gene.f $(LIBSRC)/physcnst.inc $(LIBSRC)/parcnst.inc
PDFLATEX=pdflatex


all: win32 linux win32-programs linux-programs docs

win32: 
	($(MKDIR) $(BUILD-WIN32) ; cd $(BUILD-WIN32) ; make -f ../Makefile FC=$(WFC) AR=$(WAR) RANLIB=$(WRANLIB) EXT=$(WEXT) INCDIR=$(WINCDIR) libecpack.a)
win32-programs:
	($(MKDIR) $(BUILD-WIN32) ; cd $(BUILD-WIN32) ; make -f ../Makefile FC=$(WFC) AR=$(WAR) RANLIB=$(WRANLIB) EXT=$(WEXT) INCDIR=$(WINCDIR) LIBDIR=$(WLIBDIR) ec_ncdf.exe planang.exe)

linux: 
	($(MKDIR) $(BUILD-LINUX) ; cd $(BUILD-LINUX) ; make -f ../Makefile FC=$(LFC) AR=$(LAR) RANLIB=$(LRANLIB) EXT=$(LEXT) INCDIR=$(LINCDIR) libecpack.a)
linux-programs:
	($(MKDIR) $(BUILD-LINUX) ; cd $(BUILD-LINUX) ; make -f ../Makefile FC=$(LFC) AR=$(LAR) RANLIB=$(LRANLIB) EXT=$(LEXT) INCDIR=$(LINCDIR) LIBDIR=$(LLIBDIR) ec_ncdf planang)


# By default, just the library
libecpack.a: slatec.o ec_corr.o ec_math.o ec_phys.o ec_gene.o
	$(AR) r $@ $?
	$(RANLIB) $@

# The ec_ncdf program
ec_ncdf$(EXT): ec_ncdf.o  ec_file.o ec_text.o ec_nc.o libecpack.a calibrat.o
	$(FC) -o ec_ncdf$(EXT) ec_ncdf.o ec_file.o ec_text.o ec_nc.o  calibrat.o $(LDFLAGS)

# The planar fit program
planang$(EXT): planang.o ec_file.o ec_text.o ec_nc.o calibrat.o
	$(FC) -o planang$(EXT) planang.o ec_file.o ec_text.o ec_nc.o calibrat.o $(LDFLAGS)

# This individual files of library
slatec.o: $(LIBSRC)/slatec.f
	$(FC) $(FFLAGS) -c $(LIBSRC)/slatec.f

ec_corr.o: $(LIBSRC)/ec_corr.f $(LIBSRC)/physcnst.inc $(LIBSRC)/parcnst.inc
	$(FC) $(FFLAGS) -c $(LIBSRC)/ec_corr.f

ec_text.o: $(SRC)/ec_text.f $(LIBSRC)/physcnst.inc $(LIBSRC)/parcnst.inc 
	$(FC) $(FFLAGS) -c $(SRC)/ec_text.f

ec_math.o: $(LIBSRC)/ec_math.f $(LIBSRC)/physcnst.inc $(LIBSRC)/parcnst.inc 
	$(FC) $(FFLAGS) -c $(LIBSRC)/ec_math.f

ec_phys.o: $(LIBSRC)/ec_phys.f $(LIBSRC)/physcnst.inc $(LIBSRC)/parcnst.inc 
	$(FC) $(FFLAGS) -c $(LIBSRC)/ec_phys.f

ec_gene.o: $(LIBSRC)/ec_gene.f $(LIBSRC)/physcnst.inc $(LIBSRC)/parcnst.inc 
	$(FC) $(FFLAGS) -c $(LIBSRC)/ec_gene.f

# Other individual files
ec_file.o: $(SRC)/ec_file.f $(LIBSRC)/physcnst.inc $(LIBSRC)/parcnst.inc 
	$(FC) $(FFLAGS) -c $(SRC)/ec_file.f

ec_ncdf.o: $(SRC)/ec_ncdf.f $(LIBSRC)/physcnst.inc $(LIBSRC)/parcnst.inc $(LIBSRC)/calcomm.inc $(LIBSRC)/version.inc
	$(FC) $(FFLAGS) -c $(SRC)/ec_ncdf.f

ec_nc.o: $(SRC)/ec_nc.f $(LIBSRC)/physcnst.inc $(LIBSRC)/parcnst.inc 
	$(FC) $(FFLAGS) -c $(SRC)/ec_nc.f

calibrat.o: $(SRC)/calibrat.f $(LIBSRC)/physcnst.inc $(LIBSRC)/parcnst.inc 
	$(FC) $(FFLAGS) -c $(SRC)/calibrat.f

planang.o: $(SRC)/planang.f $(LIBSRC)/physcnst.inc $(LIBSRC)/parcnst.inc  $(LIBSRC)/version.inc
	$(FC) $(FFLAGS) -c $(SRC)/planang.f

# The docs
docs: 
	$(MKDIR) $(DOCDIR)
	(cd $(DOCDIR) ; make -f ../Makefile dochtml docpdf)

dochtml: $(DOCFILES)
	$(MKDIR) $(TMPDIR)
	cp $(DOCFILES) $(TMPDIR)
	$(ROBODOC) --multidoc --doc ./  --src $(TMPDIR) --html --index
	$(RMDIR) $(TMPDIR)

docpdf: $(DOCFILES)
	$(MKDIR) $(TMPDIR)
	cp $(DOCFILES) $(TMPDIR)
	$(ROBODOC) --singledoc --doc $(PROJECT)  --src $(TMPDIR) --latex --index --sections --toc
	$(PDFLATEX) $(PROJECT).tex
	$(PDFLATEX) $(PROJECT).tex
#	$(RM) $(PROJECT).tex *.aux *.toc *.idx
	$(RMDIR) $(TMPDIR)

clean:
	$(RMDIR) $(BUILD-WIN32) $(BUILD-LINUX) $(DOCDIR)

