FC=g77
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

FFLAGS=$(INCDIR) -ff2c -O3 -Wall -Wno-unused -fexpensive-optimizations -fomit-frame-pointer
LDFLAGS=$(LIBDIR) -lnetcdf

all: ec_ncdf$(EXT) planfit$(EXT)

ec_ncdf$(EXT): ec_ncdf.o slatec.o ecpack.o getconf.o
	g77 -o ec_ncdf$(EXT) slatec.o ecpack.o getconf.o ec_ncdf.o  $(LDFLAGS)

planfit$(EXT): planfit.o slatec.o ecpack.o getconf.o
	g77 -o planfit$(EXT) slatec.o ecpack.o getconf.o planfit.o  $(LDFLAGS)


ec_ncdf.o: ec_ncdf.f physcnst.inc parcnst.inc calcomm.inc version.inc
	g77 $(FFLAGS) -c ec_ncdf.f

planfit.o: planfit.f physcnst.inc parcnst.inc calcomm.inc version.inc
	g77 $(FFLAGS) -c planfit.f

slatec.o: slatec.f
	g77 $(FFLAGS) -c slatec.f

ecpack.o: ecpack.f parcnst.inc calcomm.inc physcnst.for 
	g77 $(FFLAGS) -c ecpack.f

getconf.o: getconf.f
	g77 $(FFLAGS) -c getconf.f

	
