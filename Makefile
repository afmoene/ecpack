FC=g77
#INCDIR=-I/usr/local/include
#LIBDIR=-L/usr/local/lib
#EXT=
INCDIR=-Ic:/djgpp/include
LIBDIR=-Lc:/djgpp/lib
EXT=.exe
FFLAGS=$(INCDIR) -ff2c -O3 -Wall -Wno-unused -fexpensive-optimizations
LDFLAGS=$(LIBDIR) -lnetcdf

ec_ncdf$(EXT): ec_ncdf.o slatec.o ecpack.o getconf.o
	g77 -o ec_ncdf$(EXT) slatec.o ecpack.o getconf.o ec_ncdf.o  $(LDFLAGS)

ec_ncdf.o: ec_ncdf.f physcnst.inc parcnst.inc calcomm.inc
	g77 $(FFLAGS) -c ec_ncdf.f

slatec.o: slatec.f
	g77 $(FFLAGS) -c slatec.f

ecpack.o: ecpack.f parcnst.inc calcomm.inc physcnst.for 
	g77 $(FFLAGS) -c ecpack.f

getconf.o: getconf.f
	g77 $(FFLAGS) -c getconf.f

	
