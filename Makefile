FC=g77
FFLAGS=-Ic:/gcc-2.95.2/include -ff2c -O3 -Wall -Wno-unused
LDFLAGS=-Lc:/gcc-2.95.2/lib -lnetcdf
FFLAGS=-Ic:/djgpp/include -ff2c -O3 -Wall -Wno-unused -fexpensive-optimizations
LDFLAGS=-Lc:/djgpp/lib -lnetcdf

ec_ncdf.exe: ec_ncdf.o slatec.o ecpack.o getconf.o
	g77 -o ec_ncdf.exe slatec.o ecpack.o getconf.o ec_ncdf.o  $(LDFLAGS)

ec_ncdf.o: ec_ncdf.f
	g77 $(FFLAGS) -c ec_ncdf.f

slatec.o: slatec.f
	g77 $(FFLAGS) -c slatec.f

ecpack.o: ecpack.f
	g77 $(FFLAGS) -c ecpack.f

getconf.o: getconf.f
	g77 $(FFLAGS) -c getconf.f

	
