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

FFLAGS=$(INCDIR) -ff2c -O3 -Wall -Wno-unused -fexpensive-optimizations -fomit-frame-pointer -ffixed-line-length-none
LDFLAGS=$(LIBDIR) -lnetcdf 

all: ec_ncdf2$(EXT) planfit$(EXT)

ec_ncdf2$(EXT): ec_ncdf.o slatec.o ec_corr.o ec_file.o ec_math.o ec_phys.o ec_gene.o ec_nc.o
	$(FC) -o ec_ncdf2$(EXT) slatec.o ec_ncdf.o ec_corr.o ec_file.o ec_math.o ec_phys.o ec_gene.o ec_nc.o  $(LDFLAGS)

planfit$(EXT): planfit.o slatec.o ec_corr.o ec_file.o ec_math.o ec_phys.o ec_gene.o ec_nc.o
	$(FC) -o planfit$(EXT) slatec.o planfit.o ec_corr.o ec_file.o ec_math.o ec_phys.o ec_gene.o ec_nc.o  $(LDFLAGS)

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

	
