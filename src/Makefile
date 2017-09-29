
#SHELL = /bin/sh

#-----------------------------------------------------------------------------
# Name of the program
NAME = EMPM
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#                     HARDWARE SECTION
#-- Choose the appropriate architecture, and uncomment all lines 
#-- in that section.
#-----------------------------------------------------------------------------
#   Intel compiler 
FC   = ifort
OPTS = -I../include -O3  -xHost -ip -assume byterecl -ftz -r8
#OPTS = -I../include -g -warn all -check all -fpe0 -traceback -debug extended -check uninit
NC   = -lnetcdf -I$(NETCDF)/include -L/$(NETCDF)/lib
#-----------------------------------------------------------------------------
#   gfortran compiler
#FC   = gfortran
#OPTS =  -ffree-form -ffree-line-length-none -O2 -fbounds-check
#-----------------------------------------------------------------------------
#   g95 compiler - not tested for the moment
#FC   = g95
#OPTS = -ffree-form -ffree-line-length-huge -I../include -O2
#CPP  = cpp -C -P -traditional

#-----------------------------------------------------------------------------
#-- You shouldn't need to change anything below here
#-----------------------------------------------------------------------------

SRC   = accu.f90 \
        array.f90 \
 	assigd.f90 \
 	average.f90 \
        const.f90 \
	count_lines.f90 \
 	dgm.f90 \
 	domain.f90 \
	drop_map.f90 \
	eff_radius.f90 \
 	EMPM.f90 \
 	ew.f90 \
 	fcnkb.f90 \
 	fdtime.f90 \
 	fintrp.f90 \
	index.f90 \
	intpol.f90 \
	intspline.f90 \
	max_min.f90 \
	mean_rms.f90 \
	netcdf_check.f90 \
	netcdf_prep.f90 \
	netcdf_write.f90 \
 	odeint.f90 \
 	prob.f90 \
	psplint.f90 \
  	qt_prob.f90 \
	radius_pdf.f90 \
 	rand1.f90 \
 	rkck.f90 \
	rkqs.f90 \
 	rmax.f90 \
	size_dis.f90 \
 	triplet.f90

OBJS = $(SRC:%.f90=%.o)

FFLAGS  =  $(OPTS)


.SUFFIXES:
.SUFFIXES:      .F .f90 .o

# Inference rule
.f90.o:
			$(FC) $(FFLAGS) $(NC) $(DEBUG) -o $@ -c $<

all :empm

empm:			$(OBJS)
			$(FC)  -o $(NAME) $(OBJS) $(NC)

clean:
			rm -f *.o *.a *.mod

# DEPENDENCIES : only dependencies after this line (don't remove the word DEPENDENCIES)
dgm.o: array.o const.o
EMPM.o: const.o intspline.o
fcnkb.o: array.o const.o intspline.o
odeint.o: array.o
psplint.o: intspline.o
rkck.o: array.o
rkqs.o: array.o
size_dis.o: const.o
