MF=     Makefile

# Compilers
FC=     ftn
CC=     cc
#
#
FFLAGS= -g
CFLAGS= $(FFLAGS)
LFLAGS=
#
EXE=	benchio.x
#
FSRC= \
	mpiio.f90      \
	netcdf.f90      \
	hdf5.f90      \
	benchclock.f90 \
	benchio.f90

CSRC=

# .SUFFIXES:
.SUFFIXES: .f90 .c .o

OBJ=    $(FSRC:.f90=.o) $(CSRC:.c=.o)

MOD=    $(FSRC:.f90=.mod)

.f90.o:
	$(FC) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

all:    $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ): $(MF)

clean:
	rm -f $(OBJ) $(EXE) $(MOD) core
