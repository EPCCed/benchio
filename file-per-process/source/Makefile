MF=     Makefile

# Compilers
FC=     mpif90
CC=     mpicc
#
#
FFLAGS= -g

# Comment out a line below to disable the named tests

CFLAGS= $(FFLAGS)
LFLAGS=
#
EXE=	benchio.x
#
FSRC= \
	multifile.F90      \
	benchclock.f90 \
	benchio.F90

CSRC=

# .SUFFIXES:
.SUFFIXES: .F90 .f90 .c .o

FOBJ=	$(FSRC:.F90=.o)
OBJ=    $(FOBJ:.f90=.o) $(CSRC:.c=.o)

FMOD=	$(FSRC:.F90=.mod)
MOD=    $(FMOD:.f90=.mod)

.F90.o:
	$(FC) $(FFLAGS) -c $<

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
