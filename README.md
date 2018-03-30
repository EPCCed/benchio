# benchio
EPCC I/O benchmarking application. Tests write bandwidth to a single shared
file for a given problem size per processor (weak scaling).

Each test is performed 10 times and the minimum, maximum and average bandwidth
returned. It is recommended that the maximum bandwidth be considered in most
cases, due to variations in I/O performance from user contention.

Data layout is 3D strided - intended to more closely resemble that of a real
world application than 2D sequential. By default, array size is 256x256x256.

Supports POSIX (serial), MPI-IO, HDF5 and NetCDF backends. A run will test all
backends included at compile time.

# Compilation

Tested under all three programming environments on ARCHER: Cray, GNU, and Intel.

For an ARCHER build:

`module load cray-netcdf-hdf5parallel/4.4.0`

`module load cray-hdf5-parallel/1.8.16`  

`make clean`

`make`

For other platforms, edit the `Makefile` to point variables `FC` and `CC` to
MPI compilers (e.g. `mpif90` and `mpicc`), and edit `LFLAGS` to the location of
HDF5 and NetCDF libraries. Then `make clean && make`.

## Building with selected backends

POSIX (serial), MPI-IO, HDF5 and NetCDF backends can be disabled by commenting
out the relevant `FFLAGS` lines in the `Makefile`. e.g. Commenting out:
`FFLAGS+= -DWITH_NETCDF` will build the application with only MPI-IO and HDF5
support.

# Running

For example, to test performance with maximum striping on Lustre:

`mkdir -p benchio_files`

`lfs setstripe -c -1 benchio_files`

`aprun -n <NUMBER_OF_PROCESSORS> ./benchio.x`

Explanation of commands follows:

The application expects a working directory with subdirectory `benchio_files`.
An I/O error will be thrown if this subdirectory is not present.

Under Lustre, the appropriate striping patterns should be set on this
subdirectory using the `lfs setstripe` command.

The `benchio.x` executable should be launched with the platform job launcher,
e.g. `aprun` or `mpirun`. Results are written to STDOUT.

Complete sample run scripts are included in `source/run_scripts`.

## Adjusting volume of data

The array size (per MPI rank) is controlled through the following constants in
`benchio.F90`:

`integer, parameter :: n1 = 256`

`integer, parameter :: n2 = 256`

`integer, parameter :: n3 = 256`

