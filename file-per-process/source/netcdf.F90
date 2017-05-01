module ionetcdf

#ifdef WITH_NETCDF

  use netcdf
  use mpi
  
  implicit none

contains

subroutine netcdfwrite(filestem, iodata, n1, n2, n3, cartcomm)

  integer, parameter :: ndim = 3

  character*(*) :: filestem
  character(len=6) :: crank
  character(len=50) :: filename
  
  integer :: n1, n2, n3
  double precision, dimension(0:n1+1,0:n2+1,0:n3+1) :: iodata

  integer, dimension(ndim) :: arraysize, arraystart
  integer, dimension(ndim) :: arraysubsize

  integer :: cartcomm, ierr, rank

  integer :: ncid, varid, dimids(ndim)
  integer :: x_dimid, y_dimid, z_dimid

  call MPI_Comm_rank(cartcomm, rank, ierr)

! Each process writes to own file with rank as suffix

  write(crank, "(I0.6)") rank
  filename = trim(filestem)//crank

! Subtract halos for array subsize

  arraysize(:)      = [n1+2, n2+2, n3+2]
  arraysubsize(:)   = [n1, n2, n3]
  arraystart(:)     = [1, 1, 1]   ! Use Fortran indexing
 
  ! Create (i.e. open) the netCDF file. The NF90_NETCDF4 flag causes a
  ! HDF5/netCDF-4 type file to be created.
  call check( nf90_create(filename, nf90_netcdf4, ncid))
  ! Define the dimensions. NetCDF returns an ID for each.
  call check( nf90_def_dim(ncid, "x", n1, x_dimid) )
  call check( nf90_def_dim(ncid, "y", n2, y_dimid) )
  call check( nf90_def_dim(ncid, "z", n3, z_dimid) )
  ! The dimids array is used to pass the ID's of the dimensions of 
  ! the variables. 
  dimids = (/ x_dimid , y_dimid, z_dimid /)
  ! Define the variable. The type of the variable in this case is
  ! NF90_DOUBLE (8-byte float).
  call check( nf90_def_var(ncid, "data", NF90_DOUBLE, dimids, varid) )
  ! Make sure file it not filled with default values which doubles wrote volume
  ! Below line modified on 8th March 2016
  call check ( nf90_def_var_fill(ncid, varid, 1, 1) )
  ! End define mode. This tells netCDF we are done defining
  ! metadata.
  call check( nf90_enddef(ncid) )

  ! Write the data to file, start will equal the displacement from the 
  ! start of the file and count is the number of points each proc writes. 
  call check( nf90_put_var(ncid, varid, iodata(1:n1, 1:n2, 1:n3), &
              start = arraystart, count = arraysubsize) )
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  call check( nf90_close(ncid) )

end subroutine netcdfwrite

subroutine check(status)
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop
  end if
end subroutine check  

! WITH_NETCDF not defined. Dummy subroutine.
#else

contains
subroutine netcdfwrite(filestem, iodata, n1, n2, n3, cartcomm)
  character*(*) :: filestem
  integer :: n1, n2, n3, cartcomm
  double precision, dimension(0:n1+1,0:n2+1,0:n3+1) :: iodata
end subroutine netcdfwrite

#endif

end module ionetcdf
