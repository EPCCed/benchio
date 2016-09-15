module ionetcdf

  use netcdf
  use mpi
  
  implicit none

contains

subroutine netcdfwrite(filename, iodata, n1, n2, n3, cartcomm)

  integer, parameter :: ndim = 3

  character*(*) :: filename
  
  integer :: n1, n2, n3
  double precision, dimension(0:n1+1,0:n2+1,0:n3+1) :: iodata

  integer, dimension(ndim) :: arraysize, arraystart
  integer, dimension(ndim) :: arraygsize, arraysubsize

  integer :: cartcomm, ierr, rank, size

  integer, dimension(ndim) :: dims, coords
  logical, dimension(ndim) :: periods

  integer :: ncid, varid, oldmode, dimids(ndim)
  integer :: x_dimid, y_dimid, z_dimid

  call MPI_Comm_size(cartcomm, size, ierr)
  call MPI_Comm_rank(cartcomm, rank, ierr)

  call MPI_Cart_get(cartcomm, ndim, dims, periods, coords, ierr)

  arraysize(:) = [n1+2, n2+2, n3+2]

! Subtract halos for array subsize

  arraysubsize(:)   = [n1, n2, n3]

!
! Define filetype for this process, ie what portion of the global array
! this process owns; starting positions use C-indexing (ie counting from 0).
!

  arraygsize(:) = arraysubsize(:) * dims(:)
  arraystart(:) = arraysubsize(:) * coords(:) + 1   ! Use Fortran indexing
 
  ! Create (i.e. open) the netCDF file. The NF90_NETCDF4 flag causes a 
  ! HDF5/netCDF-4 type file to be created. The comm and info parameters 
  ! cause parallel I/O to be enabled. 
  call check( nf90_create(filename, ior(nf90_netcdf4,nf90_mpiio), ncid, comm=cartcomm, info=MPI_INFO_NULL))
  ! Define the dimensions. NetCDF returns an ID for each. Any 
  ! metadata operations must take place on ALL processors
  call check( nf90_def_dim(ncid, "x", arraygsize(1), x_dimid) )
  call check( nf90_def_dim(ncid, "y", arraygsize(2), y_dimid) )
  call check( nf90_def_dim(ncid, "z", arraygsize(3), z_dimid) )
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
  ! metadata. This operation is collective and all processors will
  ! write their metadata to disk.
  call check( nf90_enddef(ncid) )
  ! This appears also to write out all the data which is not what we want
  ! Have commented this out - appears to give exactly the same output
  ! file without this callm but at a faster rate overall. DSH 09/12/15

 ! Parallel access
  call check( nf90_var_par_access(ncid, varid, nf90_collective) )

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
 

end module ionetcdf
