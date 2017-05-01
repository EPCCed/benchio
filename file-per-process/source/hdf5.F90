module iohdf5

#ifdef WITH_HDF5

  use hdf5
  use mpi
  
  implicit none

contains

subroutine hdf5write(filestem, iodata, n1, n2, n3, cartcomm)

  integer, parameter :: ndim = 3

  character*(*) :: filestem
  character(len=6) :: crank
  character(len=50) :: filename
  
  integer :: n1, n2, n3
  double precision, dimension(0:n1+1,0:n2+1,0:n3+1) :: iodata

  integer :: info = MPI_INFO_NULL
  integer(hsize_t), dimension(ndim) :: dimsf  ! dataset dimensions.

  character(len=8), parameter :: dsetname = "IntArray" ! Dataset name

  integer(hid_t) :: file_id       ! file identifier 
  integer(hid_t) :: dset_id       ! dataset identifier 
  integer(hid_t) :: filespace     ! dataspace identifier in file 
  integer(hid_t) :: memspace      ! dataspace identifier in memory
  integer(hid_t) :: plist_id      ! property list identifier 

  integer(hsize_t), dimension(ndim) :: count  
  integer(hssize_t), dimension(ndim) :: offset 

  integer, dimension(ndim) :: arraysize
  integer, dimension(ndim) :: arraysubsize

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
  dimsf(:) = arraysubsize(:)

  ! Initialise the count and offset arrays
  count(1) = n1     ! Defines the number of values each proc dumps to 
  count(2) = n2
  count(3) = n3                                   ! the HDF5 file. 

  offset(:) = coords(:) * count(:)  ! Defines the offset used in the HDF5 file

  ! Each process writes to own file with rank as suffix
  write(crank, "(I0.6)") rank
  filename = trim(filestem)//crank

  ! Initialize FORTRAN predefined datatypes
  CALL h5open_f(ierr)
  
  ! Create a new file using default properties.
  CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr)
  
  ! Create the data space for the  dataset.
  CALL h5screate_simple_f(ndim, dimsf, filespace, ierr)
  
  ! Create the dataset with default properties.
  CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, &
                      dset_id, ierr)
  
  ! Write the dataset to this rank's own file.
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, iodata(1:n1, 1:n2, 1:n3), &
                      dimsf, ierr)

  ! Close the dataspace, dataset, and file.
  CALL h5sclose_f(filespace, ierr)
  CALL h5dclose_f(dset_id, ierr)
  CALL h5fclose_f(file_id, ierr)

  ! Close FORTRAN predefined datatypes.
  CALL h5close_f(ierr)

end subroutine hdf5write

! WITH_HDF5 not defined. Dummy subroutine.
#else

contains
subroutine hdf5write(filestem, iodata, n1, n2, n3, cartcomm)
  character*(*) :: filestem
  integer :: n1, n2, n3, cartcomm
  double precision, dimension(0:n1+1,0:n2+1,0:n3+1) :: iodata
end subroutine hdf5write

#endif

end module iohdf5
