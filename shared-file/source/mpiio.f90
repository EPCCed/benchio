module mpiio

  use mpi
  implicit none

  integer, parameter :: ndim = 3

contains

subroutine mpiiowrite(filename, iodata, n1, n2, n3, cartcomm)

  character*(*) :: filename
  
  integer :: n1, n2, n3
  double precision, dimension(0:n1+1,0:n2+1,0:n3+1) :: iodata

  integer, dimension(ndim) :: arraysize, arraystart
  integer, dimension(ndim) :: arraygsize, arraysubsize

  integer :: cartcomm, ierr, rank, size

  integer :: filetype, mpi_subarray, fh
  integer (kind=MPI_OFFSET_KIND) :: disp = 0
  integer, dimension(MPI_STATUS_SIZE) :: status

  integer, dimension(ndim) :: dims, coords
  logical, dimension(ndim) :: periods

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
  arraystart(:) = arraysubsize(:) * coords(:)

  call MPI_Type_create_subarray(ndim, arraygsize, arraysubsize, arraystart, &
                                MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                                filetype, ierr)

  call MPI_Type_commit(filetype, ierr)

!
! Define subarray for this process, ie what portion of the local array
! is to be written (excludes halos); starting positions use C-indexing.
!

  arraystart(:) = 1

  call MPI_Type_create_subarray(ndim, arraysize, arraysubsize, arraystart, &
                                MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                                mpi_subarray, ierr)

  call MPI_Type_commit(mpi_subarray, ierr)

!
!  Open the file for reading only and attach to file handle fh
!  No IO hints are passed since MPI_INFO_NULL is specified
!

  call MPI_File_open(cartcomm, filename, MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                     MPI_INFO_NULL, fh, ierr)

  if (ierr /= MPI_SUCCESS) write(*,*) 'Open error on rank ', rank

!
!  Set view for this process using appropriate datatype
!

  call MPI_File_set_view(fh, disp, MPI_DOUBLE_PRECISION, filetype, 'native', &
                         MPI_INFO_NULL, ierr)
  if (ierr /= MPI_SUCCESS) write(*,*) 'View error on rank ', rank

!
!  Write all the data for this process.
!  Remove halo data by passing an MPI subarray type
!

  call MPI_File_write_all(fh, iodata, 1, mpi_subarray, status, ierr)

  if (ierr /= MPI_SUCCESS) write(*,*) 'Write error on rank ', rank

!  Alternative: remove halo data by passing an explicit Fortran subarray
!  This should give an identical file to the above call with mpi_subaray
!
!  call MPI_File_write_all(fh, &
!       array(1:arraysubsize(1), 1:arraysubsize(2), 1:arraysubsize(3)), &
!       arraysubsize(1)*arraysubsize(2)*arraysubsize(3), &
!       MPI_DOUBLE_PRECISION, status, ierr)


!
!  Close file
!

  call MPI_File_close(fh, ierr)

  if (ierr /= MPI_SUCCESS) write(*,*) 'Close error on rank ', rank

  call MPI_Type_free(filetype, ierr)
  call MPI_Type_free(mpi_subarray, ierr)

end subroutine mpiiowrite

subroutine serialwrite(filename, iodata, n1, n2, n3, cartcomm)

  character*(*) :: filename
  
  integer :: n1, n2, n3
  double precision, dimension(0:n1+1,0:n2+1,0:n3+1) :: iodata

  integer :: cartcomm, ierr, rank, size
  integer, parameter :: iounit = 10

  integer :: i

  call MPI_Comm_size(cartcomm, size, ierr)
  call MPI_Comm_rank(cartcomm, rank, ierr)

!  Write same amount of data as the parallel write but do it all from rank 0
!  This is just to get a baseline figure for serial IO performance - note
!  that the contents of the file will be differnent from the parallel calls

  if (rank == 0) then

     open(file=filename, unit=iounit, access='stream')

     do i = 1, size
        write(unit=iounit) iodata(1:n1, 1:n2, 1:n3)
     end do

     close(iounit)

  end if

end subroutine serialwrite

end module mpiio
