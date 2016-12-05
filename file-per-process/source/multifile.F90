module mpifile

  use mpi
  implicit none

contains

! Write data to multiple files
subroutine multiwrite(filestem, iodata, n1, n2, n3, cartcomm)

  character*(*) :: filestem
  
  integer :: n1, n2, n3
  double precision, dimension(0:n1+1,0:n2+1,0:n3+1) :: iodata

  integer :: cartcomm, ierr, rank, size
  integer, parameter :: iounit = 10
  character(len=6) :: crank
  character(len=50) :: filename

  integer :: i

  call MPI_Comm_size(cartcomm, size, ierr)
  call MPI_Comm_rank(cartcomm, rank, ierr)

  write(crank, "(I0.6)") rank
  filename = trim(filestem)//crank
  open(file=filename, unit=iounit+rank, access='stream')

  write(unit=iounit+rank) iodata(1:n1, 1:n2, 1:n3)

  close(iounit+rank)

end subroutine multiwrite

subroutine multidelete(filestem, cartcomm)

  implicit none

  character *(*) :: filestem
  integer, parameter :: iounit = 15
  integer :: cartcomm, ierr, rank, size
  integer :: stat
  character(len=6) :: crank
  character(len=50) :: filename

  call MPI_Comm_size(cartcomm, size, ierr)
  call MPI_Comm_rank(cartcomm, rank, ierr)

  write(crank, "(I0.6)") rank
  filename = trim(filestem)//crank
  open(unit=iounit+rank, iostat=stat, file=filename, status='old')
  if (stat.eq.0) close(unit=iounit+rank, status='delete')

end subroutine multidelete


end module mpifile
