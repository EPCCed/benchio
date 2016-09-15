program benchio

  use benchclock
  use mpiio

  implicit none

  integer, parameter :: iolayerstart = 2
  integer, parameter :: iolayerend = 2
  integer, parameter :: numiolayer = 4
  integer, parameter :: numstriping = 3
  integer, parameter :: maxlen = 64
  integer, parameter :: numrep = 10

  character*(maxlen), dimension(numiolayer)  :: iostring, iolayername
  character*(maxlen), dimension(numstriping) :: stripestring

  character*(maxlen) :: filename

  integer :: iolayer, istriping, irep

! Set local array size - global sizes l1, l2 and l3 are scaled
! by number of processes in each dimension

  integer, parameter :: n1 = 256
  integer, parameter :: n2 = 256
  integer, parameter :: n3 = 256

  integer :: i1, i2, i3, j1, j2, j3, l1, l2, l3, p1, p2, p3

  double precision :: iodata(0:n1+1, 0:n2+1, 0:n3+1)

  integer :: rank, size, ierr, comm, cartcomm, dblesize
  integer, dimension(ndim) :: dims, coords

  integer, parameter :: iounit = 12
  integer, parameter :: mib = 1024*1024

  logical :: reorder = .false.
  logical, dimension(ndim) :: periods = [.false., .false., .false.]

  double precision :: t0, t1, time, iorate, mibdata
  double precision :: mintime, maxiorate, avgtime, avgiorate

  iostring(1) = 'Serial'
  iostring(2) = 'MPI-IO'
  iostring(3) = ' HDF5 '
  iostring(4) = 'NetCDF'

  iolayername(1) = 'serial.dat'
  iolayername(2) = 'mpiio.dat'
  iolayername(3) = 'hdf5.dat'
  iolayername(4) = 'netcdf.dat'

  stripestring(1) = 'unstriped'
  stripestring(2) = 'striped'
  stripestring(3) = 'defstriped'

  call MPI_Init(ierr)

  comm = MPI_COMM_WORLD

  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)

  dims(:) = 0

! Set 3D processor grid

  call MPI_Dims_create(size, ndim, dims, ierr)

! Reverse dimensions as MPI assumes C ordering (this is not essential)

  p1 = dims(3)
  p2 = dims(2)
  p3 = dims(1)

! Compute global sizes

  l1 = p1*n1
  l2 = p2*n2
  l3 = p3*n3

  call MPI_Type_size(MPI_DOUBLE_PRECISION, dblesize, ierr)

  mibdata = float(dblesize*n1*n2*n3)*float(p1*p2*p3)/float(mib)

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Simple Parallel IO benchmark'
     write(*,*) '----------------------------'
     write(*,*)
     write(*,*) 'Running on ', size, ' process(es)'
     write(*,*) 'Process grid is (', p1, ', ', p2, ', ', p3, ')'
     write(*,*) 'Array size is   (', n1, ', ', n2, ', ', n3, ')'
     write(*,*) 'Global size is  (', l1, ', ', l2, ', ', l3, ')'
     write(*,*)
     write(*,*) 'Total amount of data = ', mibdata, ' MiB'
     write(*,*)
     write(*,*) 'Clock resolution is ', benchtick()*1.0e6, ', usecs'
  end if
  
  dims(1) = p1
  dims(2) = p2
  dims(3) = p3

  call MPI_Cart_create(comm, ndim, dims, periods, reorder, cartcomm, ierr)

! Set halos to illegal values

  iodata(:,:,:) = -1
  
! Set iodata core to have unique values 1, 2, ..., p1*n1*p2*n2*p3*n3

  call MPI_Cart_coords(cartcomm, rank, ndim, coords, ierr)
  
  do i3 = 1, n3
     do i2 = 1, n2
        do i1 = 1, n1

           j1 = coords(1)*n1 + i1
           j2 = coords(2)*n2 + i2
           j3 = coords(3)*n3 + i3

           iodata(i1,i2,i3) = (j3-1)*l1*l2 + (j2-1)*l1 + j1

        end do
     end do
  end do

  do iolayer = iolayerstart, iolayerend

     if (rank == 0) then
        write(*,*)
        write(*,*) '------'
        write(*,*) iostring(iolayer)
        write(*,*) '------'
        write(*,*)
     end if

     do istriping = 2, numstriping

        filename = trim(stripestring(istriping))//'/'//trim(iolayername(iolayer))

        if (rank == 0) then
           write(*,*) 'Writing to ', filename
           mintime = 0
           maxiorate = 0
           avgtime = 0
           avgiorate = 0
        end if

        do irep = 1, numrep
          call MPI_Barrier(comm, ierr)
          if (rank == 0) then
            t0 = benchtime()
          end if

          select case (iolayer)

          case(1)
             call serialwrite(filename, iodata, n1, n2, n3, cartcomm)

          case(2)
             call mpiiowrite(filename, iodata, n1, n2, n3, cartcomm)

          case(3)
          !   call hdf5write(filename, iodata, n1, n2, n3, cartcomm)

          case(4)
          !   call netcdfwrite(filename, iodata, n1, n2, n3, cartcomm)

          case default
             write(*,*) 'Illegal value of iolayer = ', iolayer
             stop

          end select

          call MPI_Barrier(comm, ierr)
          if (rank == 0) then
             t1 = benchtime()

             time = t1 - t0
             iorate = mibdata/time
             avgtime = avgtime + time/numrep
             avgiorate = avgiorate + iorate/numrep

             if (maxiorate < iorate) then
               maxiorate = iorate
               mintime = time
             end if

             write(*,*) 'time = ', time, ', rate = ', iorate, ' MiB/s'
             call fdelete(filename)
          end if
        end do
        if (rank == 0) then
          write(*,*) 'mintime = ', mintime, ', maxrate = ', maxiorate, ' MiB/s'
          write(*,*) 'avgtime = ', avgtime, ', avgrate = ', avgiorate, ' MiB/s'
          write(*,*) 'Deleting: ', filename
          write(*,*)
        end if
     end do
  end do

  if (rank == 0) then
     write(*,*)
     write(*,*) '--------'
     write(*,*) 'Finished'
     write(*,*) '--------'
     write(*,*)
  end if

  call MPI_Finalize(ierr)
  
end program benchio

subroutine fdelete(filename)

  implicit none

  character *(*) :: filename
  integer, parameter :: iounit = 15
  integer :: stat

  open(unit=iounit, iostat=stat, file=filename, status='old')
  if (stat.eq.0) close(unit=iounit, status='delete')

end subroutine fdelete
