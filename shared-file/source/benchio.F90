program benchio

  use benchclock
  use mpi
  use mpiio
  use iohdf5
  use ionetcdf

  implicit none

  integer, parameter :: iolayerstart = 1
  integer, parameter :: iolayerend = 4
  integer, parameter :: numiolayer = 4
  integer, parameter :: numstriping = 3
  integer, parameter :: maxlen = 64
  integer, parameter :: numrep = 3

  character*(maxlen), dimension(numiolayer)  :: iostring, iolayername
  character*(maxlen), dimension(numstriping) :: stripestring

  character*(maxlen) :: filename

  integer :: iolayer, istriping, irep

! Set local array size - global sizes l1, l2 and l3 are scaled
! by number of processes in each dimension

  integer, parameter :: n1 = 256
  integer, parameter :: n2 = 256
  integer, parameter :: n3 = 256
  integer, parameter :: ndim = 3

  integer :: i1, i2, i3, j1, j2, j3, l1, l2, l3, p1, p2, p3

  double precision :: iodata(0:n1+1, 0:n2+1, 0:n3+1)

  integer :: rank, size, ierr, comm, cartcomm, dblesize
  integer, dimension(ndim) :: dims, coords

  integer, parameter :: iounit = 12
  integer, parameter :: mib = 1024*1024

  logical :: reorder = .false.
  logical, dimension(ndim) :: periods = [.false., .false., .false.]

  double precision :: t0, t1, time, iorate, mibdata
  double precision :: mintime_w, maxiorate_w, avgtime_w, avgiorate_w
  double precision :: mintime_r, maxiorate_r, avgtime_r, avgiorate_r

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

!  Skip layer if support is not compiled in
!  Expects iolayers in order: serial, MPI-IO, HDF5, NetCDF
#ifndef WITH_SERIAL
     if (iolayer == 1) then
       cycle
     endif
#endif
#ifndef WITH_MPIIO
     if (iolayer == 2) then
       cycle
     endif
#endif
#ifndef WITH_HDF5
     if (iolayer == 3) then
       cycle
     endif
#endif
#ifndef WITH_NETCDF
     if (iolayer == 4) then
       cycle
     endif
#endif

     if (rank == 0) then
        write(*,*)
        write(*,*) '------'
        write(*,*) iostring(iolayer)
        write(*,*) '------'
        write(*,*)
     end if

     !do istriping = 1, numstriping
     do istriping = 2, 2

        filename = trim(stripestring(istriping))//'/'//trim(iolayername(iolayer))

        if (rank == 0) then
           write(*,*) 'Writing to ', filename
           ! Write statistics
           mintime_w = 0
           maxiorate_w = 0
           avgtime_w = 0
           avgiorate_w = 0
           
           ! Read statistics
           mintime_r = 0
           maxiorate_r = 0
           avgtime_r = 0
           avgiorate_r = 0
        end if

        do irep = 1, numrep
          call MPI_Barrier(comm, ierr)
          if (rank == 0) then
            t0 = benchtime()
          end if

          ! Write test
          select case (iolayer)

          case(1)
             call serialwrite(filename, iodata, n1, n2, n3, cartcomm)

          case(2)
             call mpiiowrite(filename, iodata, n1, n2, n3, cartcomm)

          case(3)
             call hdf5write(filename, iodata, n1, n2, n3, cartcomm)

          case(4)
             call netcdfwrite(filename, iodata, n1, n2, n3, cartcomm)

          case default
             write(*,*) 'Illegal value of iolayer = ', iolayer
             stop

          end select

          call MPI_Barrier(comm, ierr)
          if (rank == 0) then
             t1 = benchtime()

             time = t1 - t0
             iorate = mibdata/time
             avgtime_w = avgtime_w + time/numrep
             avgiorate_w = avgiorate_w + iorate/numrep

             if (maxiorate_w < iorate) then
               maxiorate_w = iorate
               mintime_w = time
             end if

             write(*,*) 'write time = ', time, ', rate = ', iorate, ' MiB/s'
          end if
          
          ! Read test
          
          ! Reset data array to illegal values in preparation for values being read back in
          ! Perform write of illegal values to clear I/O cache 
          ! TODO: Extremely slow solution
          iodata(:,:,:) = -1
          call netcdfwrite('tmp.nc', iodata, n1, n2, n3, cartcomm)
          call fdelete('tmp.nc')
          
          call MPI_Barrier(comm, ierr)
          if (rank == 0) then
            t0 = benchtime()
          end if
          ! TODO: Not really a valid test.
          ! Currently just reads data back into array from file written above. 
          select case (iolayer)

          case(1)
             call serialwrite(filename, iodata, n1, n2, n3, cartcomm)

          case(2)
             call mpiiowrite(filename, iodata, n1, n2, n3, cartcomm)

          case(3)
             call hdf5write(filename, iodata, n1, n2, n3, cartcomm)

          case(4)
             call netcdfread(filename, iodata, n1, n2, n3, cartcomm)

          case default
             write(*,*) 'Illegal value of iolayer = ', iolayer
             stop

          end select
          
          call MPI_Barrier(comm, ierr)
          if (rank == 0) then
             t1 = benchtime()

             time = t1 - t0
             iorate = mibdata/time
             avgtime_r = avgtime_r + time/numrep
             avgiorate_r = avgiorate_r + iorate/numrep

             if (maxiorate_r < iorate) then
               maxiorate_r = iorate
               mintime_r = time
             end if

             write(*,*) 'read time = ', time, ', rate = ', iorate, ' MiB/s'
             call fdelete(filename)
          end if
          
        end do
        if (rank == 0) then
          write(*,*)
          write(*,*) 'Write statistics:'
          write(*,*) 'mintime = ', mintime_w, ', maxrate = ', maxiorate_w, ' MiB/s'
          write(*,*) 'avgtime = ', avgtime_w, ', avgrate = ', avgiorate_w, ' MiB/s'
          write(*,*)
          write(*,*) 'Read statistics:'
          write(*,*) 'mintime = ', mintime_r, ', maxrate = ', maxiorate_r, ' MiB/s'
          write(*,*) 'avgtime = ', avgtime_r, ', avgrate = ', avgiorate_r, ' MiB/s'
          write(*,*)
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
