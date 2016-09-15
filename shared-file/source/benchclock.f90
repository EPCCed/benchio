module benchclock

  implicit none

  logical,          save, private :: firstcall = .true.
  double precision, save, private :: ticktime = 0.0

  integer, parameter :: int32kind = selected_int_kind( 9)
  integer, parameter :: int64kind = selected_int_kind(18)

!
!  Select high resolution clock
!

  integer, parameter :: intkind = int64kind
  integer(kind = intkind) :: count,rate

contains

double precision function benchtime()

  double precision :: dummy

! Ensure clock is initialised  

  if (firstcall) dummy = benchtick()

  call system_clock(count)

  benchtime  = dble(count)*ticktime

end function benchtime


double precision function benchtick()

  if (firstcall) then

     firstcall = .false.
     call system_clock(count, rate)
     ticktime = 1.0d0/dble(rate)

  end if

  benchtick = ticktime

end function benchtick

end module benchclock



