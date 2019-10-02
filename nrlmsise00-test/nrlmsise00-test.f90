! nrlmsise00 standard atmosphere test program
! created by Daichi Hayashi 2 Oct. 2019.
module parm
  implicit none
  integer IDY,MASS
  double precision SEC,GLAT,GLONG
end module parm

module mod_nrlmsise00
  implicit none
  ! outputs for standard atmosphere
  ! rho : h -> altitude, output -> rho [kg/m^3]
  ! temp : h -> altitude, output -> temperature [K]
  ! press : h -> altitude, output -> pressure [Pa]
contains
  function rho(h)
    use parm
    implicit none
    double precision rho,h
    double precision D(9)
    double precision T(2)
    call METERC(.true.)
    call GTD7D(IDY,SEC,h,GLAT,GLONG,0.,150.,150.,4.,MASS,D,T)
    rho = D(6)
    return
  end function rho
  function temp(h)
    use parm
    implicit none
    double precision temp,h
    double precision D(9)
    double precision T(2)
    call METERC(.true.)
    call GTD7(IDY,SEC,h,GLAT,GLONG,0.,150.,150.,4.,MASS,D,T)
    temp = T(2)
    return
  end function temp
  function press(h)
    use parm
    implicit none
    double precision press,h
    press = rho(h)*286.9d0*temp(h)
    return
  end function press
end module mod_nrlmsise00

program nrlmsise00_test
  use parm
  use mod_nrlmsise00
  implicit none
  double precision h

  open(100, file='out.dat')
  write(100,*)'h     rho      temp        press'

  h = 3.75d2 ! [km]
  IDY = 20001 ! 2020/01/01
  SEC = 0.0d0
  GLAT = 0.0d0
  GLONG = 0.0d0
  MASS = 48

  do while (h > 0.0d0)
     write(100,*) h,rho(h),temp(h)
     h = h - 1.0d-1
  end do
end program nrlmsise00_test
