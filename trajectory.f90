! ******************************************************************
! trajectory calculation for artificial shooting star
! created by Daichi Hayashi 30 Sep. 2019.
! ******************************************************************
module mod_constant
  ! module for constant values
  implicit none

  double precision :: pi = dacos(-1.0d0) ! pi = 3.14
  double precision :: const_gravity = 6.67430d-11 ! constant of gravitation [m^3/kg*s^2]
  double precision :: mass_earth = 5.9736d24 ! mass of earth (NASA) [kg]
  double precision :: radius_earth = 6.371d6 ! radius of the Earth (on red road) [m]
end module mod_constant

module mod_parm
  ! module for parameters
  implicit none
  double precision :: dt,tf,drag_coeff,star_density,injection_speed,injection_angle,star_diameter,projected_area,mass_star
  double precision :: ALT,GLAT,GLONG,STL,F107A,F107,AP
  integer :: IYD,SEC,MASS
end module mod_parm

module mod_nrlmsise00
  implicit none
  ! outputs for standard atmosphere
  ! rho : h -> altitude, output -> rho [kg/m^3]
  ! temp : h -> altitude, output -> temperature [K]
  ! press : h -> altitude, output -> pressure [Pa]
contains
  function rho(h)
    use mod_parm
    implicit none
    double precision rho,h
    double precision D(9)
    double precision T(2)
    call METERC(.true.)
    call GTD7D(IYD,SEC,h/1.0d3,GLAT,GLONG,0.,150.,150.,4.,MASS,D,T)
    rho = D(6)
    return
  end function rho
  function temp(h)
    use mod_parm
    implicit none
    double precision temp,h
    double precision D(9)
    double precision T(2)
    call METERC(.true.)
    call GTD7(IYD,SEC,h/1.0d3,GLAT,GLONG,0.,150.,150.,4.,MASS,D,T)
    temp = T(2)
    return
  end function temp
  function press(h)
    use mod_parm
    implicit none
    double precision press,h
    press = rho(h)*286.9d0*temp(h)
    return
  end function press
end module mod_nrlmsise00

program trajectory
  use mod_constant
  use mod_parm
  use mod_nrlmsise00
  implicit none
  double precision :: t0,t,r0,r,dr0,dr,th0,th,dth0,dth
  double precision :: k1(4),k2(4),k3(4),k4(4)
  double precision :: func1, func2, func3, func4, velo

  ! read run_time parameter.
  call parm

  ! initialize
  t0 = 0.0d0 ! initial time [s]
  t = t0
  r0 = ALT + radius_earth ! initial r [m]
  r = r0
  th0 = 0.0d0 ! initial theta [rad.]
  th = th0
  dr0 = injection_speed*dsin(injection_angle)
  dr = dr0
  dth0 = injection_speed*dcos(injection_angle)/r ! V_th = r*dth
  dth = dth0

  ! output data
  open(100, file='./out.dat')
  write(100,*) 'Time      Altitude      Velosity'

  do while ((t < tf).and.((r-radius_earth) > 0.0d0))
    ! outputs
    write(100,*) t,(r-radius_earth),velo(r,dr,dth)

    ! 4th-order runge-ketta
    k1(1) = dt*func1(t,r,dr,th,dth)
    k1(2) = dt*func2(t,r,dr,th,dth)
    k1(3) = dt*func3(t,r,dr,th,dth)
    k1(4) = dt*func4(t,r,dr,th,dth)

    k2(1) = dt*func1(t+dt/2.0d0,r+k1(1)/2.0d0,dr+k1(2)/2.0d0,th+k1(3)/2.0d0,dth+k1(4)/2.0d0)
    k2(2) = dt*func2(t+dt/2.0d0,r+k1(1)/2.0d0,dr+k1(2)/2.0d0,th+k1(3)/2.0d0,dth+k1(4)/2.0d0)
    k2(3) = dt*func3(t+dt/2.0d0,r+k1(1)/2.0d0,dr+k1(2)/2.0d0,th+k1(3)/2.0d0,dth+k1(4)/2.0d0)
    k2(4) = dt*func4(t+dt/2.0d0,r+k1(1)/2.0d0,dr+k1(2)/2.0d0,th+k1(3)/2.0d0,dth+k1(4)/2.0d0)

    k3(1) = dt*func1(t+dt/2.0d0,r+k2(1)/2.0d0,dr+k2(2)/2.0d0,th+k2(3)/2.0d0,dth+k2(4)/2.0d0)
    k3(2) = dt*func2(t+dt/2.0d0,r+k2(1)/2.0d0,dr+k2(2)/2.0d0,th+k2(3)/2.0d0,dth+k2(4)/2.0d0)
    k3(3) = dt*func3(t+dt/2.0d0,r+k2(1)/2.0d0,dr+k2(2)/2.0d0,th+k2(3)/2.0d0,dth+k2(4)/2.0d0)
    k3(4) = dt*func4(t+dt/2.0d0,r+k2(1)/2.0d0,dr+k2(2)/2.0d0,th+k2(3)/2.0d0,dth+k2(4)/2.0d0)

    k4(1) = dt*func1(t+dt,r+k3(1),dr+k3(2),th+k3(3),dth+k3(4))
    k4(2) = dt*func2(t+dt,r+k3(1),dr+k3(2),th+k3(3),dth+k3(4))
    k4(3) = dt*func3(t+dt,r+k3(1),dr+k3(2),th+k3(3),dth+k3(4))
    k4(4) = dt*func4(t+dt,r+k3(1),dr+k3(2),th+k3(3),dth+k3(4))

    ! calc. next step
    r   = r   + (k1(1) + 2.0*k2(1) + 2.0*k3(1) + k4(1))/6.0d0
    dr  = dr  + (k1(2) + 2.0*k2(2) + 2.0*k3(2) + k4(2))/6.0d0
    th  = th  + (k1(3) + 2.0*k2(3) + 2.0*k3(3) + k4(3))/6.0d0
    dth = dth + (k1(4) + 2.0*k2(4) + 2.0*k3(4) + k4(4))/6.0d0

    ! integration
    t = t + dt
  end do

end program trajectory

! ******************************************************************
! run_time parameters
! ******************************************************************
subroutine parm
  use mod_constant
  use mod_parm
  implicit none

  write(6,*) 'input values'
  read(5,*) dt ! delta t [s]
  read(5,*) tf ! finish time [s]
  read(5,*) drag_coeff
  read(5,*) star_density
  read(5,*) injection_speed
  read(5,*) injection_angle
  read(5,*) star_diameter
  read(5,*) IYD ! Year and Days
  read(5,*) SEC ! Universal Time
  read(5,*) ALT ! altitude
  read(5,*) GLAT ! geodetic latitude
  read(5,*) GLONG ! geodetic longtitude
  GLAT = GLAT*pi/180.0d0
  GLONG = GLONG*pi/180.0d0
  STL = SEC/3600.0d0 + GLONG/15.0d0
  read(5,*) MASS ! Mass Number (0:temp,48:all)
  write(6,*) 'input complete'

  ! variables
  projected_area = pi/4*star_diameter**2 ! projected area [m^2]
  mass_star = star_density*pi/6.0d0*(star_diameter**3) ! mass = density * volume [kg]
end subroutine parm

function func1(t,r,dr,th,dth)
  ! function for dr/dt
  use mod_constant
  use mod_parm
  use mod_nrlmsise00
  implicit none
  double precision t,r,dr,th,dth,func1
  func1 = dr
  return
end function func1

function func2(t,r,dr,th,dth)
  ! function for d/dt(dr/dt)
  use mod_constant
  use mod_parm
  use mod_nrlmsise00
  implicit none
  double precision t,r,dr,th,dth,func2,velo
  func2 = -const_gravity*mass_earth/(r*r) - 1/(2*mass_star)*drag_coeff*projected_area*rho(r-radius_earth)*velo(r,dr,dth)*dr + r*dth**2

  return
end function func2

function func3(t,r,dr,th,dth)
  ! function for dth/dt
  use mod_constant
  use mod_parm
  use mod_nrlmsise00
  implicit none
  double precision t,r,dr,th,dth,func3
  func3 = dth
  return
end function func3

function func4(t,r,dr,th,dth)
  ! function for d/dt(dth/dt)
  use mod_constant
  use mod_parm
  use mod_nrlmsise00
  implicit none
  double precision t,r,dr,th,dth,func4,velo
  func4 = -1/(2*mass_star*r)*drag_coeff*projected_area*rho(r-radius_earth)*velo(r,dr,dth)*r*dth - 2/r*dr*dth
  return
end function func4

function velo(r,dr,dth)
  ! function for calculate Velosity
  use mod_constant
  use mod_parm
  implicit none
  double precision r,dr,dth,velo
  velo = dsqrt(dr**2 + r**2 * dth**2)
  return
end function velo

