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
  double precision :: shape = 2.0d0/3.0d0 ! for sphere (0.66666...)
  double precision :: abl_heat = 1.0d6 ! heat for ablation [J/kg]
  double precision :: bright_coeff = 1.0d0
end module mod_constant

module mod_parm
  ! module for parameters
  implicit none
  double precision dt,tf,m_f,star_density,injection_speed,injection_angle,d_init,area_init,m_init
  double precision drag_coeff_inp,heat_coeff_inp
  double precision SEC,ALT,GLAT,GLONG,STL
  double precision :: F107A = 150.0d0
  double precision :: F107  = 150.0d0
  double precision :: AP = 4.0d0
  integer IYD,MASS
end module mod_parm

module mod_nrlmsise00
  implicit none
  ! outputs for standard atmosphere
  ! rho : h -> altitude, output -> rho [kg/m^3]
  ! temp : h -> altitude, output -> temperature [K]
  ! press : h -> altitude, output -> pressure [Pa]
contains
  function rho(h)
    use mod_parm, only: IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,MASS
    implicit none
    double precision rho,h
    double precision D(9),T(2)
    call METERS(.true.)
    call GTD7(IYD,SEC,h*1d-3,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
    rho = D(6)
  end function rho
  function temp(h)
    use mod_parm, only: IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,MASS
    implicit none
    double precision temp,h
    double precision D(9),T(2)
    call METERS(.true.)
    call GTD7(IYD,SEC,h*1d-3,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
    temp = T(2)
  end function temp
  function press(h)
    use mod_parm
    implicit none
    double precision press,h
    press = rho(h)*286.9d0*temp(h)
    return
  end function press
end module mod_nrlmsise00

! module for coeffiecients (C_d, C_h)
module mod_coeff
  implicit none
contains
  ! C_d model (Henderson)
  function drag_coeff(h,velo,m)
    use mod_nrlmsise00
    implicit none
    double precision h,velo,m,gamma,temp_w,a1,b1,c1,sa,drag_coeff
    gamma = 1.4d0 ! gas constant
    temp_w = 3.0d3 ! wall temperature [K]
    a1 = 0.9d0 + 0.34d0/mach_num(h,velo)**2
    b1 = 1.86d0*dsqrt(mach_num(h,velo)/reynolds(h,velo,m))
    sa = mach_num(h,velo)*dsqrt(gamma/2.0d0)
    c1 = 2.0d0 + 2.0d0/sa**2 + 1.058d0/sa*dsqrt(temp_w/temp(h)) - 2.0d0/sa**4

    drag_coeff = (a1 + b1*c1)/(1.0d0 + b1)
  end function drag_coeff
  ! C_h model (Y. Prevereaud, 2014.)
  function heat_coeff(h,velo,m)
    use mod_constant, only: pi
    use mod_nrlmsise00
    implicit none
    double precision h,velo,m,gamma,gm_mc,q,a2,b2,c2,heat_coeff
    gamma = 1.4d0
    gm_mc = gamma*mach_num(h,velo)**2 ! use just for makeing simpler
    a2 = q_conv(h,velo,m)/(0.5d0*rho(h)*velo**3)
    b2 = pi*(3.0d0*gm_mc+8.0d0)/(64.0d0*(gm_mc+1.0d0))
    c2 = (pi**2*(gm_mc+4.0d0)+16.0d0)/(64.0d0*(gm_mc+1.0d0))
    heat_coeff = a2*b2/dsqrt(c2)
  end function heat_coeff

  ! conductive heating q_conv
  function q_conv(h,velo,m)
    use mod_nrlmsise00
    implicit none
    double precision h,velo,m,q_conv,diameter
    ! 7.9248 [km/s] is defined sqrt(radius_earth * gravity(at h = 0))
    q_conv = 110.35d6/dsqrt(diameter(m)*0.5d0)*dsqrt(rho(h)/rho(0.0d0))*(velo/7.9248d3)**3.15d0
  end function q_conv

  ! Mach number
  function mach_num(h,velo)
    use mod_nrlmsise00
    implicit none
    double precision h,velo,gamma,gass_const,ss,mach_num
    gamma = 1.4d0
    gass_const = 286.9d0 ! gass const of Air (M = 28.966[g/mol])
    ss = dsqrt(gamma*gass_const*temp(h)) ! sound speed
    mach_num = velo/ss
  end function mach_num

  ! Reynolds number
  function reynolds(h,velo,m)
    use mod_nrlmsise00
    implicit none
    double precision h,reynolds,velo,m,diameter
    reynolds = rho(h)*velo*diameter(m)/mu(h)
  end function reynolds

  ! viscosity coefficient
  function mu(h)
    use mod_nrlmsise00
    implicit none
    double precision h,mu,mu_0,temp_0,sutherland
    ! constants (from U.S. stand. atmos. 1976.)
    mu_0 = 1.7894d-5 ! [Pa*s]
    temp_0 = 2.8815d2 ! [K]
    sutherland = 1.104d2 ! [K]
    ! calc. mu
    mu = mu_0*((temp(h)/temp_0)**1.5d0)*(temp_0+sutherland)/(temp(h)+sutherland)
  end function mu
end module mod_coeff

! main program
program trajectory
  use mod_constant
  use mod_parm, only: injection_speed,injection_angle,dt,tf,m_f,ALT,m_init
  use mod_nrlmsise00
  use mod_coeff
  implicit none
  double precision t,r,dr,th,dth,m,brightness
  double precision k1(5),k2(5),k3(5),k4(5)
  double precision func1,func2,func3,func4,func5,velo,area,diameter
  double precision altitude,air_dens,air_temp,air_visc,velocity,mach,re_num,c_d,c_h

  ! read run_time parameter.
  call parm

  ! initialize
  t = 0.0d0 ! initial time [s]
  r = ALT + radius_earth ! initial r [m]
  th = 0.0d0 ! initial theta [rad.]
  dr = injection_speed*dsin(injection_angle)
  dth = injection_speed*dcos(injection_angle)/r ! V_th = r*dth
  m = m_init

  ! output data
  open(100, file='./out.dat')
  open(101, file='./coeff.dat')
  open(102, file='./mach_drag.dat')
  open(103, file='./atmos_model.dat')
  open(104, file='./brightness.dat')
  write(100,*) '       Time  Altitude  Velosity      mass'
  write(101,*) '  Altitude   mach    reynolds     Cd      Ch'
  write(102,*) '  mach     Cd'
  write(103,*) ' Altitude    density        temp        mu'
  write(104,*) ' Altitude brightness'

  do while ((t < tf).and.((r-radius_earth) > 0.0d0).and.(m > m_f))
    ! calc. parameters
    altitude = r-radius_earth
    air_dens = rho(altitude)
    air_temp = temp(altitude)
    air_visc = mu(altitude)
    velocity = velo(r,dr,dth)
    mach = mach_num(altitude,velocity)
    re_num = reynolds(altitude,velocity,m)
    c_d = drag_coeff(altitude,velocity,m)
    c_h = heat_coeff(altitude,velocity,m)

    ! brightness
    brightness = -bright_coeff*(0.5d0*velocity**2*func5(t,r,dr,th,dth,m) + m*velocity*dsqrt((func2(t,r,dr,th,dth,m)-r*dth**2)**2 + (r*func4(t,r,dr,th,dth,m)+2.0d0*dr*dth)**2))

    ! outputs
    write(100,200) t,altitude,velocity,m
    write(101,201) altitude,mach,re_num,c_d,c_h
    write(102,202) mach,c_d
    write(103,203) altitude,air_dens,air_temp,air_visc
    write(104,204) altitude,brightness
200 format(e12.4,2(f10.2),e12.4)
201 format(f10.2,f7.3,e12.4,f7.3,e12.4)
202 format(2(f7.3))
203 format(f10.2,3(e12.4))
204 format(f10.2,e12.4)

    ! 4th-order runge-ketta
    k1(1) = dt*func1(t,r,dr,th,dth,m)
    k1(2) = dt*func2(t,r,dr,th,dth,m)
    k1(3) = dt*func3(t,r,dr,th,dth,m)
    k1(4) = dt*func4(t,r,dr,th,dth,m)
    k1(5) = dt*func5(t,r,dr,th,dth,m)

    k2(1) = dt*func1(t+dt/2.0d0,r+k1(1)/2.0d0,dr+k1(2)/2.0d0,th+k1(3)/2.0d0,dth+k1(4)/2.0d0,m+k1(5)/2.0d0)
    k2(2) = dt*func2(t+dt/2.0d0,r+k1(1)/2.0d0,dr+k1(2)/2.0d0,th+k1(3)/2.0d0,dth+k1(4)/2.0d0,m+k1(5)/2.0d0)
    k2(3) = dt*func3(t+dt/2.0d0,r+k1(1)/2.0d0,dr+k1(2)/2.0d0,th+k1(3)/2.0d0,dth+k1(4)/2.0d0,m+k1(5)/2.0d0)
    k2(4) = dt*func4(t+dt/2.0d0,r+k1(1)/2.0d0,dr+k1(2)/2.0d0,th+k1(3)/2.0d0,dth+k1(4)/2.0d0,m+k1(5)/2.0d0)
    k2(5) = dt*func5(t+dt/2.0d0,r+k1(1)/2.0d0,dr+k1(2)/2.0d0,th+k1(3)/2.0d0,dth+k1(4)/2.0d0,m+k1(5)/2.0d0)

    k3(1) = dt*func1(t+dt/2.0d0,r+k2(1)/2.0d0,dr+k2(2)/2.0d0,th+k2(3)/2.0d0,dth+k2(4)/2.0d0,m+k2(5)/2.0d0)
    k3(2) = dt*func2(t+dt/2.0d0,r+k2(1)/2.0d0,dr+k2(2)/2.0d0,th+k2(3)/2.0d0,dth+k2(4)/2.0d0,m+k2(5)/2.0d0)
    k3(3) = dt*func3(t+dt/2.0d0,r+k2(1)/2.0d0,dr+k2(2)/2.0d0,th+k2(3)/2.0d0,dth+k2(4)/2.0d0,m+k2(5)/2.0d0)
    k3(4) = dt*func4(t+dt/2.0d0,r+k2(1)/2.0d0,dr+k2(2)/2.0d0,th+k2(3)/2.0d0,dth+k2(4)/2.0d0,m+k2(5)/2.0d0)
    k3(5) = dt*func5(t+dt/2.0d0,r+k2(1)/2.0d0,dr+k2(2)/2.0d0,th+k2(3)/2.0d0,dth+k2(4)/2.0d0,m+k2(5)/2.0d0)

    k4(1) = dt*func1(t+dt,r+k3(1),dr+k3(2),th+k3(3),dth+k3(4),m+k3(5))
    k4(2) = dt*func2(t+dt,r+k3(1),dr+k3(2),th+k3(3),dth+k3(4),m+k3(5))
    k4(3) = dt*func3(t+dt,r+k3(1),dr+k3(2),th+k3(3),dth+k3(4),m+k3(5))
    k4(4) = dt*func4(t+dt,r+k3(1),dr+k3(2),th+k3(3),dth+k3(4),m+k3(5))
    k4(5) = dt*func5(t+dt,r+k3(1),dr+k3(2),th+k3(3),dth+k3(4),m+k3(5))

    ! calc. next step
    r   = r   + (k1(1) + 2.0d0*k2(1) + 2.0d0*k3(1) + k4(1))/6.0d0
    dr  = dr  + (k1(2) + 2.0d0*k2(2) + 2.0d0*k3(2) + k4(2))/6.0d0
    th  = th  + (k1(3) + 2.0d0*k2(3) + 2.0d0*k3(3) + k4(3))/6.0d0
    dth = dth + (k1(4) + 2.0d0*k2(4) + 2.0d0*k3(4) + k4(4))/6.0d0
    m   = m   + (k1(5) + 2.0d0*k2(5) + 2.0d0*k3(5) + k4(5))/6.0d0

    ! integration
    t = t + dt
  end do

  ! close files
  close(100)
  close(101)
  close(102)
  close(103)
  close(104)

end program trajectory

! ******************************************************************
! run_time parameters
! ******************************************************************
subroutine parm
  use mod_constant
  use mod_parm, only: dt,tf,m_f,drag_coeff_inp,star_density,injection_speed,injection_angle,d_init,m_init,area_init,IYD,SEC,ALT,GLAT,GLONG,STL,MASS
  implicit none

  write(6,*) 'dt :'
  read(5,*) dt ! delta t [s]
  write(6,*) 't_f :'
  read(5,*) tf ! finish time [s]
  write(6,*) 'm_f :'
  read(5,*) m_f ! finish mass [kg]
  write(6,*) 'Cd'
  read(5,*) drag_coeff_inp
  write(6,*) 'dens_star'
  read(5,*) star_density
  write(6,*) 'inj_spee'
  read(5,*) injection_speed
  write(6,*) 'inj_ang'
  read(5,*) injection_angle
  write(6,*) 'star_d'
  read(5,*) d_init
  write(6,*) 'IYD'
  read(5,*) IYD ! Year and Days
  write(6,*) 'UTC'
  read(5,*) SEC ! Universal Time
  write(6,*) 'ALT'
  read(5,*) ALT ! altitude
  write(6,*) 'GLAT'
  read(5,*) GLAT ! geodetic latitude
  write(6,*) 'GLONG'
  read(5,*) GLONG ! geodetic longtitude
  write(6,*) 'MASS'
  read(5,*) MASS ! Mass Number (0:temp,48:all)

  ! variables
  STL = SEC/3600.0d0 + GLONG/15.0d0
  area_init = pi/4.0d0*d_init**2 ! initial projected area [m^2]
  m_init = star_density*pi/6.0d0*(d_init**3) ! mass = density * volume [kg]
end subroutine parm

function func1(t,r,dr,th,dth,m)
  ! function for dr/dt
  implicit none
  double precision t,r,dr,th,dth,m,func1
  func1 = dr
end function func1

function func2(t,r,dr,th,dth,m)
  ! function for d/dt(dr/dt)
  use mod_constant, only: const_gravity,mass_earth,radius_earth
  use mod_nrlmsise00
  use mod_coeff
  implicit none
  double precision t,r,dr,th,dth,m,func2,velo,area
  func2 = -const_gravity*mass_earth/(r*r) - 1.0d0/(2.0d0*m)*drag_coeff(r-radius_earth,velo(r,dr,dth),m)*area(m)*rho(r-radius_earth)*velo(r,dr,dth)*dr + r*dth**2.0d0
end function func2

function func3(t,r,dr,th,dth,m)
  ! function for dth/dt
  implicit none
  double precision t,r,dr,th,dth,m,func3
  func3 = dth
end function func3

function func4(t,r,dr,th,dth,m)
  ! function for d/dt(dth/dt)
  use mod_constant, only: radius_earth
  use mod_nrlmsise00
  use mod_coeff
  implicit none
  double precision t,r,dr,th,dth,m,func4,velo,area
  func4 = -1.0d0/(2.0d0*m*r)*drag_coeff(r-radius_earth,velo(r,dr,dth),m)*area(m)*rho(r-radius_earth)*velo(r,dr,dth)*r*dth - 2.0d0/r*dr*dth
end function func4

! function for dm/dt
function func5(t,r,dr,th,dth,m)
  use mod_constant, only: radius_earth,shape,abl_heat
  use mod_nrlmsise00
  use mod_coeff
  double precision t,r,dr,th,dth,m,func5,area
  func5 = -0.5d0*heat_coeff(r-radius_earth,velo(r,dr,dth),m)*area(m)*rho(r-radius_earth)*(velo(r,dr,dth)**3)/abl_heat
end function func5

function velo(r,dr,dth)
  ! function for calculate Velosity
  implicit none
  double precision r,dr,dth,velo
  velo = dsqrt(dr**2 + r**2*dth**2)
end function velo

! function for calc. Projected Area
function area(m)
  use mod_constant, only: shape
  use mod_parm, only: area_init, m_init
  double precision m,area
  area = area_init*((m/m_init)**shape)
end function area

! function for calc. diameter
function diameter(m)
  use mod_constant, only: pi
  implicit none
  double precision diameter,area,m
  diameter = 2.0d0*dsqrt(area(m)/pi)
end function diameter

