! nrlmsise00 standard atmosphere test program
! created by Daichi Hayashi 2 Oct. 2019.
module mod_parm
  implicit none
  integer IDY,MASS
  double precision SEC,GLAT,GLONG
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
    call GTD7D(IDY,SEC,h,GLAT,GLONG,0.,150.,150.,4.,MASS,D,T)
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
    call GTD7(IDY,SEC,h,GLAT,GLONG,0.,150.,150.,4.,MASS,D,T)
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

!! module for coeffiecients (C_d, C_h)
!module mod_coeff
!  implicit none
!contains
!  ! C_d model (Henderson)
!  function drag_coeff(h,velo)
!    use mod_parm
!    use mod_nrlmsise00
!    implicit none
!    double precision h,velo,gamma,temp_w,a1,b1,c1,sa,drag_coeff
!    gamma = 1.4d0 ! gas constant
!    temp_w = 3.0d3 ! wall temperature [K]
!    a1 = 0.9d0 + 0.34d0/mach_num(h,velo)**2
!    b1 = 1.86d0*dsqrt(mach_num(h,velo)/reynolds(h,velo))
!    sa = mach_num(h,velo)*dsqrt(gamma/2.0d0)
!    c1 = 2.0d0 + 2.0d0/sa**2 + 1.058d0/sa*dsqrt(temp_w/temp(h)) - 2.0d0/sa**4
!
!    drag_coeff = (a1 + b1*c1)/(1.0d0 + b1)
!    return
!  end function drag_coeff
!  ! Mach number
!  function mach_num(h,velo)
!    use mod_parm
!    use mod_nrlmsise00
!    implicit none
!    double precision h,velo,gamma,gass_const,ss,mach_num
!    gamma = 1.4d0
!    gass_const = 286.9d0 ! gass const of Air (M = 28.966[g/mol])
!    ss = dsqrt(gamma*gass_const*temp(h)) ! sound speed
!    mach_num = velo/ss
!    return
!  end function mach_num
!  ! Reynolds number
!  function reynolds(h,velo)
!    use mod_parm
!    use mod_nrlmsise00
!    implicit none
!    double precision h,reynolds,velo,star_diameter
!    reynolds = rho(h)*velo*star_diameter/mu(h)
!    return
!  end function reynolds
!  ! viscosity coefficient
!  function mu(h)
!    use mod_nrlmsise00
!    implicit none
!    double precision h,mu,mu_0,temp_0,sutherland
!    ! constants (from U.S. stand. atmos. 1976.)
!    mu_0 = 1.7894d-5 ! [Pa*s]
!    temp_0 = 2.8815d2 ! [K]
!    sutherland = 1.104d2 ! [K]
!    ! calc. mu
!    mu = mu_0*(temp(h)/temp_0)**1.5d0*(temp_0+sutherland)/(temp(h)+sutherland)
!    return
!  end function mu
!end module mod_coeff

program nrlmsise00_test
  use mod_parm
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
     h = h - 1.0d-2
  end do
end program nrlmsise00_test
