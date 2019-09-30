! ******************************************************************
! trajectory calculation for artificial shooting star
! created by Daichi Hayashi 30 Sep. 2019.
! ******************************************************************

program trajectory
  double precision pi,const_gravity,mass_earth,diameter,projected_area,star_density,mass_star,drag_coeff
  ! constant values
  pi = dacos(-1.0d0) ! pi = 3.14...
  const_gravity = 6.67430d-11 ! constant of gravitation [m^3/kg*s^2]
  mass_earth = 5.9736d24 ! mass of earth (NASA) [kg]

  ! variables (とりあえずテスト用)
  diameter = 1.0d-2 ! star diameter [m]
  projected_area = pi/4*diameter**2 ! projected area [m^2]
  star_density = 5.0d3 ! mass density [kg/m^3]
  mass_star = star_density*pi/6.0d0*(diameter**3) ! mass = density * volume [kg]
  drag_coeff = 0.1d0 ! drag coefficient C_d [-]

  ! read run_time parameter.
  call parm
end program trajectory

subroutine parm
end subroutine parm
