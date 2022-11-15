subroutine printout_plot
!
  use phys_constant, only : c, g, solmas, pi, nnrg
  use def_metric_1D
  use CB_fR_param_iomis
  use grid_parameter_1D, only : nrf
  use def_matter_1D		!ber
  use coordinate_grav_r_1D, only : rg
  use weight_grav_1D, only : drf
  implicit none
!
  real(8) :: rgcgs, hh, pre, rho, ene, &
  &          rhocgs, precgs, epsiloncgs, &
  &          proper_radius, sou(0:nnrg)
  integer :: irg
!
  sou(0:nrf) = psif(0:nrf)**2
  call halfsou(sou)
  proper_radius = 0.0d0
!
  open(3,file='phys_plot.dat',status='unknown')
!
  do irg = 0, nrf
    call peos_q2hprho(emdg(irg), hh, pre, rho, ene)
    rhocgs = rho*c**6/(g**3*solmas**2)
    precgs = pre*c**8/(g**3*solmas**2)
    epsiloncgs = ene*c**6/(g**3*solmas**2)
    rgcgs = rg(irg)*radi*g*solmas*1.0d-5/c**2
!
    if (irg > 0) then
      proper_radius = proper_radius &
      &             + sou(irg)*drf(irg)*radi*g*solmas*1.0d-5/c**2
    end if
!
    write(3,'(8es20.12)') rgcgs, psi(irg), alph(irg), emdg(irg), &
    &                     precgs, epsiloncgs, rhocgs, proper_radius
  end do
!
  close(3)
!
end subroutine printout_plot
