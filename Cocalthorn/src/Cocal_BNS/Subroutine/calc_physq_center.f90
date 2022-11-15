subroutine calc_physq_center
  use def_matter, only : emd
  use def_matter_parameter, only : pinx
  use def_quantities, only : rho_max, pre_max, epsi_max, q_max
  implicit none
!
  rho_max  = emd(0,0,0)**pinx
  pre_max  = emd(0,0,0)**(pinx+1.0d0)
  epsi_max = emd(0,0,0)**pinx*(1.0d0+pinx*emd(0,0,0))
  q_max    = emd(0,0,0)
!
end subroutine calc_physq_center
