subroutine artificial_deformation
!
  use grid_parameter, only : ntg, npg
  use trigonometry_grav_phi, only : cosmpg
  use def_matter, only : rs
  implicit none
  integer :: it, ip
!
do it = 0, ntg
do ip = 0, npg
  rs(it,ip) = rs(it,ip)*(0.8+0.2d0*cosmpg(2,ip))
end do
end do
!
end subroutine artificial_deformation
