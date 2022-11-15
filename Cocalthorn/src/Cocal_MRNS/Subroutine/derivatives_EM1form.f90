subroutine derivatives_EM1form
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_emfield, only : vaxd, vayd, vazd
  use def_emfield_derivatives, only : pdvaxd, pdvayd, pdvazd, &
  &                                   cdvaxd, cdvayd, cdvazd
  use interface_grgrad_vec_partial
  use interface_grgrad_vec_covariant
  implicit none
!
  character(len=1) :: index
!
  call grgrad_vec_partial(vaxd,vayd,vazd,pdvaxd,pdvayd,pdvazd)
  call grgrad_vec_covariant('d',vaxd,vayd,vazd,pdvaxd,pdvayd,pdvazd, &
       &                                       cdvaxd,cdvayd,cdvazd)
!
end subroutine derivatives_EM1form
