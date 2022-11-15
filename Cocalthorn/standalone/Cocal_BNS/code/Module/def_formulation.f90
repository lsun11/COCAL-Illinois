module def_formulation
  use phys_constant, only : long
  use grid_parameter, only : chrot, chgra, chope
  implicit none
  real(long) :: swlp(1:4), swls(1:4), swflu
  integer :: iswl
  integer :: sw_hij, sw_hxx, sw_hxy, sw_hxz, sw_hyy, sw_hyz, sw_hzz
  integer :: sw_psi, sw_alps, sw_bi, sw_disk
  character(1) :: char
  real(8) :: detgt
!  character(1) :: chrot, chgra, chope, char
!
contains
subroutine choose_formulation
  use grid_parameter, only : chrot, chgra, chope
  implicit none
  chrot = 'c'
  chgra = 'w'
  chope = 'L'
!
! --- Switch for hij source term.
!
  swlp(1) = 1.0d0
  swlp(2) = 1.0d0
  swlp(3) = 0.0d0
  swlp(4) = 0.0d0
  swls(1) = 0.0d0
  swls(2) = 1.0d0
  swls(3) = 0.0d0
  swls(4) = 1.0d0
!
! --- Switch for fluid source term and the gravity.
!
  if (chrot == 'c') then
    write(6,*) ' ### Co-rotating solutions. ###'
  else if (chrot == 'i') then 
    write(6,*) ' ### Irrotating solutions. ###'
  else
    write(6,*) ' ### INVARID PARAMETER --- chrot. ###'
    stop
  end if
  swflu = 1.0d0
  if (chrot == 'c') swflu = 0.0d0
!
  if (chgra == 'i') then 
    write(6,*) ' ### IWM formalism. ###'
  else if (chgra == 'w') then 
    write(6,*) ' ### SUF-Waveless formalism. ###'
  else if (chgra == 'k') then 
    write(6,*) ' ### SUF-Waveless formalism, with Kerr-Schild gauge and CTT decomposition. ###'
  else if (chgra == 'c') then 
    write(6,*) ' ### Waveless-Helical cut off formalism. ###'
  else if (chgra == 'C') then 
    write(6,*) ' ### Waveless-Helical cut off formalism. ###'
  else if (chgra == 'h') then 
    write(6,*) ' ### Helical formalism. ###'
  else if (chgra == 'H') then 
    write(6,*) ' ### Helical cut off formalism (simple cutoff) ###'
  else if (chgra == 'W') then 
    write(6,*) ' ### Helical -> Waveless source ###'
  else 
    write(6,*) ' ### INVARID PARAMETER --- chgra. ###'
    stop
  end if
!
  if (chope == 'L') then 
    write(6,*) ' ### Invert Laplacian for hij. ###'
  else if (chope == 'H') then 
    write(6,*) ' ### Invert Helmholtz operator for hij. ###'
  else
    write(6,*) ' ### INVARID PARAMETER --- chope. ###'
    stop
  end if
!
  char = chrot//chgra//chope
!
end subroutine choose_formulation
!
subroutine choose_formulation_v1(ca,cb,cc)
  use grid_parameter, only : chrot, chgra, chope
  implicit none
  character(1), intent(in) :: ca,cb,cc
!
  chrot = ca
  chgra = cb
  chope = cc
!
  detgt = 1.0d0

  sw_disk = 1    ! 0 = iteration for CO spacetime only, no self gravity disk.
  sw_psi  = 1    ! Turn on/off iteration over psi
  sw_alps = 1    ! Turn on/off iteration over alps
  sw_bi   = 1    ! Turn on/off iteration over beta_i
  sw_hxx  = 1    ! Turn on/off iteration over h_{ij}
  sw_hxy  = 1    ! Turn on/off iteration over h_{ij}
  sw_hxz  = 1    ! Turn on/off iteration over h_{ij}
  sw_hyy  = 1    ! Turn on/off iteration over h_{ij}
  sw_hyz  = 1    ! Turn on/off iteration over h_{ij}
  sw_hzz  = 1    ! Turn on/off iteration over h_{ij}
  sw_hij  = sw_hxx+sw_hxy+sw_hxz+sw_hyy+sw_hyz+sw_hzz

!
! --- Switch for hij source term.
!
  swlp(1) = 1.0d0
  swlp(2) = 1.0d0
  swlp(3) = 0.0d0
  swlp(4) = 0.0d0
  swls(1) = 0.0d0
  swls(2) = 1.0d0
  swls(3) = 0.0d0
  swls(4) = 1.0d0
!
  write(6,'(a56,a56)')   "--------------------------------------------------------", &
                    &    "--------------------------------------------------------"

  if (chrot == 'c') then
    write(6,*) ' ### Co-rotating solutions. ###'
  else if (chrot == 'i') then 
    write(6,*) ' ### Irrotating solutions. ###'
  else
    write(6,*) ' ### INVARID PARAMETER --- chrot. ###'
    stop
  end if
  swflu = 1.0d0
  if (chrot == 'c') swflu = 0.0d0
!
  if (chgra == 'i') then 
    write(6,*) ' ### IWM formalism. ###'
  else if (chgra == 'w') then 
    write(6,*) ' ### SUF-Waveless formalism. ###'
  else if (chgra == 'k') then 
    write(6,*) ' ### SUF-Waveless formalism, with Kerr-Schild gauge and CTT decomposition. ###'
    write(6,'(a80,10i2,a5)') ' ### (sw_disk,sw_psi,sw_alps,sw_bi,sw_hxx,sw_hxy,sw_hxz,sw_hyy,sw_hyz,sw_hzz)=(', &
       &                            sw_disk,sw_psi,sw_alps,sw_bi,sw_hxx,sw_hxy,sw_hxz,sw_hyy,sw_hyz,sw_hzz,') ###'
  else if (chgra == 'b') then 
    write(6,*) ' ### SUF-Waveless formalism, with Kerr-Schild type coordinates. ###'
  else if (chgra == 'c') then 
    write(6,*) ' ### Waveless-Helical cut off formalism. ###'
  else if (chgra == 'C') then 
    write(6,*) ' ### Waveless-Helical cut off formalism. ###'
  else if (chgra == 'h') then 
    write(6,*) ' ### Helical formalism. ###'
  else if (chgra == 'H') then 
    write(6,*) ' ### Helical cut off formalism (simple cutoff) ###'
  else if (chgra == 'W') then 
    write(6,*) ' ### Helical -> Waveless source ###'
  else 
    write(6,*) ' ### INVARID PARAMETER --- chgra. ###'
    stop
  end if
!
  if (chope == 'L') then 
    write(6,*) ' ### Invert Laplacian for hij. ###'
  else if (chope == 'H') then 
    write(6,*) ' ### Invert Helmholtz operator for hij. ###'
  else
    write(6,*) ' ### INVARID PARAMETER --- chope. ###'
    stop
  end if
  write(6,'(a56,a56)')   "--------------------------------------------------------", & 
                    &    "--------------------------------------------------------"
!
  char = chrot//chgra//chope
end subroutine choose_formulation_v1
!
end module def_formulation
