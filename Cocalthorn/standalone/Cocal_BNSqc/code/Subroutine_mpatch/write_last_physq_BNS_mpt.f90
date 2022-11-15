subroutine write_last_physq_BNS_mpt
  use phys_constant, only  : long, nmpt
  use def_quantities_mpt
  use def_matter_parameter_mpt
  use grid_parameter_mpt
  implicit none 
  real(long) :: emdenc, omega, radius, m_rest, j_star, m_adm, m_kom, j_adm
  real(long) :: m_grav_sph, m_rest_sph, M_R_sph, rhoc, qc, rhomax, qmax
  real(long) :: chicusp, rsurf, rgmid
  character(40) :: char1, char2, char3, char4, char5
  character(100) :: dircommand

  emdenc     = def_matter_param_real_(2,1)  ! same as q_c
  omega      = def_matter_param_real_(3,1)
  radius     = def_matter_param_real_(5,1)
  m_rest     = def_quantities_real_(4,1)  
  j_star     = def_quantities_real_(6,1)  
  m_adm      = def_quantities_real_(7,nmpt)  
  m_kom      = def_quantities_real_(8,nmpt)  
  j_adm      = def_quantities_real_(9,nmpt)  
  m_grav_sph = def_quantities_real_(28,1) 
  m_rest_sph = def_quantities_real_(29,1)
  M_R_sph    = def_quantities_real_(31,1)
  rhoc       = def_quantities_real_(40,1)  
  qc         = def_quantities_real_(42,1)  
  rhomax     = def_quantities_real_(43,1)  
  qmax       = def_quantities_real_(46,1)  
  chicusp    = def_quantities_real_(121,1)
  rsurf      = surf_param_real_(1,1)
  rgmid      = grid_param_real_(2,1)

  open(15,file='physical_quantities.txt',status='unknown')
  write(15,'(1p,20e20.12)') omega, radius, m_rest, j_star, m_adm, m_kom, j_adm, &
  &  m_grav_sph, m_rest_sph, M_R_sph, rhoc, qc, rhomax, qmax, chicusp, rsurf, rgmid 

  write (15,'(a30,1p,e20.12)') '#Omega*Radius               = ', omega
  write (15,'(a30,1p,e20.12)') '#Radius                     = ', radius
  write (15,'(a30,1p,e20.12)') '#Rest mass                  = ', m_rest
  write (15,'(a30,1p,e20.12)') '#NS angular momentum        = ', j_star
  write (15,'(a30,1p,e20.12)') '#ADM mass                   = ', m_adm
  write (15,'(a30,1p,e20.12)') '#Komar mass                 = ', m_kom
  write (15,'(a30,1p,e20.12)') '#Angular momentum           = ', j_adm
  write (15,'(a30,1p,e20.12)') '#Gravitational spher. mass  = ', m_grav_sph
  write (15,'(a30,1p,e20.12)') '#Rest spher. mass           = ', m_rest_sph
  write (15,'(a30,1p,e20.12)') '#Spherical compactness M/R  = ', M_R_sph
  write (15,'(a30,1p,e20.12)') '#Central restmass density   = ', rhoc
  write (15,'(a30,1p,e20.12)') '#Central Emden function     = ', qc
  write (15,'(a30,1p,e20.12)') '#Maximum restmass density   = ', rhomax
  write (15,'(a30,1p,e20.12)') '#Maximum Emden function     = ', qmax
  write (15,'(a30,1p,e20.12)') '#chi at cusp                = ', chicusp
  write (15,'(a30,1p,e20.12)') '#r_surf                     = ', rsurf
  write (15,'(a30,1p,e20.12)') '#rgmid                      = ', rgmid
  close(15)

end subroutine write_last_physq_BNS_mpt
