!______________________________________________
include 'include_modulefiles_1D_qeos.f90'
include 'include_subroutines_1D_qeos.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Import_isotr_qeos
!
  use phys_constant, only : nnrg
  use def_metric_1D, only : alphf, psif, psi, alph
  use def_matter_1D
  use coordinate_grav_r_1D, only : rg
  use grid_parameter_1D, only : iter_max, rgmid, nrgin, eps, emdc, nrg, nrf, rhoc_qs
  use def_matter_1D, only : ber, emd, emdg, radi, rhog, rhof
  implicit none
!
  real(8) :: m_iso(0:nnrg), h_iso(0:nnrg), prm_iso(0:nnrg), pm_iso(0:nnrg), rs_iso(0:nnrg), &
    &        rho0_iso(0:nnrg), pre_iso(0:nnrg), ene_iso(0:nnrg), rg_iso(0:nnrg)
  real(8) :: xecgs, dummy, rhocgs, precgs, epsiloncgs, xe
  real(8) ::  mo2r, m_normal, rs_normal, compa, rho0, pre, ene, &
  &       xi, xicgs, restmass_sph, propermass_sph, adm, gravmass_sph
  real(8) :: c_sound(0:nnrg)
!
  character(LEN=100) :: cmd
  integer :: irg, i, j, k, nlines
!
  write(6,*)' coordinates etc. '
  call coordinate_patch_kit_grav_1D
!
  nlines = 0
  call system("cat ovlas.dat | grep '[^ ]' | wc -l > nlines.txt")
  open(1,file='nlines.txt')
  read(1,*) nlines
  close(1)
  write(6,'(a21,i3)') "#Number of lines are ", nlines
  call system('rm nlines.txt')
  if (nlines <= 0) then
    write(6,*) " Number of lines in ovlas.dat =", nlines, "....exiting"
    stop
  endif

  ber = 0.0d0
  open(3,file='ovlas.dat',status='unknown')
  do i=0, nlines-1
!    read(3,'(5es15.7)') rg_iso(i), psi(i), alph(i), emdg(i), h_iso(i)
    read(3, '(18(es15.7))') xecgs, m_iso(i), rho0_iso(i), prm_iso(i), pm_iso(i), &
    &     psi(i), dummy, alph(i), rs_iso(i), h_iso(i), rhocgs,   &
    &     pre_iso(i), precgs, ene_iso(i), epsiloncgs, xe, rg_iso(i), &
    &     c_sound(i)
    ber = ber + alph(i)*h_iso(i)
    rhog(i) = rho0_iso(i)
  end do
  ber = ber/nlines
  close(3)

  open(4,file='ovphy_plot.dat',status='unknown')
  read(4 , '(13es14.6)') compa, emdc, rhoc_qs, rhocgs, pre, ene, &
  &       xi, xicgs, restmass_sph, propermass_sph, adm, gravmass_sph, dummy
  close(4)
  radi = xi
  rs_normal = 2.0d0/(1.0d0 + (1.0d0-2.0d0*compa)**0.5 - compa)
  m_normal = compa*rs_normal

  call interpg_ini_iso(nlines-1,rg_iso,psi)
  call interpg_ini_iso(nlines-1,rg_iso,alph)
  call interpg_ini_iso(nlines-1,rg_iso,rhog)
!
  open(3,file='rnsgra_1D.las',status='unknown')
  write(3,'(5i5)') nrg
  do irg = 0, nrf
    write(3,'(6es20.12)') rg(irg), psi(irg), alph(irg), rhog(irg)
  end do
  do irg =nrf+1, nrg
    mo2r = m_normal/2.0d0/rg(irg)
    psi(irg)  = 1.0d0 + mo2r
    alph(irg) = (1.0d0 - mo2r)/(1.0d0 + mo2r)
    rhog(irg) = 0.0d0
    write(3,'(6es20.12)') rg(irg), psi(irg), alph(irg), rhog(irg)
  end do
  write(3,'(6es20.12)') ber, radi
  close(3)
!
  open(4,file='rnsflu_1D.las',status='unknown')
  write(4,'(5i5)') nrf
  do irg = 0, nrf
    write(4,'(6es20.12)') rg(irg), rhog(irg)
  end do
  write(4,'(6es20.12)') ber, radi
  close(4)
!
! --  For plotter
!
  open(3,file='frggraplot.dat',status='unknown')
  do irg = 0, nrg
    write(3,'(6es20.12)') rg(irg), psi(irg), alph(irg), rhog(irg)
  end do
  close(3)
!
!  call printout_plot
!
END PROGRAM Import_isotr_qeos
!______________________________________________
!
function fn_lagint(x,y,v)
  implicit none
  real(8), intent(in) :: x(4),y(4), v
  real(8) :: dx12, dx13, dx14, dx21, dx23, dx24, dx31, dx32, dx34, &
  &          dx41, dx42, dx43, wex1, wex2, wex3, wex4, &
  &          xv1, xv2, xv3, xv4, fn_lagint
!
  dx12 = x(1) - x(2)
  dx13 = x(1) - x(3)
  dx14 = x(1) - x(4)
  dx23 = x(2) - x(3)
  dx24 = x(2) - x(4)
  dx34 = x(3) - x(4)
  dx21 = - dx12
  dx31 = - dx13
  dx32 = - dx23
  dx41 = - dx14
  dx42 = - dx24
  dx43 = - dx34
  xv1 = v - x(1)
  xv2 = v - x(2)
  xv3 = v - x(3)
  xv4 = v - x(4)
  wex1 = xv2*xv3*xv4/(dx12*dx13*dx14)
  wex2 = xv1*xv3*xv4/(dx21*dx23*dx24)
  wex3 = xv1*xv2*xv4/(dx31*dx32*dx34)
  wex4 = xv1*xv2*xv3/(dx41*dx42*dx43)
!
  fn_lagint = wex1*y(1) + wex2*y(2) + wex3*y(3) + wex4*y(4)
!
end function fn_lagint
!
function fn_lagint_2nd(x,y,v)
  implicit none
  real(8), intent(in) :: x(2),y(2), v
  real(8) :: dx12, dx21, wex1, wex2, xv1, xv2, fn_lagint_2nd
!
  dx12 = x(1) - x(2)
  dx21 = - dx12
  xv1 = v - x(1)
  xv2 = v - x(2)
  wex1 = xv2/dx12
  wex2 = xv1/dx21
!
  fn_lagint_2nd = wex1*y(1) + wex2*y(2)
!
end function fn_lagint_2nd

