subroutine IO_output_solution_IDall_3D
  use def_matter, only    : rhog, utg, uxg, uyg, uzg
  use def_metric, only    : psi, alph, bvxd, bvyd, bvzd
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig
  use grid_parameter, only : nrg, ntg, npg
  implicit none
  integer :: ir, it, ip
!
! --- 
  open(12,file='ID/rho.dat',status='unknown')
  do ip = 0, npg; do it = 0, ntg; do ir = 0, nrg
    write(12,'(1p,7e20.12)') rhog(ir,it,ip)
  end do; end do; end do
  close(12)
!
  open(12,file='ID/ut.dat',status='unknown')
  do ip = 0, npg; do it = 0, ntg; do ir = 0, nrg
    write(12,'(1p,7e20.12)') utg(ir,it,ip)
  end do; end do; end do
  close(12)
!
  open(12,file='ID/ux.dat',status='unknown')
  do ip = 0, npg; do it = 0, ntg; do ir = 0, nrg
    write(12,'(1p,7e20.12)') uxg(ir,it,ip)
  end do; end do; end do
  close(12)
!
  open(12,file='ID/uy.dat',status='unknown')
  do ip = 0, npg; do it = 0, ntg; do ir = 0, nrg
    write(12,'(1p,7e20.12)') uyg(ir,it,ip)
  end do; end do; end do
  close(12)
!
  open(12,file='ID/uz.dat',status='unknown')
  do ip = 0, npg; do it = 0, ntg; do ir = 0, nrg
    write(12,'(1p,7e20.12)') uzg(ir,it,ip)
  end do; end do; end do
  close(12)
!
  open(12,file='ID/psi.dat',status='unknown')
  do ip = 0, npg; do it = 0, ntg; do ir = 0, nrg
    write(12,'(1p,7e20.12)') psi(ir,it,ip)
  end do; end do; end do
  close(12)
!
  open(12,file='ID/alpha.dat',status='unknown')
  do ip = 0, npg; do it = 0, ntg; do ir = 0, nrg
    write(12,'(1p,7e20.12)') alph(ir,it,ip)
  end do; end do; end do
  close(12)
!
  open(12,file='ID/betax.dat',status='unknown')
  do ip = 0, npg; do it = 0, ntg; do ir = 0, nrg
    write(12,'(1p,7e20.12)') bvxd(ir,it,ip)
  end do; end do; end do
  close(12)
!
  open(12,file='ID/betay.dat',status='unknown')
  do ip = 0, npg; do it = 0, ntg; do ir = 0, nrg
    write(12,'(1p,7e20.12)') bvyd(ir,it,ip)
  end do; end do; end do
  close(12)
!
  open(12,file='ID/betaz.dat',status='unknown')
  do ip = 0, npg; do it = 0, ntg; do ir = 0, nrg
    write(12,'(1p,7e20.12)') bvzd(ir,it,ip)
  end do; end do; end do
  close(12)
!
  open(12,file='ID/SphericalCoordinates.dat',status='unknown')
  do ip = 0, npg; do it = 0, ntg; do ir = 0, nrg
    write(12,'(1p,7e20.12)') rg(ir), thg(it), phig(ip)
  end do; end do; end do
  close(12)
!
end subroutine IO_output_solution_IDall_3D
