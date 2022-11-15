subroutine IO_output_cartesian_contour_EMF
  use phys_constant, only : long
  use def_matter_parameter, only  : radi
  use def_metric_cartesian, only : alphca, bvxuca, bvyuca, bvzuca
  use def_emfield_cartesian, only : vaca, vaxdca, vaydca, vazdca
  use def_faraday_tensor_cartesian, only  : fxd_gridca, fyd_gridca, & 
  &                                         fzd_gridca, fijd_gridca
  use grid_parameter_cartesian, only  : nx, ny, nz, nx_mid
  use coordinate_grav_xyz, only  : x, y, z
  implicit none
  real(long) :: small = 1.0d-20
  real(long) :: va, vaxdc, vaydc, vazdc, Exca, Eyca, Ezca, Bxca, Byca, Bzca
  real(long) :: alph, bxu, byu, bzu, At, Aphi
  integer :: ix, iy, iz
!
  open(14,file='rns_contour_xy_EMF.dat',status='unknown')
  iz = nx_mid
  do iy = 1, ny
    write(14,*) ' '
    do ix = 1, nx
!
      va  =  vaca(ix,iy,iz)
      vaxdc = vaxdca(ix,iy,iz)
      vaydc = vaydca(ix,iy,iz)
      vazdc = vazdca(ix,iy,iz)
      Exca = fxd_gridca(ix,iy,iz)/radi
      Eyca = fyd_gridca(ix,iy,iz)/radi
      Ezca = fzd_gridca(ix,iy,iz)/radi
      Bxca =   fijd_gridca(ix,iy,iz,3)/radi
      Byca = - fijd_gridca(ix,iy,iz,2)/radi
      Bzca =   fijd_gridca(ix,iy,iz,1)/radi
      if (dabs(va).le.small) va = 0.0d0
      if (dabs(vaxdc).le.small) vaxdc = 0.0d0
      if (dabs(vaydc).le.small) vaydc = 0.0d0
      if (dabs(vazdc).le.small) vazdc = 0.0d0
      if (dabs(Exca).le.small) Exca = 0.0d0
      if (dabs(Eyca).le.small) Eyca = 0.0d0
      if (dabs(Ezca).le.small) Ezca = 0.0d0
      if (dabs(Bxca).le.small) Bxca = 0.0d0
      if (dabs(Byca).le.small) Byca = 0.0d0
      if (dabs(Bzca).le.small) Bzca = 0.0d0
      alph = alphca(ix,iy,iz)
      bxu  = bvxuca(ix,iy,iz)
      byu  = bvyuca(ix,iy,iz)
      bzu  = bvzuca(ix,iy,iz)
      At   = - alph*va + vaxdc*bxu + vaydc*byu + vazdc*bzu
      Aphi = - y(iy)*vaxdc + x(ix)*vaydc
!
      write(14,'(20es14.6)') x(ix), y(iy), &
     &  va, vaxdc, vaydc, vazdc, Exca, Eyca, Ezca, Bxca, Byca, Bzca, At, Aphi
!
    end do
  end do
  close(14)
!
  open(14,file='rns_contour_xz_EMF.dat',status='unknown')
  iy = nx_mid
  do iz = 1, nz
    write(14,*) ' '
    do ix = 1, nx
!
      va  =  vaca(ix,iy,iz)
      vaxdc = vaxdca(ix,iy,iz)
      vaydc = vaydca(ix,iy,iz)
      vazdc = vazdca(ix,iy,iz)
      Exca = fxd_gridca(ix,iy,iz)/radi
      Eyca = fyd_gridca(ix,iy,iz)/radi
      Ezca = fzd_gridca(ix,iy,iz)/radi
      Bxca =   fijd_gridca(ix,iy,iz,3)/radi
      Byca = - fijd_gridca(ix,iy,iz,2)/radi
      Bzca =   fijd_gridca(ix,iy,iz,1)/radi
      if (dabs(va).le.small) va = 0.0d0
      if (dabs(vaxdc).le.small) vaxdc = 0.0d0
      if (dabs(vaydc).le.small) vaydc = 0.0d0
      if (dabs(vazdc).le.small) vazdc = 0.0d0
      if (dabs(Exca).le.small) Exca = 0.0d0
      if (dabs(Eyca).le.small) Eyca = 0.0d0
      if (dabs(Ezca).le.small) Ezca = 0.0d0
      if (dabs(Bxca).le.small) Bxca = 0.0d0
      if (dabs(Byca).le.small) Byca = 0.0d0
      if (dabs(Bzca).le.small) Bzca = 0.0d0
      alph = alphca(ix,iy,iz)
      bxu  = bvxuca(ix,iy,iz)
      byu  = bvyuca(ix,iy,iz)
      bzu  = bvzuca(ix,iy,iz)
      At   = - alph*va + vaxdc*bxu + vaydc*byu + vazdc*bzu
      Aphi = - y(iy)*vaxdc + x(ix)*vaydc
!
      write(14,'(20es14.6)') x(ix), z(iz), &
     &  va, vaxdc, vaydc, vazdc, Exca, Eyca, Ezca, Bxca, Byca, Bzca, At, Aphi
!
    end do
  end do
  close(14)
!
  open(14,file='rns_contour_yz_EMF.dat',status='unknown')
  ix = nx_mid
  do iz = 1, nz
    write(14,*) ' '
    do iy = 1, ny
!
      va  =  vaca(ix,iy,iz)
      vaxdc = vaxdca(ix,iy,iz)
      vaydc = vaydca(ix,iy,iz)
      vazdc = vazdca(ix,iy,iz)
      Exca = fxd_gridca(ix,iy,iz)/radi
      Eyca = fyd_gridca(ix,iy,iz)/radi
      Ezca = fzd_gridca(ix,iy,iz)/radi
      Bxca =   fijd_gridca(ix,iy,iz,3)/radi
      Byca = - fijd_gridca(ix,iy,iz,2)/radi
      Bzca =   fijd_gridca(ix,iy,iz,1)/radi
      if (dabs(va).le.small) va = 0.0d0
      if (dabs(vaxdc).le.small) vaxdc = 0.0d0
      if (dabs(vaydc).le.small) vaydc = 0.0d0
      if (dabs(vazdc).le.small) vazdc = 0.0d0
      if (dabs(Exca).le.small) Exca = 0.0d0
      if (dabs(Eyca).le.small) Eyca = 0.0d0
      if (dabs(Ezca).le.small) Ezca = 0.0d0
      if (dabs(Bxca).le.small) Bxca = 0.0d0
      if (dabs(Byca).le.small) Byca = 0.0d0
      if (dabs(Bzca).le.small) Bzca = 0.0d0
      alph = alphca(ix,iy,iz)
      bxu  = bvxuca(ix,iy,iz)
      byu  = bvyuca(ix,iy,iz)
      bzu  = bvzuca(ix,iy,iz)
      At   = - alph*va + vaxdc*bxu + vaydc*byu + vazdc*bzu
      Aphi = - y(iy)*vaxdc + x(ix)*vaydc
!
      write(14,'(20es14.6)') y(iy), z(iz), &
     &  va, vaxdc, vaydc, vazdc, Exca, Eyca, Ezca, Bxca, Byca, Bzca, At, Aphi
!
    end do
  end do
  close(14)
!
end subroutine IO_output_cartesian_contour_EMF
