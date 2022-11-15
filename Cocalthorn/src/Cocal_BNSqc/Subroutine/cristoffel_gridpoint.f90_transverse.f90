subroutine cristoffel_gridpoint
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_cristoffel_grid, only : cri_grid, crid_grid
  use def_gamma_crist_grid, only : gmcrix_grid, gmcriy_grid, gmcriz_grid
  use def_vector_x, only : vec_xg
  use interface_grgrad_4th_gridpoint
  use make_array_3d
  implicit none
  real(8) :: dhyxdx, dhyxdy, dhyxdz, dhzxdx, dhzxdy, dhzxdz, &
  &          dhzydx, dhzydy, dhzydz, &
  &          dhxxdx, dhxxdy, dhxxdz, dhxydx, dhxydy, dhxydz, &
  &          dhxzdx, dhxzdy, dhxzdz, dhyydx, dhyydy, dhyydz, &
  &          dhyzdx, dhyzdy, dhyzdz, dhzzdx, dhzzdy, dhzzdz, &
  &          gmxxu, gmxyu, gmxzu, gmyyu, gmyzu, gmzzu, &
  &          gmyxu, gmzxu, gmzyu
  integer :: ipg, itg, irg
!
! --- Compute Cristoffel symbols,   
!     whose value is assigned on the grid points. 
!     C^a_bc=1/2 gm^ad(pa_b h_dc+...)
! --- for interporation.
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        gmxxu=hxxu(irg,itg,ipg)
        gmxyu=hxyu(irg,itg,ipg)
        gmxzu=hxzu(irg,itg,ipg)
        gmyyu=hyyu(irg,itg,ipg)
        gmyzu=hyzu(irg,itg,ipg)
        gmzzu=hzzu(irg,itg,ipg)
        gmyyu = gmyyu + 1.0d0
        gmxxu = gmxxu + 1.0d0
        gmzzu = gmzzu + 1.0d0
        gmyxu = gmxyu
        gmzxu = gmxzu
        gmzyu = gmyzu
!
        call grgrad_4th_gridpoint(hxxd,dhxxdx,dhxxdy,dhxxdz,irg,itg,ipg)
        call grgrad_4th_gridpoint(hxyd,dhxydx,dhxydy,dhxydz,irg,itg,ipg)
        call grgrad_4th_gridpoint(hxzd,dhxzdx,dhxzdy,dhxzdz,irg,itg,ipg)
        call grgrad_4th_gridpoint(hyyd,dhyydx,dhyydy,dhyydz,irg,itg,ipg)
        call grgrad_4th_gridpoint(hyzd,dhyzdx,dhyzdy,dhyzdz,irg,itg,ipg)
        call grgrad_4th_gridpoint(hzzd,dhzzdx,dhzzdy,dhzzdz,irg,itg,ipg)
!
        dhyxdx = dhxydx
        dhyxdy = dhxydy
        dhyxdz = dhxydz
!
        dhzxdx = dhxzdx
        dhzxdy = dhxzdy
        dhzxdz = dhxzdz
!
        dhzydx = dhyzdx
        dhzydy = dhyzdy
        dhzydz = dhyzdz
!
! --- Impose Dirac gauge.
!
!Dirac      dhx = gmxxu*dhxxdx + gmxyu*dhxxdy + gmxzu*dhxxdz
!Dirac     &    + gmyxu*dhyxdx + gmyyu*dhyxdy + gmyzu*dhyxdz
!Dirac     &    + gmzxu*dhzxdx + gmzyu*dhzxdy + gmzzu*dhzxdz
!Dirac      dhy = gmxxu*dhxydx + gmxyu*dhxydy + gmxzu*dhxydz
!Dirac     &    + gmyxu*dhyydx + gmyyu*dhyydy + gmyzu*dhyydz
!Dirac     &    + gmzxu*dhzydx + gmzyu*dhzydy + gmzzu*dhzydz
!Dirac      dhz = gmxxu*dhxzdx + gmxyu*dhxzdy + gmxzu*dhxzdz
!Dirac     &    + gmyxu*dhyzdx + gmyyu*dhyzdy + gmyzu*dhyzdz
!Dirac     &    + gmzxu*dhzzdx + gmzyu*dhzzdy + gmzzu*dhzzdz
!Dirac      gzzinv = 1.0d0/gmzzu
!Dirac      dhxzdz = - gzzinv*dhx + dhxzdz
!Dirac      dhyzdz = - gzzinv*dhy + dhyzdz
!Dirac      dhzzdz = - gzzinv*dhz + dhzzdz
!Dirac      dhzxdz = dhxzdz
!Dirac      dhzydz = dhyzdz
!
! --- end
!
        crid_grid(irg,itg,ipg,1,1) = 0.5d0*(dhxxdx + dhxxdx - dhxxdx)
        crid_grid(irg,itg,ipg,1,2) = 0.5d0*(dhxxdy + dhxydx - dhxydx)
        crid_grid(irg,itg,ipg,1,3) = 0.5d0*(dhxxdz + dhxzdx - dhxzdx)
        crid_grid(irg,itg,ipg,1,4) = 0.5d0*(dhxydy + dhxydy - dhyydx)
        crid_grid(irg,itg,ipg,1,5) = 0.5d0*(dhxydz + dhxzdy - dhyzdx)
        crid_grid(irg,itg,ipg,1,6) = 0.5d0*(dhxzdz + dhxzdz - dhzzdx)
        crid_grid(irg,itg,ipg,2,1) = 0.5d0*(dhxydx + dhxydx - dhxxdy)
        crid_grid(irg,itg,ipg,2,2) = 0.5d0*(dhxydy + dhyydx - dhxydy)
        crid_grid(irg,itg,ipg,2,3) = 0.5d0*(dhxydz + dhyzdx - dhxzdy)
        crid_grid(irg,itg,ipg,2,4) = 0.5d0*(dhyydy + dhyydy - dhyydy)
        crid_grid(irg,itg,ipg,2,5) = 0.5d0*(dhyydz + dhyzdy - dhyzdy)
        crid_grid(irg,itg,ipg,2,6) = 0.5d0*(dhyzdz + dhyzdz - dhzzdy)
        crid_grid(irg,itg,ipg,3,1) = 0.5d0*(dhxzdx + dhxzdx - dhxxdz)
        crid_grid(irg,itg,ipg,3,2) = 0.5d0*(dhxzdy + dhyzdx - dhxydz)
        crid_grid(irg,itg,ipg,3,3) = 0.5d0*(dhxzdz + dhzzdx - dhxzdz)
        crid_grid(irg,itg,ipg,3,4) = 0.5d0*(dhyzdy + dhyzdy - dhyydz)
        crid_grid(irg,itg,ipg,3,5) = 0.5d0*(dhyzdz + dhzzdy - dhyzdz)
        crid_grid(irg,itg,ipg,3,6) = 0.5d0*(dhzzdz + dhzzdz - dhzzdz)
!
        cri_grid(irg,itg,ipg,1,1) = &
           &      (-gmxzu*dhxxdz - gmxyu*dhxxdy + gmxxu*dhxxdx + &
           &  2.0d0*gmxyu*dhxydx + 2.0d0*gmxzu*dhxzdx)*0.5d0
        cri_grid(irg,itg,ipg,1,2) = &
           &      (-gmxzu*dhxydz + gmxxu*dhxxdy + gmxzu*dhxzdy + &
           &        gmxyu*dhyydx + gmxzu*dhyzdx)*0.5d0
        cri_grid(irg,itg,ipg,1,3) = &
           &       (gmxxu*dhxxdz + gmxyu*dhxydz - gmxyu*dhxzdy + &
           &        gmxyu*dhyzdx + gmxzu*dhzzdx)*0.5d0
        cri_grid(irg,itg,ipg,1,4) = &
           &      (-gmxzu*dhyydz + 2.0d0*gmxxu*dhxydy + gmxyu*dhyydy + &
           &  2.0d0*gmxzu*dhyzdy - gmxxu*dhyydx)*0.5d0
        cri_grid(irg,itg,ipg,1,5) = &
           &       (gmxxu*dhxydz + gmxyu*dhyydz + gmxxu*dhxzdy + &
           &        gmxzu*dhzzdy - gmxxu*dhyzdx)*0.5d0
        cri_grid(irg,itg,ipg,1,6) = &
           & (2.0d0*gmxxu*dhxzdz + 2.0d0*gmxyu*dhyzdz + &
           &        gmxzu*dhzzdz - gmxyu*dhzzdy - gmxxu*dhzzdx)*0.5d0
!
        cri_grid(irg,itg,ipg,2,1) = &
           &      (-gmyzu*dhxxdz - gmyyu*dhxxdy + gmxyu*dhxxdx + &
           &  2.0d0*gmyyu*dhxydx + 2.0d0*gmyzu*dhxzdx)*0.5d0
        cri_grid(irg,itg,ipg,2,2) = &
           &      (-gmyzu*dhxydz + gmxyu*dhxxdy + gmyzu*dhxzdy + &
           &        gmyyu*dhyydx + gmyzu*dhyzdx)*0.5d0
        cri_grid(irg,itg,ipg,2,3) = &
           &       (gmxyu*dhxxdz + gmyyu*dhxydz - gmyyu*dhxzdy + &
           &        gmyyu*dhyzdx + gmyzu*dhzzdx)*0.5d0
        cri_grid(irg,itg,ipg,2,4) = &
           &      (-gmyzu*dhyydz + 2.0d0*gmxyu*dhxydy + gmyyu*dhyydy + &
           &  2.0d0*gmyzu*dhyzdy - gmxyu*dhyydx)*0.5d0
        cri_grid(irg,itg,ipg,2,5) = &
           &       (gmxyu*dhxydz + gmyyu*dhyydz + gmxyu*dhxzdy + &
           &        gmyzu*dhzzdy - gmxyu*dhyzdx)*0.5d0
        cri_grid(irg,itg,ipg,2,6) = &
           & (2.0d0*gmxyu*dhxzdz + 2.0d0*gmyyu*dhyzdz + &
           &        gmyzu*dhzzdz - gmyyu*dhzzdy - gmxyu*dhzzdx)*0.5d0
!
        cri_grid(irg,itg,ipg,3,1) = &
           &      (-gmzzu*dhxxdz - gmyzu*dhxxdy + gmxzu*dhxxdx + &
           &  2.0d0*gmyzu*dhxydx + 2.0d0*gmzzu*dhxzdx)*0.5d0
        cri_grid(irg,itg,ipg,3,2) = &
           &      (-gmzzu*dhxydz + gmxzu*dhxxdy + gmzzu*dhxzdy + &
           &        gmyzu*dhyydx + gmzzu*dhyzdx)*0.5d0
        cri_grid(irg,itg,ipg,3,3) = &
           &       (gmxzu*dhxxdz + gmyzu*dhxydz - gmyzu*dhxzdy + &
           &        gmyzu*dhyzdx + gmzzu*dhzzdx)*0.5d0
        cri_grid(irg,itg,ipg,3,4) = &
           &      (-gmzzu*dhyydz + 2.0d0*gmxzu*dhxydy + gmyzu*dhyydy + &
           &  2.0d0*gmzzu*dhyzdy - gmxzu*dhyydx)*0.5d0
        cri_grid(irg,itg,ipg,3,5) = &
           &       (gmxzu*dhxydz + gmyzu*dhyydz + gmxzu*dhxzdy + &
           &        gmzzu*dhzzdy - gmxzu*dhyzdx)*0.5d0
        cri_grid(irg,itg,ipg,3,6) = &
           & (2.0d0*gmxzu*dhxzdz + 2.0d0*gmyzu*dhyzdz + &
           &        gmzzu*dhzzdz - gmyzu*dhzzdy - gmxzu*dhzzdx)*0.5d0
!
! --  gamma^ab C^c_ab
!
        gmcrix_grid(irg,itg,ipg) = gmxxu*cri_grid(irg,itg,ipg,1,1)       &
           &                     + gmxyu*cri_grid(irg,itg,ipg,1,2)*2.0d0 &
           &                     + gmxzu*cri_grid(irg,itg,ipg,1,3)*2.0d0 &
           &                     + gmyyu*cri_grid(irg,itg,ipg,1,4)       &
           &                     + gmyzu*cri_grid(irg,itg,ipg,1,5)*2.0d0 &
           &                     + gmzzu*cri_grid(irg,itg,ipg,1,6)
        gmcriy_grid(irg,itg,ipg) = gmxxu*cri_grid(irg,itg,ipg,2,1)       &
           &                     + gmxyu*cri_grid(irg,itg,ipg,2,2)*2.0d0 &
           &                     + gmxzu*cri_grid(irg,itg,ipg,2,3)*2.0d0 &
           &                     + gmyyu*cri_grid(irg,itg,ipg,2,4)       &
           &                     + gmyzu*cri_grid(irg,itg,ipg,2,5)*2.0d0 &
           &                     + gmzzu*cri_grid(irg,itg,ipg,2,6)
        gmcriz_grid(irg,itg,ipg) = gmxxu*cri_grid(irg,itg,ipg,3,1)       &
           &                     + gmxyu*cri_grid(irg,itg,ipg,3,2)*2.0d0 &
           &                     + gmxzu*cri_grid(irg,itg,ipg,3,3)*2.0d0 &
           &                     + gmyyu*cri_grid(irg,itg,ipg,3,4)       &
           &                     + gmyzu*cri_grid(irg,itg,ipg,3,5)*2.0d0 &
           &                     + gmzzu*cri_grid(irg,itg,ipg,3,6)
!
! --- Impose Dirac gauge.
!Dirac        gmcrix_grid(irg,itg,ipg) = 0.0d0
!Dirac        gmcriy_grid(irg,itg,ipg) = 0.0d0
!Dirac        gmcriz_grid(irg,itg,ipg) = 0.0d0
! --- end
      end do
    end do
  end do
!
!irg = 2; itg = ntg/2-2; ipg = 2
!write(6,'(1p,3e16.8)') vec_xg(irg,itg,ipg,1), vec_xg(irg,itg,ipg,2), &
!&                      vec_xg(irg,itg,ipg,3)
!write(6,'(1p,3e16.8)') gmcrix_grid(irg,itg,ipg), gmcriy_grid(irg,itg,ipg), &
!&                      gmcriz_grid(irg,itg,ipg)
!write(6,'(1p,3e16.8)') cri_grid(irg,itg,ipg,1,1), cri_grid(irg,itg,ipg,2,1), &
!&                      cri_grid(irg,itg,ipg,3,1)
!write(6,'(1p,3e16.8)') cri_grid(irg,itg,ipg,1,2), cri_grid(irg,itg,ipg,2,2), &
!&                      cri_grid(irg,itg,ipg,3,2)
!write(6,'(1p,3e16.8)') cri_grid(irg,itg,ipg,1,3), cri_grid(irg,itg,ipg,2,3), &
!&                      cri_grid(irg,itg,ipg,3,3)
!write(6,'(1p,3e16.8)') cri_grid(irg,itg,ipg,1,4), cri_grid(irg,itg,ipg,2,4), &
!&                      cri_grid(irg,itg,ipg,3,4)
!write(6,'(1p,3e16.8)') cri_grid(irg,itg,ipg,1,5), cri_grid(irg,itg,ipg,2,5), &
!&                      cri_grid(irg,itg,ipg,3,5)
!write(6,'(1p,3e16.8)') cri_grid(irg,itg,ipg,1,6), cri_grid(irg,itg,ipg,2,6), &
!&                      cri_grid(irg,itg,ipg,3,6)
!stop                                                                          
!gauge      write(6,*) 'gmcri',gmcrix(8,7,12), gmcriy(8,7,12), gmcriz(8,7,12)
!
end subroutine cristoffel_gridpoint
