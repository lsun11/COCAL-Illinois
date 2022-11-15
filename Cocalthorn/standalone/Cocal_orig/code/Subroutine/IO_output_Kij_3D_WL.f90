subroutine IO_output_Kij_3D_WL
  use phys_constant, only : long
  use grid_parameter, only  :   nrg, ntg, npg
  use def_metric_excurve_grid, only : tfkij_grid, trk_grid
  use def_metric_hij, only : hxxd, hxyd, hxzd, &
     &                       hyyd, hyzd, hzzd
  implicit none
  real(long) :: axx,axy,axz,ayy,ayz,azz, &
      &         kxx,kxy,kxz,kyy,kyz,kzz, &
      &         gxx,gxy,gxz,gyy,gyz,gzz
  real(long) :: oo3, tk
  integer :: irg, itg, ipg
!
! --- Metric potentials.
  oo3 = 1.0d0/3.0d0

  open(13,file='rnsgra_Kij_3D.las',status='unknown')
  write(13,'(5i5)') nrg, ntg, npg
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        axx = tfkij_grid(irg,itg,ipg,1,1)
        axy = tfkij_grid(irg,itg,ipg,1,2)
        axz = tfkij_grid(irg,itg,ipg,1,3)
        ayy = tfkij_grid(irg,itg,ipg,2,2)
        ayz = tfkij_grid(irg,itg,ipg,2,3)
        azz = tfkij_grid(irg,itg,ipg,3,3)

        gxx = 1.0d0 + hxxd(irg,itg,ipg)
        gxy =         hxyd(irg,itg,ipg)
        gxz =         hxzd(irg,itg,ipg)
        gyy = 1.0d0 + hyyd(irg,itg,ipg)
        gyz =         hyzd(irg,itg,ipg)
        gzz = 1.0d0 + hzzd(irg,itg,ipg)

        tk  = trk_grid(irg,itg,ipg)

        kxx = axx + oo3*gxx*tk        
        kxy = axy + oo3*gxy*tk        
        kxz = axz + oo3*gxz*tk        
        kyy = ayy + oo3*gyy*tk        
        kyz = ayz + oo3*gyz*tk        
        kzz = azz + oo3*gzz*tk        

        write(13,'(1p,6e23.15)') kxx,kxy,kxz,kyy,kyz,kzz
      end do
    end do
  end do
  close(13)
!
end subroutine IO_output_Kij_3D_WL
