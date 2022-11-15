subroutine interpolation_cartesian_BHT_WL
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : rg
  use trigonometry_grav_theta, only : sinthg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use def_metric
  use def_metric_hij
  use def_matter, only : emdg, omeg
  use def_metric_cartesian
  use def_metric_hij_cartesian
  use def_matter_cartesian
  use make_array_3d
  use interface_modules_cartesian
  implicit none
  real(long) :: xxxx, yyyy, omew, qq
  real(long), pointer :: vxu(:,:,:), vyu(:,:,:), vzu(:,:,:)
  integer :: ir, it, ip
!
  call alloc_array3d(vxu,0,nrg,0,ntg,0,npg)
  call alloc_array3d(vyu,0,nrg,0,ntg,0,npg)
  call alloc_array3d(vzu,0,nrg,0,ntg,0,npg)

  write(6,*) "Interpolation psi..."
  call interpolation_metric_bh(psi,psica)
  write(6,*) "Interpolation alpha..."
  call interpolation_metric_bh(alph,alphca)
  write(6,*) "Interpolation bx..."
  call interpolation_metric_bh(bvxd,bvxdca)
  write(6,*) "Interpolation by..."
  call interpolation_metric_bh(bvyd,bvydca)
  write(6,*) "Interpolation bz..."
  call interpolation_metric_bh(bvzd,bvzdca)
  write(6,*) "Interpolation hxx..."
  call interpolation_metric_bh(hxxd,hxxdca)
  write(6,*) "Interpolation hxy..."
  call interpolation_metric_bh(hxyd,hxydca)
  write(6,*) "Interpolation hxz..."
  call interpolation_metric_bh(hxzd,hxzdca)
  write(6,*) "Interpolation hyy..."
  call interpolation_metric_bh(hyyd,hyydca)
  write(6,*) "Interpolation hyz..."
  call interpolation_metric_bh(hyzd,hyzdca)
  write(6,*) "Interpolation hzz..."
  call interpolation_metric_bh(hzzd,hzzdca)

  write(6,*) "Interpolation emd..."
  call interpolation_metric_bh(emdg,emdca)
  write(6,*) "Interpolation ome..."
  call interpolation_metric_bh(omeg,omeca)
!
!  call invhij
!  call calc_shift_down2up
!  call interpolation_metric_bh(bvxu,bvxuca)
!  call interpolation_metric_bh(bvyu,bvyuca)
!  call interpolation_metric_bh(bvzu,bvzuca)
!
  do ip = 0, npg
    do it = 0, ntg
      do ir = 0, nrg
        xxxx = rg(ir)*sinthg(it)*cosphig(ip)
        yyyy = rg(ir)*sinthg(it)*sinphig(ip)
        omew = omeg(ir,it,ip)
        qq   = emdg(ir,it,ip)
        if (qq>1.0d-13 .or. omew.ne.0.0d0)  then
          vxu(ir,it,ip) = - omew*yyyy
          vyu(ir,it,ip) =   omew*xxxx
          vzu(ir,it,ip) = 0.0d0
        end if
      end do
    end do
  end do
!
  write(6,*) "Interpolation vx..."
  call interpolation_metric_bh(vxu,vxca)
  write(6,*) "Interpolation vy..."
  call interpolation_metric_bh(vyu,vyca)
  write(6,*) "Interpolation vz..."
  call interpolation_metric_bh(vzu,vzca)
!
  deallocate(vxu)
  deallocate(vyu)
  deallocate(vzu)
!
end subroutine interpolation_cartesian_BHT_WL
