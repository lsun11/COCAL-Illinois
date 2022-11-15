subroutine gauge_hiju(conv_gra)
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : rg
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu, &
  &                          gaugex, gaugey, gaugez
  use def_vector_x, only : hvec_xg
  use interface_poisson_solver_1bh_homosol
  use interface_grgrad_gridpoint_type0
  use interface_grgrad_midpoint_type0
  use interface_update_grfield
  use make_array_2d
  use make_array_3d
  implicit none
!
  real(long), pointer :: gx(:,:,:), gy(:,:,:), gz(:,:,:)
  real(long), pointer :: sou_bhsurf(:,:), dsou_bhsurf(:,:)
  real(long), pointer :: sou_outsurf(:,:), dsou_outsurf(:,:)
  real(long), pointer :: divg(:,:,:), &
  &                     sougx(:,:,:), sougy(:,:,:), sougz(:,:,:)
  real(long) :: dgxdx, dgxdy, dgxdz, dgydx, dgydy, dgydz, &
  &             dgzdx, dgzdy, dgzdz, &
  &             dhxxux, dhxxuy, dhxxuz, dhxyux, dhxyuy, dhxyuz, &
  &             dhxzux, dhxzuy, dhxzuz, &
  &             dhyxux, dhyxuy, dhyxuz, dhyyux, dhyyuy, dhyyuz, &
  &             dhyzux, dhyzuy, dhyzuz, &
  &             dhzxux, dhzxuy, dhzxuz, dhzyux, dhzyuy, dhzyuz, &
  &             dhzzux, dhzzuy, dhzzuz, &
  &             divgg, dlgxdx, dlgxdy, dlgxdz, dlgydy, dlgydz, dlgzdz, &
  &             hxidiv, hxxuc, hxyuc, hxzuc, &
  &             hyidiv, hyxuc, hyyuc, hyzuc, hzidiv, hzxuc, hzyuc, hzzuc
  real(long) :: dxdivg, dydivg, dzdivg, fac13, fac23, conv_gra
  real(long) :: dgamma(3), x, y, z
  integer :: ipg, irg, itg
!
  call alloc_array3d(gx,0,nrg,0,ntg,0,npg)
  call alloc_array3d(gy,0,nrg,0,ntg,0,npg)
  call alloc_array3d(gz,0,nrg,0,ntg,0,npg)
  call alloc_array3d(divg,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sougx,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sougy,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sougz,0,nrg,0,ntg,0,npg)
  call alloc_array2d(sou_bhsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_bhsurf,0,ntg,0,npg)
  call alloc_array2d(sou_outsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_outsurf,0,ntg,0,npg)
!
! --- Impose extended Dirac gauge condition.
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        call grgrad_gridpoint_type0(gaugex,dgxdx,dgxdy,dgxdz,irg,itg,ipg)
        call grgrad_gridpoint_type0(gaugey,dgydx,dgydy,dgydz,irg,itg,ipg)
        call grgrad_gridpoint_type0(gaugez,dgzdx,dgzdy,dgzdz,irg,itg,ipg)
        divg(irg,itg,ipg) = dgxdx + dgydy + dgzdz
      end do
    end do
  end do
!
  fac13 = 1.0d0/3.0d0
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call grgrad_midpoint_type0(divg,dxdivg,dydivg,dzdivg,irg,itg,ipg)
        call grgrad_midpoint_type0(hxxu,dhxxux,dhxxuy,dhxxuz,irg,itg,ipg)
        call grgrad_midpoint_type0(hxyu,dhxyux,dhxyuy,dhxyuz,irg,itg,ipg)
        call grgrad_midpoint_type0(hxzu,dhxzux,dhxzuy,dhxzuz,irg,itg,ipg)
        call grgrad_midpoint_type0(hyyu,dhyyux,dhyyuy,dhyyuz,irg,itg,ipg)
        call grgrad_midpoint_type0(hyzu,dhyzux,dhyzuy,dhyzuz,irg,itg,ipg)
        call grgrad_midpoint_type0(hzzu,dhzzux,dhzzuy,dhzzuz,irg,itg,ipg)
        hxidiv = dhxxux + dhxyuy + dhxzuz
        hyidiv = dhxyux + dhyyuy + dhyzuz
        hzidiv = dhxzux + dhyzuy + dhzzuz
        x = hvec_xg(irg,itg,ipg,1)
        y = hvec_xg(irg,itg,ipg,2)
        z = hvec_xg(irg,itg,ipg,3)
        call kerr_schild_transverse_part(x,y,z,dgamma)
!
        sougx(irg,itg,ipg) = hxidiv - dgamma(1) - fac13*dxdivg
        sougy(irg,itg,ipg) = hyidiv - dgamma(2) - fac13*dydivg
        sougz(irg,itg,ipg) = hzidiv - dgamma(3) - fac13*dzdivg
!
!if (itg.eq.1.and.ipg.eq.2) then
!write(6,*) hxidiv, - dgamma(1)
!write(6,*) hyidiv, - dgamma(2)
!write(6,*) hzidiv, - dgamma(3)
!end if
      end do
    end do
  end do
!
   sou_bhsurf( 0:ntg,0:npg) = 0.0d0
  dsou_bhsurf( 0:ntg,0:npg) = 0.0d0
   sou_outsurf(0:ntg,0:npg) = 0.0d0
  dsou_outsurf(0:ntg,0:npg) = 0.0d0
!$omp parallel sections num_threads(3)
!$omp section
  call poisson_solver_1bh_homosol('nd',sougx, &
  &                         sou_bhsurf,dsou_bhsurf, &
  &                        sou_outsurf,dsou_outsurf,gx)
!$omp section
  call poisson_solver_1bh_homosol('nd',sougy, &
  &                         sou_bhsurf,dsou_bhsurf, &
  &                        sou_outsurf,dsou_outsurf,gy)
!$omp section
  call poisson_solver_1bh_homosol('nd',sougz, &
  &                         sou_bhsurf,dsou_bhsurf, &
  &                        sou_outsurf,dsou_outsurf,gz)
!$omp end parallel sections
!
  call update_grfield(gx,gaugex,conv_gra)
  call update_grfield(gy,gaugey,conv_gra)
  call update_grfield(gz,gaugez,conv_gra)
!!  call error_metric_type2(gaugex,potx_bak,error_all(1),flag_all(1),'ns')
!!  call error_metric_type2(gaugey,poty_bak,error_all(2),flag_all(2),'ns')
!!  call error_metric_type2(gaugez,potz_bak,error_all(3),flag_all(3),'ns')
!!  call printout_error_metric_combined_3(iter_count,error_all(1), &
!!  &                                   error_all(2),error_all(3))
!
  open(16,file='test_gauge',status='unknown')

  fac23 = 2.0d0/3.0d0
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        call grgrad_gridpoint_type0(gaugex,dgxdx,dgxdy,dgxdz,irg,itg,ipg)
        call grgrad_gridpoint_type0(gaugey,dgydx,dgydy,dgydz,irg,itg,ipg)
        call grgrad_gridpoint_type0(gaugez,dgzdx,dgzdy,dgzdz,irg,itg,ipg)
!
        divgg  = dgxdx + dgydy + dgzdz
        dlgxdx = dgxdx + dgxdx - fac23*divgg
        dlgxdy = dgxdy + dgydx
        dlgxdz = dgxdz + dgzdx
        dlgydy = dgydy + dgydy - fac23*divgg
        dlgydz = dgydz + dgzdy
        dlgzdz = dgzdz + dgzdz - fac23*divgg
!
if(itg.le.2.and.ipg.eq.2) write(16,'(1p,4e15.6)') &
& rg(irg),hxxu(irg,itg,ipg),dlgxdx
!
        hxxu(irg,itg,ipg) = hxxu(irg,itg,ipg) - dlgxdx
        hxyu(irg,itg,ipg) = hxyu(irg,itg,ipg) - dlgxdy
        hxzu(irg,itg,ipg) = hxzu(irg,itg,ipg) - dlgxdz
        hyyu(irg,itg,ipg) = hyyu(irg,itg,ipg) - dlgydy
        hyzu(irg,itg,ipg) = hyzu(irg,itg,ipg) - dlgydz
        hzzu(irg,itg,ipg) = hzzu(irg,itg,ipg) - dlgzdz
!
      end do
    end do
  end do
!
close(16)
!
  deallocate(gx)
  deallocate(gy)
  deallocate(gz)
  deallocate(divg)
  deallocate(sougx)
  deallocate(sougy)
  deallocate(sougz)
  deallocate(sou_bhsurf)
  deallocate(dsou_bhsurf)
  deallocate(sou_outsurf)
  deallocate(dsou_outsurf)
end subroutine gauge_hiju
