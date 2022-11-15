subroutine gauge_hiju(conv_gra)
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu, &
  &                          gaugex, gaugey, gaugez
  use interface_poisson_solver
  use interface_grgrad_gridpoint_type0
  use interface_grgrad_midpoint_type0
  use interface_update_grfield
  use make_array_3d
  implicit none
!
  real(8), pointer :: gx(:,:,:), gy(:,:,:), gz(:,:,:)
  real(8), pointer :: divg(:,:,:), &
              &       sougx(:,:,:), sougy(:,:,:), sougz(:,:,:)
  real(8) :: dgxdx, dgxdy, dgxdz, dgydx, dgydy, dgydz, &
  &          dgzdx, dgzdy, dgzdz, &
  &          dhxxux, dhxxuy, dhxxuz, dhxyux, dhxyuy, dhxyuz, &
  &          dhxzux, dhxzuy, dhxzuz, &
  &          dhyxux, dhyxuy, dhyxuz, dhyyux, dhyyuy, dhyyuz, &
  &          dhyzux, dhyzuy, dhyzuz, &
  &          dhzxux, dhzxuy, dhzxuz, dhzyux, dhzyuy, dhzyuz, &
  &          dhzzux, dhzzuy, dhzzuz, &
  &          divgg, dlgxdx, dlgxdy, dlgxdz, dlgydy, dlgydz, dlgzdz, &
  &          hxidiv, hxxuc, hxyuc, hxzuc, &
  &          hyidiv, hyxuc, hyyuc, hyzuc, hzidiv, hzxuc, hzyuc, hzzuc
  real(8) :: dxdivg, dydivg, dzdivg, fac13, fac23, conv_gra
  integer :: ipg, irg, itg
!
  call alloc_array3d(gx,0,nrg,0,ntg,0,npg)
  call alloc_array3d(gy,0,nrg,0,ntg,0,npg)
  call alloc_array3d(gz,0,nrg,0,ntg,0,npg)
  call alloc_array3d(divg,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sougx,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sougy,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sougz,0,nrg,0,ntg,0,npg)
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
!
        sougx(irg,itg,ipg) = hxidiv - fac13*dxdivg
        sougy(irg,itg,ipg) = hyidiv - fac13*dydivg
        sougz(irg,itg,ipg) = hzidiv - fac13*dzdivg
!
      end do
    end do
  end do
!
!$omp parallel sections num_threads(3)
!$omp section
  call poisson_solver(sougx,gx)
!$omp section
  call poisson_solver(sougy,gy)
!$omp section
  call poisson_solver(sougz,gz)
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
! -- For get subrotine
! call gauge_hiju.f90_kerr_schild
!
  deallocate(gx)
  deallocate(gy)
  deallocate(gz)
  deallocate(divg)
  deallocate(sougx)
  deallocate(sougy)
  deallocate(sougz)
end subroutine gauge_hiju
