subroutine gauge_hiju_bhex(conv_gra, iter_count)
  use phys_constant, only :  long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu, &
  &                          gaugex, gaugey, gaugez
  use def_bh_parameter, only : bh_bctype, bh_soltype, ome_bh
  use def_metric_dihiju,  only : dihixu, dihiyu, dihizu
  use def_formulation, only : chgra
  use interface_poisson_solver_1bh
  use interface_bh_boundary_WL
  use interface_outer_boundary_WL
!  use interface_bc_kerr_schild
  use interface_sourceterm_surface_int
  use interface_grgrad_4th_gridpoint_bhex
  use interface_grgrad_gridpoint_type0
  use interface_grgrad_midpoint_type0
  use interface_update_grfield_inAH
  use interface_update_grfield
  use interface_IO_output_1D_general
  use make_array_3d
  use make_array_2d
  implicit none
!
  real(long), pointer :: gx(:,:,:), gy(:,:,:), gz(:,:,:)
  real(long), pointer :: divg(:,:,:), &
              &       sougx(:,:,:), sougy(:,:,:), sougz(:,:,:)
  real(long), pointer :: sou_bhsurf(:,:), dsou_bhsurf(:,:)
  real(long), pointer :: sou_outsurf(:,:), dsou_outsurf(:,:)
  real(long) :: dgxdx, dgxdy, dgxdz, dgydx, dgydy, dgydz, &
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
  real(long) :: dxdivg, dydivg, dzdivg, fac13, fac23, conv_gra
  integer :: ipg, irg, itg, iter_count
  character(len=2) :: chgreen, chpa, chpB
  character(30) :: char1, char2, char3, char4, char5
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
        call grgrad_4th_gridpoint_bhex(gaugex,dgxdx,dgxdy,dgxdz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(gaugey,dgydx,dgydy,dgydz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(gaugez,dgzdx,dgzdy,dgzdz,irg,itg,ipg)
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
!        if (chgra == 'k') then
          sougx(irg,itg,ipg) = sougx(irg,itg,ipg) - dihixu(irg,itg,ipg)    
          sougy(irg,itg,ipg) = sougy(irg,itg,ipg) - dihiyu(irg,itg,ipg)    
          sougz(irg,itg,ipg) = sougz(irg,itg,ipg) - dihizu(irg,itg,ipg)    
!          sougx(irg,itg,ipg) = dihixu(irg,itg,ipg) - fac13*dxdivg
!          sougy(irg,itg,ipg) = dihiyu(irg,itg,ipg) - fac13*dydivg    
!          sougz(irg,itg,ipg) = dihizu(irg,itg,ipg) - fac13*dzdivg   
!        end if
!
      end do
    end do
  end do
!
!  char2 = adjustl("sougx.txt")
!  call IO_output_1D_general(char2, 'g', 'm', sougx, -1, ntg/2, 1)
!  char2 = adjustl("sougy.txt")
!  call IO_output_1D_general(char2, 'g', 'm', sougy, -1, ntg/2, 1)
!  char2 = adjustl("sougz.txt")
!  call IO_output_1D_general(char2, 'g', 'm', sougz, -1, ntg/2, 1)
!
!  char2 = adjustl("dihixu.txt")
!  call IO_output_1D_general(char2, 'g', 'm', dihixu, -1, ntg/2, 1)
!  char2 = adjustl("dihiyu.txt")
!  call IO_output_1D_general(char2, 'g', 'm', dihiyu, -1, ntg/2, 1)
!  char2 = adjustl("dihizu.txt")
!  call IO_output_1D_general(char2, 'g', 'm', dihizu, -1, ntg/2, 1)


!$omp parallel sections num_threads(3)
!$omp section
  sou_bhsurf(0:ntg,0:npg) = 0.0d0 ; dsou_bhsurf(0:ntg,0:npg) = 0.0d0
  sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
  chgreen = 'dd'

  call poisson_solver_1bh(chgreen,sougx, &
  &                       sou_bhsurf, dsou_bhsurf, &
  &                       sou_outsurf,dsou_outsurf,gx)

  call poisson_solver_1bh(chgreen,sougy, &
  &                       sou_bhsurf, dsou_bhsurf, &
  &                       sou_outsurf,dsou_outsurf,gy)

  call poisson_solver_1bh(chgreen,sougz, &
  &                       sou_bhsurf, dsou_bhsurf, &
  &                       sou_outsurf,dsou_outsurf,gz)

  if (mod(iter_count,100)==0.or.iter_count==0) then
    write(char1, '(i5)') iter_count
    char2 = adjustl(char1)
    char3 = 'gx' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', gx, -1, ntg/2, 0)
    char3 = 'gy' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', gy, -1, ntg/2, 0)
    char3 = 'gz' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', gz, -1, ntg/2, 0)
  endif

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
        call grgrad_4th_gridpoint_bhex(gaugex,dgxdx,dgxdy,dgxdz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(gaugey,dgydx,dgydy,dgydz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(gaugez,dgzdx,dgzdy,dgzdz,irg,itg,ipg)
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
!
end subroutine gauge_hiju_bhex
