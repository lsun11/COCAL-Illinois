subroutine gauge
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use interface_poisson_solver
  use make_array_3d
  use interface_grgrad_4th_gridpoint
  use interface_grgrad_midpoint_type0
  implicit none
!
  real(8), pointer :: gx(:,:,:), gy(:,:,:), gz(:,:,:)
  real(8), pointer :: divg(:,:,:), &
              &       sougx(:,:,:), sougy(:,:,:), sougz(:,:,:), &
              &       bfn(:,:,:)
  real(8) :: dgxdx, dgxdy, dgxdz, dgydx, dgydy, dgydz, &
  &          dgzdx, dgzdy, dgzdz, &
  &          dhxxdx, dhxxdy, dhxxdz, dhxydx, dhxydy, dhxydz, &
  &          dhxzdx, dhxzdy, dhxzdz, &
  &          dhyxdx, dhyxdy, dhyxdz, dhyydx, dhyydy, dhyydz, &
  &          dhyzdx, dhyzdy, dhyzdz, &
  &          dhzxdx, dhzxdy, dhzxdz, dhzydx, dhzydy, dhzydz, &
  &          dhzzdx, dhzzdy, dhzzdz, &
  &          divgg, dlgxdx, dlgxdy, dlgxdz, dlgydy, dlgydz, dlgzdz, &
  &          fac23, hdhx, hdhy, hdhz, hxidiv, hxxuc, hxyuc, hxzuc, &
  &          hyidiv, hyxuc, hyyuc, hyzuc, hzidiv, hzxuc, hzyuc, hzzuc
  real(8) :: dbdx, dbdy, dbdz
  integer :: ipg, irg, itg
!
  call alloc_array3d(gx,0,nrg,0,ntg,0,npg)
  call alloc_array3d(gy,0,nrg,0,ntg,0,npg)
  call alloc_array3d(gz,0,nrg,0,ntg,0,npg)
  call alloc_array3d(divg,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sougx,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sougy,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sougz,0,nrg,0,ntg,0,npg)
  call alloc_array3d(bfn,0,nrg,0,ntg,0,npg)

! --- Impose extended Dirac gauge condition.
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call grgrad_midpoint_type0(hxxd,dhxxdx,dhxxdy,dhxxdz,irg,itg,ipg)
        call grgrad_midpoint_type0(hxyd,dhxydx,dhxydy,dhxydz,irg,itg,ipg)
        call grgrad_midpoint_type0(hxzd,dhxzdx,dhxzdy,dhxzdz,irg,itg,ipg)
        call grgrad_midpoint_type0(hyyd,dhyydx,dhyydy,dhyydz,irg,itg,ipg)
        call grgrad_midpoint_type0(hyzd,dhyzdx,dhyzdy,dhyzdz,irg,itg,ipg)
        call grgrad_midpoint_type0(hzzd,dhzzdx,dhzzdy,dhzzdz,irg,itg,ipg)
        hxidiv = dhxxdx + dhxydy + dhxzdz
        hyidiv = dhxydx + dhyydy + dhyzdz
        hzidiv = dhxzdx + dhyzdy + dhzzdz
!
! --  higher order
! --  NOT USED
!        call interpo_linear_type0(hxxuc,hxxu,irg,itg,ipg)
!        call interpo_linear_type0(hxyuc,hxyu,irg,itg,ipg)
!        call interpo_linear_type0(hxzuc,hxzu,irg,itg,ipg)
!        call interpo_linear_type0(hyxuc,hyxu,irg,itg,ipg)
!        call interpo_linear_type0(hyyuc,hyyu,irg,itg,ipg)
!        call interpo_linear_type0(hyzuc,hyzu,irg,itg,ipg)
!        call interpo_linear_type0(hzxuc,hzxu,irg,itg,ipg)
!        call interpo_linear_type0(hzyuc,hzyu,irg,itg,ipg)
!        call interpo_linear_type0(hzzuc,hzzu,irg,itg,ipg)
!        hdhx = hxxuc*dhxxdx + hxyuc*dhxxdy + hxzuc*dhxxdz &
!           &     + hyxuc*dhxydx + hyyuc*dhxydy + hyzuc*dhxydz &
!           &     + hzxuc*dhxzdx + hzyuc*dhxzdy + hzzuc*dhxzdz
!        hdhy = hxxuc*dhyxdx + hxyuc*dhyxdy + hxzuc*dhyxdz &
!           &     + hyxuc*dhyydx + hyyuc*dhyydy + hyzuc*dhyydz &
!           &     + hzxuc*dhyzdx + hzyuc*dhyzdy + hzzuc*dhyzdz
!        hdhz = hxxuc*dhzxdx + hxyuc*dhzxdy + hxzuc*dhzxdz &
!           &     + hyxuc*dhzydx + hyyuc*dhzydy + hyzuc*dhzydz &
!           &     + hzxuc*dhzzdx + hzyuc*dhzzdy + hzzuc*dhzzdz
! --
! -- original
        sougx(irg,itg,ipg) = hxidiv
        sougy(irg,itg,ipg) = hyidiv
        sougz(irg,itg,ipg) = hzidiv
!
      end do
    end do
  end do
!
  call poisson_solver(sougx,gx)
  call poisson_solver(sougy,gy)
  call poisson_solver(sougz,gz)
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call grgrad_midpoint_type0(gx,dgxdx,dgxdy,dgxdz,irg,itg,ipg)
        call grgrad_midpoint_type0(gy,dgydx,dgydy,dgydz,irg,itg,ipg)
        call grgrad_midpoint_type0(gz,dgzdx,dgzdy,dgzdz,irg,itg,ipg)
        divg(irg,itg,ipg) = dgxdx + dgydy + dgzdz
!shi      xxx = rg(irg)*sinthg(itg)*cosphig(ipg)
!shi      yyy = rg(irg)*sinthg(itg)*sinphig(ipg)
!shi      zzz = rg(irg)*costhg(itg)
!shi      divg(irg,itg,ipg) = xxx*gx(irg,itg,ipg)
!shi     &                  + yyy*gy(irg,itg,ipg)
!shi     &                  + zzz*gz(irg,itg,ipg)
      end do
    end do
  end do
!
  call poisson_solver(divg,bfn)
!
!shi      do 50 ipg = 0, npg
!shi      do 50 itg = 0, ntg
!shi      do 50 irg = 0, nrg
!shi      xxx = rg(irg)*sinthg(itg)*cosphig(ipg)
!shi      yyy = rg(irg)*sinthg(itg)*sinphig(ipg)
!shi      zzz = rg(irg)*costhg(itg)
!shi      bfn(irg,itg,ipg) = bfn(irg,itg,ipg) 
!shi     &                 - (xxx*gx(irg,itg,ipg)
!shi     &                  + yyy*gy(irg,itg,ipg)
!shi     &                  + zzz*gz(irg,itg,ipg))
!shi   50 continue
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        call grgrad_4th_gridpoint(bfn,dbdx,dbdy,dbdz,irg,itg,ipg)
        gx(irg,itg,ipg) = gx(irg,itg,ipg) -0.25d0*dbdx
        gy(irg,itg,ipg) = gy(irg,itg,ipg) -0.25d0*dbdy
        gz(irg,itg,ipg) = gz(irg,itg,ipg) -0.25d0*dbdz
!shi      gx(irg,itg,ipg) = gx(irg,itg,ipg) + 0.125d0*dbdx
!shi      gy(irg,itg,ipg) = gy(irg,itg,ipg) + 0.125d0*dbdy
!shi      gz(irg,itg,ipg) = gz(irg,itg,ipg) + 0.125d0*dbdz
      end do
    end do
  end do
!
  fac23 = 2.0d0/3.0d0
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        call grgrad_4th_gridpoint(gx,dgxdx,dgxdy,dgxdz,irg,itg,ipg)
        call grgrad_4th_gridpoint(gy,dgydx,dgydy,dgydz,irg,itg,ipg)
        call grgrad_4th_gridpoint(gz,dgzdx,dgzdy,dgzdz,irg,itg,ipg)
!
        divgg  = dgxdx + dgydy + dgzdz
        dlgxdx = dgxdx + dgxdx - fac23*divgg
        dlgxdy = dgxdy + dgydx
        dlgxdz = dgxdz + dgzdx
        dlgydy = dgydy + dgydy - fac23*divgg
        dlgydz = dgydz + dgzdy
        dlgzdz = dgzdz + dgzdz - fac23*divgg
!
        hxxd(irg,itg,ipg) = hxxd(irg,itg,ipg) - dlgxdx
        hxyd(irg,itg,ipg) = hxyd(irg,itg,ipg) - dlgxdy
        hxzd(irg,itg,ipg) = hxzd(irg,itg,ipg) - dlgxdz
        hyyd(irg,itg,ipg) = hyyd(irg,itg,ipg) - dlgydy
        hyzd(irg,itg,ipg) = hyzd(irg,itg,ipg) - dlgydz
        hzzd(irg,itg,ipg) = hzzd(irg,itg,ipg) - dlgzdz
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
  deallocate(bfn)
end subroutine gauge
