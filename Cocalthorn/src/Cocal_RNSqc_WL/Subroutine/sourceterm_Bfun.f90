subroutine sourceterm_Bfun(sou,Bfun,potx,poty,potz)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric
  use make_array_3d
  use interface_grgrad_midpoint
  use interface_grgrad_midpoint_r3rd_type0
  use interface_grgrad_midpoint_r4th_type0
  implicit none
  integer :: irg, itg, ipg
  real(long) :: gradB2, divG
  real(long) :: dfdx, dfdy, dfdz
  real(long), pointer :: sou(:,:,:), Bfun(:,:,:), potx(:,:,:), poty(:,:,:), potz(:,:,:) 
  real(long),pointer :: dBdx(:,:,:), dBdy(:,:,:), dBdz(:,:,:)
  real(long),pointer :: dpotxdx(:,:,:), dpotxdy(:,:,:), dpotxdz(:,:,:)
  real(long),pointer :: dpotydx(:,:,:), dpotydy(:,:,:), dpotydz(:,:,:)
  real(long),pointer :: dpotzdx(:,:,:), dpotzdy(:,:,:), dpotzdz(:,:,:)
!
  call alloc_array3d(dpotxdx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dpotxdy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dpotxdz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dpotydx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dpotydy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dpotydz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dpotzdx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dpotzdy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dpotzdz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dBdx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dBdy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dBdz, 0, nrg, 0, ntg, 0, npg)
!
!  call grgrad_midpoint(potx,dpotxdx,dpotxdy,dpotxdz)
!  call grgrad_midpoint(poty,dpotydx,dpotydy,dpotydz)
!  call grgrad_midpoint(potz,dpotzdx,dpotzdy,dpotzdz)
!  call grgrad_midpoint(Bfun,dBdx,dBdy,dBdz)
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call grgrad_midpoint_r4th_type0(potx,dfdx,dfdy,dfdz,irg,itg,ipg,'bh')
        dpotxdx(irg,itg,ipg) = dfdx
        dpotxdy(irg,itg,ipg) = dfdy
        dpotxdz(irg,itg,ipg) = dfdz
        call grgrad_midpoint_r4th_type0(poty,dfdx,dfdy,dfdz,irg,itg,ipg,'bh')
        dpotydx(irg,itg,ipg) = dfdx
        dpotydy(irg,itg,ipg) = dfdy
        dpotydz(irg,itg,ipg) = dfdz
        call grgrad_midpoint_r4th_type0(potz,dfdx,dfdy,dfdz,irg,itg,ipg,'bh')
        dpotzdx(irg,itg,ipg) = dfdx
        dpotzdy(irg,itg,ipg) = dfdy
        dpotzdz(irg,itg,ipg) = dfdz
        call grgrad_midpoint_r4th_type0(Bfun,dfdx,dfdy,dfdz,irg,itg,ipg,'bh')
        dBdx(irg,itg,ipg) = dfdx
        dBdy(irg,itg,ipg) = dfdy
        dBdz(irg,itg,ipg) = dfdz
!
        sou(irg,itg,ipg) = dpotxdx(irg,itg,ipg) + dpotydy(irg,itg,ipg) + dpotzdz(irg,itg,ipg)
!
!       divG = dpotxdx(irg,itg,ipg) + dpotydy(irg,itg,ipg) + dpotzdz(irg,itg,ipg)
!       gradB2 = dBdx(irg,itg,ipg)**2 + dBdy(irg,itg,ipg)**2 + dBdz(irg,itg,ipg)**2
!       sou(irg,itg,ipg) = gradB2/Bfun(irg,itg,ipg) + Bfun(irg,itg,ipg)*divG
      end do
    end do
  end do 
!
  deallocate(dpotxdx)
  deallocate(dpotxdy)
  deallocate(dpotxdz)
  deallocate(dpotydx)
  deallocate(dpotydy)
  deallocate(dpotydz)
  deallocate(dpotzdx)
  deallocate(dpotzdy)
  deallocate(dpotzdz)
  deallocate(dBdx)
  deallocate(dBdy)
  deallocate(dBdz)
end subroutine sourceterm_Bfun
