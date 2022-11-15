subroutine dadbscalar_type2(fnc,dadbfnc,cobj)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use interface_grgrad_4th_gridpoint
  use interface_grgrad_4th_gridpoint_bhex
  use interface_grgrad_midpoint_r2nd
  use interface_grgrad_midpoint_r3rd
  use make_array_3d
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long), pointer :: dadbfnc(:,:,:,:,:)
  real(long), pointer :: dfdx(:,:,:), dfdy(:,:,:), dfdz(:,:,:)
  real(long), pointer :: d2fdxdx(:,:,:), d2fdxdy(:,:,:), d2fdxdz(:,:,:), &
                       & d2fdydx(:,:,:), d2fdydy(:,:,:), d2fdydz(:,:,:), &
                       & d2fdzdx(:,:,:), d2fdzdy(:,:,:), d2fdzdz(:,:,:)
  real(long) :: dfncdx, dfncdy, dfncdz
  integer    :: irg, itg, ipg
  character(len=2), intent(in) :: cobj
!
  call alloc_array3d(dfdx,0,nrg,0,ntg,0,npg)
  call alloc_array3d(dfdy,0,nrg,0,ntg,0,npg)
  call alloc_array3d(dfdz,0,nrg,0,ntg,0,npg)
  call alloc_array3d(d2fdxdx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(d2fdydx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(d2fdzdx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(d2fdxdy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(d2fdydy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(d2fdzdy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(d2fdxdz,1,nrg,1,ntg,1,npg)
  call alloc_array3d(d2fdydz,1,nrg,1,ntg,1,npg)
  call alloc_array3d(d2fdzdz,1,nrg,1,ntg,1,npg)
!
! --- Compute the dadbf
! --- The gradient is evaluated at the mid points
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        if (cobj.eq.'bh') &
        &  call grgrad_4th_gridpoint_bhex(fnc,dfncdx,dfncdy,dfncdz,irg,itg,ipg)
        if (cobj.eq.'ns') &
        &  call grgrad_4th_gridpoint(fnc,dfncdx,dfncdy,dfncdz,irg,itg,ipg)
        dfdx(irg,itg,ipg) = dfncdx
        dfdy(irg,itg,ipg) = dfncdy
        dfdz(irg,itg,ipg) = dfncdz
      end do
    end do
  end do
!
  call grgrad_midpoint_r2nd(dfdx,d2fdxdx,d2fdxdy,d2fdxdz)
  call grgrad_midpoint_r2nd(dfdy,d2fdydx,d2fdydy,d2fdydz)
  call grgrad_midpoint_r2nd(dfdz,d2fdzdx,d2fdzdy,d2fdzdz)
!p  call grgrad_midpoint_r3rd(dfdx,d2fdxdx,d2fdxdy,d2fdxdz,cobj)
!p  call grgrad_midpoint_r3rd(dfdy,d2fdydx,d2fdydy,d2fdydz,cobj)
!p  call grgrad_midpoint_r3rd(dfdz,d2fdzdx,d2fdzdy,d2fdzdz,cobj)
!
  dadbfnc(1:nrg,1:ntg,1:npg,1,1) = d2fdxdx(1:nrg,1:ntg,1:npg)
  dadbfnc(1:nrg,1:ntg,1:npg,1,2) =(d2fdxdy(1:nrg,1:ntg,1:npg) &
  &                              + d2fdydx(1:nrg,1:ntg,1:npg))*0.5d0
  dadbfnc(1:nrg,1:ntg,1:npg,1,3) =(d2fdxdz(1:nrg,1:ntg,1:npg) &
  &                              + d2fdzdx(1:nrg,1:ntg,1:npg))*0.5d0
  dadbfnc(1:nrg,1:ntg,1:npg,2,1) =(d2fdxdy(1:nrg,1:ntg,1:npg) &
  &                              + d2fdydx(1:nrg,1:ntg,1:npg))*0.5d0
  dadbfnc(1:nrg,1:ntg,1:npg,2,2) = d2fdydy(1:nrg,1:ntg,1:npg)
  dadbfnc(1:nrg,1:ntg,1:npg,2,3) =(d2fdydz(1:nrg,1:ntg,1:npg) &
  &                              + d2fdzdy(1:nrg,1:ntg,1:npg))*0.5d0
  dadbfnc(1:nrg,1:ntg,1:npg,3,1) =(d2fdxdz(1:nrg,1:ntg,1:npg) &
  &                              + d2fdzdx(1:nrg,1:ntg,1:npg))*0.5d0
  dadbfnc(1:nrg,1:ntg,1:npg,3,2) =(d2fdydz(1:nrg,1:ntg,1:npg) &
  &                              + d2fdzdy(1:nrg,1:ntg,1:npg))*0.5d0
  dadbfnc(1:nrg,1:ntg,1:npg,3,3) = d2fdzdz(1:nrg,1:ntg,1:npg)
!
  deallocate(dfdx,dfdy,dfdz)
  deallocate(d2fdxdx,d2fdxdy,d2fdxdz)
  deallocate(d2fdydx,d2fdydy,d2fdydz)
  deallocate(d2fdzdx,d2fdzdy,d2fdzdz)
!
end subroutine dadbscalar_type2
