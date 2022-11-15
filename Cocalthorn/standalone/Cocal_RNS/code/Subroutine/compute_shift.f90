subroutine compute_shift(potx, poty, potz, gvec, bfnc)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use trigonometry_grav_theta, only : sinthg,  costhg
  use trigonometry_grav_phi,   only : sinphig, cosphig
  use coordinate_grav_r, only : rg
  use def_vector_x, only : vec_xg
  use make_array_3d
  use interface_grgrad_4th_gridpoint
  implicit none
  real(long), pointer :: potx(:,:,:), poty(:,:,:), potz(:,:,:)
  real(long), pointer :: bfnc(:,:,:), gvec(:,:,:,:)
  real(long), pointer :: xdotg(:,:,:), bminusxdotg(:,:,:)
  real(long) :: gradx, grady, gradz
  real(long) :: xxx, yyy, zzz 
  integer :: irg, itg, ipg
!
  call alloc_array3d(xdotg,0,nrg,0,ntg,0,npg)
  call alloc_array3d(bminusxdotg,0,nrg,0,ntg,0,npg)
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        xxx = vec_xg(irg,itg,ipg,1)
        yyy = vec_xg(irg,itg,ipg,2)
        zzz = vec_xg(irg,itg,ipg,3)
        xdotg(irg,itg,ipg) = xxx*gvec(irg,itg,ipg,1)  &
      &                    + yyy*gvec(irg,itg,ipg,2)  &
      &                    + zzz*gvec(irg,itg,ipg,3)
        bminusxdotg(irg,itg,ipg) = bfnc(irg,itg,ipg) - xdotg(irg,itg,ipg)
      end do
    end do
  end do
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        call grgrad_4th_gridpoint(bminusxdotg,gradx,grady,gradz,irg,itg,ipg)
        potx(irg,itg,ipg) = gvec(irg,itg,ipg,1) + 0.125d0*gradx
        poty(irg,itg,ipg) = gvec(irg,itg,ipg,2) + 0.125d0*grady
        potz(irg,itg,ipg) = gvec(irg,itg,ipg,3) + 0.125d0*gradz
      end do
    end do
  end do
!
deallocate(xdotg)
deallocate(bminusxdotg)
!
end subroutine compute_shift
